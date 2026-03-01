"""Stages 5, 6, 7: Collinearity detection, PC-metadata associations, and ICC."""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal, spearmanr
from sklearn.linear_model import LinearRegression

from .exceptions import BatchDetectiveDataError

logger = logging.getLogger(__name__)


# ─── Stage 5: Collinearity ───────────────────────────────────────────────────

def cramers_v(x: pd.Series, y: pd.Series) -> float:
    """Compute Cramér's V for two categorical variables.

    Args:
        x: First categorical series.
        y: Second categorical series.

    Returns:
        Cramér's V [0, 1].
    """
    contingency = pd.crosstab(x, y)
    r, c = contingency.shape
    if r < 2 or c < 2:
        return 0.0
    chi2_stat, _, _, _ = chi2_contingency(contingency, correction=False)
    n = len(x)
    denom = n * (min(r, c) - 1)
    if denom <= 0:
        return 0.0
    return float(np.sqrt(max(0.0, chi2_stat / denom)))


def eta_squared_anova(continuous: pd.Series, categorical: pd.Series) -> float:
    """Compute eta-squared from one-way ANOVA.

    Args:
        continuous: Continuous variable.
        categorical: Categorical grouping variable.

    Returns:
        Eta-squared [0, 1].
    """
    groups = [continuous[categorical == g].dropna() for g in categorical.unique()]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return 0.0

    overall_mean = continuous.mean()
    ss_between = sum(
        len(g) * (g.mean() - overall_mean) ** 2 for g in groups
    )
    ss_total = ((continuous - overall_mean) ** 2).sum()
    if ss_total == 0:
        return 0.0
    return float(np.clip(ss_between / ss_total, 0, 1))


def detect_collinearity(
    metadata: pd.DataFrame,
    covariate_info: Dict,
) -> List[Dict]:
    """Detect collinearity between all covariate pairs.

    Args:
        metadata: Metadata DataFrame.
        covariate_info: Covariate classification dict.

    Returns:
        List of collinearity warning dicts.
    """
    warnings = []
    valid_covs = [
        name for name, info in covariate_info.items()
        if not info.get("skip")
    ]

    for i, cov_a in enumerate(valid_covs):
        for cov_b in valid_covs[i + 1:]:
            info_a = covariate_info[cov_a]
            info_b = covariate_info[cov_b]
            type_a = info_a.get("cov_type", "unknown")
            type_b = info_b.get("cov_type", "unknown")

            # Common valid mask
            mask = (
                info_a.get("valid_mask", pd.Series(True, index=metadata.index))
                & info_b.get("valid_mask", pd.Series(True, index=metadata.index))
            )
            col_a = metadata.loc[mask, cov_a]
            col_b = metadata.loc[mask, cov_b]

            if len(col_a) < 5:
                continue

            metric = None
            metric_val = 0.0

            if type_a == "categorical" and type_b == "categorical":
                metric = "Cramér's V"
                metric_val = cramers_v(col_a, col_b)
                threshold = 0.70

            elif type_a == "continuous" and type_b == "continuous":
                metric = "|Pearson r|"
                r = col_a.corr(pd.to_numeric(col_b, errors="coerce"))
                metric_val = abs(r) if not np.isnan(r) else 0.0
                threshold = 0.70

            elif type_a == "categorical" and type_b == "continuous":
                metric = "eta²"
                col_b_num = pd.to_numeric(col_b, errors="coerce").dropna()
                col_a_aligned = col_a[col_b_num.index]
                metric_val = eta_squared_anova(col_b_num, col_a_aligned)
                threshold = 0.50

            elif type_a == "continuous" and type_b == "categorical":
                metric = "eta²"
                col_a_num = pd.to_numeric(col_a, errors="coerce").dropna()
                col_b_aligned = col_b[col_a_num.index]
                metric_val = eta_squared_anova(col_a_num, col_b_aligned)
                threshold = 0.50
            else:
                continue

            if metric_val >= threshold:
                warning = {
                    "cov_a": cov_a,
                    "cov_b": cov_b,
                    "metric": metric,
                    "value": metric_val,
                    "mixed": (type_a != type_b),
                }

                # Extra guidance for mixed collinearity
                if type_a == "continuous" and type_b == "categorical":
                    warning["guidance"] = (
                        f"When a continuous variable ({cov_a}) and a categorical "
                        f"variable ({cov_b}) are collinear, the continuous variable "
                        f"often captures more variance because it preserves within-batch "
                        f"temporal gradients. Consider including {cov_a} as a continuous "
                        f"covariate in your DE model (~{cov_a} + treatment) rather than "
                        f"using {cov_b} as a categorical covariate."
                    )
                elif type_a == "categorical" and type_b == "continuous":
                    warning["guidance"] = (
                        f"When a continuous variable ({cov_b}) and a categorical "
                        f"variable ({cov_a}) are collinear, the continuous variable "
                        f"often captures more variance because it preserves within-batch "
                        f"temporal gradients. Consider including {cov_b} as a continuous "
                        f"covariate in your DE model (~{cov_b} + treatment) rather than "
                        f"using {cov_a} as a categorical covariate."
                    )

                warnings.append(warning)
                logger.warning(
                    f"⚠️ COLLINEARITY WARNING: {cov_a} and {cov_b} are strongly "
                    f"collinear ({metric} = {metric_val:.2f})."
                )

    return warnings


# ─── Stage 6: PC-Metadata Associations ───────────────────────────────────────

def kruskal_eta_squared(H: float, k: int, n: int) -> float:
    """Rank-based eta-squared for Kruskal-Wallis test.

    Args:
        H: Kruskal-Wallis H statistic.
        k: Number of groups.
        n: Total sample size.

    Returns:
        Eta-squared clamped to [0, 1].
    """
    eta_sq = (H - k + 1) / (n - k)
    return float(np.clip(eta_sq, 0.0, 1.0))


def bootstrap_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_bootstrap: int = 1000,
    alpha: float = 0.05,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[float, float]:
    """Compute bootstrap 95% CI for Spearman rho.

    Args:
        x: First variable.
        y: Second variable.
        n_bootstrap: Number of bootstrap samples.
        alpha: Significance level.
        rng: Optional numpy random generator.

    Returns:
        Tuple (ci_lower, ci_upper).
    """
    if rng is None:
        rng = np.random.default_rng()
    n = len(x)
    boot_rhos = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        r, _ = spearmanr(x[idx], y[idx])
        boot_rhos.append(r)
    ci_lower = float(np.percentile(boot_rhos, 100 * alpha / 2))
    ci_upper = float(np.percentile(boot_rhos, 100 * (1 - alpha / 2)))
    return ci_lower, ci_upper


def bh_correct(pvals: List[float]) -> List[float]:
    """Apply Benjamini-Hochberg FDR correction.

    Args:
        pvals: List of p-values.

    Returns:
        List of adjusted p-values (q-values).
    """
    n = len(pvals)
    if n == 0:
        return []
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]
    qvals = np.minimum(1.0, sorted_pvals * n / (np.arange(n) + 1))

    # Enforce monotonicity (take cumulative minimum from right)
    for i in range(n - 2, -1, -1):
        qvals[i] = min(qvals[i], qvals[i + 1])

    result = np.empty(n)
    result[sorted_idx] = qvals
    return result.tolist()


def build_conditioning_design_matrix(
    metadata: pd.DataFrame,
    covariate_names: List[str],
    covariate_info: Dict,
) -> np.ndarray:
    """Build combined design matrix for conditioning covariates.

    Args:
        metadata: Metadata DataFrame.
        covariate_names: Names of covariates to condition on.
        covariate_info: Covariate type info dict.

    Returns:
        Design matrix (n_samples, total_df).
    """
    matrices = []
    for cov in covariate_names:
        if cov not in metadata.columns:
            logger.warning(f"Conditioning covariate {cov} not found in metadata.")
            continue
        info = covariate_info.get(cov, {})
        cov_type = info.get("cov_type", "continuous")
        if cov_type == "categorical":
            dm = pd.get_dummies(metadata[cov], drop_first=True).values.astype(float)
        else:
            dm = metadata[[cov]].values.astype(float)
        matrices.append(dm)

    if not matrices:
        return np.zeros((len(metadata), 1))

    return np.hstack(matrices)


def test_pc_associations(
    pc_scores: np.ndarray,
    explained_variance_ratio: np.ndarray,
    metadata: pd.DataFrame,
    covariate_info: Dict,
    n_pcs_test: int = 5,
    condition_on: Optional[List[str]] = None,
    rng: Optional[np.random.Generator] = None,
) -> pd.DataFrame:
    """Test association between PCs and metadata covariates.

    Args:
        pc_scores: PC scores (n_samples, n_pcs).
        explained_variance_ratio: Variance explained per PC.
        metadata: Metadata DataFrame.
        covariate_info: Covariate classification dict.
        n_pcs_test: Number of PCs to test.
        condition_on: List of covariates to condition on.
        rng: Optional numpy random generator.

    Returns:
        Association results DataFrame.
    """
    n_pcs_avail = pc_scores.shape[1]
    n_pcs_test = min(n_pcs_test, n_pcs_avail)

    # Apply conditioning
    pc_residuals = pc_scores.copy()
    if condition_on:
        conditioning_matrix = build_conditioning_design_matrix(
            metadata, condition_on, covariate_info
        )
        for pc_idx in range(n_pcs_avail):
            lr = LinearRegression().fit(
                conditioning_matrix, pc_scores[:, pc_idx]
            )
            pc_residuals[:, pc_idx] = (
                pc_scores[:, pc_idx] - lr.predict(conditioning_matrix)
            )

    if rng is None:
        rng = np.random.default_rng(42)

    rows = []
    for cov_name, info in covariate_info.items():
        if info.get("skip"):
            continue

        cov_type = info.get("cov_type", "unknown")
        valid_mask = info.get(
            "valid_mask", pd.Series(True, index=metadata.index)
        )
        valid_idx = np.where(valid_mask.values)[0]

        if len(valid_idx) < 5:
            continue

        pc_vals_all = pc_residuals[valid_idx, :]
        cov_values = metadata.iloc[valid_idx][cov_name]

        pvals = []

        for pc_idx in range(n_pcs_test):
            pc_vals = pc_vals_all[:, pc_idx]
            var_explained = (
                explained_variance_ratio[pc_idx]
                if pc_idx < len(explained_variance_ratio)
                else 0.0
            )

            if cov_type == "categorical":
                # Check small groups
                if "small_groups" in info:
                    rows.append({
                        "covariate": cov_name,
                        "covariate_type": cov_type,
                        "label": info.get("label", "unknown"),
                        "pc": pc_idx + 1,
                        "pc_variance_explained": float(var_explained),
                        "effect_size": np.nan,
                        "effect_size_type": "eta2_kruskal",
                        "pval_raw": np.nan,
                        "pval_adjusted_bh": np.nan,
                        "significant_q05": False,
                        "ci_lower": np.nan,
                        "ci_upper": np.nan,
                        "skip_reason": "small group",
                    })
                    pvals.append(np.nan)
                    continue

                groups_data = [
                    pc_vals[cov_values == g]
                    for g in info.get("groups", [])
                    if (cov_values == g).sum() > 0
                ]
                if len(groups_data) < 2:
                    pvals.append(np.nan)
                    continue

                try:
                    H_stat, pval = kruskal(*groups_data)
                except ValueError:
                    pvals.append(np.nan)
                    continue

                k = len(groups_data)
                n = len(pc_vals)
                eta2 = kruskal_eta_squared(H_stat, k, n)

                rows.append({
                    "covariate": cov_name,
                    "covariate_type": cov_type,
                    "label": info.get("label", "unknown"),
                    "pc": pc_idx + 1,
                    "pc_variance_explained": float(var_explained),
                    "effect_size": eta2,
                    "effect_size_type": "eta2_kruskal",
                    "pval_raw": float(pval),
                    "pval_adjusted_bh": np.nan,  # filled below
                    "significant_q05": False,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                })
                pvals.append(float(pval))

            else:
                # Continuous
                try:
                    cov_num = pd.to_numeric(cov_values, errors="coerce")
                    valid2 = ~cov_num.isna()
                    if valid2.sum() < 5:
                        pvals.append(np.nan)
                        continue
                    rho, pval = spearmanr(
                        pc_vals[valid2.values], cov_num[valid2].values
                    )
                    ci_lower, ci_upper = bootstrap_spearman_ci(
                        pc_vals[valid2.values],
                        cov_num[valid2].values,
                        rng=rng,
                    )
                except Exception:
                    pvals.append(np.nan)
                    continue

                rows.append({
                    "covariate": cov_name,
                    "covariate_type": cov_type,
                    "label": info.get("label", "unknown"),
                    "pc": pc_idx + 1,
                    "pc_variance_explained": float(var_explained),
                    "effect_size": float(rho),
                    "effect_size_type": "spearman_rho",
                    "pval_raw": float(pval),
                    "pval_adjusted_bh": np.nan,
                    "significant_q05": False,
                    "ci_lower": ci_lower,
                    "ci_upper": ci_upper,
                })
                pvals.append(float(pval))

        # BH correction within covariate
        valid_pvals = [(i, p) for i, p in enumerate(pvals) if not np.isnan(p)]
        if valid_pvals:
            idxs, pvs = zip(*valid_pvals)
            qvals = bh_correct(list(pvs))
            # Assign back
            row_offset = len(rows) - n_pcs_test
            pc_row_idx = 0
            for row in rows[row_offset:]:
                if row["covariate"] == cov_name and not np.isnan(row.get("pval_raw", np.nan)):
                    if pc_row_idx < len(qvals):
                        row["pval_adjusted_bh"] = qvals[pc_row_idx]
                        row["significant_q05"] = qvals[pc_row_idx] < 0.05
                    pc_row_idx += 1

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "covariate", "covariate_type", "label", "pc", "pc_variance_explained",
        "effect_size", "effect_size_type", "pval_raw", "pval_adjusted_bh",
        "significant_q05", "ci_lower", "ci_upper",
    ])
    return df


# ─── Stage 7: ICC ─────────────────────────────────────────────────────────────

def compute_icc11_vectorized(
    X: np.ndarray,
    group_labels: np.ndarray,
) -> Tuple[np.ndarray, float, Tuple[float, float]]:
    """Vectorized ICC(1,1) for multiple genes simultaneously.

    Uses one-way random effects ANOVA decomposition with harmonic mean
    for unbalanced designs.

    Args:
        X: Expression matrix (n_genes, n_samples).
        group_labels: Group assignments (n_samples,).

    Returns:
        Tuple of (icc_values, icc_median, (ci_lower, ci_upper)).
    """
    unique_groups = np.unique(group_labels)
    n_groups = len(unique_groups)
    n_samples = X.shape[1]

    if n_groups < 2:
        return np.zeros(X.shape[0]), 0.0, (0.0, 0.0)

    group_sizes = np.array([np.sum(group_labels == g) for g in unique_groups])

    # Harmonic mean for unbalanced designs
    k_harmonic = n_groups / np.sum(1.0 / group_sizes)

    # Compute group means (vectorized)
    group_means = np.zeros((X.shape[0], n_groups))
    for j, g in enumerate(unique_groups):
        mask = group_labels == g
        group_means[:, j] = X[:, mask].mean(axis=1)

    # Overall mean
    overall_mean = X.mean(axis=1, keepdims=True)  # (n_genes, 1)

    # SS_between
    ss_between = np.sum(
        group_sizes[np.newaxis, :] * (group_means - overall_mean) ** 2,
        axis=1,
    )

    # SS_within
    ss_within = np.zeros(X.shape[0])
    for j, g in enumerate(unique_groups):
        mask = group_labels == g
        ss_within += np.sum(
            (X[:, mask] - group_means[:, j:j+1]) ** 2, axis=1
        )

    df_between = n_groups - 1
    df_within = n_samples - n_groups

    if df_between <= 0 or df_within <= 0:
        return np.zeros(X.shape[0]), 0.0, (0.0, 0.0)

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within

    denom = ms_between + (k_harmonic - 1) * ms_within
    icc_values = np.where(
        denom > 0,
        (ms_between - ms_within) / denom,
        0.0,
    )
    icc_values = np.clip(icc_values, 0.0, 1.0)

    icc_median = float(np.median(icc_values))
    return icc_values, icc_median, (0.0, 0.0)


def compute_icc_with_bootstrap(
    X: np.ndarray,
    group_labels: np.ndarray,
    n_bootstrap: int = 1000,
    rng: Optional[np.random.Generator] = None,
) -> Dict:
    """Compute ICC with bootstrap CI over samples.

    Args:
        X: Expression matrix (n_genes, n_samples).
        group_labels: Group assignments (n_samples,).
        n_bootstrap: Number of bootstrap iterations.
        rng: Optional random generator.

    Returns:
        Dict with median_icc, iqr_lower, iqr_upper, ci_lower_95, ci_upper_95,
        prop_genes_moderate_plus, icc_tier, icc_values.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    n_samples = X.shape[1]
    icc_values, icc_median, _ = compute_icc11_vectorized(X, group_labels)

    iqr_lower = float(np.percentile(icc_values, 25))
    iqr_upper = float(np.percentile(icc_values, 75))

    # Bootstrap CI over samples
    bootstrap_medians = []
    for _ in range(n_bootstrap):
        sample_idx = rng.integers(0, n_samples, size=n_samples)
        X_boot = X[:, sample_idx]
        labels_boot = group_labels[sample_idx]
        if len(np.unique(labels_boot)) < 2:
            continue
        _, med, _ = compute_icc11_vectorized(X_boot, labels_boot)
        bootstrap_medians.append(med)

    if len(bootstrap_medians) >= 10:
        ci_lower = float(np.percentile(bootstrap_medians, 2.5))
        ci_upper = float(np.percentile(bootstrap_medians, 97.5))
    else:
        ci_lower = max(0.0, icc_median - 0.1)
        ci_upper = min(1.0, icc_median + 0.1)

    # Proportion of genes with moderate or higher ICC
    prop_moderate = float(np.mean(icc_values >= 0.30))

    # Tier
    if icc_median < 0.10:
        tier = "negligible"
    elif icc_median < 0.30:
        tier = "mild"
    elif icc_median <= 0.60:
        tier = "moderate"
    else:
        tier = "strong"

    return {
        "median_icc": icc_median,
        "iqr_lower": iqr_lower,
        "iqr_upper": iqr_upper,
        "ci_lower_95": ci_lower,
        "ci_upper_95": ci_upper,
        "prop_genes_moderate_plus": prop_moderate,
        "icc_tier": tier,
        "icc_values": icc_values,
    }


def compute_all_icc(
    log1p_cpm_selected: pd.DataFrame,
    metadata: pd.DataFrame,
    covariate_info: Dict,
    n_bootstrap: int = 1000,
    rng: Optional[np.random.Generator] = None,
) -> pd.DataFrame:
    """Compute ICC for all categorical covariates.

    Args:
        log1p_cpm_selected: Selected gene expression matrix (genes × samples).
        metadata: Metadata DataFrame.
        covariate_info: Covariate classification dict.
        n_bootstrap: Bootstrap iterations.
        rng: Optional random generator.

    Returns:
        ICC table DataFrame.
    """
    X = log1p_cpm_selected.values.astype(float)
    rows = []

    for cov_name, info in covariate_info.items():
        if info.get("skip") or info.get("cov_type") != "categorical":
            continue
        if "small_groups" in info:
            logger.warning(
                f"Covariate {cov_name}: small groups present, skipping ICC."
            )
            continue

        valid_mask = info.get(
            "valid_mask", pd.Series(True, index=metadata.index)
        )
        valid_idx = np.where(valid_mask.values)[0]

        if len(valid_idx) < 5:
            continue

        X_valid = X[:, valid_idx]
        labels = metadata.iloc[valid_idx][cov_name].values.astype(str)

        result = compute_icc_with_bootstrap(X_valid, labels, n_bootstrap, rng)

        rows.append({
            "covariate": cov_name,
            "label": info.get("label", "unknown"),
            "n_samples_used": len(valid_idx),
            "n_groups": info.get("n_groups", 2),
            "median_icc": result["median_icc"],
            "iqr_lower": result["iqr_lower"],
            "iqr_upper": result["iqr_upper"],
            "ci_lower_95": result["ci_lower_95"],
            "ci_upper_95": result["ci_upper_95"],
            "prop_genes_moderate_plus": result["prop_genes_moderate_plus"],
            "icc_tier": result["icc_tier"],
            "_icc_values": result["icc_values"],
        })

    if rows:
        df = pd.DataFrame(rows)
    else:
        df = pd.DataFrame(columns=[
            "covariate", "label", "n_samples_used", "n_groups",
            "median_icc", "iqr_lower", "iqr_upper", "ci_lower_95", "ci_upper_95",
            "prop_genes_moderate_plus", "icc_tier",
        ])

    return df
