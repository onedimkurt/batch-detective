"""Stage 2: Quality control for batch-detective."""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2, ncx2

from .exceptions import BatchDetectiveDataError

logger = logging.getLogger(__name__)


def estimate_kruskal_power(
    n_samples: int,
    n_groups: int,
    effect_size_eta2: float,
    alpha: float = 0.05,
) -> float:
    """Approximate power for KW test to detect eta2 effect size.

    Args:
        n_samples: Total number of samples.
        n_groups: Number of groups.
        effect_size_eta2: Expected eta-squared effect size.
        alpha: Significance level.

    Returns:
        Power as float [0, 1].
    """
    df = n_groups - 1
    if df <= 0:
        return 0.0
    critical_value = chi2.ppf(1 - alpha, df=df)
    denom = 1 - effect_size_eta2
    if denom <= 0:
        return 1.0
    ncp = n_samples * effect_size_eta2 / denom
    power = 1 - ncx2.cdf(critical_value, df=df, nc=ncp)
    return float(power)


def find_min_detectable_icc(
    n_samples: int,
    n_groups: int,
    target_power: float = 0.8,
    alpha: float = 0.05,
) -> float:
    """Find minimum detectable ICC at given power.

    Args:
        n_samples: Number of samples.
        n_groups: Number of groups.
        target_power: Desired power (default 0.80).
        alpha: Significance level.

    Returns:
        Minimum detectable ICC as float.
    """
    for icc in np.arange(0.05, 1.0, 0.05):
        power = estimate_kruskal_power(n_samples, n_groups, icc, alpha)
        if power >= target_power:
            return float(icc)
    return 1.0


class QualityController:
    """Performs sample-level and gene-level QC.

    Args:
        counts_working: Working copy of count matrix.
        metadata_working: Working copy of metadata.
        min_cpm: Minimum CPM threshold for gene filtering.
        min_samples_expressing: Minimum samples expressing a gene.
        technical_covariates: List of technical covariate names.
        biological_covariates: List of biological covariate names.
    """

    def __init__(
        self,
        counts_working: pd.DataFrame,
        metadata_working: pd.DataFrame,
        min_cpm: float = 1.0,
        min_samples_expressing: Optional[int] = None,
        technical_covariates: Optional[List[str]] = None,
        biological_covariates: Optional[List[str]] = None,
    ):
        self.counts = counts_working.copy()
        self.metadata = metadata_working.copy()
        self.min_cpm = min_cpm
        self.min_samples_expressing = min_samples_expressing
        self.technical_covariates = technical_covariates or []
        self.biological_covariates = biological_covariates or []
        self.qc_summary: Optional[pd.DataFrame] = None
        self.covariate_info: Dict = {}
        self.power_warnings: List[str] = []
        self.repeated_measures_covariates: List[str] = []

    def _classify_covariate(self, series: pd.Series, name: str) -> Dict:
        """Classify a metadata column and check for issues.

        Args:
            series: Metadata column.
            name: Column name.

        Returns:
            Dict with type, n_unique, groups, etc.
        """
        info = {"name": name, "skip": False, "skip_reason": None}

        n_unique = series.nunique(dropna=True)
        n_samples = len(series.dropna())

        # Single value check
        if n_unique < 2:
            info["skip"] = True
            info["skip_reason"] = f"only 1 unique value"
            logger.warning(
                f"Covariate {name} has only 1 unique value — skipped."
            )
            info["cov_type"] = "constant"
            return info

        # Repeated measures detection — only meaningful for categorical-looking columns
        # Skip this check if the column is numeric with many unique values (it's continuous)
        numeric_check = pd.to_numeric(series, errors="coerce")
        is_likely_numeric = numeric_check.notna().mean() > 0.9
        avg_group_size = n_samples / n_unique if n_unique > 0 else 0
        if (not is_likely_numeric
                and n_unique > n_samples * 0.4
                and avg_group_size < 3):
            info["skip"] = True
            info["skip_reason"] = "likely subject/patient ID (repeated measures)"
            info["repeated_measures"] = True
            logger.warning(
                f"⚠️ REPEATED MEASURES WARNING: Covariate {name} appears to be a "
                f"subject identifier (n_unique={n_unique}, avg group size={avg_group_size:.1f}). "
                "batch-detective assumes independent samples. Statistical tests "
                "(Kruskal-Wallis, ICC) are INVALID for repeated measures designs "
                "and will produce misleading results. This covariate will be skipped."
            )
            self.repeated_measures_covariates.append(name)
            return info

        # Type detection
        numeric_series = pd.to_numeric(series, errors="coerce")
        has_non_numeric = numeric_series.isna().any() and not series.isna().any()

        if has_non_numeric or series.dtype == object or series.dtype.name == "category":
            cov_type = "categorical"
        elif n_unique < 10:
            cov_type = "categorical"
        else:
            cov_type = "continuous"

        # High cardinality categorical with non-numeric
        if cov_type == "continuous" and has_non_numeric:
            info["skip"] = True
            info["skip_reason"] = "non-numeric values with too many levels"
            info["cov_type"] = "high_cardinality_categorical"
            logger.warning(
                f"Covariate {name} has non-numeric values and too many unique "
                "levels — skipped. Consider binning this variable."
            )
            return info

        info["cov_type"] = cov_type
        info["n_unique"] = n_unique

        if cov_type == "categorical":
            groups = series.dropna().unique()
            group_sizes = {g: (series == g).sum() for g in groups}
            info["groups"] = sorted(groups, key=str)
            info["group_sizes"] = group_sizes
            info["n_groups"] = len(groups)

            # Small group safeguards
            small_groups = {g: s for g, s in group_sizes.items() if s < 5}
            if small_groups:
                info["small_groups"] = small_groups
                for g, s in small_groups.items():
                    logger.warning(
                        f"Covariate {name}: group {g} has n={s} (<5 required). "
                        "Kruskal-Wallis and ICC will be skipped."
                    )
        else:
            info["range"] = (
                float(numeric_series.min()),
                float(numeric_series.max()),
            )

        # Label type
        if name in self.technical_covariates:
            info["label"] = "technical"
        elif name in self.biological_covariates:
            info["label"] = "biological"
        else:
            info["label"] = "unknown"

        # Missing values
        n_missing = series.isna().sum()
        if n_missing > 0:
            info["n_missing"] = int(n_missing)
            info["valid_mask"] = ~series.isna()
            logger.info(
                f"Covariate {name}: {n_missing} samples excluded due to missing values"
            )
        else:
            info["valid_mask"] = pd.Series(True, index=series.index)

        return info

    def run_covariate_prescreening(self) -> Dict:
        """Classify all metadata columns.

        Returns:
            Dict mapping covariate name to info dict.
        """
        for col in self.metadata.columns:
            self.covariate_info[col] = self._classify_covariate(
                self.metadata[col], col
            )
        return self.covariate_info

    def run_sample_qc(self) -> pd.DataFrame:
        """Run sample-level QC checks.

        Returns:
            QC summary DataFrame.
        """
        n_samples = self.counts.shape[1]
        samples = self.counts.columns.tolist()

        library_sizes = self.counts.sum(axis=0)
        lib_median = float(np.median(library_sizes))
        lib_mad = float(np.median(np.abs(library_sizes - lib_median)))

        lower = lib_median - 3 * lib_mad
        upper = lib_median + 3 * lib_mad
        lib_outlier_flag = (library_sizes < lower) | (library_sizes > upper)

        n_outliers = lib_outlier_flag.sum()
        if n_outliers > n_samples * 0.25:
            logger.warning(
                f"⚠️ A high proportion of samples ({n_outliers/n_samples:.0%}) are "
                "flagged as library size outliers. This may indicate systematic depth "
                "differences between sample groups (e.g., different platforms, sequencing "
                "runs, or batches) rather than individual problematic samples."
            )

        # CPM for dominant gene check
        cpm = self.counts.div(library_sizes, axis=1) * 1e6
        max_gene_cpm = cpm.max(axis=0)
        max_gene_name = cpm.idxmax(axis=0)
        dominant_gene_pct = max_gene_cpm / library_sizes * 100

        dominant_gene_flag = dominant_gene_pct > 30

        # Zero count saturation
        zero_pct = (self.counts == 0).sum(axis=0) / self.counts.shape[0] * 100
        zero_flag = zero_pct > 80

        qc_df = pd.DataFrame({
            "sample_id": samples,
            "library_size": library_sizes.values,
            "lib_outlier_flag": lib_outlier_flag.values,
            "dominant_gene": max_gene_name.values,
            "dominant_gene_pct": dominant_gene_pct.values,
            "zero_count_pct": zero_pct.values,
            "zero_flag": zero_flag.values,
            "mahal_distance": np.nan,
            "mahal_outlier_flag": False,
            "iqr_outlier_flag": False,
        })

        self.qc_summary = qc_df
        return qc_df

    def run_power_assessment(self) -> List[str]:
        """Assess statistical power for each categorical covariate.

        Returns:
            List of power warning strings.
        """
        n_samples = self.counts.shape[1]
        warnings = []

        any_low_power = False
        max_min_detectable = 0.0

        for name, info in self.covariate_info.items():
            if info.get("skip") or info.get("cov_type") != "categorical":
                continue
            n_groups = info.get("n_groups", 2)
            min_icc = find_min_detectable_icc(n_samples, n_groups)
            info["min_detectable_icc"] = min_icc

            msg = (
                f"Statistical power: With n={n_samples} samples and "
                f"{n_groups} groups, this analysis can detect ICC ≥ "
                f"{min_icc:.2f} with 80% power at α=0.05."
            )
            warnings.append(msg)

            if min_icc > 0.40:
                any_low_power = True
                max_min_detectable = max(max_min_detectable, min_icc)

        if any_low_power:
            self.low_power = True
            self.min_detectable_icc = max_min_detectable
            logger.warning(
                f"⚠️ LOW STATISTICAL POWER: With {n_samples} samples, this "
                f"analysis cannot reliably detect moderate batch effects "
                f"(ICC < {max_min_detectable:.2f}). A negative result does NOT "
                "rule out batch effects."
            )
        else:
            self.low_power = False
            self.min_detectable_icc = max_min_detectable

        self.power_warnings = warnings
        return warnings

    def filter_genes(
        self,
        primary_covariate: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, Dict]:
        """Apply gene-level filtering.

        Args:
            primary_covariate: Primary categorical covariate for min_samples logic.

        Returns:
            Tuple of (filtered_counts, filter_stats).
        """
        stats = {
            "n_input": self.counts.shape[0],
            "n_after_zero_removal": 0,
            "n_after_cpm_filter": 0,
        }

        # Step 1: Remove all-zero genes
        nonzero_mask = self.counts.sum(axis=1) > 0
        counts_filtered = self.counts[nonzero_mask]
        stats["n_after_zero_removal"] = counts_filtered.shape[0]

        # Step 2: CPM-based filter
        cpm_matrix = counts_filtered.div(counts_filtered.sum(axis=0), axis=1) * 1e6

        # Determine min_samples_expressing
        if self.min_samples_expressing is not None:
            min_samples = self.min_samples_expressing
        else:
            # Auto-compute
            cat_covs = [
                name for name, info in self.covariate_info.items()
                if not info.get("skip") and info.get("cov_type") == "categorical"
            ]
            if primary_covariate and primary_covariate in self.covariate_info:
                info = self.covariate_info[primary_covariate]
                if not info.get("skip") and info.get("cov_type") == "categorical":
                    min_samples = min(info["group_sizes"].values())
                else:
                    min_samples = self._auto_min_samples(cat_covs)
            elif cat_covs:
                min_samples = self._auto_min_samples(cat_covs)
            else:
                min_samples = max(3, counts_filtered.shape[1] // 4)

        self.min_samples_expressing_resolved = min_samples

        keep_mask = (cpm_matrix > self.min_cpm).sum(axis=1) >= min_samples
        counts_filtered = counts_filtered[keep_mask]
        stats["n_after_cpm_filter"] = counts_filtered.shape[0]

        n_zero_removed = stats["n_input"] - stats["n_after_zero_removal"]
        n_low_removed = stats["n_after_zero_removal"] - stats["n_after_cpm_filter"]
        logger.info(
            f"Gene filtering: {stats['n_input']} input → "
            f"{stats['n_after_cpm_filter']} retained "
            f"({n_zero_removed} all-zero removed, {n_low_removed} low-expression removed)"
        )

        self.counts = counts_filtered
        return counts_filtered, stats

    def _auto_min_samples(self, cat_covs: List[str]) -> int:
        """Compute auto min_samples from group sizes."""
        all_sizes = []
        for name in cat_covs:
            info = self.covariate_info[name]
            if "group_sizes" in info:
                all_sizes.extend(info["group_sizes"].values())
        if all_sizes:
            return min(all_sizes)
        return max(3, self.counts.shape[1] // 4)
