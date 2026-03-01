"""Stage 10+11: Report assembly for batch-detective."""

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from . import __version__

logger = logging.getLogger(__name__)


def _determine_traffic_light(
    icc_table: pd.DataFrame,
    assoc_df: pd.DataFrame,
    low_power: bool = False,
    technical_covariates: Optional[List[str]] = None,
) -> Dict[str, str]:
    """Determine traffic light color and text.

    Args:
        icc_table: ICC results.
        assoc_df: Association results.
        low_power: Whether low power warning is active.
        technical_covariates: Names of technical covariates.

    Returns:
        Dict with emoji, label, text, class.
    """
    if technical_covariates is None:
        technical_covariates = []

    if low_power and icc_table.empty:
        return {
            "emoji": "⚪",
            "label": "INSUFFICIENT POWER TO ASSESS",
            "text": (
                "Sample size is too small to reliably assess batch effects. "
                "Results are inconclusive. Consider this analysis exploratory only."
            ),
            "class": "banner-gray",
        }

    if not technical_covariates:
        return {
            "emoji": "🔵",
            "label": "ANALYSIS COMPLETE — INTERPRETATION REQUIRES COVARIATE LABELING",
            "text": (
                "Run with --technical-covariates and --biological-covariates "
                "for contextual interpretation. See ICC and association results below."
            ),
            "class": "banner-blue",
        }

    # Filter to technical covariates
    tech_icc = icc_table[icc_table["covariate"].isin(technical_covariates)]
    tech_assoc = assoc_df[assoc_df["covariate"].isin(technical_covariates)] if not assoc_df.empty else pd.DataFrame()

    max_tech_icc = tech_icc["median_icc"].max() if not tech_icc.empty else 0.0
    any_sig = (
        tech_assoc["significant_q05"].any()
        if not tech_assoc.empty and "significant_q05" in tech_assoc.columns
        else False
    )

    # Determine top offender
    top_cov = ""
    if not tech_icc.empty:
        top_row = tech_icc.loc[tech_icc["median_icc"].idxmax()]
        top_cov = top_row["covariate"]
        top_icc = top_row["median_icc"]

    if max_tech_icc >= 0.30 or (any_sig and max_tech_icc >= 0.15):
        return {
            "emoji": "🔴",
            "label": "SIGNIFICANT BATCH EFFECT DETECTED",
            "text": (
                f"[{top_cov}] accounts for {top_icc:.0%} of transcriptomic variance "
                "and is significantly associated with the major axis of variation. "
                "Batch correction or covariate modeling is recommended before "
                "differential expression analysis."
            ),
            "class": "banner-red",
        }
    elif max_tech_icc >= 0.15 or any_sig:
        return {
            "emoji": "🟡",
            "label": "POSSIBLE BATCH EFFECT — REVIEW RECOMMENDED",
            "text": (
                f"A potential batch effect was detected for [{top_cov}] "
                f"(ICC={max_tech_icc:.2f}). Review the PCA plots and discuss "
                "with your bioinformatician before proceeding."
            ),
            "class": "banner-yellow",
        }
    else:
        return {
            "emoji": "🟢",
            "label": "NO SIGNIFICANT BATCH EFFECT DETECTED",
            "text": (
                "No significant technical batch effects were detected. "
                "Your data appears suitable for standard differential expression analysis."
            ),
            "class": "banner-green",
        }


def _build_executive_summary(
    icc_table: pd.DataFrame,
    assoc_df: pd.DataFrame,
    covariate_info: Dict,
    traffic_light: Dict,
    condition_on: Optional[List[str]] = None,
    technical_covariates: Optional[List[str]] = None,
    biological_covariates: Optional[List[str]] = None,
    n_samples: int = 0,
) -> str:
    """Build executive summary paragraphs.

    Args:
        icc_table: ICC results.
        assoc_df: Association results.
        covariate_info: Covariate info.
        traffic_light: Traffic light dict.
        condition_on: Conditioning covariates.
        technical_covariates: Technical covariate names.
        biological_covariates: Biological covariate names.
        n_samples: Number of samples.

    Returns:
        HTML string.
    """
    paras = []

    if icc_table.empty and assoc_df.empty:
        return "<p>No covariates could be analyzed. Check data overview for details.</p>"

    # Main finding paragraph
    if not icc_table.empty:
        sorted_icc = icc_table.sort_values("median_icc", ascending=False)
        for _, row in sorted_icc.iterrows():
            cov = row["covariate"]
            icc = row["median_icc"]
            ci_lo = row["ci_lower_95"]
            ci_hi = row["ci_upper_95"]
            tier = row["icc_tier"]
            label = row.get("label", "unknown")

            # Find significant PC associations for this covariate
            cov_assoc = (
                assoc_df[
                    (assoc_df["covariate"] == cov) & (assoc_df.get("significant_q05", False))
                ]
                if not assoc_df.empty and "significant_q05" in assoc_df.columns
                else pd.DataFrame()
            )

            if label == "biological":
                interp = (
                    f"<strong>{cov}</strong> shows {tier} ICC (median ICC = {icc:.2f}, "
                    f"95% CI: {ci_lo:.2f}–{ci_hi:.2f}) — this is expected for a biological "
                    f"variable and suggests treatment groups have distinct transcriptomic profiles."
                )
            elif label == "technical":
                interp = (
                    f"<strong>{cov}</strong> shows {tier} ICC (median ICC = {icc:.2f}, "
                    f"95% CI: {ci_lo:.2f}–{ci_hi:.2f}) — this may represent a technical batch effect."
                )
            else:
                interp = (
                    f"<strong>{cov}</strong> shows {tier} ICC (median ICC = {icc:.2f}, "
                    f"95% CI: {ci_lo:.2f}–{ci_hi:.2f}) — label as biological or technical "
                    "for contextual interpretation."
                )

            paras.append(f"<p>{interp}</p>")

    # Recommendation based on traffic light
    tl_emoji = traffic_light.get("emoji", "")
    if tl_emoji == "🔴":
        if not assoc_df.empty and not icc_table.empty:
            # Check collinearity between batch and treatment
            top_tech = ""
            if technical_covariates and not icc_table.empty:
                tech_icc = icc_table[icc_table["covariate"].isin(technical_covariates)]
                if not tech_icc.empty:
                    top_tech = tech_icc.loc[tech_icc["median_icc"].idxmax(), "covariate"]
            paras.append(
                f"<p>If batch and treatment groups are fully crossed in your design, "
                f"consider including <code>{top_tech or 'batch'}</code> as a covariate in "
                "your DE model (e.g., <code>~batch + treatment</code> in DESeq2) or "
                "applying ComBat correction. If batch and treatment are collinear "
                "(see collinearity warning above), correction is not recommended — "
                "it would likely remove biological signal.</p>"
            )
    elif tl_emoji == "🟢":
        paras.append(
            "<p>No metadata variable shows significant association with the top 5 PCs after "
            "FDR correction (all q &gt; 0.05), and ICC values are negligible or mild for all "
            "covariates. Your data appears free of major technical confounders.</p>"
        )

    paras.append(
        "<p><em>Note: These are marginal (univariate) associations. Covariates correlated "
        "with biology may show elevated associations reflecting biological rather than technical "
        "variation. Use --primary-covariate to condition on known biological variables. "
        "High ICC values for known biological variables (e.g., treatment, disease status) "
        "are expected and do not indicate a technical problem — they confirm your experiment worked. "
        "Only unexpected associations (e.g., operator, plate, date) with high ICC suggest "
        "technical batch effects.</em></p>"
    )

    if condition_on:
        paras.append(
            f"<p><em>Results above are partial, conditioning on {', '.join(condition_on)}. "
            f"They show associations beyond what {', '.join(condition_on)} already explain.</em></p>"
        )

    if n_samples < 15:
        paras.append(
            f"<p><em>Note: Bootstrap confidence intervals are wide at small n (n={n_samples}). "
            "ICC estimates should be interpreted cautiously. Collect more samples if possible.</em></p>"
        )

    return "\n".join(paras)


def _build_next_steps(
    traffic_light: Dict,
    collinearity_warnings: List[Dict],
    technical_covariates: Optional[List[str]] = None,
    icc_table: Optional[pd.DataFrame] = None,
) -> str:
    """Build Next Steps section HTML.

    Args:
        traffic_light: Traffic light result.
        collinearity_warnings: Collinearity warning list.
        technical_covariates: Technical covariate names.
        icc_table: ICC table.

    Returns:
        HTML string.
    """
    emoji = traffic_light.get("emoji", "")

    # Check if batch/treatment are collinear
    batch_treatment_collinear = any(
        (w.get("cov_a") in (technical_covariates or [])
         or w.get("cov_b") in (technical_covariates or []))
        for w in collinearity_warnings
    )

    if emoji == "🔴":
        if batch_treatment_collinear:
            return """
<div class="banner banner-red">
  <strong>⚠️ IMPORTANT:</strong> Batch and treatment are collinear in your design.
  Batch correction IS NOT RECOMMENDED — it will remove biological signal.
  Proceed with standard analysis but acknowledge this limitation in your methods.
</div>"""
        else:
            return """
<p><strong>What to do next:</strong></p>
<pre>
If you use DESeq2:
  Add batch to your design formula: design = ~ batch + treatment
  Example: dds &lt;- DESeqDataSetFromMatrix(countData, colData, design = ~ batch + treatment)

If you use edgeR:
  Include batch in your model matrix: design &lt;- model.matrix(~ batch + treatment)

If you want to remove batch effects for visualization:
  Use ComBat from the sva package in R:
  library(sva); adjusted &lt;- ComBat(dat=log_counts, batch=batch_vector)
  ⚠️ Only use ComBat-corrected data for visualization — use the model-based
  approach above for differential expression to avoid double-correction.

Useful resources:
  DESeq2 vignette: https://bioconductor.org/packages/DESeq2
  sva/ComBat tutorial: https://bioconductor.org/packages/sva
</pre>"""
    elif emoji == "🟡":
        return """
<p>A potential batch effect was detected. Review the PCA plots in this report.
If the suspected batch variable clearly separates samples in PCA space, consider
modeling it in your DE analysis (see 🔴 guidance above). Discuss with your
bioinformatician before applying batch correction.</p>"""
    elif emoji == "🟢":
        return """
<p>Your data appears free of major batch effects.
Proceed with standard differential expression analysis using your preferred tool
(DESeq2, edgeR, limma-voom). No batch correction is indicated based on this analysis.</p>"""
    else:
        return """
<p>Label covariates using --technical-covariates and --biological-covariates
to receive contextual recommendations. See the ICC and association results above
for quantitative evidence.</p>"""


def _build_methods_text(
    version: str,
    timestamp: str,
    n_variable_genes: int,
    outlier_pval: float,
    n_pcs_mahal: int,
    skip_mahalanobis: bool,
    condition_on: Optional[List[str]] = None,
) -> str:
    """Build methods blurb text.

    Args:
        version: Tool version.
        timestamp: Run timestamp.
        n_variable_genes: Number of variable genes.
        outlier_pval: Outlier p-value threshold.
        n_pcs_mahal: PCs used for Mahalanobis.
        skip_mahalanobis: Whether Mahalanobis was skipped.
        condition_on: Conditioning covariates.

    Returns:
        HTML string.
    """
    try:
        import sklearn
        sklearn_version = sklearn.__version__
    except ImportError:
        sklearn_version = "unknown"

    mahal_text = (
        f"Sample outliers were detected using Mahalanobis distance with "
        f"LedoitWolf-regularized covariance (Ledoit &amp; Wolf, 2004), "
        f"thresholded at p&lt;{outlier_pval} (chi-squared, df={n_pcs_mahal})."
        if not skip_mahalanobis
        else "Sample outliers were detected using IQR-based method on PC1/PC2 "
             "(Mahalanobis skipped due to insufficient samples)."
    )

    return f"""
<p>Raw count data was processed using batch-detective v{version} (run ID: {timestamp}).
Counts were normalized to counts per million (CPM) and log1p-transformed. The top
{n_variable_genes} most variable genes (by median absolute deviation, MAD) were
selected for principal component analysis (PCA) using mean-centering without
variance scaling, implemented via scikit-learn {sklearn_version}
(Pedregosa et al., 2011). Note: log1p(CPM) normalization does not achieve full
variance stabilization (cf. DESeq2 VST); high-expression genes may contribute
disproportionately to variance. Association between PC scores and metadata covariates
was assessed using the Kruskal-Wallis test with rank-based eta-squared as effect size
for categorical covariates (eta-squared = max(0, (H-k+1)/(n-k))), and Spearman rank
correlation with bootstrap 95% CI (n=1000 resamplings) for continuous covariates
P-values were adjusted for multiple testing within each covariate using
the Benjamini-Hochberg FDR procedure (Benjamini &amp; Hochberg, 1995). Intraclass
correlation coefficients [ICC(1,1), one-way random effects model] were computed across
variable genes using a vectorized ANOVA decomposition; 95% confidence intervals were
estimated by bootstrap resampling of samples (n=1000). ICC values were interpreted per
Koo &amp; Li (2016). Batch labels were treated as randomly sampled from a population
(ICC(1,1) assumption); if batches are fixed, ICC values may underestimate true effect
magnitude. {mahal_text}</p>

<p><strong>LIMITATIONS:</strong></p>
<ol>
  <li>Assumes independent samples. Repeated measures and paired designs are not
      supported; apply longitudinal-specific methods for such data.</li>
  <li>log1p(CPM) normalization is not variance-stabilized; results may differ from
      analyses using DESeq2 VST or similar transformations.</li>
  <li>Spearman correlation detects monotonic associations only; non-monotonic
      relationships between continuous covariates and PC scores are not detected.</li>
  <li>Statistical power depends on sample size. Results marked ⚠️ have low power.</li>
  <li>Association tests are marginal (univariate); use --condition-on for partial analysis.</li>
</ol>

<p>This report was generated from read-only analysis. Original data files were not
modified. See run_manifest.json for full reproducibility record.</p>
"""


def assemble_report(
    output_path: Path,
    qc_summary: pd.DataFrame,
    icc_table: pd.DataFrame,
    assoc_df: pd.DataFrame,
    outlier_df: pd.DataFrame,
    plots: Dict[str, str],
    top_genes: Dict[str, Any],
    covariate_info: Dict,
    collinearity_warnings: List[Dict],
    repeated_measures_covariates: List[str],
    config: Dict[str, Any],
    analysis_meta: Dict[str, Any],
    detection_info: Dict,
    timestamp: str,
) -> None:
    """Assemble and write the HTML report.

    Args:
        output_path: Path for the HTML file.
        qc_summary: QC summary DataFrame.
        icc_table: ICC table.
        assoc_df: Association results.
        outlier_df: Outlier report.
        plots: Dict of plot name → base64.
        top_genes: Top batch-associated genes per covariate.
        covariate_info: Covariate info dict.
        collinearity_warnings: Collinearity warnings.
        repeated_measures_covariates: Repeated measures covariate names.
        config: Resolved config dict.
        analysis_meta: Analysis metadata.
        detection_info: Outlier detection info.
        timestamp: Run timestamp string.
    """
    technical_covariates = config.get("technical_covariates", [])
    biological_covariates = config.get("biological_covariates", [])
    condition_on = config.get("condition_on_list", [])
    low_power = analysis_meta.get("low_power", False)
    min_detectable_icc = analysis_meta.get("min_detectable_icc", 0.0)

    traffic_light = _determine_traffic_light(
        icc_table, assoc_df, low_power, technical_covariates
    )

    executive_summary = _build_executive_summary(
        icc_table, assoc_df, covariate_info, traffic_light,
        condition_on, technical_covariates, biological_covariates,
        n_samples=analysis_meta.get("n_samples", 0),
    )

    next_steps = _build_next_steps(
        traffic_light, collinearity_warnings, technical_covariates, icc_table
    )

    methods_text = _build_methods_text(
        version=__version__,
        timestamp=timestamp,
        n_variable_genes=config.get("n_variable_genes", 2000),
        outlier_pval=config.get("outlier_pval", 0.001),
        n_pcs_mahal=detection_info.get("n_pcs_mahal", 0),
        skip_mahalanobis=detection_info.get("skip_mahalanobis", False),
        condition_on=condition_on,
    )

    # QC flagged samples
    qc_flagged = (
        qc_summary[
            qc_summary["lib_outlier_flag"] |
            qc_summary["zero_flag"] |
            qc_summary["mahal_outlier_flag"]
        ]
        if not qc_summary.empty
        else pd.DataFrame()
    )

    # Outlier metadata columns
    outlier_meta_cols = [
        c for c in (outlier_df.columns if not outlier_df.empty else [])
        if c not in ("sample_id", "detection_method", "mahal_distance")
    ]

    # Top genes for template
    top_genes_data = {}
    for cov_name, df in top_genes.items():
        top_genes_data[cov_name] = df.to_dict("records")

    # Covariates tested/skipped
    covariates_tested = [
        n for n, i in covariate_info.items() if not i.get("skip")
    ]
    covariates_skipped = [
        f"{n} ({i.get('skip_reason', '')})"
        for n, i in covariate_info.items() if i.get("skip")
    ]

    template_dir = Path(__file__).parent / "templates"
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report_template.html")

    html = template.render(
        timestamp=timestamp,
        version=__version__,
        n_samples=analysis_meta.get("n_samples", 0),
        n_genes_input=analysis_meta.get("n_genes_input", 0),
        n_genes_analyzed=analysis_meta.get("n_genes_analyzed", 0),
        normalization_method=analysis_meta.get("normalization_method", "CPM + log1p"),
        n_variable_genes=config.get("n_variable_genes", 2000),
        n_pcs=config.get("n_pcs", 10),
        min_cpm=config.get("min_cpm", 1.0),
        covariates_tested=covariates_tested,
        covariates_skipped=covariates_skipped,
        condition_on=condition_on,
        anonymized=config.get("anonymize_samples", False),
        # Banners
        repeated_measures_covariates=repeated_measures_covariates,
        low_power=low_power,
        min_detectable_icc=min_detectable_icc,
        collinearity_warnings=collinearity_warnings,
        # Traffic light
        traffic_light_emoji=traffic_light["emoji"],
        traffic_light_label=traffic_light["label"],
        traffic_light_text=traffic_light["text"],
        traffic_light_class=traffic_light["class"],
        # Content
        executive_summary=executive_summary,
        next_steps=next_steps,
        qc_flagged_samples=qc_flagged.to_dict("records") if not qc_flagged.empty else [],
        icc_table_rows=(
            icc_table.drop(columns=["_icc_values"], errors="ignore").to_dict("records")
            if not icc_table.empty else []
        ),
        assoc_table_rows=(
            assoc_df.to_dict("records") if not assoc_df.empty else []
        ),
        n_covariates_tested=len(covariates_tested),
        outlier_rows=(
            outlier_df.to_dict("records") if not outlier_df.empty else []
        ),
        outlier_meta_cols=outlier_meta_cols,
        top_genes_data=top_genes_data,
        plots=plots,
        methods_text=methods_text,
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info(f"Report written to {output_path}")
