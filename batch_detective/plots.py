"""Stage 9: Plot generation for batch-detective."""

import base64
import io
import logging
import sys
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from .compat import tqdm

logger = logging.getLogger(__name__)

COLORBLIND_PALETTE = "tab10"


def _fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64


def _save_figure(fig: plt.Figure, output_dir, name: str, dpi: int = 300) -> None:
    """Save figure as PNG and SVG."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(output_dir / f"{name}.svg", bbox_inches="tight")


def plot_library_sizes(
    qc_summary: pd.DataFrame,
    dpi: int = 150,
    figures_dir=None,
) -> str:
    """Plot library size bar chart.

    Args:
        qc_summary: QC summary with library_size and lib_outlier_flag.
        dpi: Figure resolution.
        figures_dir: Optional output directory for PNG/SVG.

    Returns:
        Base64-encoded PNG string.
    """
    fig, ax = plt.subplots(figsize=(max(8, len(qc_summary) * 0.3), 4))

    colors = [
        "#e74c3c" if flag else "#3498db"
        for flag in qc_summary["lib_outlier_flag"]
    ]
    ax.bar(range(len(qc_summary)), qc_summary["library_size"], color=colors)
    ax.set_xticks(range(len(qc_summary)))
    ax.set_xticklabels(
        qc_summary["sample_id"], rotation=90, fontsize=7
    )
    ax.set_ylabel("Library size (total counts)")
    ax.set_title("Library Sizes by Sample")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#3498db", label="Normal"),
        Patch(facecolor="#e74c3c", label="Outlier (±3 MAD)"),
    ]
    ax.legend(handles=legend_elements)
    plt.tight_layout()

    if figures_dir:
        _save_figure(fig, figures_dir, "library_sizes")
    return _fig_to_base64(fig, dpi)


def plot_scree(
    explained_variance_ratio: np.ndarray,
    dpi: int = 150,
    figures_dir=None,
) -> str:
    """Plot PCA scree plot.

    Args:
        explained_variance_ratio: Variance explained per PC.
        dpi: Figure resolution.
        figures_dir: Optional output directory.

    Returns:
        Base64-encoded PNG string.
    """
    n_pcs = len(explained_variance_ratio)
    fig, ax = plt.subplots(figsize=(7, 4))
    pct = explained_variance_ratio * 100
    ax.bar(range(1, n_pcs + 1), pct, color="#3498db")
    ax.plot(range(1, n_pcs + 1), np.cumsum(pct), "ro-", markersize=4, label="Cumulative %")
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance Explained (%)")
    ax.set_title("PCA Scree Plot")
    ax.legend()
    plt.tight_layout()

    if figures_dir:
        _save_figure(fig, figures_dir, "scree_plot")
    return _fig_to_base64(fig, dpi)


def plot_icc_barplot(
    icc_table: pd.DataFrame,
    low_power: bool = False,
    dpi: int = 150,
    figures_dir=None,
) -> str:
    """Plot ICC bar chart with 95% CI error bars.

    Args:
        icc_table: ICC results.
        low_power: If True, add warning indicator.
        dpi: Figure resolution.
        figures_dir: Optional output directory.

    Returns:
        Base64-encoded PNG string.
    """
    if icc_table.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No ICC data", ha="center", va="center",
                transform=ax.transAxes)
        return _fig_to_base64(fig, dpi)

    tier_colors = {
        "negligible": "#2ecc71",
        "mild": "#f1c40f",
        "moderate": "#e67e22",
        "strong": "#e74c3c",
    }

    fig, ax = plt.subplots(figsize=(max(6, len(icc_table) * 1.2), 5))

    x = range(len(icc_table))
    colors = [tier_colors.get(t, "#95a5a6") for t in icc_table["icc_tier"]]

    bars = ax.bar(x, icc_table["median_icc"], color=colors, alpha=0.8)

    # Error bars (95% CI)
    yerr_lower = icc_table["median_icc"] - icc_table["ci_lower_95"]
    yerr_upper = icc_table["ci_upper_95"] - icc_table["median_icc"]
    ax.errorbar(
        x, icc_table["median_icc"],
        yerr=[yerr_lower.clip(0), yerr_upper.clip(0)],
        fmt="none", color="black", capsize=4, linewidth=1.5,
    )

    ax.axhline(0.30, color="orange", linestyle="--", alpha=0.7, label="Moderate threshold (0.30)")
    ax.axhline(0.60, color="red", linestyle="--", alpha=0.5, label="Strong threshold (0.60)")

    ax.set_xticks(list(x))
    labels = icc_table["covariate"].tolist()
    if low_power:
        labels = [f"⚠️ {l}" for l in labels]
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Median ICC")
    ax.set_title("ICC per Covariate (Headline Metric)" + (" — ⚠️ LOW POWER" if low_power else ""))
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8)
    plt.tight_layout()

    if figures_dir:
        _save_figure(fig, figures_dir, "icc_barplot")
    return _fig_to_base64(fig, dpi)


def plot_association_heatmap(
    assoc_df: pd.DataFrame,
    dpi: int = 150,
    figures_dir=None,
) -> str:
    """Plot PC-metadata association heatmap.

    Args:
        assoc_df: Association results DataFrame.
        dpi: Figure resolution.
        figures_dir: Optional output directory.

    Returns:
        Base64-encoded PNG string.
    """
    if assoc_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No association data", ha="center", va="center",
                transform=ax.transAxes)
        return _fig_to_base64(fig, dpi)

    # Pivot: covariates × PCs, values = absolute effect size
    pivot_effect = assoc_df.pivot_table(
        index="covariate", columns="pc",
        values="effect_size", aggfunc=lambda x: abs(x.mean())
    )
    pivot_q = assoc_df.pivot_table(
        index="covariate", columns="pc",
        values="pval_adjusted_bh", aggfunc="min"
    )

    fig, ax = plt.subplots(
        figsize=(max(5, pivot_effect.shape[1] * 1.2), max(4, pivot_effect.shape[0] * 0.8))
    )

    sns.heatmap(
        pivot_effect.fillna(0),
        ax=ax,
        cmap="YlOrRd",
        vmin=0, vmax=1,
        annot=False,
        linewidths=0.5,
        cbar_kws={"label": "Effect size (|η² or |ρ|)"},
    )

    # Add asterisks
    for i, cov in enumerate(pivot_effect.index):
        for j, pc in enumerate(pivot_effect.columns):
            if cov in pivot_q.index and pc in pivot_q.columns:
                q = pivot_q.loc[cov, pc]
                if pd.notna(q):
                    if q < 0.001:
                        marker = "***"
                    elif q < 0.01:
                        marker = "**"
                    elif q < 0.05:
                        marker = "*"
                    else:
                        marker = ""
                    if marker:
                        ax.text(
                            j + 0.5, i + 0.5, marker,
                            ha="center", va="center",
                            fontsize=10, color="white", fontweight="bold",
                        )

    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Covariate")
    ax.set_title("PC-Metadata Association Heatmap\n(* q<0.05, ** q<0.01, *** q<0.001)")
    plt.tight_layout()

    if figures_dir:
        _save_figure(fig, figures_dir, "association_heatmap")
    return _fig_to_base64(fig, dpi)


def plot_pca_scatter(
    pc_scores: np.ndarray,
    sample_ids: List[str],
    explained_variance_ratio: np.ndarray,
    metadata: pd.DataFrame,
    covariate: str,
    cov_type: str,
    outlier_samples: List[str],
    pc_x: int = 0,
    pc_y: int = 1,
    dpi: int = 150,
    figures_dir=None,
    fig_name: str = "",
) -> str:
    """Plot PCA scatter colored by covariate.

    Args:
        pc_scores: PC scores array.
        sample_ids: Sample ID list.
        explained_variance_ratio: Variance per PC.
        metadata: Metadata DataFrame.
        covariate: Column to color by.
        cov_type: 'categorical' or 'continuous'.
        outlier_samples: List of outlier sample IDs.
        pc_x: Index of X PC (0-based).
        pc_y: Index of Y PC (0-based).
        dpi: Resolution.
        figures_dir: Optional output directory.
        fig_name: File name for saving.

    Returns:
        Base64-encoded PNG string.
    """
    fig, ax = plt.subplots(figsize=(7, 6))

    pct_x = explained_variance_ratio[pc_x] * 100 if pc_x < len(explained_variance_ratio) else 0
    pct_y = explained_variance_ratio[pc_y] * 100 if pc_y < len(explained_variance_ratio) else 0

    x = pc_scores[:, pc_x]
    y = pc_scores[:, pc_y]

    cov_vals = metadata.reindex(sample_ids)[covariate]

    if cov_type == "categorical":
        groups = sorted(cov_vals.dropna().unique(), key=str)
        palette = plt.cm.get_cmap(COLORBLIND_PALETTE, len(groups))
        for i, g in enumerate(groups):
            mask = cov_vals == g
            ax.scatter(x[mask], y[mask], c=[palette(i)], label=str(g), s=60, alpha=0.8)
        ax.legend(title=covariate, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    else:
        cov_num = pd.to_numeric(cov_vals, errors="coerce")
        sc = ax.scatter(x, y, c=cov_num, cmap="viridis", s=60, alpha=0.8)
        plt.colorbar(sc, ax=ax, label=covariate)

    # Annotate outliers
    for i, sid in enumerate(sample_ids):
        if sid in outlier_samples:
            ax.annotate(sid, (x[i], y[i]), fontsize=6, ha="left", va="bottom")

    ax.set_xlabel(f"PC{pc_x+1} ({pct_x:.1f}% variance)")
    ax.set_ylabel(f"PC{pc_y+1} ({pct_y:.1f}% variance)")
    ax.set_title(f"PCA — colored by {covariate}")
    plt.tight_layout()

    if figures_dir and fig_name:
        _save_figure(fig, figures_dir, fig_name)
    return _fig_to_base64(fig, dpi)


def plot_sample_distance_heatmap(
    pc_scores: np.ndarray,
    sample_ids: List[str],
    metadata: pd.DataFrame,
    covariate_info: Dict,
    n_pcs: int = 10,
    dpi: int = 150,
    figures_dir=None,
) -> str:
    """Plot sample distance heatmap with metadata strips.

    Args:
        pc_scores: PC scores (n_samples × n_pcs).
        sample_ids: Sample IDs.
        metadata: Metadata DataFrame.
        covariate_info: Covariate info dict.
        n_pcs: Number of PCs to use for distances.
        dpi: Resolution.
        figures_dir: Optional output dir.

    Returns:
        Base64-encoded PNG string.
    """
    n_use = min(n_pcs, pc_scores.shape[1])
    pc_use = pc_scores[:, :n_use]

    dist_mat = squareform(pdist(pc_use, metric="euclidean"))
    dist_df = pd.DataFrame(dist_mat, index=sample_ids, columns=sample_ids)

    # Clustering
    linkage = hierarchy.linkage(pdist(pc_use), method="average")
    order = hierarchy.leaves_list(linkage)

    dist_ordered = dist_df.iloc[order, order]

    # Build metadata annotation colors
    valid_covs = [
        name for name, info in covariate_info.items()
        if not info.get("skip") and info.get("cov_type") == "categorical"
    ]

    n_annotation_rows = len(valid_covs)
    fig_height = max(6, len(sample_ids) * 0.2 + n_annotation_rows * 0.5)

    fig, axes = plt.subplots(
        1 + n_annotation_rows, 1,
        figsize=(max(8, len(sample_ids) * 0.25), fig_height),
        gridspec_kw={"height_ratios": [4] + [0.3] * n_annotation_rows},
    ) if n_annotation_rows > 0 else plt.subplots(figsize=(max(8, len(sample_ids) * 0.25), 6))

    if n_annotation_rows == 0:
        ax_main = axes
        ax_annots = []
    elif n_annotation_rows == 1:
        ax_main = axes[0]
        ax_annots = [axes[1]]
    else:
        ax_main = axes[0]
        ax_annots = axes[1:]

    sns.heatmap(
        dist_ordered,
        ax=ax_main,
        cmap="Blues",
        xticklabels=False,
        yticklabels=True,
        cbar_kws={"label": "Euclidean distance (PC space)"},
        linewidths=0,
    )
    ax_main.set_yticklabels(ax_main.get_yticklabels(), fontsize=6)
    ax_main.set_title("Sample Distance Heatmap (PC space)")

    # Metadata annotation strips
    for idx, cov_name in enumerate(valid_covs):
        ax = ax_annots[idx]
        cov_series = metadata.reindex(dist_ordered.columns)[cov_name]
        groups = sorted(cov_series.dropna().unique(), key=str)
        palette = plt.cm.get_cmap(COLORBLIND_PALETTE, len(groups))
        color_map = {g: palette(i) for i, g in enumerate(groups)}
        colors = [color_map.get(cov_series.get(s, None), (0.7, 0.7, 0.7)) for s in dist_ordered.columns]
        ax.imshow([colors], aspect="auto")
        ax.set_yticks([0])
        ax.set_yticklabels([cov_name], fontsize=7)
        ax.set_xticks([])

    plt.tight_layout()

    if figures_dir:
        _save_figure(fig, figures_dir, "sample_distance_heatmap")
    return _fig_to_base64(fig, dpi)


def plot_gene_icc_scatter(
    top_genes_df: pd.DataFrame,
    covariate: str,
    dpi: int = 150,
) -> str:
    """Plot gene ICC vs mean expression scatter for a covariate.

    Args:
        top_genes_df: Top genes DataFrame.
        covariate: Covariate name.
        dpi: Resolution.

    Returns:
        Base64-encoded PNG string.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(
        top_genes_df["mean_expression_log_cpm"],
        top_genes_df["gene_icc"],
        c=["#e74c3c" if t else "#3498db" for t in top_genes_df.get("is_housekeeping", [False] * len(top_genes_df))],
        s=40, alpha=0.7,
    )
    ax.set_xlabel("Mean expression (log CPM)")
    ax.set_ylabel("ICC")
    ax.set_title(f"Gene ICC vs Expression — {covariate}")
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color="#e74c3c", label="Housekeeping"),
        Patch(color="#3498db", label="Other"),
    ])
    plt.tight_layout()
    return _fig_to_base64(fig, dpi)


def generate_all_plots(
    qc_summary: pd.DataFrame,
    pc_scores: np.ndarray,
    explained_variance_ratio: np.ndarray,
    metadata: pd.DataFrame,
    covariate_info: Dict,
    icc_table: pd.DataFrame,
    assoc_df: pd.DataFrame,
    outlier_samples: List[str],
    sample_ids: List[str],
    top_genes: Dict[str, pd.DataFrame],
    dpi: int = 150,
    figures_dir=None,
    low_power: bool = False,
) -> Dict[str, str]:
    """Generate all plots and return as base64 strings.

    Args:
        qc_summary: QC summary DataFrame.
        pc_scores: PC scores array.
        explained_variance_ratio: Variance per PC.
        metadata: Metadata DataFrame.
        covariate_info: Covariate info dict.
        icc_table: ICC table.
        assoc_df: Association results.
        outlier_samples: Outlier sample IDs.
        sample_ids: All sample IDs.
        top_genes: Top genes per covariate.
        dpi: Resolution.
        figures_dir: Optional figures output directory.
        low_power: Whether low power warning is active.

    Returns:
        Dict of plot name → base64 string.
    """
    plots = {}
    disable_tqdm = not sys.stdout.isatty()

    plot_tasks = [
        ("library_sizes", lambda: plot_library_sizes(qc_summary, dpi, figures_dir)),
        ("scree_plot", lambda: plot_scree(explained_variance_ratio, dpi, figures_dir)),
        ("icc_barplot", lambda: plot_icc_barplot(icc_table, low_power, dpi, figures_dir)),
        ("association_heatmap", lambda: plot_association_heatmap(assoc_df, dpi, figures_dir)),
        ("sample_distance_heatmap", lambda: plot_sample_distance_heatmap(
            pc_scores, sample_ids, metadata, covariate_info, dpi=dpi, figures_dir=figures_dir
        )),
    ]

    for name, func in tqdm(plot_tasks, desc="Generating plots", disable=disable_tqdm):
        try:
            plots[name] = func()
        except Exception as e:
            logger.warning(f"Failed to generate {name}: {e}")
            plots[name] = ""

    # PCA scatter plots per covariate
    valid_covs = [
        (name, info) for name, info in covariate_info.items()
        if not info.get("skip")
    ]

    for cov_name, info in tqdm(valid_covs, desc="PCA plots", disable=disable_tqdm):
        cov_type = info.get("cov_type", "categorical")

        for pc_pair, suffix in [((0, 1), "pc1v2"), ((0, 2), "pc1v3")]:
            pc_x, pc_y = pc_pair
            if pc_scores.shape[1] <= max(pc_x, pc_y):
                continue
            fig_name = f"pca_{cov_name}_{suffix}" if figures_dir else ""
            try:
                img = plot_pca_scatter(
                    pc_scores, sample_ids, explained_variance_ratio,
                    metadata, cov_name, cov_type, outlier_samples,
                    pc_x=pc_x, pc_y=pc_y, dpi=dpi,
                    figures_dir=figures_dir, fig_name=fig_name,
                )
                plots[f"pca_{cov_name}_{suffix}"] = img
            except Exception as e:
                logger.warning(f"Failed PCA plot for {cov_name} {suffix}: {e}")

    # Gene ICC scatter plots
    for cov_name, df in top_genes.items():
        try:
            plots[f"gene_icc_scatter_{cov_name}"] = plot_gene_icc_scatter(df, cov_name, dpi)
        except Exception as e:
            logger.warning(f"Failed gene ICC scatter for {cov_name}: {e}")

    return plots
