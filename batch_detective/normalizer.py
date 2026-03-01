"""Stage 3: Normalization and variable gene selection for batch-detective."""

import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def normalize_counts(
    counts: pd.DataFrame,
    n_variable_genes: int = 2000,
    normalized: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Normalize counts and select variable genes.

    Args:
        counts: Filtered count matrix (genes × samples).
        n_variable_genes: Number of top variable genes to select.
        normalized: If True, data is pre-normalized; skip CPM step.

    Returns:
        Tuple of (log1p_cpm_all, log1p_cpm_selected, gene_stats_df).
    """
    n_genes_input = counts.shape[0]

    if normalized:
        # Use as-is, just log1p-transform if values look like CPM
        log1p_cpm = np.log1p(counts.values.astype(float))
        log1p_cpm_df = pd.DataFrame(
            log1p_cpm, index=counts.index, columns=counts.columns
        )
        logger.info("Pre-normalized data: skipping CPM normalization.")
    else:
        # CPM normalization
        library_sizes = counts.sum(axis=0)
        cpm = counts.div(library_sizes, axis=1) * 1e6
        log1p_cpm = np.log1p(cpm.values.astype(float))
        log1p_cpm_df = pd.DataFrame(
            log1p_cpm, index=counts.index, columns=counts.columns
        )

    # Step 3: Remove zero-variance genes
    gene_vars = log1p_cpm_df.var(axis=1)
    nonzero_var = gene_vars > 0
    log1p_cpm_df = log1p_cpm_df[nonzero_var]
    n_after_var = log1p_cpm_df.shape[0]

    # Step 4: Remove bottom 1st percentile variance genes
    var_threshold = np.percentile(log1p_cpm_df.var(axis=1), 1)
    above_threshold = log1p_cpm_df.var(axis=1) > var_threshold
    log1p_cpm_df = log1p_cpm_df[above_threshold]
    n_after_pct = log1p_cpm_df.shape[0]

    # MAD-based variable gene selection
    gene_medians = log1p_cpm_df.median(axis=1)
    gene_mad = (log1p_cpm_df.sub(gene_medians, axis=0)).abs().median(axis=1)

    n_select = min(n_variable_genes, log1p_cpm_df.shape[0])
    top_genes = gene_mad.nlargest(n_select).index
    log1p_cpm_selected = log1p_cpm_df.loc[top_genes]

    logger.info(
        f"Normalization: {n_genes_input} input → "
        f"{n_after_var} after zero-var removal → "
        f"{n_after_pct} after 1st-pct filter → "
        f"{n_select} after MAD selection"
    )

    # Gene stats
    gene_stats = pd.DataFrame({
        "gene_id": log1p_cpm_selected.index,
        "mean_log_cpm": log1p_cpm_selected.mean(axis=1).values,
        "mad": gene_mad.loc[top_genes].values,
    })

    return log1p_cpm_df, log1p_cpm_selected, gene_stats
