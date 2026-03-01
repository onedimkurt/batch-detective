"""Stage 11: Top batch-associated genes analysis."""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .association import compute_icc11_vectorized

logger = logging.getLogger(__name__)

HOUSEKEEPING_GENES = {
    "GAPDH", "ACTB", "B2M", "HPRT1", "TBP",
    "RPLP0", "SDHA", "HMBS", "UBC", "YWHAZ",
}


def get_top_batch_genes(
    log1p_cpm_selected: pd.DataFrame,
    metadata: pd.DataFrame,
    covariate_info: Dict,
    icc_table: pd.DataFrame,
    icc_threshold: float = 0.20,
    n_top: int = 20,
) -> Dict[str, pd.DataFrame]:
    """Get top batch-associated genes per technical covariate.

    Args:
        log1p_cpm_selected: Selected gene expression (genes × samples).
        metadata: Metadata DataFrame.
        covariate_info: Covariate info dict.
        icc_table: ICC results table.
        icc_threshold: Minimum median ICC to analyze genes.
        n_top: Number of top genes to return.

    Returns:
        Dict mapping covariate name to top genes DataFrame.
    """
    X = log1p_cpm_selected.values.astype(float)
    gene_ids = log1p_cpm_selected.index.tolist()
    mean_expr = log1p_cpm_selected.mean(axis=1).values
    expr_rank = pd.Series(mean_expr).rank(ascending=False).values

    results = {}

    for _, row in icc_table.iterrows():
        cov_name = row["covariate"]
        median_icc = row["median_icc"]

        if median_icc < icc_threshold:
            continue

        info = covariate_info.get(cov_name, {})
        if info.get("skip"):
            continue

        valid_mask = info.get(
            "valid_mask", pd.Series(True, index=metadata.index)
        )
        valid_idx = np.where(valid_mask.values)[0]

        if len(valid_idx) < 5:
            continue

        X_valid = X[:, valid_idx]
        labels = metadata.iloc[valid_idx][cov_name].values.astype(str)

        gene_iccs, _, _ = compute_icc11_vectorized(X_valid, labels)

        top_idx = np.argsort(gene_iccs)[-n_top:][::-1]

        top_df = pd.DataFrame({
            "gene_id": [gene_ids[i] for i in top_idx],
            "gene_icc": gene_iccs[top_idx],
            "mean_expression_log_cpm": mean_expr[top_idx],
            "expression_rank_by_mean": expr_rank[top_idx].astype(int),
            "is_housekeeping": [
                str(gene_ids[i]).upper() in HOUSEKEEPING_GENES for i in top_idx
            ],
        })

        results[cov_name] = top_df
        logger.info(
            f"Top {n_top} genes for {cov_name}: "
            f"highest ICC = {gene_iccs[top_idx[0]]:.3f}"
        )

    return results
