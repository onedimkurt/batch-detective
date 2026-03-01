"""Stage 4: PCA analysis for batch-detective."""

import logging
from typing import Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)


def run_pca(
    log1p_cpm_selected: pd.DataFrame,
    n_pcs: int = 10,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run PCA on log1p CPM selected genes.

    Mean-center ONLY (do NOT scale to unit variance).

    Args:
        log1p_cpm_selected: Log1p-CPM matrix (genes × samples).
        n_pcs: Number of PCs to compute.

    Returns:
        Tuple of (pc_scores, explained_variance_ratio, components, n_components_used).
    """
    # Transpose: samples × genes
    X = log1p_cpm_selected.T.values.astype(float)
    n_samples, n_genes = X.shape

    # Mean-center only
    X_centered = X - X.mean(axis=0)

    n_components = min(n_samples - 1, n_pcs, n_genes)
    if n_components < 2:
        logger.warning(
            f"Too few samples ({n_samples}) for meaningful PCA. "
            "Need at least 3 samples."
        )
        n_components = max(2, n_components)

    pca = PCA(n_components=n_components)
    pc_scores = pca.fit_transform(X_centered)
    explained_variance_ratio = pca.explained_variance_ratio_
    components = pca.components_

    sample_ids = log1p_cpm_selected.columns.tolist()

    logger.info(
        f"PCA: {n_components} PCs computed. "
        f"PC1 explains {explained_variance_ratio[0]:.1%} of variance."
    )

    return pc_scores, explained_variance_ratio, components, sample_ids
