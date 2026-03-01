"""Stage 8: Outlier detection for batch-detective."""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.covariance import LedoitWolf

logger = logging.getLogger(__name__)


def detect_outliers(
    pc_scores: np.ndarray,
    sample_ids: List[str],
    metadata: pd.DataFrame,
    n_pcs: int = 10,
    outlier_pval: float = 0.001,
) -> Tuple[pd.DataFrame, Dict]:
    """Detect outlier samples using Mahalanobis distance and IQR.

    Args:
        pc_scores: PC scores (n_samples, n_pcs).
        sample_ids: Sample ID list.
        metadata: Metadata DataFrame.
        n_pcs: Max PCs to use.
        outlier_pval: P-value threshold for Mahalanobis.

    Returns:
        Tuple of (outlier_df, detection_info).
    """
    n_samples = pc_scores.shape[1] if pc_scores.ndim > 1 else len(pc_scores)
    n_samples = pc_scores.shape[0]

    detection_info = {
        "method": [],
        "n_pcs_mahal": 0,
        "skip_mahalanobis": False,
    }

    mahal_distances = np.full(n_samples, np.nan)
    mahal_outliers = np.zeros(n_samples, dtype=bool)
    iqr_outliers = np.zeros(n_samples, dtype=bool)

    # Mahalanobis with LedoitWolf
    n_pcs_mahal = min(n_pcs, n_samples // 3)
    detection_info["n_pcs_mahal"] = n_pcs_mahal

    if n_pcs_mahal < 2:
        logger.warning(
            f"n_samples={n_samples} too small for Mahalanobis outlier detection "
            f"(requires n >= 6). Using IQR-based detection only."
        )
        detection_info["skip_mahalanobis"] = True
    else:
        pc_mahal = pc_scores[:, :n_pcs_mahal]
        try:
            lw = LedoitWolf().fit(pc_mahal)
            mahal_sq = lw.mahalanobis(pc_mahal)
            threshold = chi2.ppf(1 - outlier_pval, df=n_pcs_mahal)
            mahal_outliers = mahal_sq > threshold
            mahal_distances = np.sqrt(np.maximum(mahal_sq, 0))
        except Exception as e:
            logger.warning(f"Mahalanobis computation failed: {e}. Using IQR only.")
            detection_info["skip_mahalanobis"] = True

    # IQR-based on PC1 and PC2
    for pc_idx in [0, 1]:
        if pc_scores.shape[1] <= pc_idx:
            break
        scores = pc_scores[:, pc_idx]
        q75 = np.percentile(scores, 75)
        q25 = np.percentile(scores, 25)
        iqr_val = q75 - q25
        median = np.median(scores)
        threshold_iqr = 3 * iqr_val
        outlier_mask = np.abs(scores - median) > threshold_iqr
        iqr_outliers |= outlier_mask

    # Build outlier report
    is_outlier = mahal_outliers | iqr_outliers
    outlier_indices = np.where(is_outlier)[0]

    rows = []
    for idx in outlier_indices:
        method = []
        if mahal_outliers[idx]:
            method.append("Mahalanobis")
        if iqr_outliers[idx]:
            method.append("IQR-PC1/PC2")

        row = {
            "sample_id": sample_ids[idx],
            "detection_method": "/".join(method),
            "mahal_distance": mahal_distances[idx],
        }
        # Add metadata
        if sample_ids[idx] in metadata.index:
            for col in metadata.columns:
                row[col] = metadata.loc[sample_ids[idx], col]
        rows.append(row)

    outlier_df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["sample_id", "detection_method", "mahal_distance"]
    )

    return outlier_df, {
        "mahal_distances": mahal_distances,
        "mahal_outliers": mahal_outliers,
        "iqr_outliers": iqr_outliers,
        "n_pcs_mahal": n_pcs_mahal,
        "skip_mahalanobis": detection_info["skip_mahalanobis"],
    }
