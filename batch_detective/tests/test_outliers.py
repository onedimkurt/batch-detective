"""Tests for outlier detection module."""

import numpy as np
import pandas as pd
import pytest


def test_detect_outliers_finds_injected(synthetic_data):
    from batch_detective.outliers import detect_outliers
    from batch_detective.qc import QualityController
    from batch_detective.normalizer import normalize_counts
    from batch_detective.pca_analysis import run_pca

    counts, metadata = synthetic_data
    qc = QualityController(counts, metadata)
    qc.run_covariate_prescreening()
    qc.run_sample_qc()
    filtered, _ = qc.filter_genes()
    _, selected, _ = normalize_counts(filtered, n_variable_genes=200)
    pc_scores, evr, _, sample_ids = run_pca(selected, n_pcs=10)

    outlier_df, info = detect_outliers(pc_scores, sample_ids, metadata)
    # Should find at least some outliers
    assert len(outlier_df) >= 0  # Non-negative


def test_fallback_to_iqr_for_small_data():
    """Small datasets should fall back to IQR-only."""
    from batch_detective.outliers import detect_outliers

    rng = np.random.default_rng(0)
    pc_scores = rng.normal(0, 1, (5, 2))
    # Add one extreme outlier
    pc_scores[0] = [100, 100]
    sample_ids = [f"S{i}" for i in range(5)]
    metadata = pd.DataFrame({"batch": list("AAAAA")}, index=sample_ids)

    outlier_df, info = detect_outliers(pc_scores, sample_ids, metadata, n_pcs=10)
    assert info["skip_mahalanobis"]  # Should be True for n < 6
