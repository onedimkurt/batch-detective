"""Tests for association module."""

import numpy as np
import pandas as pd
import pytest
from batch_detective.association import (
    cramers_v,
    kruskal_eta_squared,
    compute_icc11_vectorized,
    bh_correct,
    detect_collinearity,
    run_pc_associations,
    compute_all_icc,
)


def test_cramers_v_identical():
    x = pd.Series(["A", "A", "B", "B"])
    y = pd.Series(["A", "A", "B", "B"])
    v = cramers_v(x, y)
    assert abs(v - 1.0) < 0.01


def test_cramers_v_independent():
    x = pd.Series(["A", "B", "A", "B"])
    y = pd.Series(["C", "C", "D", "D"])
    v = cramers_v(x, y)
    assert v < 0.3


def test_kruskal_eta_squared_range():
    eta = kruskal_eta_squared(H=10, k=3, n=30)
    assert 0.0 <= eta <= 1.0


def test_icc_vectorized_high():
    """ICC should be high when batch signal is strong."""
    rng = np.random.default_rng(0)
    n_genes, n_samples = 100, 24
    groups = np.array(["B1"] * 6 + ["B2"] * 6 + ["B3"] * 6 + ["B4"] * 6)
    X = rng.normal(0, 0.1, (n_genes, n_samples))
    # Add large group effect
    for j, g in enumerate(["B1", "B2", "B3", "B4"]):
        mask = groups == g
        X[:, mask] += j * 2.0
    icc_vals, median, _ = compute_icc11_vectorized(X, groups)
    assert median > 0.3


def test_icc_vectorized_low():
    """ICC should be low for pure noise."""
    rng = np.random.default_rng(1)
    n_genes, n_samples = 100, 24
    groups = np.array(["B1"] * 6 + ["B2"] * 6 + ["B3"] * 6 + ["B4"] * 6)
    X = rng.normal(0, 1.0, (n_genes, n_samples))
    _, median, _ = compute_icc11_vectorized(X, groups)
    assert median < 0.3


def test_bh_correct():
    pvals = [0.001, 0.01, 0.05, 0.2, 0.5]
    qvals = bh_correct(pvals)
    assert len(qvals) == len(pvals)
    assert all(q >= p for q, p in zip(qvals, sorted(pvals)))


def test_detect_collinearity_warns(synthetic_data):
    counts, metadata = synthetic_data
    from batch_detective.qc import QualityController
    qc = QualityController(counts, metadata)
    covariate_info = qc.run_covariate_prescreening()
    warnings = detect_collinearity(metadata, covariate_info)
    # any collinearity warning should be detected (batch/treatment are perfectly collinear)
    assert len(warnings) > 0


def test_full_icc_pipeline(synthetic_data):
    """End-to-end ICC on synthetic data should show batch ICC ~ 0.30-0.60."""
    counts, metadata = synthetic_data
    from batch_detective.qc import QualityController
    from batch_detective.normalizer import normalize_counts

    qc = QualityController(counts, metadata)
    covariate_info = qc.run_covariate_prescreening()
    qc.run_sample_qc()
    filtered, _ = qc.filter_genes()
    _, selected, _ = normalize_counts(filtered, n_variable_genes=300)

    icc_table = compute_all_icc(selected, metadata, covariate_info, n_bootstrap=50)
    assert not icc_table.empty
    batch_row = icc_table[icc_table["covariate"] == "batch"]
    assert not batch_row.empty
    batch_icc = float(batch_row["median_icc"].iloc[0])
    assert 0.20 < batch_icc < 0.80
