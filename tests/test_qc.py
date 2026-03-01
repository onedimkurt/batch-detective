"""Tests for QC module."""

import numpy as np
import pandas as pd
import pytest

from batch_detective.qc import QualityController, estimate_kruskal_power


def test_kruskal_power_basic():
    """Power should increase with sample size."""
    p1 = estimate_kruskal_power(20, 2, 0.2)
    p2 = estimate_kruskal_power(100, 2, 0.2)
    assert p2 > p1


def test_covariate_prescreening(synthetic_data):
    counts, metadata = synthetic_data
    qc = QualityController(counts, metadata)
    info = qc.run_covariate_prescreening()
    assert "batch" in info
    assert info["batch"]["cov_type"] == "categorical"
    assert info["rin_score"]["cov_type"] == "continuous"


def test_constant_covariate_skipped():
    counts = pd.DataFrame({"S1": [1, 2], "S2": [3, 4]}, index=["g1", "g2"])
    metadata = pd.DataFrame({"species": ["human", "human"]}, index=["S1", "S2"])
    qc = QualityController(counts, metadata)
    info = qc.run_covariate_prescreening()
    assert info["species"]["skip"]


def test_sample_qc_flags_lib_outliers(synthetic_data):
    counts, metadata = synthetic_data
    qc = QualityController(counts, metadata)
    qc.run_covariate_prescreening()
    qc_summary = qc.run_sample_qc()
    # Should have some lib_outlier_flag = True
    assert qc_summary["lib_outlier_flag"].any()


def test_gene_filter(synthetic_data):
    counts, metadata = synthetic_data
    # Use very high min_cpm AND high min_samples_expressing to force filtering
    qc = QualityController(counts, metadata, min_cpm=500.0, min_samples_expressing=20)
    qc.run_covariate_prescreening()
    qc.run_sample_qc()
    filtered, stats = qc.filter_genes()
    assert filtered.shape[0] < counts.shape[0]


def test_repeated_measures_detection():
    """Subject ID with avg group size < 3 should be flagged."""
    n = 12
    counts = pd.DataFrame(
        np.random.randint(1, 100, (20, n)),
        index=[f"g{i}" for i in range(20)],
        columns=[f"S{i}" for i in range(n)],
    )
    # Each patient has 2 samples
    metadata = pd.DataFrame({
        "patient_id": [f"P{i//2}" for i in range(n)],
        "treatment": ["ctrl"] * 6 + ["treat"] * 6,
    }, index=[f"S{i}" for i in range(n)])
    qc = QualityController(counts, metadata)
    info = qc.run_covariate_prescreening()
    assert info["patient_id"]["skip"]
    assert "patient_id" in qc.repeated_measures_covariates
