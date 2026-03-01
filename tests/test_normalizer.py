"""Tests for normalizer module."""

import numpy as np
import pandas as pd
import pytest
from batch_detective.normalizer import normalize_counts


def test_normalize_counts_shape(synthetic_data):
    counts, _ = synthetic_data
    _, selected, _ = normalize_counts(counts, n_variable_genes=200)
    assert selected.shape[0] <= 200
    assert selected.shape[1] == counts.shape[1]


def test_normalize_selects_variable_genes():
    """MAD selection should prefer high-variance genes."""
    rng = np.random.default_rng(0)
    n = 20
    # Base low-variance counts
    counts = pd.DataFrame(
        rng.negative_binomial(50, 0.5, (50, n)),
        index=[f"g{i}" for i in range(50)],
        columns=[f"S{i}" for i in range(n)],
    )
    # Make first 5 genes high variance by adding large group-specific signal
    # Half samples get large addition, half do not - creates real MAD signal
    counts.iloc[:5, :10] += 500
    _, selected, _ = normalize_counts(counts, n_variable_genes=5)
    top_genes = set(selected.index)
    high_var = {f"g{i}" for i in range(5)}
    overlap = len(top_genes & high_var)
    assert overlap >= 3  # Most high-var genes should be selected


def test_prenormalized_skips_cpm():
    rng = np.random.default_rng(0)
    n = 6
    data = pd.DataFrame(
        rng.random((20, n)) * 10,
        index=[f"g{i}" for i in range(20)],
        columns=[f"S{i}" for i in range(n)],
    )
    _, selected, _ = normalize_counts(data, n_variable_genes=10, normalized=True)
    assert selected.shape[0] <= 10
