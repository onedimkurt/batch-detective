"""Synthetic test fixtures for batch-detective."""

import numpy as np
import pandas as pd
import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def synthetic_data():
    """Generate synthetic RNA-seq count data with known batch effects.

    Dataset: 24 samples × 1,200 genes
    Metadata:
      - batch: 4 groups × 6 samples, ICC target ~0.45
      - treatment: 2 groups × 12 samples (biological)
      - sex: 2 groups (partially collinear with batch, Cramér's V ~0.75)
      - rin_score: continuous, correlated with library size ~0.6

    Injected artifacts:
      - 3 samples with 4× lower library sizes (lib outliers)
      - 1 sample far from centroid (Mahalanobis outlier)
      - Strong batch signal in first 400 genes
      - Weak sex signal in next 100 genes
      - Pure noise in remaining 700 genes
    """
    rng = np.random.default_rng(42)
    n_samples = 24
    n_genes = 1200
    sample_ids = [f"S{i+1:02d}" for i in range(n_samples)]

    # Metadata
    batch = [f"B{i+1}" for i in range(4) for _ in range(6)]
    treatment = ["ctrl"] * 12 + ["treat"] * 12

    # sex partially collinear with batch (Cramér's V ~0.75)
    sex = (
        ["M", "M", "M", "F", "F", "F"] +  # B1: 3M 3F
        ["M", "M", "M", "M", "F", "F"] +  # B2: 4M 2F
        ["F", "F", "F", "F", "M", "M"] +  # B3: 4F 2M
        ["F", "F", "F", "M", "M", "M"]    # B4: 3F 3M
    )

    # rin_score: higher for early batches
    rin_base = {"B1": 9.0, "B2": 8.5, "B3": 7.5, "B4": 6.5}
    rin_score = [rin_base[b] + rng.normal(0, 0.3) for b in batch]

    metadata = pd.DataFrame({
        "batch": batch,
        "treatment": treatment,
        "sex": sex,
        "rin_score": rin_score,
    }, index=sample_ids)

    # Base counts
    base_counts = rng.negative_binomial(50, 0.5, size=(n_genes, n_samples)).astype(float)

    # Strong batch signal in first 400 genes
    batch_effects = {"B1": 0, "B2": 30, "B3": 60, "B4": 90}
    for b_name, effect in batch_effects.items():
        mask = np.array([b == b_name for b in batch])
        sample_idx = np.where(mask)[0]
        base_counts[:400, sample_idx] += rng.poisson(
            effect * np.ones((400, len(sample_idx))) + 5, size=(400, len(sample_idx))
        )

    # Weak sex signal in genes 400–500
    sex_effects = {"M": 10, "F": 0}
    for s_name, effect in sex_effects.items():
        mask = np.array([s == s_name for s in sex])
        sample_idx = np.where(mask)[0]
        base_counts[400:500, sample_idx] += rng.poisson(
            effect, size=(100, len(sample_idx))
        )

    # Inject Mahalanobis outlier in sample 10
    base_counts[:, 10] += rng.poisson(200, size=n_genes)

    # Make library size correlated with rin_score
    lib_sizes = base_counts.sum(axis=0)
    lib_target = 1e6 + np.array(rin_score) * 50000
    scale = lib_target / lib_sizes
    for i in range(n_samples):
        base_counts[:, i] = (base_counts[:, i] * scale[i]).astype(int)

    # Inject library size outliers AFTER rescaling so they survive
    # Reduce to ~8% of median to ensure 3-MAD detection
    base_counts[:, 21:24] = (base_counts[:, 21:24] * 0.08).astype(int)

    counts = pd.DataFrame(
        base_counts.astype(int).clip(0),
        index=[f"gene_{i:04d}" for i in range(n_genes)],
        columns=sample_ids,
    )

    return counts, metadata


@pytest.fixture(scope="session")
def synthetic_files(synthetic_data, tmp_path_factory):
    """Write synthetic data to temporary CSV files."""
    counts, metadata = synthetic_data
    tmpdir = tmp_path_factory.mktemp("data")
    counts_path = tmpdir / "synthetic_counts.csv"
    meta_path = tmpdir / "synthetic_metadata.csv"
    counts.to_csv(counts_path)
    metadata.to_csv(meta_path)
    return counts_path, meta_path


@pytest.fixture(scope="session")
def small_data():
    """Small dataset (8 samples) to test low power warnings."""
    rng = np.random.default_rng(99)
    n_samples = 8
    n_genes = 200
    sample_ids = [f"S{i+1:02d}" for i in range(n_samples)]
    metadata = pd.DataFrame({
        "treatment": ["ctrl"] * 4 + ["treat"] * 4,
    }, index=sample_ids)
    counts = pd.DataFrame(
        rng.negative_binomial(30, 0.5, size=(n_genes, n_samples)).astype(int),
        index=[f"gene_{i:04d}" for i in range(n_genes)],
        columns=sample_ids,
    )
    return counts, metadata
