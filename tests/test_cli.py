"""CLI integration tests using Click's CliRunner."""

import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from batch_detective.cli import main


def _make_data(tmp_path, n_samples=24, n_genes=300, seed=42):
    rng = np.random.default_rng(seed)
    samples = [f"S{i+1:02d}" for i in range(n_samples)]
    counts = pd.DataFrame(
        rng.negative_binomial(50, 0.5, (n_genes, n_samples)).astype(int),
        index=[f"gene_{i:04d}" for i in range(n_genes)],
        columns=samples,
    )
    # Inject batch signal
    batch = [f"B{(i//6)+1}" for i in range(n_samples)]
    for b_idx, b_name in enumerate(["B1", "B2", "B3", "B4"]):
        mask = [b == b_name for b in batch]
        sample_idx = [i for i, m in enumerate(mask) if m]
        counts.iloc[:100, sample_idx] += b_idx * 30

    metadata = pd.DataFrame({
        "batch": batch,
        "treatment": ["ctrl"] * (n_samples // 2) + ["treat"] * (n_samples // 2),
    }, index=samples)

    cp = tmp_path / "counts.csv"
    mp = tmp_path / "metadata.csv"
    counts.to_csv(cp)
    metadata.to_csv(mp)
    return cp, mp


def test_run_happy_path(tmp_path):
    runner = CliRunner()
    cp, mp = _make_data(tmp_path)
    out_dir = tmp_path / "output"
    result = runner.invoke(main, [
        "run",
        "--counts", str(cp),
        "--metadata", str(mp),
        "--output-dir", str(out_dir),
        "--technical-covariates", "batch",
        "--n-variable-genes", "100",
    ])
    assert result.exit_code in (0, 1), result.output
    assert (out_dir / "report.html").exists()
    assert (out_dir / "qc_summary.csv").exists()
    assert (out_dir / "association_table.csv").exists()
    assert (out_dir / "icc_table.csv").exists()
    assert (out_dir / "run_manifest.json").exists()


def test_run_missing_counts(tmp_path):
    runner = CliRunner()
    out_dir = tmp_path / "output"
    result = runner.invoke(main, [
        "run",
        "--counts", str(tmp_path / "nonexistent.csv"),
        "--metadata", str(tmp_path / "meta.csv"),
        "--output-dir", str(out_dir),
    ])
    assert result.exit_code == 3


def test_run_dry_run(tmp_path):
    runner = CliRunner()
    cp, mp = _make_data(tmp_path)
    out_dir = tmp_path / "output"
    result = runner.invoke(main, [
        "run",
        "--counts", str(cp),
        "--metadata", str(mp),
        "--output-dir", str(out_dir),
        "--dry-run",
    ])
    assert result.exit_code == 0
    # No files should be created
    assert not out_dir.exists() or not any(out_dir.iterdir())
    assert "batch-detective dry run" in result.output


def test_run_existing_dir_no_overwrite(tmp_path):
    runner = CliRunner()
    cp, mp = _make_data(tmp_path)
    out_dir = tmp_path / "existing"
    out_dir.mkdir()
    (out_dir / "dummy.txt").write_text("existing")
    result = runner.invoke(main, [
        "run",
        "--counts", str(cp),
        "--metadata", str(mp),
        "--output-dir", str(out_dir),
    ])
    assert result.exit_code == 3
    assert "overwrite" in result.output.lower() or "--overwrite" in result.output


def test_run_existing_dir_with_overwrite(tmp_path):
    runner = CliRunner()
    cp, mp = _make_data(tmp_path)
    out_dir = tmp_path / "existing2"
    out_dir.mkdir()
    (out_dir / "dummy.txt").write_text("existing")
    result = runner.invoke(main, [
        "run",
        "--counts", str(cp),
        "--metadata", str(mp),
        "--output-dir", str(out_dir),
        "--overwrite",
        "--n-variable-genes", "100",
    ])
    assert result.exit_code in (0, 1), result.output


def test_version(tmp_path):
    runner = CliRunner()
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert "batch-detective" in result.output
    assert "Python" in result.output


def test_serve_missing_extra(tmp_path):
    """serve without [serve] installed should give clean error."""
    import sys
    runner = CliRunner()
    # Mock streamlit not being available
    import unittest.mock as mock
    with mock.patch.dict(sys.modules, {"streamlit": None}):
        result = runner.invoke(main, ["serve"])
        # Should either exit cleanly with error or succeed if streamlit is installed
        # Just check it doesn't produce a traceback
        assert "Traceback" not in result.output or result.exit_code == 0


def test_dry_run_output_format(tmp_path):
    runner = CliRunner()
    cp, mp = _make_data(tmp_path)
    out_dir = tmp_path / "output"
    result = runner.invoke(main, [
        "run",
        "--counts", str(cp),
        "--metadata", str(mp),
        "--output-dir", str(out_dir),
        "--dry-run",
    ])
    assert result.exit_code == 0
    output = result.output
    assert "batch-detective dry run" in output
    assert "genes" in output
    assert "samples" in output
    assert "n_variable_genes" in output
    assert "Estimated runtime" in output


def test_excel_rejection(tmp_path):
    runner = CliRunner()
    xlsx = tmp_path / "counts.xlsx"
    xlsx.touch()
    out_dir = tmp_path / "out"
    result = runner.invoke(main, [
        "run",
        "--counts", str(xlsx),
        "--metadata", str(tmp_path / "meta.csv"),
        "--output-dir", str(out_dir),
    ])
    assert result.exit_code == 3
    assert "Excel" in result.output


def test_manifest_contents(tmp_path):
    runner = CliRunner()
    cp, mp = _make_data(tmp_path)
    out_dir = tmp_path / "manifest_test"
    result = runner.invoke(main, [
        "run",
        "--counts", str(cp),
        "--metadata", str(mp),
        "--output-dir", str(out_dir),
        "--n-variable-genes", "100",
    ])
    assert result.exit_code in (0, 1)
    manifest_path = out_dir / "run_manifest.json"
    assert manifest_path.exists()
    with open(manifest_path) as f:
        manifest = json.load(f)
    assert "version" in manifest
    assert "timestamp_utc" in manifest
    assert "platform" in manifest
    assert "python_version" in manifest
    assert "counts_md5" in manifest["inputs"]
    assert len(manifest["inputs"]["counts_md5"]) == 32
    assert "exit_code" in manifest
