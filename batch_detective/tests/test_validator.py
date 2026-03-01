"""Tests for validator module."""

import pandas as pd
import numpy as np
import pytest
from pathlib import Path
import tempfile

from batch_detective.exceptions import BatchDetectiveValidationError


def write_csv(df, path):
    df.to_csv(path)


def test_validate_happy_path(synthetic_files):
    from batch_detective.validator import validate_inputs
    counts_path, meta_path = synthetic_files
    counts_raw, metadata_raw, counts_working, metadata_working = validate_inputs(
        counts_path, meta_path
    )
    assert counts_raw.shape[1] > 0
    assert metadata_raw.shape[0] > 0
    # Working copies are distinct objects
    assert counts_working is not counts_raw


def test_excel_rejection(tmp_path):
    from batch_detective.validator import validate_inputs
    xlsx_path = tmp_path / "counts.xlsx"
    xlsx_path.touch()
    with pytest.raises(BatchDetectiveValidationError, match="Excel"):
        validate_inputs(xlsx_path, tmp_path / "meta.csv")


def test_missing_file(tmp_path):
    from batch_detective.validator import validate_inputs
    with pytest.raises(BatchDetectiveValidationError, match="not found"):
        validate_inputs(tmp_path / "nonexistent.csv", tmp_path / "meta.csv")


def test_nan_in_counts(tmp_path):
    from batch_detective.validator import validate_inputs
    counts = pd.DataFrame({"S1": [1.0, np.nan], "S2": [3.0, 4.0]}, index=["g1", "g2"])
    meta = pd.DataFrame({"batch": ["A", "A"]}, index=["S1", "S2"])
    cp = tmp_path / "counts.csv"
    mp = tmp_path / "meta.csv"
    write_csv(counts, cp)
    write_csv(meta, mp)
    with pytest.raises(BatchDetectiveValidationError, match="NaN"):
        validate_inputs(cp, mp)


def test_duplicate_samples(tmp_path):
    from batch_detective.validator import validate_inputs
    # Create CSV manually with duplicate column
    cp = tmp_path / "counts.csv"
    cp.write_text("gene,S1,S1\ng1,1,2\ng2,3,4\n")
    mp = tmp_path / "meta.csv"
    meta = pd.DataFrame({"batch": ["A", "B"]}, index=["S1", "S2"])
    write_csv(meta, mp)
    with pytest.raises(BatchDetectiveValidationError, match="Duplicate"):
        validate_inputs(cp, mp)


def test_negative_values(tmp_path):
    from batch_detective.validator import validate_inputs
    counts = pd.DataFrame({"S1": [1, -1], "S2": [3, 4]}, index=["g1", "g2"])
    meta = pd.DataFrame({"batch": ["A", "A"]}, index=["S1", "S2"])
    cp, mp = tmp_path / "c.csv", tmp_path / "m.csv"
    write_csv(counts, cp)
    write_csv(meta, mp)
    with pytest.raises(BatchDetectiveValidationError, match="negative"):
        validate_inputs(cp, mp)


def test_star_summary_rows(tmp_path):
    from batch_detective.validator import validate_inputs
    cp = tmp_path / "counts.csv"
    cp.write_text(",S1,S2\nN_unmapped,1,2\ngene1,10,20\n")
    mp = tmp_path / "meta.csv"
    meta = pd.DataFrame({"batch": ["A", "B"]}, index=["S1", "S2"])
    write_csv(meta, mp)
    with pytest.raises(BatchDetectiveValidationError, match="N_unmapped"):
        validate_inputs(cp, mp)


def test_star_force_passes(tmp_path):
    from batch_detective.validator import validate_inputs
    cp = tmp_path / "counts.csv"
    cp.write_text(",S1,S2\nN_unmapped,1,2\ngene1,10,20\n")
    mp = tmp_path / "meta.csv"
    meta = pd.DataFrame({"batch": ["A", "B"]}, index=["S1", "S2"])
    write_csv(meta, mp)
    # With force, should not raise
    validate_inputs(cp, mp, force=True)


def test_featurecounts_annotation_removed(tmp_path):
    from batch_detective.validator import validate_inputs
    cp = tmp_path / "counts.csv"
    cp.write_text(",Chr,Start,End,Strand,Length,S1,S2\ngene1,chr1,100,200,+,100,10,20\n")
    mp = tmp_path / "meta.csv"
    meta = pd.DataFrame({"batch": ["A", "B"]}, index=["S1", "S2"])
    write_csv(meta, mp)
    counts_raw, *_ = validate_inputs(cp, mp)
    assert "Chr" not in counts_raw.columns
    assert "S1" in counts_raw.columns
