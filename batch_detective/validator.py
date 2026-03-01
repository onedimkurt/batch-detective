"""Stage 1: Structural input validation for batch-detective."""

import csv
import io
import logging
import re
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd

from .exceptions import BatchDetectiveValidationError

logger = logging.getLogger(__name__)

FEATURECOUNTS_COLS = {"chr", "start", "end", "strand", "length"}
STAR_PREFIXES = ("N_",)
HTSEQ_PREFIXES = ("__",)

GENE_SYMBOL_PATTERN = re.compile(
    r"^(ENSG\d|ENSMUSG\d|[A-Z][A-Z0-9]{1,9}$)", re.IGNORECASE
)


def _detect_delimiter(filepath: Path) -> str:
    """Auto-detect CSV or TSV delimiter."""
    with open(filepath, encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        return ","


def _load_counts(filepath: Path) -> pd.DataFrame:
    """Load count matrix from CSV or TSV."""
    delim = _detect_delimiter(filepath)
    df = pd.read_csv(
        filepath,
        index_col=0,
        sep=delim,
        encoding="utf-8",
        comment="#",
    )
    return df


def _load_metadata(filepath: Path) -> pd.DataFrame:
    """Load metadata from CSV or TSV."""
    delim = _detect_delimiter(filepath)
    df = pd.read_csv(
        filepath,
        index_col=0,
        sep=delim,
        encoding="utf-8",
    )
    return df


def validate_inputs(
    counts_path: Path,
    metadata_path: Path,
    force: bool = False,
    normalized: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Validate and load input files.

    Args:
        counts_path: Path to count matrix CSV/TSV.
        metadata_path: Path to metadata CSV/TSV.
        force: If True, proceed past non-fatal warnings.
        normalized: If True, allow non-integer counts.

    Returns:
        Tuple of (counts_raw, metadata_raw, counts_working, metadata_working).

    Raises:
        BatchDetectiveValidationError: On structural issues.
    """
    # Check extensions
    if counts_path.suffix.lower() in (".xlsx", ".xls"):
        raise BatchDetectiveValidationError(
            "Excel files are not supported. Please convert to CSV:\n"
            "  In Excel: File → Save As → CSV (Comma delimited).\n"
            "  In Python: pd.read_excel('file.xlsx').to_csv('file.csv')"
        )
    if metadata_path.suffix.lower() in (".xlsx", ".xls"):
        raise BatchDetectiveValidationError(
            "Excel files are not supported. Please convert to CSV:\n"
            "  In Excel: File → Save As → CSV (Comma delimited).\n"
            "  In Python: pd.read_excel('file.xlsx').to_csv('file.csv')"
        )

    # Check files exist
    if not counts_path.exists():
        raise BatchDetectiveValidationError(
            f"Counts file not found: {counts_path}"
        )
    if not metadata_path.exists():
        raise BatchDetectiveValidationError(
            f"Metadata file not found: {metadata_path}"
        )

    # Load
    counts_raw = _load_counts(counts_path)
    metadata_raw = _load_metadata(metadata_path)

    # Detect featureCounts annotation columns
    fc_cols = [
        c for c in counts_raw.columns
        if c.lower() in FEATURECOUNTS_COLS
    ]
    if fc_cols:
        logger.warning(
            f"Detected featureCounts format. Annotation columns "
            f"({', '.join(fc_cols)}) removed automatically. "
            f"Proceeding with sample columns only."
        )
        counts_raw = counts_raw.drop(columns=fc_cols)

    # Detect STAR/HTSeq special summary rows
    bad_rows = [
        idx for idx in counts_raw.index
        if (
            isinstance(idx, str)
            and (
                any(idx.startswith(p) for p in STAR_PREFIXES)
                or any(idx.startswith(p) for p in HTSEQ_PREFIXES)
            )
        )
    ]
    if bad_rows:
        msg = (
            f"Count matrix contains aligner summary rows: {bad_rows}\n"
            "These are NOT genes and will corrupt library size estimates and normalization.\n"
            "Please remove these rows before analysis.\n"
            f"Affected row indices: {bad_rows}"
        )
        if force:
            logger.warning(msg + "\n--force passed: proceeding (results may be unreliable).")
        else:
            raise BatchDetectiveValidationError(msg)

    # Check for NaN
    if counts_raw.isna().any().any():
        nan_genes = counts_raw.index[counts_raw.isna().any(axis=1)].tolist()
        nan_samples = counts_raw.columns[counts_raw.isna().any(axis=0)].tolist()
        raise BatchDetectiveValidationError(
            f"Count matrix contains NaN values. Please remove or impute before analysis.\n"
            f"Affected genes: {nan_genes[:10]}{'...' if len(nan_genes) > 10 else ''}\n"
            f"Affected samples: {nan_samples[:10]}{'...' if len(nan_samples) > 10 else ''}"
        )

    # Check duplicate column/index names
    # Read raw header to catch duplicates before pandas auto-renames them
    with open(counts_path, encoding="utf-8", errors="replace") as _f:
        _header = _f.readline().rstrip("\n").split(_detect_delimiter(counts_path))
    _sample_cols = _header[1:]  # skip index column
    _seen, _dups = set(), []
    for _c in _sample_cols:
        if _c in _seen:
            _dups.append(_c)
        _seen.add(_c)
    if _dups:
        raise BatchDetectiveValidationError(
            f"Duplicate sample names in count matrix: {_dups}"
        )
    dup_samples = counts_raw.columns[counts_raw.columns.duplicated()].tolist()
    if dup_samples:
        raise BatchDetectiveValidationError(
            f"Duplicate sample names in count matrix: {dup_samples}"
        )
    dup_genes = counts_raw.index[counts_raw.index.duplicated()].tolist()
    if dup_genes:
        raise BatchDetectiveValidationError(
            f"Duplicate gene names in count matrix: {dup_genes[:20]}"
        )
    dup_meta = metadata_raw.index[metadata_raw.index.duplicated()].tolist()
    if dup_meta:
        raise BatchDetectiveValidationError(
            f"Duplicate sample names in metadata: {dup_meta}"
        )

    # Check for negative values
    if (counts_raw.values < 0).any():
        raise BatchDetectiveValidationError(
            "Count matrix contains negative values. Counts must be >= 0."
        )

    # Check non-integer values
    try:
        numeric_vals = counts_raw.values.astype(float)
        non_int_frac = (numeric_vals % 1 != 0).mean()
        if non_int_frac > 0.05 and not normalized and not force:
            raise BatchDetectiveValidationError(
                f"Count matrix contains non-integer values ({non_int_frac:.1%} of entries).\n"
                "If you have pre-normalized data, use --normalized flag.\n"
                "To proceed anyway: use --force."
            )
        elif non_int_frac > 0 and normalized:
            logger.warning(
                "Pre-normalized data detected (--normalized flag). "
                "CPM normalization will be skipped."
            )
    except (TypeError, ValueError):
        raise BatchDetectiveValidationError(
            "Count matrix contains non-numeric values."
        )

    # Auto-detect transposed matrix
    n_rows, n_cols = counts_raw.shape
    if n_rows < n_cols:
        col_names = list(counts_raw.columns)
        gene_like = sum(
            1 for c in col_names
            if (
                str(c).startswith("ENSG")
                or str(c).startswith("ENSMUSG")
                or GENE_SYMBOL_PATTERN.match(str(c))
            )
        )
        if gene_like > len(col_names) * 0.3:
            logger.warning(
                f"Auto-detected transposed matrix (rows={n_rows} < cols={n_cols} "
                f"and columns look like gene IDs). Transposing."
            )
            counts_raw = counts_raw.T

    # Sample overlap
    counts_samples = set(counts_raw.columns)
    meta_samples = set(metadata_raw.index)
    intersection = counts_samples & meta_samples

    in_counts_not_meta = sorted(counts_samples - meta_samples)
    in_meta_not_counts = sorted(meta_samples - counts_samples)

    if in_counts_not_meta:
        logger.warning(f"In counts but not metadata: {in_counts_not_meta}")
    if in_meta_not_counts:
        logger.warning(f"In metadata but not counts: {in_meta_not_counts}")

    min_size = min(len(counts_samples), len(meta_samples))
    if len(intersection) < min_size * 0.5:
        msg = (
            f"Intersection of count samples and metadata samples "
            f"({len(intersection)}) is less than 50% of either file "
            f"({len(counts_samples)} counts, {len(meta_samples)} metadata)."
        )
        if force:
            logger.warning(msg + " --force passed: proceeding.")
        else:
            raise BatchDetectiveValidationError(msg)

    if len(intersection) < min_size * 0.8:
        pct_lost = (1 - len(intersection) / min_size) * 100
        logger.warning(
            f"Sample intersection loses {pct_lost:.0f}% of samples. "
            "Check that sample IDs match between files."
        )

    # Subset to intersection
    common = sorted(intersection)
    counts_raw = counts_raw[common]
    metadata_raw = metadata_raw.loc[common]

    counts_working = counts_raw.copy()
    metadata_working = metadata_raw.copy()

    logger.info(
        f"Validated inputs: {counts_raw.shape[0]} genes × "
        f"{counts_raw.shape[1]} samples"
    )
    return counts_raw, metadata_raw, counts_working, metadata_working
