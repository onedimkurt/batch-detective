"""Run manifest generation for batch-detective."""

import hashlib
import json
import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import __version__

logger = logging.getLogger(__name__)


def compute_md5(filepath: Path) -> str:
    """Compute MD5 checksum of a file.

    Args:
        filepath: Path to file.

    Returns:
        Hex MD5 string.
    """
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _get_dep_versions() -> Dict[str, str]:
    """Get installed versions of key dependencies."""
    deps = {}
    for pkg, import_name in [
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("scikit-learn", "sklearn"),
        ("pingouin", "pingouin"),
    ]:
        try:
            mod = __import__(import_name)
            deps[pkg] = getattr(mod, "__version__", "unknown")
        except ImportError:
            deps[pkg] = "not installed"
    return deps


def write_manifest(
    output_dir: Path,
    counts_path: Path,
    metadata_path: Path,
    resolved_params: Dict[str, Any],
    analysis_summary: Dict[str, Any],
    exit_code: int,
    anonymized: bool = False,
) -> Path:
    """Write run_manifest.json to output directory.

    Args:
        output_dir: Output directory path.
        counts_path: Absolute path to counts file.
        metadata_path: Absolute path to metadata file.
        resolved_params: Final merged parameters.
        analysis_summary: Analysis summary dict.
        exit_code: Final exit code.
        anonymized: Whether --anonymize-samples was used.

    Returns:
        Path to written manifest file.
    """
    manifest = {
        "tool": "batch-detective",
        "version": __version__,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "platform": f"{platform.system()}-{platform.release()}",
        "python_version": sys.version.split()[0],
        "dependency_versions": _get_dep_versions(),
        "inputs": {
            "counts_path": str(counts_path.resolve()),
            "metadata_path": str(metadata_path.resolve()),
            "counts_md5": compute_md5(counts_path),
            "metadata_md5": compute_md5(metadata_path),
        },
        "resolved_parameters": resolved_params,
        "analysis_summary": analysis_summary,
        "anonymization_applied": anonymized,
        "exit_code": exit_code,
    }

    manifest_path = output_dir / "run_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, default=str)

    logger.info(f"Run manifest written to {manifest_path}")
    return manifest_path
