"""Configuration loading and merging for batch-detective."""

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from .exceptions import BatchDetectiveValidationError

logger = logging.getLogger(__name__)

DEFAULTS = {
    "n_variable_genes": 2000,
    "n_pcs": 10,
    "min_cpm": 1.0,
    "min_samples_expressing": None,
    "outlier_pval": 0.001,
    "normalized": False,
    "primary_covariate": None,
    "condition_on": None,
    "technical_covariates": [],
    "biological_covariates": [],
    "dry_run": False,
    "overwrite": False,
    "verbose": False,
    "quiet": False,
    "pdf": False,
    "dpi": 150,
    "export_figures": False,
    "anonymize_samples": False,
    "force": False,
    "log_file": None,
}


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load configuration from YAML file.

    Args:
        config_path: Path to the YAML config file.

    Returns:
        Dict of configuration parameters.

    Raises:
        BatchDetectiveValidationError: If config file not found.
    """
    if not config_path.exists():
        raise BatchDetectiveValidationError(
            f"Config file not found: {config_path}"
        )

    with open(config_path, encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}

    config_dir = config_path.parent.resolve()

    # Resolve relative paths in config relative to config file directory
    for key in ("counts", "metadata", "output_dir"):
        if key in raw and raw[key] is not None:
            raw[key] = (config_dir / raw[key]).resolve()

    # Warn on unknown keys
    known_keys = set(DEFAULTS.keys()) | {
        "counts", "metadata", "output_dir", "config", "log_file"
    }
    for key in raw:
        if key not in known_keys:
            logger.warning(f"Unknown config key: {key} — ignored.")

    return raw


def merge_config(
    cli_params: Dict[str, Any],
    config_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Merge CLI parameters with config file and defaults.

    Priority: CLI > config file > defaults.

    Args:
        cli_params: Parameters from CLI (None values = not set by user).
        config_path: Optional path to YAML config file.

    Returns:
        Merged configuration dict.
    """
    merged = dict(DEFAULTS)

    if config_path is not None:
        file_config = load_config(config_path)
        for key, val in file_config.items():
            if key in merged or key in ("counts", "metadata", "output_dir"):
                merged[key] = val

    # CLI overrides (only if not None / not default sentinel)
    for key, val in cli_params.items():
        if val is not None:
            merged[key] = val

    return merged
