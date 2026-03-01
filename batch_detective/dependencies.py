"""Dependency checking for batch-detective."""

import logging
from .exceptions import BatchDetectiveDependencyError

logger = logging.getLogger(__name__)

# Hard requirements (must be present)
REQUIRED_DEPS = [
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
    ("matplotlib", "matplotlib"),
    ("seaborn", "seaborn"),
    ("click", "click"),
    ("jinja2", "jinja2"),
    ("yaml", "pyyaml"),
]

# Optional but recommended
SOFT_DEPS = [
    ("tqdm", "tqdm"),
    ("pingouin", "pingouin"),
]


def check_dependencies() -> None:
    """Check that all required dependencies are importable.

    Raises:
        BatchDetectiveDependencyError: If any required dependency is missing.
    """
    missing = []
    for import_name, pkg_name in REQUIRED_DEPS:
        try:
            __import__(import_name)
        except ImportError:
            missing.append(pkg_name)

    if missing:
        raise BatchDetectiveDependencyError(
            f"Missing required dependencies: {', '.join(missing)}\n"
            f"Install with: pip install batch-detective"
        )

    for import_name, pkg_name in SOFT_DEPS:
        try:
            __import__(import_name)
        except ImportError:
            logger.debug(f"Optional dependency not available: {pkg_name}")

    logger.debug("Dependency check passed.")
