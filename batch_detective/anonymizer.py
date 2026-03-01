"""Sample anonymization for batch-detective."""

import logging
from typing import Dict, List

import pandas as pd

logger = logging.getLogger(__name__)


def anonymize_samples(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple:
    """Replace sample IDs with Sample_001, Sample_002, etc.

    Args:
        counts: Count matrix with samples as columns.
        metadata: Metadata with samples as index.

    Returns:
        Tuple of (anonymized_counts, anonymized_metadata, id_mapping).
    """
    original_ids = list(counts.columns)
    anon_ids = [f"Sample_{i+1:03d}" for i in range(len(original_ids))]
    id_mapping: Dict[str, str] = dict(zip(original_ids, anon_ids))

    anon_counts = counts.copy()
    anon_counts.columns = anon_ids

    anon_metadata = metadata.copy()
    anon_metadata.index = [id_mapping.get(idx, idx) for idx in metadata.index]

    logger.info(f"Anonymized {len(original_ids)} sample IDs.")
    return anon_counts, anon_metadata, id_mapping
