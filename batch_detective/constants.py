"""
batch_detective/constants.py
============================
Central definition of all detection thresholds used by batch-detective.

Empirical basis
---------------
These values were calibrated against the 6-dataset validation panel
(see tests/validation/README.md).  Key reference points:

  Dataset 1  (GSE83083 + GSE59765, strong batch)   → ICC = 0.878
  Dataset 3  (GSE271332, large, moderate batch)     → ICC = 0.194
  Dataset 4  (GSE120099, multi-batch, weak)         → ICC = 0.075
  Dataset 6  (GSE185263, clinical, weak)            → ICC = 0.074
  Dataset 2  (GSE48035, true negative)              → ICC = N/A

The ICC_STRONG threshold of 0.10 sits above datasets 4 and 6 (weak batch)
and well below dataset 1 (strong batch), providing a clean separation.

Cramér's V collinearity threshold of 0.70 is the conventional "strong
association" cut-off in the social-science literature and was validated
against Dataset 5 (GFRN confounded, V ≈ 1.0) which correctly triggers
the warning, while Dataset 3 (V = 0.289) and Dataset 4 (V = 0.171)
correctly do not.
"""

# ---------------------------------------------------------------------------
# ICC thresholds
# ---------------------------------------------------------------------------

#: Batch effect strong enough to warrant correction before any downstream
#: analysis.  Corresponds to >10 % of within-gene variance explained by batch.
ICC_STRONG: float = 0.10

#: Batch effect moderate — correction recommended, especially for
#: differential expression.
ICC_MODERATE: float = 0.05

#: Weak batch signal.  Monitor; may not require correction for all analyses.
ICC_WEAK: float = 0.01

# ---------------------------------------------------------------------------
# Collinearity / confounding threshold
# ---------------------------------------------------------------------------

#: Cramér's V between a technical covariate and a biological covariate above
#: this value triggers a collinearity warning.  At this level the two variables
#: are strongly associated and batch correction may remove biological signal.
CRAMERS_V_COLLINEARITY: float = 0.70

# ---------------------------------------------------------------------------
# Exit codes
# ---------------------------------------------------------------------------

EXIT_OK: int = 0          # No batch effect detected
EXIT_BATCH: int = 1       # Batch effect detected, correction recommended
EXIT_CONFOUNDED: int = 2  # Batch effect detected but design is confounded

# ---------------------------------------------------------------------------
# Minimum samples per batch for reliable ICC estimation
# ---------------------------------------------------------------------------

MIN_SAMPLES_PER_BATCH: int = 3

#: Minimum number of batches required to compute ICC
MIN_BATCHES: int = 2

#: Cramér's V threshold for moderate association warning (continuous vs categorical)
CRAMERS_V_MODERATE: float = 0.50

#: ICC upper bounds for traffic-light classification
ICC_MODERATE_UPPER: float = 0.30
ICC_STRONG_UPPER: float = 0.60
