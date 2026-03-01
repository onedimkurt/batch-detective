# batch-detective Validation Panel

Six real-world GEO datasets that cover the core detection scenarios batch-detective
is designed to handle. All datasets use **NCBI-generated count matrices** (uniform
alignment pipeline) so differences between datasets reflect biology/batch, not
pipeline variation.

Run the full panel from the project root:

```bash
python tests/validation/run_validation_panel.py
```

Results are written to `tests/validation/results/validation_panel_results.json`
and a human-readable HTML report to `tests/validation/results/validation_report.html`.

---

## Dataset Inventory

### Dataset 1 — Strong known batch effect
| Field | Value |
|-------|-------|
| Accessions | GSE83083 (batch 1) + GSE59765 (batch 2) |
| Samples | 29 (17 HER2/GFP + 12 EGFR/GFP) |
| Batch variable | study of origin (two independent labs) |
| Bio variable | oncogene treatment (HER2 vs EGFR vs GFP control) |
| Expected | exit 1, ICC ≥ 0.5 |
| Observed | exit 1, ICC = 0.878 ✅ |

**Rationale:** The canonical ComBat-seq benchmark (Zhang et al. 2020, NAR Genomics &
Bioinformatics). Two completely independent GEO submissions of the same experiment
type, merged post-hoc. Lab-of-origin is the batch. ICC of 0.878 confirms strong
inter-batch variance — the tool's primary detection use case.

---

### Dataset 2 — True negative (no batch effect)
| Field | Value |
|-------|-------|
| Accession | GSE48035 |
| Samples | 20 |
| Batch variable | none (single sequencing run) |
| Bio variable | cell type |
| Expected | exit 0, no batch variable detected |
| Observed | exit 0, ICC = N/A ✅ |

**Rationale:** Single-batch study. Verifies the tool does not produce false positives
and exits cleanly with code 0 when no technical batch variable is present.

---

### Dataset 3 — Large dataset (>100 samples)
| Field | Value |
|-------|-------|
| Accession | GSE271332 |
| Samples | 286 (173 with NCBI counts) |
| Batch variable | `batch` (1 / 2 / 3) — author-labeled |
| Bio variable | chemical treatment (11 bisphenol analogues × 5 doses) |
| Cramér's V | 0.289 (balanced, not confounded) |
| Expected | exit 1, ICC > 0, no crash, no false collinearity warning |
| Observed | exit 1, ICC = 0.194 ✅ |

**Rationale:** High-throughput toxicology screen (MCF-7 cells). 286 samples with
explicit author-provided batch labels across 3 sequencing runs. Chemicals are
distributed across all batches (Cramér's V = 0.289), so this tests large-matrix
handling and true batch detection without confounding. Confirms the tool scales
beyond 100 samples without crash or memory error.

---

### Dataset 4 — Multi-batch (3+ batches)
| Field | Value |
|-------|-------|
| Accession | GSE120099 |
| Samples | 92 |
| Batch variable | `sequencing batch` (Set1 / Set2 / Set3) — author-labeled |
| Bio variable | genotype (KO / WT) |
| Cramér's V | 0.171 (balanced) |
| Expected | exit 1, 3 batches detected, ICC > 0 |
| Observed | exit 1, 3 batches, ICC = 0.075 ✅ |

**Rationale:** Hematopoiesis study with the batch variable explicitly named
`sequencing batch` by the authors in the GEO metadata. Confirms the tool correctly
handles ≥3 batches and reports per-batch ICC. The low ICC (0.075) is expected —
the authors designed a balanced experiment.

---

### Dataset 5 — Confounded design
| Field | Value |
|-------|-------|
| Accessions | GSE83083 (batches 1 + 3) + GSE59765 (batch 2) |
| Samples | 47 (HER2×5, GFP18×12, KRAS×9, GFP30×9, EGFR×6, GFP59×6) |
| Batch variable | study/sequencing batch (batch1_HER2 / batch2_EGFR / batch3_KRAS) |
| Bio variable | oncogene treatment |
| Cramér's V | ~1.0 (each treatment nested in exactly one batch) |
| Expected | exit 1 or 2, collinearity warning fired |
| Observed | exit 1, ICC = 0.878, collinearity warning = True ✅ |

**Rationale:** The full 3-batch GFRN dataset from the original ComBat-seq paper.
Each oncogene treatment (HER2, EGFR, KRAS) is completely nested within a single
sequencing batch — a textbook confounded design. GFP controls appear in all three
batches providing the only shared reference. Cramér's V (batch vs treatment) = 1.0
by construction. Verifies the collinearity detection logic fires correctly and the
tool warns the user rather than silently proceeding.

**Sample selection note:** GSE83083 contains additional oncogenes (IGF1R, RAF1,
BAD, AKT1) that are NOT part of the ComBat-seq 3-batch subset. These are excluded
by title-keyword matching (GFP18/GFP30 distinguish which GFP control group belongs
to which batch).

---

### Dataset 6 — Real-world multi-batch clinical
| Field | Value |
|-------|-------|
| Accession | GSE185263 |
| Samples | varies |
| Batch variable | sequencing batch |
| Bio variable | sepsis status |
| Expected | exit 1, ICC > 0 |
| Observed | exit 1, ICC = 0.074 ✅ |

**Rationale:** Clinical sepsis study from a published multi-centre cohort.
Represents the most common real-world use case: a clinical researcher who has
combined samples from multiple sequencing runs and wants to know whether batch
correction is needed before differential expression.

---

## ICC Thresholds

These thresholds are defined in `batch_detective/constants.py` and used by the
detection and reporting logic:

| Threshold | Value | Meaning |
|-----------|-------|---------|
| `ICC_STRONG` | 0.10 | Batch effect present, correction strongly recommended |
| `ICC_MODERATE` | 0.05 | Batch effect present, correction recommended |
| `ICC_WEAK` | 0.01 | Minor batch signal, monitor |
| `CRAMERS_V_COLLINEARITY` | 0.70 | Warn that batch and biology are confounded |

**Empirical basis:** Dataset 1 (strong batch, ICC=0.878) and Dataset 6 (weak batch,
ICC=0.074) bracket the realistic range. Dataset 4 (ICC=0.075) sits just above the
`ICC_STRONG` threshold. Dataset 2 (true negative, ICC=N/A) anchors the lower bound.

---

## Adding New Datasets

To add a dataset to the panel:

1. Confirm NCBI-generated counts exist:
   ```bash
   curl -sI "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSEXXXXXX&format=file&file=GSEXXXXXX_raw_counts_GRCh38.p13_NCBI.tsv.gz" | grep HTTP
   ```
2. Check that `Sample_characteristics_ch1` in the SOFT file contains a parseable
   batch variable and a biological covariate.
3. Cross-tabulate batch vs biology and compute Cramér's V — document the value.
4. Add an entry to `run_validation_panel.py` following the existing pattern.
5. Update this README with the dataset details and rationale.

---

## References

- Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect adjustment for
  RNA-seq count data. *NAR Genomics and Bioinformatics* 2020;2(3):lqaa078.
  https://doi.org/10.1093/nargab/lqaa078
- NCBI GEO RNA-seq counts pipeline:
  https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html
