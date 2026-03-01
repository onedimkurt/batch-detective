# batch-detective

**Statistical diagnosis of batch effects in bulk RNA-seq count matrices.**

[![PyPI version](https://img.shields.io/pypi/v/batch-detective)](https://pypi.org/project/batch-detective/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/onedimkurt/batch-detective/actions/workflows/tests.yml/badge.svg)](https://github.com/onedimkurt/batch-detective/actions)

---

`batch-detective` is a zero-configuration, read-only CLI tool that diagnoses technical batch effects in bulk RNA-seq count matrices. It accepts a raw count matrix and a sample metadata table, and produces a self-contained HTML report with a quantitative assessment of variance attributable to technical vs. biological covariates.

**It does not modify your data.** It only reads and reports.

---

## Table of Contents

- [What it does](#what-it-does)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input format](#input-format)
- [CLI reference](#cli-reference)
- [Interpreting results](#interpreting-results)
- [Output files](#output-files)
- [Config file](#config-file)
- [Pipeline integration](#pipeline-integration)
- [HPC / cluster usage](#hpc--cluster-usage)
- [Troubleshooting](#troubleshooting)
- [Known limitations](#known-limitations)
- [Citation](#citation)
- [License](#license)

---

## What it does

`batch-detective` runs a seven-stage analysis pipeline:

1. **Validates** input files — detects featureCounts annotation columns, STAR/HTSeq summary rows, mismatched sample IDs, non-integer counts, and transposed matrices
2. **Quality control** — library size outliers, dominant gene flags, zero saturation
3. **Normalizes** — CPM + log1p, MAD-based variable gene selection
4. **PCA** — mean-centered, top variable genes
5. **Collinearity detection** — flags when batch and biology are confounded before you spend time on correction
6. **PC–covariate associations** — Kruskal-Wallis η² (categorical) and Spearman ρ (continuous), BH FDR-corrected
7. **ICC** — Intraclass Correlation Coefficient ICC(1,1) per covariate, with bootstrap 95% CIs

**Key outputs in the HTML report:**
- 🔴🟡🟢 **Traffic light** — one-line summary for quick interpretation
- **ICC table** — headline batch effect magnitude per covariate
- **PC–covariate association heatmap** — which principal components are driven by which variables
- **PCA scatter plots** — colored by each covariate
- **Collinearity warnings** — with specific guidance when batch and treatment are confounded
- **Top batch-associated genes** — gene-level ICC for the highest-variance genes
- **Sample outlier report** — Mahalanobis distance + IQR detection
- **Next-steps recommendations** — DESeq2/edgeR/ComBat code snippets, context-aware

---

## Installation

### Base install (CLI only)

```bash
pip install batch-detective
```

### With progress bars

```bash
pip install "batch-detective[progress]"
```

### With local web UI

```bash
pip install "batch-detective[serve]"
```

### With PDF report export

```bash
pip install "batch-detective[pdf]"

# Additional system libraries required:
# Linux:   sudo apt install libpango-1.0-0 libcairo2 libgdk-pixbuf2.0-0
# macOS:   brew install pango cairo gdk-pixbuf
# Windows: https://doc.courtbouillon.org/weasyprint/stable/first_steps.html
```

### Full install (progress + stats extras)

```bash
pip install "batch-detective[full]"
```

### Developer install

```bash
git clone https://github.com/onedimkurt/batch-detective
cd batch-detective
pip install -e ".[dev]"
pytest  # run test suite
```

---

## Quick start

```bash
batch-detective run \
  --counts counts.csv \
  --metadata metadata.csv \
  --output-dir ./results/ \
  --technical-covariates batch \
  --biological-covariates treatment

# Open the report
open results/report.html        # macOS
xdg-open results/report.html    # Linux
```

**That's it.** The report is a single self-contained HTML file — no server needed.

### With multiple covariates

```bash
batch-detective run \
  --counts counts.csv \
  --metadata metadata.csv \
  --output-dir ./results/ \
  --technical-covariates "batch,plate,operator" \
  --biological-covariates "treatment,sex,age"
```

### Preview without running (dry run)

```bash
batch-detective run \
  --counts counts.csv \
  --metadata metadata.csv \
  --output-dir ./results/ \
  --dry-run
```

This prints detected covariates, sample counts, and resolved parameters — without running the analysis.

### Check tool is working

```bash
batch-detective validate
```

Runs on built-in synthetic data and confirms the tool is functioning correctly.

---

## Input format

### Count matrix (`--counts`)

A gene × sample count matrix in CSV or TSV format.

- **Rows** = genes (ENSEMBL IDs or gene symbols)
- **Columns** = samples
- **Values** = raw integer counts (not normalized)
- First column = gene IDs (used as row index)
- First row = sample IDs (used as column headers)

```
gene_id,Sample_A,Sample_B,Sample_C,Sample_D
ENSG00000000003,423,387,512,449
ENSG00000000005,0,1,0,2
ENSG00000000419,891,1023,876,934
```

**Accepted formats:** CSV, TSV (delimiter auto-detected)

**Not accepted:** Excel (`.xlsx`/`.xls`), pre-normalized matrices (use `--normalized`)

**Automatic handling:**
- featureCounts annotation columns (Chr, Start, End, Strand, Length) are detected and stripped automatically
- Transposed matrices (samples as rows) are detected and corrected automatically
- STAR (`N_*`) and HTSeq (`__*`) summary rows trigger a fatal error with instructions to remove them

### Metadata table (`--metadata`)

A sample × covariate table in CSV or TSV format.

- **Rows** = samples (must match sample IDs in counts)
- **Columns** = covariate names
- First column = sample IDs (used as row index)

```
sample_id,batch,treatment,sex
Sample_A,B1,control,M
Sample_B,B1,treated,F
Sample_C,B2,control,F
Sample_D,B2,treated,M
```

**Tips:**
- Sample IDs must match between files (order does not matter)
- Missing values are allowed — affected samples are excluded from that covariate's analysis
- Categorical covariates are detected automatically (< 10 unique values, or non-numeric)
- Continuous covariates are tested with Spearman correlation instead of Kruskal-Wallis

---

## CLI reference

### `batch-detective run`

```
Usage: batch-detective run [OPTIONS]

  Run batch effect diagnosis.

Options:
  --counts PATH                  Count matrix CSV/TSV (required)
  --metadata PATH                Sample metadata CSV/TSV (required)
  --output-dir PATH              Output directory (required)

  --normalized                   Input is pre-normalized; skip CPM step
  --n-variable-genes INT         Top variable genes for PCA [default: 2000]
  --n-pcs INT                    Number of principal components [default: 10]
  --min-cpm FLOAT                Minimum CPM for gene filtering [default: 1.0]
  --min-samples-expressing INT   Min samples expressing a gene [default: auto]
  --outlier-pval FLOAT           Mahalanobis outlier p-value [default: 0.001]

  --technical-covariates TEXT    Comma-separated list of technical covariate names
                                 (e.g. "batch,plate,operator")
  --biological-covariates TEXT   Comma-separated list of biological covariate names
                                 (e.g. "treatment,sex,age")
  --primary-covariate TEXT       Single covariate to condition on
  --condition-on TEXT            Comma-separated covariates to regress out of PCs
                                 before association testing

  --config PATH                  YAML config file (see Config file section)
  --dry-run                      Validate inputs and preview, no analysis run
  --overwrite                    Overwrite existing output directory
  --force                        Proceed past non-fatal validation warnings
  --anonymize-samples            Replace sample IDs with Sample_001, etc.

  --pdf                          Export PDF report (requires [pdf] extra)
  --dpi INT                      Figure resolution [default: 150 HTML, 300 PDF]
  --export-figures               Export all figures as PNG + SVG files

  --log-file PATH                Write full log to file
  --verbose                      Enable debug-level logging
  --quiet                        Suppress all non-error output

  --help                         Show this message and exit.
```

### `batch-detective serve`

```bash
batch-detective serve [--port INT]
# Default port: 8501
# Requires: pip install "batch-detective[serve]"
```

Launches a local Streamlit web UI for interactive analysis (upload files, explore results).

### `batch-detective validate`

```bash
batch-detective validate
```

Runs a built-in self-test on synthetic data (80 samples, 4 batches, known ICC ≈ 0.45). Exits 0 on success.

### Exit codes

| Code | Meaning |
|------|---------|
| `0` | Analysis complete — no significant batch effects detected |
| `1` | Analysis complete — **significant batch effects detected** (not an error) |
| `2` | Analysis complete with data quality warnings |
| `3` | Fatal error — validation failed or unrecoverable error |

> **Important:** Exit code `1` means the tool ran successfully and found batch effects. It is **not** an error. See [Pipeline integration](#pipeline-integration) for correct usage.

---

## Interpreting results

### ICC (Intraclass Correlation Coefficient)

ICC quantifies what fraction of transcriptome-wide variance is explained by a grouping variable, computed using ICC(1,1) — a one-way random effects model across all selected variable genes.

| ICC | Interpretation |
|-----|---------------|
| < 0.10 | Negligible |
| 0.10 – 0.30 | Mild |
| 0.30 – 0.60 | Moderate |
| > 0.60 | Strong |

*(Thresholds after Koo & Li, 2016)*

**High ICC for a biological variable** (treatment, disease status) means your experiment worked — it is expected and not a problem.

**High ICC for a technical variable** (batch, plate, operator, sequencing date) indicates a batch effect that may confound downstream analysis.

### Traffic light

The traffic light summarizes the highest technical-covariate ICC and PC associations:

| Signal | Meaning |
|--------|---------|
| 🟢 Green | No significant batch effects detected |
| 🟡 Yellow | Possible batch effect — review recommended |
| 🔴 Red | Significant batch effect — correction recommended |

The traffic light is only informative when `--technical-covariates` are specified. Without labeling, all covariates are reported but not classified.

### PC–covariate association heatmap

Shows effect size (η² for categorical, \|ρ\| for continuous) between each metadata covariate and the top 5 PCs. Asterisks show FDR significance: `*` q < 0.05, `**` q < 0.01, `***` q < 0.001.

### Collinearity warnings

When batch and treatment are strongly collinear (Cramér's V ≥ 0.70), the report warns that batch correction is **not recommended** — it would remove biological signal along with the batch effect.

### What to do when batch effects are found

**If batch and treatment are not collinear:**

```r
# DESeq2
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ batch + treatment)

# edgeR
design <- model.matrix(~ batch + treatment)

# Visualization only (do NOT use for DE):
library(sva)
adjusted <- ComBat(dat = log_counts, batch = batch_vector)
```

**If batch and treatment are collinear:** do not apply correction. Proceed with standard analysis and acknowledge the limitation.

---

## Output files

| File | Description |
|------|-------------|
| `report.html` | Self-contained HTML report (main output) |
| `qc_summary.csv` | Per-sample QC metrics (library size, outlier flags) |
| `icc_table.csv` | ICC results per covariate with bootstrap CIs |
| `association_table.csv` | PC–covariate association test results (p-values, effect sizes) |
| `run_manifest.json` | Full reproducibility record (parameters, MD5 checksums, versions) |
| `figures/` | PNG + SVG exports (if `--export-figures` used) |
| `report.pdf` | PDF version (if `--pdf` used) |

---

## Config file

All CLI parameters can be specified in a YAML config file:

```yaml
# batch_detective_config.yaml
counts: counts.csv
metadata: metadata.csv
output_dir: ./results/

n_variable_genes: 2000
n_pcs: 10
min_cpm: 1.0
outlier_pval: 0.001

technical_covariates: [batch, plate, operator]
biological_covariates: [treatment, sex, age]
```

```bash
batch-detective run --config batch_detective_config.yaml
```

Paths in the config file are resolved relative to the config file's location. CLI flags override config file values. Config file values override defaults.

---

## Pipeline integration

`batch-detective` is designed for use in automated pipelines. Exit code `1` signals detected batch effects — **this is a normal result, not an error.**

```bash
# Correct pipeline pattern:
batch-detective run \
  --counts counts.csv \
  --metadata metadata.csv \
  --output-dir ./results/ \
  --technical-covariates batch \
  --quiet

EXIT=$?

if [ $EXIT -eq 0 ]; then
    echo "No batch effects. Proceeding with standard DE."
elif [ $EXIT -eq 1 ]; then
    echo "Batch effects detected. Adding batch to DE model."
elif [ $EXIT -eq 3 ]; then
    echo "Validation failed. Check inputs." >&2
    exit 1
fi
```

### Nextflow / Snakemake

```groovy
// Nextflow example
process BATCH_DETECTIVE {
    input:
    path counts
    path metadata

    output:
    path "results/report.html"
    path "results/run_manifest.json"

    script:
    """
    batch-detective run \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --output-dir results/ \\
        --technical-covariates batch \\
        --biological-covariates treatment || true
    """
}
```

```python
# Snakemake example
rule batch_detective:
    input:
        counts="counts.csv",
        metadata="metadata.csv"
    output:
        report="results/report.html",
        manifest="results/run_manifest.json"
    shell:
        """
        batch-detective run \
            --counts {input.counts} \
            --metadata {input.metadata} \
            --output-dir results/ \
            --technical-covariates batch \
            --biological-covariates treatment || true
        """
```

> The `|| true` prevents the pipeline from aborting on exit code 1 (batch detected). Check the manifest for the actual exit code.

---

## HPC / cluster usage

`batch-detective` runs entirely on CPU. Memory usage is typically < 2 GB for datasets up to 20,000 genes × 500 samples.

```bash
# Batch submission (SLURM example)
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

batch-detective run \
  --counts counts.csv \
  --metadata metadata.csv \
  --output-dir ./results/ \
  --log-file ./results/run.log
```

**Web UI on HPC (SSH port forwarding):**

```bash
# On HPC
batch-detective serve --port 8501

# On your local machine (in a separate terminal)
ssh -L 8501:localhost:8501 user@hpc.server.com

# Then open http://localhost:8501 in your local browser
```

---

## Troubleshooting

### `ERROR: Count matrix contains aligner summary rows`

Your count matrix contains STAR or HTSeq summary rows (`N_unmapped`, `N_multimapping`, `__no_feature`, etc.). Remove them before running:

```python
import pandas as pd
df = pd.read_csv("counts.csv", index_col=0)
df = df[~df.index.str.startswith(("N_", "__"))]
df.to_csv("counts_clean.csv")
```

### `ERROR: Count matrix contains non-integer values`

If your data is already normalized (RPKM, FPKM, TPM, log-transformed), use `--normalized`. If it contains floats due to rounding, use `--force`.

### `ERROR: Output directory already contains files`

Use `--overwrite` to replace existing outputs:

```bash
batch-detective run ... --overwrite
```

### weasyprint PDF dependencies missing

```bash
# Linux
sudo apt install libpango-1.0-0 libcairo2 libgdk-pixbuf2.0-0

# macOS
brew install pango cairo gdk-pixbuf

# Windows
# See: https://doc.courtbouillon.org/weasyprint/stable/first_steps.html
```

### Port conflict with `serve`

batch-detective automatically finds the next available port if 8501 is in use. The port used is printed on startup.

### Sample IDs do not match between files

Use `--dry-run` to inspect detected sample IDs before running the full analysis. Check for whitespace, case differences, or trailing characters in your CSV headers.

### Large datasets (> 50,000 genes)

Gene filtering (CPM threshold) reduces the working set before PCA. For very large matrices, increase `--min-cpm` or decrease `--n-variable-genes` to reduce memory usage.

---

## Known limitations

1. **Repeated measures / longitudinal designs** — Not supported. The tool detects likely subject-ID columns and warns, but statistical tests (Kruskal-Wallis, ICC) assume independent samples. Use DESeq2 LRT, lme4, or dream for repeated measures.

2. **ICC underestimation under log1p(CPM)** — log1p-CPM normalization is not variance-stabilized (cf. DESeq2 VST). High-expression genes may dominate variance. ICC values derived from log1p-CPM may underestimate true batch effect magnitude by approximately 20–35%.

3. **Non-monotonic continuous relationships** — Spearman correlation detects monotonic associations only. U-shaped or non-monotonic relationships between continuous covariates and PC scores are not detected.

4. **Single-cell, spatial, and ATAC-seq** — Not supported. The tool is designed for bulk RNA-seq count matrices.

5. **Collinearity suppression** — When batch and treatment are strongly collinear (Cramér's V ≥ 0.70), the tool reports ICC = 0 and suppresses PC association tests for the affected covariate. This is the correct behavior — associations cannot be interpreted when covariates are confounded.

6. **Low statistical power at small n** — With < 15 samples per group, the tool may miss moderate batch effects (ICC < 0.40). Power warnings are displayed in the report when this applies.

7. **Fixed vs. random batch effects** — ICC(1,1) treats batch labels as randomly sampled from a population. If your batches are fixed (i.e., you care specifically about these batches and would not replicate the study with different ones), ICC values may underestimate true effect magnitude.

---

## Citation

If you use batch-detective in published research, please cite:

> **batch-detective v0.1.0** — Statistical diagnosis of batch effects in bulk RNA-seq count matrices.

And the underlying statistical methods:

- **ICC interpretation:** Koo, T.K. & Li, M.Y. (2016). A guideline of selecting and reporting intraclass correlation coefficients for reliability research. *Journal of Chiropractic Medicine*, 15(2), 155–163.
- **FDR correction:** Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.
- **PCA implementation:** Pedregosa, F. et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.
- **Outlier detection (LedoitWolf):** Ledoit, O. & Wolf, M. (2004). A well-conditioned estimator for large-dimensional covariance matrices. *Journal of Multivariate Analysis*, 88(2), 365–411.

---

## Contributing

Pull requests are welcome. Please:

1. Run `pytest` before submitting
2. Run `flake8 batch_detective/` to check style
3. Add tests for any new functionality
4. Update this README if CLI options change

---

## License

MIT — see [LICENSE](LICENSE) for details.
