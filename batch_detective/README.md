# batch-detective

A statistically rigorous, zero-configuration, read-only diagnostic CLI tool for diagnosing technical batch effects in bulk RNA-seq count matrices.

## What it does

`batch-detective` accepts a raw count matrix and a sample metadata table, then produces a self-contained HTML report with statistical diagnosis of potential batch effects. It quantifies how much transcriptomic variance is attributable to technical variables (batch, plate, operator, run date) vs biological variables (treatment, disease, sex).

Key outputs:
- **Traffic light summary** (🔴/🟡/🟢) for quick interpretation
- **ICC (Intraclass Correlation Coefficient)** — the headline metric for batch effect magnitude
- **PC-covariate association heatmap** — which PCs are driven by which covariates
- **Collinearity detection** — warns when batch and biology are confounded
- **PCA scatter plots** colored by each covariate

## Installation

### Base (CLI)
```bash
pip install batch-detective
```

### With local web UI [serve]
```bash
pip install "batch-detective[serve]"
```

### With PDF export [pdf]
```bash
pip install "batch-detective[pdf]"
# Linux: sudo apt install libpango-1.0-0 libcairo2 libgdk-pixbuf2.0-0
# macOS: brew install pango cairo gdk-pixbuf
```

### Developer install [dev]
```bash
git clone https://github.com/example/batch-detective
cd batch-detective
pip install -e ".[dev]"
```

## Quick start

```bash
# Run analysis
batch-detective run \
  --counts counts.csv \
  --metadata metadata.csv \
  --output-dir ./results/ \
  --technical-covariates batch,plate \
  --biological-covariates treatment,sex

# Open report
open results/report.html

# Local web UI
batch-detective serve
```

## CLI Reference

### batch-detective run

```
Options:
  --counts PATH                 Counts CSV/TSV file (required)
  --metadata PATH               Metadata CSV/TSV file (required)
  --output-dir PATH             Output directory (required)
  --normalized                  Input is pre-normalized (skip CPM)
  --n-variable-genes INT        Number of variable genes [2000]
  --n-pcs INT                   Number of PCs [10]
  --min-cpm FLOAT               Min CPM threshold [1.0]
  --min-samples-expressing INT  Min samples expressing a gene (auto)
  --outlier-pval FLOAT          Outlier p-value threshold [0.001]
  --primary-covariate NAME      Single covariate to condition on
  --condition-on "c1,c2"        Multiple covariates to regress out
  --technical-covariates "b,p"  Label covariates as technical
  --biological-covariates "t,s" Label covariates as biological
  --config PATH                 YAML config file
  --dry-run                     Preview only, no analysis
  --overwrite                   Overwrite existing output directory
  --export-figures              Export PNG+SVG figures
  --anonymize-samples           Replace sample IDs with Sample_001, etc.
  --pdf                         Export PDF report (requires [pdf])
  --force                       Proceed past non-fatal warnings
  --verbose                     Debug logging
  --quiet                       Suppress non-error output
  --log-file PATH               Write log to file
```

### batch-detective serve

```bash
batch-detective serve [--port INT]   # default: 8501
```

### Exit codes

| Code | Meaning |
|------|---------|
| 0 | Analysis complete — no significant batch effects |
| 1 | Analysis complete — significant batch effects detected (**NOT an error**) |
| 2 | Analysis complete with data quality warnings |
| 3 | Fatal error — validation failed or unrecoverable error |

> **Note:** Exit code 1 means the tool ran successfully and found significant batch effects. Use in pipelines:
> ```bash
> batch-detective run ... && echo "clean" || echo "batch_detected"
> ```

## Interpreting your report

### Understanding ICC

ICC (Intraclass Correlation Coefficient) quantifies what fraction of gene expression variance is attributable to a grouping variable. Values are computed using ICC(1,1) — a one-way random effects model.

| ICC range | Interpretation |
|-----------|----------------|
| < 0.10 | Negligible |
| 0.10 – 0.30 | Mild |
| 0.30 – 0.60 | Moderate |
| > 0.60 | Strong |

For **biological variables** (treatment, disease), high ICC is expected — it means your experiment worked. Only unexpected high ICC for **technical variables** (batch, plate, operator, date) indicates a batch effect.

### Understanding the association heatmap

Color intensity shows effect size (η² for categorical, |ρ| for continuous covariates). Asterisks indicate FDR-corrected significance (* q<0.05, ** q<0.01, *** q<0.001).

### What to do if batch effects are found

If batch and treatment are **not** collinear:
- **DESeq2:** `design = ~ batch + treatment`
- **edgeR:** `design <- model.matrix(~ batch + treatment)`
- **Visualization only:** ComBat from the `sva` R package

If batch and treatment **are** collinear (confounded design): do not apply batch correction — it will remove biological signal.

## Config file

```yaml
# batch_detective_config.yaml
counts: counts.csv
metadata: metadata.csv
output_dir: ./results/
n_variable_genes: 2000
n_pcs: 10
min_cpm: 1.0
technical_covariates: [batch, plate]
biological_covariates: [treatment, sex]
```

```bash
batch-detective run --config batch_detective_config.yaml
```

Paths in the config file are resolved relative to the config file's directory.

## HPC / cluster usage

```bash
# On HPC, run analysis
batch-detective run --counts ... --output-dir ./results/

# To use the web UI on a cluster, use SSH port forwarding:
# On your local machine:
ssh -L 8501:localhost:8501 user@hpc.server.com
# On the HPC:
batch-detective serve
# Then open http://localhost:8501 on your local browser
```

## Troubleshooting

### weasyprint system dependencies
```bash
# Linux
sudo apt install libpango-1.0-0 libcairo2 libgdk-pixbuf2.0-0

# macOS
brew install pango cairo gdk-pixbuf

# Windows
# See: https://doc.courtbouillon.org/weasyprint/stable/first_steps.html
```

### Port conflicts with --serve
If port 8501 is in use, batch-detective will automatically find the next available port.

### Non-integer count matrix
Use `--normalized` if your data is already normalized (RPKM, FPKM, TPM, log-transformed).
Use `--force` to proceed despite non-integer values with a warning.

### featureCounts output
featureCounts adds annotation columns (Chr, Start, End, Strand, Length). batch-detective
detects and removes these automatically.

### STAR/HTSeq summary rows
Remove rows like `N_unmapped`, `__no_feature` before running batch-detective.

## Known Limitations

1. **Repeated measures / longitudinal designs:** Not supported. Tool warns when a likely subject ID column is detected but cannot correctly analyze paired data. Use DESeq2 LRT, lme4, or dream.

2. **Multi-platform merged datasets:** Library size outlier detection may flag entire platforms. Interpret QC flags in context of platform metadata.

3. **Non-bulk RNA-seq:** Single-cell, spatial, and ATAC-seq data are not supported.

4. **FFPE vs fresh-frozen:** Tissue preservation method differences create technical variation not distinguishable from batch effects.

5. **Variance stabilization:** log1p(CPM) is used for speed. Results may differ from DESeq2 VST analyses.

6. **Non-monotonic continuous covariate relationships:** Spearman correlation will miss U-shaped associations.

7. **Batch mode for multiple projects:** Not supported in v1. Run tool separately per project.

## Contributing

Pull requests welcome. Please run `pytest` and `flake8` before submitting.

## License

MIT

## Citation

If you use batch-detective in published research, please cite:

> batch-detective v0.1.0. Statistical diagnosis of batch effects in bulk RNA-seq.
> References: Koo & Mae (2016) for ICC interpretation;
> Benjamini & Hochberg (1995) for FDR correction;
> Pedregosa et al. (2011) for scikit-learn PCA.
