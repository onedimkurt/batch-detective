"""CLI entry point for batch-detective."""

import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import click
from .compat import tqdm

from . import __version__
from .exceptions import BatchDetectiveError, BatchDetectiveDependencyError
from .config import merge_config

logger = logging.getLogger(__name__)

VERSION_STRING = (
    f"batch-detective {__version__} "
    f"(Python {sys.version.split()[0]}, {platform.system()})"
)


def setup_logging(verbose: bool, quiet: bool, log_file: Optional[Path]) -> None:
    """Configure logging based on flags."""
    level = (
        logging.DEBUG if verbose
        else logging.WARNING if quiet
        else logging.INFO
    )
    handlers = [logging.StreamHandler()]
    if log_file:
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        handlers.append(fh)
    logging.basicConfig(
        level=level,
        handlers=handlers,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )


@click.group()
@click.version_option(version=__version__, message=VERSION_STRING)
def main():
    """batch-detective — Batch effect diagnosis for bulk RNA-seq count matrices.

    Exit codes:
      0 — Analysis complete, no significant batch effects detected
      1 — Analysis complete, significant batch effects detected
      2 — Analysis complete with data quality warnings
      3 — Fatal error: validation failed or unrecoverable error

    Note: Exit code 1 means the tool ran successfully and found significant
    batch effects. This is NOT an error. Use in pipelines:
      batch-detective run ... && echo clean || echo batch_detected
    """
    pass


@main.command()
@click.option("--counts", required=False, type=click.Path(), help="Counts CSV/TSV file")
@click.option("--metadata", required=False, type=click.Path(), help="Metadata CSV/TSV file")
@click.option("--output-dir", required=False, type=click.Path(), help="Output directory")
@click.option("--normalized", is_flag=True, default=False, help="Input is pre-normalized")
@click.option("--n-variable-genes", type=int, default=None, help="Number of variable genes [2000]")
@click.option("--n-pcs", type=int, default=None, help="Number of PCs [10]")
@click.option("--min-cpm", type=float, default=None, help="Min CPM threshold [1.0]")
@click.option("--min-samples-expressing", type=int, default=None, help="Min samples expressing a gene")
@click.option("--outlier-pval", type=float, default=None, help="Outlier p-value [0.001]")
@click.option("--primary-covariate", type=str, default=None, help="Single covariate to condition on")
@click.option("--condition-on", type=str, default=None, help="Covariates to regress out (comma-sep)")
@click.option("--technical-covariates", type=str, default=None, help="Technical covariate names (comma-sep)")
@click.option("--biological-covariates", type=str, default=None, help="Biological covariate names (comma-sep)")
@click.option("--config", "config_path", type=click.Path(), default=None, help="YAML config file")
@click.option("--dry-run", is_flag=True, default=False, help="Validate and preview, no analysis")
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory")
@click.option("--log-file", type=click.Path(), default=None, help="Write log to file")
@click.option("--verbose", is_flag=True, default=False, help="Debug-level logging")
@click.option("--quiet", is_flag=True, default=False, help="Suppress non-error output")
@click.option("--pdf", is_flag=True, default=False, help="Export PDF report (requires [pdf])")
@click.option("--dpi", type=int, default=None, help="Figure DPI [150 HTML, 300 PDF]")
@click.option("--export-figures", is_flag=True, default=False, help="Export figures as PNG+SVG")
@click.option("--anonymize-samples", is_flag=True, default=False, help="Replace sample IDs")
@click.option("--force", is_flag=True, default=False, help="Proceed past non-fatal warnings")
def run(
    counts, metadata, output_dir,
    normalized, n_variable_genes, n_pcs, min_cpm,
    min_samples_expressing, outlier_pval, primary_covariate,
    condition_on, technical_covariates, biological_covariates,
    config_path, dry_run, overwrite, log_file, verbose, quiet,
    pdf, dpi, export_figures, anonymize_samples, force,
):
    """Run batch effect diagnosis.

    Accepts a raw count matrix and sample metadata table and produces a
    self-contained HTML report diagnosing potential batch effects.
    """
    setup_logging(verbose, quiet, Path(log_file) if log_file else None)

    from .dependencies import check_dependencies
    try:
        check_dependencies()
    except BatchDetectiveDependencyError as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(3)

    # Parse list params
    def _split(s):
        return [x.strip() for x in s.split(",") if x.strip()] if s else []

    tech_covs = _split(technical_covariates)
    bio_covs = _split(biological_covariates)
    condition_on_list = _split(condition_on)
    if primary_covariate and primary_covariate not in condition_on_list:
        condition_on_list = [primary_covariate] + condition_on_list

    cli_params = {
        k: v for k, v in {
            "counts": counts,
            "metadata": metadata,
            "output_dir": output_dir,
            "normalized": normalized or None,
            "n_variable_genes": n_variable_genes,
            "n_pcs": n_pcs,
            "min_cpm": min_cpm,
            "min_samples_expressing": min_samples_expressing,
            "outlier_pval": outlier_pval,
            "primary_covariate": primary_covariate,
            "technical_covariates": tech_covs or None,
            "biological_covariates": bio_covs or None,
            "dry_run": dry_run or None,
            "overwrite": overwrite or None,
            "verbose": verbose or None,
            "quiet": quiet or None,
            "pdf": pdf or None,
            "dpi": dpi,
            "export_figures": export_figures or None,
            "anonymize_samples": anonymize_samples or None,
            "force": force or None,
        }.items() if v is not None
    }

    config_file_path = Path(config_path) if config_path else None

    try:
        config = merge_config(cli_params, config_file_path)
    except BatchDetectiveError as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(3)

    # Store parsed list params
    config["condition_on_list"] = condition_on_list
    config["technical_covariates"] = tech_covs or config.get("technical_covariates", [])
    config["biological_covariates"] = bio_covs or config.get("biological_covariates", [])

    # Validate required params
    if not config.get("counts"):
        click.echo("ERROR: --counts is required", err=True)
        sys.exit(3)
    if not config.get("metadata"):
        click.echo("ERROR: --metadata is required", err=True)
        sys.exit(3)
    if not config.get("output_dir"):
        click.echo("ERROR: --output-dir is required", err=True)
        sys.exit(3)

    counts_path = Path(config["counts"])
    metadata_path = Path(config["metadata"])
    output_path = Path(config["output_dir"])

    dpi_val = config.get("dpi") or (300 if pdf else 150)

    # Output directory management
    if not dry_run:
        if output_path.exists():
            existing_files = list(output_path.iterdir())
            if existing_files and not config.get("overwrite"):
                click.echo(
                    f"ERROR: Output directory {output_path} already contains files.\n"
                    "Use --overwrite to replace existing outputs.",
                    err=True,
                )
                sys.exit(3)
        else:
            output_path.mkdir(parents=True, exist_ok=False)

    # Always write manifest (even on failure)
    manifest_data = {}

    try:
        exit_code = _run_analysis(
            counts_path=counts_path,
            metadata_path=metadata_path,
            output_path=output_path,
            config=config,
            dpi=dpi_val,
            dry_run=dry_run,
            condition_on_list=condition_on_list,
        )
    except BatchDetectiveError as e:
        if not verbose:
            click.echo(f"ERROR: {e}", err=True)
        else:
            import traceback
            click.echo(traceback.format_exc(), err=True)

        # Write failure manifest
        if not dry_run and output_path.exists():
            try:
                from .manifest import write_manifest
                write_manifest(
                    output_path, counts_path, metadata_path,
                    config, {"status": "failed", "error": str(e)}, 3,
                    anonymized=config.get("anonymize_samples", False),
                )
            except Exception:
                pass
        sys.exit(3)

    sys.exit(exit_code)


def _run_analysis(
    counts_path: Path,
    metadata_path: Path,
    output_path: Path,
    config: dict,
    dpi: int,
    dry_run: bool,
    condition_on_list: List[str],
) -> int:
    """Execute the full analysis pipeline.

    Args:
        counts_path: Path to counts file.
        metadata_path: Path to metadata file.
        output_path: Output directory.
        config: Resolved config dict.
        dpi: Figure resolution.
        dry_run: If True, only validate.
        condition_on_list: Conditioning covariate names.

    Returns:
        Exit code (0, 1, or 2).
    """
    import numpy as np
    from .validator import validate_inputs
    from .qc import QualityController
    from .normalizer import normalize_counts
    from .pca_analysis import run_pca
    from .association import detect_collinearity, run_pc_associations, compute_all_icc
    from .outliers import detect_outliers
    from .gene_impact import get_top_batch_genes
    from .plots import generate_all_plots
    from .report import assemble_report
    from .manifest import write_manifest
    from .anonymizer import anonymize_samples

    rng = np.random.default_rng(42)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    disable_tqdm = not sys.stdout.isatty()

    steps = [
        "[1/7] Validating inputs...",
        "[2/7] Running quality control...",
        "[3/7] Normalizing and selecting variable genes...",
        "[4/7] Running PCA...",
        "[5/7] Testing PC-covariate associations...",
        "[6/7] Detecting outliers...",
        "[7/7] Generating report...",
    ]

    # STAGE 1: Validate
    click.echo(steps[0])
    counts_raw, metadata_raw, counts_working, metadata_working = validate_inputs(
        counts_path=counts_path,
        metadata_path=metadata_path,
        force=config.get("force", False),
        normalized=config.get("normalized", False),
    )

    if dry_run:
        _print_dry_run(
            counts_raw, metadata_working, config, condition_on_list, output_path
        )
        return 0

    # Apply anonymization
    if config.get("anonymize_samples"):
        counts_working, metadata_working, _ = anonymize_samples(counts_working, metadata_working)

    n_genes_input = counts_working.shape[0]

    # STAGE 2: QC
    click.echo(steps[1])
    qc = QualityController(
        counts_working=counts_working,
        metadata_working=metadata_working,
        min_cpm=config.get("min_cpm", 1.0),
        min_samples_expressing=config.get("min_samples_expressing"),
        technical_covariates=config.get("technical_covariates", []),
        biological_covariates=config.get("biological_covariates", []),
    )

    covariate_info = qc.run_covariate_prescreening()
    qc_summary = qc.run_sample_qc()
    qc.run_power_assessment()

    primary_covariate = config.get("primary_covariate")
    counts_filtered, filter_stats = qc.filter_genes(primary_covariate)

    # STAGE 3: Normalize
    click.echo(steps[2])
    log1p_cpm_all, log1p_cpm_selected, gene_stats = normalize_counts(
        counts_filtered,
        n_variable_genes=config.get("n_variable_genes", 2000),
        normalized=config.get("normalized", False),
    )

    n_genes_analyzed = log1p_cpm_selected.shape[0]

    # STAGE 4: PCA
    click.echo(steps[3])
    pc_scores, explained_variance_ratio, components, sample_ids = run_pca(
        log1p_cpm_selected,
        n_pcs=config.get("n_pcs", 10),
    )

    # STAGE 5: Collinearity
    click.echo(steps[4])
    collinearity_warnings = detect_collinearity(metadata_working, covariate_info)

    # STAGE 6: Associations
    assoc_df = run_pc_associations(
        pc_scores=pc_scores,
        explained_variance_ratio=explained_variance_ratio,
        metadata=metadata_working,
        covariate_info=covariate_info,
        n_pcs_test=min(5, pc_scores.shape[1]),
        condition_on=condition_on_list if condition_on_list else None,
        rng=rng,
    )

    # STAGE 7: ICC
    icc_table = compute_all_icc(
        log1p_cpm_selected,
        metadata_working,
        covariate_info,
        n_bootstrap=1000,
        rng=rng,
    )

    # STAGE 8: Outliers
    click.echo(steps[5])
    outlier_df, detection_info = detect_outliers(
        pc_scores=pc_scores,
        sample_ids=sample_ids,
        metadata=metadata_working,
        n_pcs=config.get("n_pcs", 10),
        outlier_pval=config.get("outlier_pval", 0.001),
    )

    # Update qc_summary with Mahalanobis
    mahal_distances = detection_info.get("mahal_distances", [])
    mahal_outliers = detection_info.get("mahal_outliers", [])
    if len(mahal_distances) == len(sample_ids):
        qc_summary["mahal_distance"] = mahal_distances
        qc_summary["mahal_outlier_flag"] = mahal_outliers

    # Top batch-associated genes
    top_genes = get_top_batch_genes(
        log1p_cpm_selected, metadata_working, covariate_info, icc_table
    ) if not icc_table.empty else {}

    # Figures dir
    figures_dir = None
    if config.get("export_figures"):
        figures_dir = output_path / "figures"
        figures_dir.mkdir(exist_ok=True)

    # STAGE 9: Plots
    outlier_sample_ids = outlier_df["sample_id"].tolist() if not outlier_df.empty else []

    click.echo(steps[6])
    plots = generate_all_plots(
        qc_summary=qc_summary,
        pc_scores=pc_scores,
        explained_variance_ratio=explained_variance_ratio,
        metadata=metadata_working,
        covariate_info=covariate_info,
        icc_table=icc_table,
        assoc_df=assoc_df,
        outlier_samples=outlier_sample_ids,
        sample_ids=sample_ids,
        top_genes=top_genes,
        dpi=dpi,
        figures_dir=figures_dir,
        low_power=getattr(qc, "low_power", False),
    )

    # Assemble report
    analysis_meta = {
        "n_samples": len(sample_ids),
        "n_genes_input": n_genes_input,
        "n_genes_analyzed": n_genes_analyzed,
        "normalization_method": "Pre-normalized (user-supplied)" if config.get("normalized") else "CPM + log1p",
        "low_power": getattr(qc, "low_power", False),
        "min_detectable_icc": getattr(qc, "min_detectable_icc", 0.0),
    }

    assemble_report(
        output_path=output_path / "report.html",
        qc_summary=qc_summary,
        icc_table=icc_table,
        assoc_df=assoc_df,
        outlier_df=outlier_df,
        plots=plots,
        top_genes=top_genes,
        covariate_info=covariate_info,
        collinearity_warnings=collinearity_warnings,
        repeated_measures_covariates=getattr(qc, "repeated_measures_covariates", []),
        config=config,
        analysis_meta=analysis_meta,
        detection_info=detection_info,
        timestamp=timestamp,
    )

    # Write CSV outputs
    qc_summary.to_csv(output_path / "qc_summary.csv", index=False)

    if not assoc_df.empty:
        assoc_df.to_csv(output_path / "association_table.csv", index=False)
    else:
        import pandas as pd
        pd.DataFrame().to_csv(output_path / "association_table.csv", index=False)

    if not icc_table.empty:
        icc_table.drop(columns=["_icc_values"], errors="ignore").to_csv(
            output_path / "icc_table.csv", index=False
        )
    else:
        import pandas as pd
        pd.DataFrame().to_csv(output_path / "icc_table.csv", index=False)

    # PDF export
    if config.get("pdf"):
        try:
            import weasyprint
            pdf_path = output_path / "report.pdf"
            with open(output_path / "report.html", encoding="utf-8") as f:
                html_content = f.read()
            weasyprint.HTML(string=html_content).write_pdf(str(pdf_path))
            click.echo(f"PDF report written to {pdf_path}")
        except ImportError:
            from .exceptions import BatchDetectiveDependencyError
            raise BatchDetectiveDependencyError(
                "PDF export requires additional dependencies.\n"
                "Install: pip install 'batch-detective[pdf]'\n"
                "Linux:   sudo apt install libpango-1.0-0 libcairo2 libgdk-pixbuf2.0-0\n"
                "macOS:   brew install pango cairo gdk-pixbuf\n"
                "Windows: https://doc.courtbouillon.org/weasyprint/stable/first_steps.html"
            )

    # Determine exit code
    tech_covs = config.get("technical_covariates", [])
    exit_code = _determine_exit_code(icc_table, assoc_df, tech_covs, getattr(qc, "low_power", False))

    # Write manifest
    sig_assoc = []
    if not assoc_df.empty and "significant_q05" in assoc_df.columns:
        sig_rows = assoc_df[assoc_df["significant_q05"] == True]
        sig_assoc = [f"{r.covariate}:PC{r.pc}" for _, r in sig_rows.iterrows()]

    max_icc_info = {}
    if not icc_table.empty:
        top_row = icc_table.loc[icc_table["median_icc"].idxmax()]
        max_icc_info = {
            "covariate": top_row["covariate"],
            "median_icc": float(top_row["median_icc"]),
        }

    analysis_summary = {
        "n_samples": len(sample_ids),
        "n_genes_input": n_genes_input,
        "n_genes_analyzed": n_genes_analyzed,
        "covariates_tested": [n for n, i in covariate_info.items() if not i.get("skip")],
        "covariates_skipped": [n for n, i in covariate_info.items() if i.get("skip")],
        "significant_associations_q05": sig_assoc,
        "max_icc": max_icc_info,
        "n_outliers_detected": len(outlier_df) if not outlier_df.empty else 0,
    }

    resolved_params = {
        "n_variable_genes": config.get("n_variable_genes", 2000),
        "n_pcs": config.get("n_pcs", 10),
        "min_cpm": config.get("min_cpm", 1.0),
        "normalized": config.get("normalized", False),
        "technical_covariates": tech_covs,
        "biological_covariates": config.get("biological_covariates", []),
    }

    write_manifest(
        output_dir=output_path,
        counts_path=counts_path,
        metadata_path=metadata_path,
        resolved_params=resolved_params,
        analysis_summary=analysis_summary,
        exit_code=exit_code,
        anonymized=config.get("anonymize_samples", False),
    )

    click.echo(f"\n✅ Analysis complete. Report: {output_path / 'report.html'}")
    return exit_code


def _determine_exit_code(
    icc_table,
    assoc_df,
    technical_covariates,
    low_power: bool,
) -> int:
    """Determine exit code from results."""
    import pandas as pd

    if icc_table.empty and assoc_df.empty:
        return 0

    tech_icc = icc_table[icc_table["covariate"].isin(technical_covariates)] if not icc_table.empty else icc_table

    max_tech_icc = tech_icc["median_icc"].max() if not tech_icc.empty else 0.0
    any_sig = (
        assoc_df[assoc_df["covariate"].isin(technical_covariates)]["significant_q05"].any()
        if not assoc_df.empty and "significant_q05" in assoc_df.columns and technical_covariates
        else False
    )

    if max_tech_icc >= 0.30 or any_sig:
        return 1
    return 0


def _print_dry_run(counts_raw, metadata_working, config, condition_on_list, output_path):
    """Print dry-run output."""
    from .qc import QualityController
    n_genes, n_samples = counts_raw.shape

    click.echo("batch-detective dry run")
    click.echo("=======================")
    click.echo("Input files")
    click.echo(f"  Counts:   {config.get('counts')} ({n_genes:,} genes × {n_samples} samples, CSV)")
    click.echo(f"  Metadata: {config.get('metadata')} ({len(metadata_working)} samples × {len(metadata_working.columns)} covariates)")
    click.echo("")
    click.echo("Sample overlap")
    click.echo(f"  Counts samples:   {n_samples}")
    click.echo(f"  Metadata samples: {len(metadata_working)}")
    click.echo(f"  Intersection:     {n_samples} (100%)")
    click.echo("")
    click.echo("Detected covariates")

    import pandas as pd
    import numpy as np

    for col in metadata_working.columns:
        series = metadata_working[col]
        n_unique = series.nunique(dropna=True)
        numeric_series = pd.to_numeric(series, errors="coerce")
        is_numeric = not numeric_series.isna().any()

        if n_unique < 10 or not is_numeric:
            groups = sorted(series.dropna().unique(), key=str)
            sizes = [int((series == g).sum()) for g in groups]
            click.echo(
                f"  {col:<16} categorical  {n_unique} groups  "
                f"[sizes: {', '.join(str(s) for s in sizes)}]"
            )
        else:
            lo = numeric_series.min()
            hi = numeric_series.max()
            click.echo(
                f"  {col:<16} continuous   range: {lo:.1f} – {hi:.1f}"
            )

    click.echo("")
    click.echo("Resolved parameters")

    qc = QualityController(
        counts_working=counts_raw,
        metadata_working=metadata_working,
        min_cpm=config.get("min_cpm", 1.0),
    )
    covariate_info = qc.run_covariate_prescreening()
    _, filter_stats = qc.filter_genes(config.get("primary_covariate"))

    min_samples = getattr(qc, "min_samples_expressing_resolved", "auto")

    click.echo(f"  n_variable_genes:       {config.get('n_variable_genes', 2000)}")
    click.echo(f"  n_pcs:                  {config.get('n_pcs', 10)}")
    click.echo(f"  min_cpm:                {config.get('min_cpm', 1.0)}")
    click.echo(f"  min_samples_expressing: {min_samples} (auto)")
    click.echo(f"  outlier_pval:           {config.get('outlier_pval', 0.001)}")
    click.echo("")
    click.echo("Estimated runtime: < 1 minute")
    click.echo("")
    click.echo(f"Output files that WOULD be created in {output_path}")
    click.echo("  report.html")
    click.echo("  qc_summary.csv")
    click.echo("  association_table.csv")
    click.echo("  icc_table.csv")
    click.echo("  run_manifest.json")
    click.echo("")
    click.echo("Run without --dry-run to execute analysis.")


@main.command()
@click.option("--port", type=int, default=8501, help="Port for Streamlit UI [8501]")
def serve(port):
    """Launch local Streamlit web UI.

    Requires: pip install 'batch-detective[serve]'
    """
    try:
        import streamlit
    except ImportError:
        raise BatchDetectiveDependencyError(
            "Streamlit web UI requires additional dependencies.\n"
            "Install: pip install 'batch-detective[serve]'"
        )

    from .utils import find_available_port
    available_port = find_available_port(port)
    if available_port != port:
        click.echo(f"Port {port} in use. Starting on port {available_port}.")

    import subprocess
    app_path = Path(__file__).parent / "streamlit_app.py"
    subprocess.run(
        [sys.executable, "-m", "streamlit", "run", str(app_path), "--server.port", str(available_port)],
        check=True,
    )


@main.command()
def validate():
    """Run validation on built-in synthetic data.

    Confirms the tool works correctly by running on a known dataset
    with expected ICC ~0.45 for the batch covariate.
    """
    # Note: don't call setup_logging here — can cause I/O issues in test runners

    click.echo("Running batch-detective validation on built-in synthetic data...")

    import tempfile
    import numpy as np
    import pandas as pd

    from .dependencies import check_dependencies
    from .validator import validate_inputs
    from .qc import QualityController
    from .normalizer import normalize_counts
    from .pca_analysis import run_pca
    from .association import compute_all_icc

    try:
        check_dependencies()
    except BatchDetectiveDependencyError as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(3)

    # Generate synthetic data
    rng = np.random.default_rng(42)
    n_samples = 24
    n_genes = 1200

    metadata = pd.DataFrame({
        "batch": [f"B{i+1}" for i in range(4) for _ in range(6)],
        "treatment": ["ctrl"] * 12 + ["treat"] * 12,
    }, index=[f"S{i+1:02d}" for i in range(n_samples)])

    counts = rng.negative_binomial(50, 0.5, size=(n_genes, n_samples)).astype(float)

    # Inject batch signal in first 400 genes
    batch_effect = np.array([0, 2, 4, 6])
    for b_idx, batch_name in enumerate(["B1", "B2", "B3", "B4"]):
        mask = metadata["batch"] == batch_name
        sample_idx = np.where(mask)[0]
        counts[:400, sample_idx] += rng.poisson(batch_effect[b_idx] * 20, size=(400, len(sample_idx)))

    counts = pd.DataFrame(
        counts.astype(int),
        index=[f"gene_{i:04d}" for i in range(n_genes)],
        columns=metadata.index,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        counts_path = tmpdir / "counts.csv"
        meta_path = tmpdir / "metadata.csv"
        counts.to_csv(counts_path)
        metadata.to_csv(meta_path)

        counts_raw, metadata_raw, counts_working, metadata_working = validate_inputs(
            counts_path, meta_path
        )
        qc = QualityController(counts_working, metadata_working)
        covariate_info = qc.run_covariate_prescreening()
        qc.run_sample_qc()
        counts_filtered, _ = qc.filter_genes()
        _, log1p_selected, _ = normalize_counts(counts_filtered, n_variable_genes=500)
        icc_table = compute_all_icc(log1p_selected, metadata_working, covariate_info, n_bootstrap=100)

    if icc_table.empty or "batch" not in icc_table["covariate"].values:
        click.echo("ERROR: Validation failed — no ICC computed for batch covariate.", err=True)
        sys.exit(3)

    batch_icc = float(icc_table.loc[icc_table["covariate"] == "batch", "median_icc"].iloc[0])

    if 0.20 <= batch_icc <= 0.80:
        click.echo(f"✅ Validation passed: batch ICC = {batch_icc:.2f} (expected ~0.20–0.80)")
        sys.exit(0)
    else:
        click.echo(
            f"❌ Validation FAILED: batch ICC = {batch_icc:.2f} "
            "(expected 0.20–0.80). Tool may be broken.",
            err=True,
        )
        sys.exit(3)
