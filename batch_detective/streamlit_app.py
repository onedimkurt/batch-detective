"""Local Streamlit web UI for batch-detective.

Only imported when 'batch-detective serve' is called.
Requires: pip install 'batch-detective[serve]'
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import streamlit as st


def main():
    """Run the Streamlit app."""
    import streamlit as st
    import tempfile
    import pandas as pd
    from pathlib import Path

    st.set_page_config(
        page_title="batch-detective",
        page_icon="🔬",
        layout="wide",
    )

    st.title("🔬 batch-detective")
    st.info("🔒 All analysis runs locally on your machine. No data leaves your computer.")

    col1, col2 = st.columns(2)
    with col1:
        counts_file = st.file_uploader("Upload counts matrix (CSV/TSV)", type=["csv", "tsv", "txt"])
    with col2:
        meta_file = st.file_uploader("Upload metadata (CSV/TSV)", type=["csv", "tsv", "txt"])

    st.subheader("Parameters")
    col3, col4, col5 = st.columns(3)
    with col3:
        n_variable_genes = st.slider("Variable genes", 500, 5000, 2000, 100)
    with col4:
        n_pcs = st.slider("Number of PCs", 5, 20, 10)
    with col5:
        min_cpm = st.slider("Min CPM", 0.1, 5.0, 1.0, 0.1)

    if counts_file and meta_file:
        if st.button("Run Analysis", type="primary"):
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                counts_path = tmpdir / "counts.csv"
                meta_path = tmpdir / "metadata.csv"
                output_dir = tmpdir / "output"
                output_dir.mkdir()

                counts_path.write_bytes(counts_file.getvalue())
                meta_path.write_bytes(meta_file.getvalue())

                progress = st.progress(0)
                status = st.empty()

                stages = [
                    "Validating inputs",
                    "Running quality control",
                    "Normalizing",
                    "Running PCA",
                    "Testing associations",
                    "Detecting outliers",
                    "Generating report",
                ]

                try:
                    from batch_detective.validator import validate_inputs
                    status.text(f"[1/7] {stages[0]}...")
                    progress.progress(10)
                    counts_raw, metadata_raw, counts_working, metadata_working = validate_inputs(
                        counts_path, meta_path
                    )

                    from batch_detective.qc import QualityController
                    status.text(f"[2/7] {stages[1]}...")
                    progress.progress(25)
                    qc = QualityController(counts_working, metadata_working, min_cpm=min_cpm)
                    covariate_info = qc.run_covariate_prescreening()
                    qc_summary = qc.run_sample_qc()
                    qc.run_power_assessment()
                    counts_filtered, _ = qc.filter_genes()

                    from batch_detective.normalizer import normalize_counts
                    status.text(f"[3/7] {stages[2]}...")
                    progress.progress(40)
                    _, log1p_selected, _ = normalize_counts(counts_filtered, n_variable_genes)

                    from batch_detective.pca_analysis import run_pca
                    status.text(f"[4/7] {stages[3]}...")
                    progress.progress(55)
                    pc_scores, evr, components, sample_ids = run_pca(log1p_selected, n_pcs)

                    from batch_detective.association import (
                        detect_collinearity, run_pc_associations, compute_all_icc
                    )
                    status.text(f"[5/7] {stages[4]}...")
                    progress.progress(70)
                    coll_warnings = detect_collinearity(metadata_working, covariate_info)
                    assoc_df = run_pc_associations(pc_scores, evr, metadata_working, covariate_info)
                    icc_table = compute_all_icc(log1p_selected, metadata_working, covariate_info)

                    from batch_detective.outliers import detect_outliers
                    status.text(f"[6/7] {stages[5]}...")
                    progress.progress(83)
                    outlier_df, detection_info = detect_outliers(pc_scores, sample_ids, metadata_working)

                    from batch_detective.report import assemble_report
                    from batch_detective.plots import generate_all_plots
                    status.text(f"[7/7] {stages[6]}...")
                    progress.progress(95)

                    plots = generate_all_plots(
                        qc_summary, pc_scores, evr, metadata_working,
                        covariate_info, icc_table, assoc_df,
                        outlier_df["sample_id"].tolist() if not outlier_df.empty else [],
                        sample_ids, {}, dpi=100,
                        low_power=getattr(qc, "low_power", False),
                    )

                    from datetime import datetime, timezone
                    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                    config = {
                        "n_variable_genes": n_variable_genes,
                        "n_pcs": n_pcs,
                        "min_cpm": min_cpm,
                        "technical_covariates": [],
                        "biological_covariates": [],
                        "condition_on_list": [],
                    }
                    analysis_meta = {
                        "n_samples": len(sample_ids),
                        "n_genes_input": counts_working.shape[0],
                        "n_genes_analyzed": log1p_selected.shape[0],
                        "normalization_method": "CPM + log1p",
                        "low_power": getattr(qc, "low_power", False),
                        "min_detectable_icc": getattr(qc, "min_detectable_icc", 0.0),
                    }
                    report_path = output_dir / "report.html"
                    assemble_report(
                        report_path, qc_summary, icc_table, assoc_df, outlier_df,
                        plots, {}, covariate_info, coll_warnings,
                        getattr(qc, "repeated_measures_covariates", []),
                        config, analysis_meta, detection_info, ts,
                    )

                    progress.progress(100)
                    status.text("✅ Analysis complete!")

                    with open(report_path, encoding="utf-8") as f:
                        html_content = f.read()

                    st.components.v1.html(html_content, height=1000, scrolling=True)

                    # Download buttons
                    st.download_button(
                        "📥 Download report.html",
                        data=html_content.encode("utf-8"),
                        file_name="batch_detective_report.html",
                        mime="text/html",
                    )
                    if not qc_summary.empty:
                        st.download_button(
                            "📥 Download qc_summary.csv",
                            data=qc_summary.to_csv(index=False).encode("utf-8"),
                            file_name="qc_summary.csv",
                        )
                    if not assoc_df.empty:
                        st.download_button(
                            "📥 Download association_table.csv",
                            data=assoc_df.to_csv(index=False).encode("utf-8"),
                            file_name="association_table.csv",
                        )
                    if not icc_table.empty:
                        export_icc = icc_table.drop(columns=["_icc_values"], errors="ignore")
                        st.download_button(
                            "📥 Download icc_table.csv",
                            data=export_icc.to_csv(index=False).encode("utf-8"),
                            file_name="icc_table.csv",
                        )

                except Exception as e:
                    st.error(f"Analysis failed: {e}")
                    raise

    st.markdown(
        "---\n*Running on a compute cluster? See README for SSH port forwarding.*"
    )


if __name__ == "__main__":
    main()
