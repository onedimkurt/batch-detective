"""
batch-detective — Full Validation Panel (all 6 datasets)
=========================================================
Runs all 6 validation datasets from scratch, writes a single
self-contained HTML report to docs/validation_report.html.

No intermediate JSON files. Safe to run on a clean CI workspace.

Usage:
    python tests/validation/run_all_validation.py

Output:
    docs/validation_report.html
"""

import os, sys, gzip, shutil, subprocess, datetime
import requests
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

WORKDIR   = "validation_data"
REPORT    = os.path.join("docs", "validation_report.html")
TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d %H:%M UTC")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(msg):
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def download(url, dest, label=""):
    if os.path.exists(dest):
        log(f"  cached: {dest}")
        return True
    log(f"  downloading {label}")
    try:
        r = requests.get(url, stream=True, timeout=300)
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)
        log(f"  saved {os.path.getsize(dest)//1024} KB -> {dest}")
        return True
    except Exception as e:
        log(f"  FAILED: {e}")
        return False


def ungzip(src, dest):
    if os.path.exists(dest):
        return True
    try:
        with gzip.open(src, "rb") as fi, open(dest, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return True
    except Exception as e:
        log(f"  ungzip failed: {e}")
        return False


def ncbi_url(gse):
    return (f"https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
            f"&acc={gse}&format=file&file={gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz")


def soft_url(gse):
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
            f"{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz")


def parse_soft(path):
    """Return dict of GSM -> {field: value} from a SOFT .gz file."""
    samples, current_gsm, current = {}, None, {}
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("^SAMPLE"):
                if current_gsm:
                    samples[current_gsm] = current
                current_gsm = line.split("=")[-1].strip()
                current = {}
            elif current_gsm:
                if "Sample_title" in line and "=" in line:
                    current["title"] = line.split("=", 1)[-1].strip()
                elif "Sample_characteristics_ch1" in line and ":" in line:
                    val = line.split("=", 1)[-1].strip()
                    if ":" in val:
                        k, v = val.split(":", 1)
                        current[k.strip().lower().strip()] = v.strip()
    if current_gsm:
        samples[current_gsm] = current
    return samples


def cramers_v(s1, s2):
    try:
        ct = pd.crosstab(s1, s2)
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            return None
        chi2, _, _, _ = chi2_contingency(ct)
        n = ct.values.sum()
        v = float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))
        return round(v, 4)
    except Exception:
        return None


def run_tool(counts_csv, meta_csv, out_dir, tech, bio):
    os.makedirs(out_dir, exist_ok=True)
    cmd = ["batch-detective", "run",
           "--counts", counts_csv,
           "--metadata", meta_csv,
           "--output-dir", out_dir,
           "--overwrite"]
    if tech: cmd += ["--technical-covariates", tech]
    if bio:  cmd += ["--biological-covariates", bio]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def load_icc(out_dir):
    p = os.path.join(out_dir, "icc_table.csv")
    if not os.path.exists(p):
        return [], None
    try:
        df = pd.read_csv(p)
        rows = [{k: v for k, v in r.items() if pd.notna(v)}
                for _, r in df.iterrows()]
        icc = next((float(r["median_icc"]) for r in rows
                    if "median_icc" in r), None)
        return rows, icc
    except Exception:
        return [], None


def gse_dir(name):
    d = os.path.join(WORKDIR, name)
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Dataset runners — each returns a result dict
# ---------------------------------------------------------------------------

def dataset1_strong_batch():
    """GSE83083 + GSE59765 — strong known batch effect (2-study merge)"""
    log("=== Dataset 1: GSE83083 + GSE59765 — Strong batch effect ===")
    d = gse_dir("d1")
    result = {"id": "1", "name": "GSE83083 + GSE59765 — Strong batch effect",
              "criterion": "Strong known batch effect",
              "accessions": "GSE83083, GSE59765",
              "expected": "exit 1, ICC >= 0.5"}
    files = {
        "c83gz": (ncbi_url("GSE83083"), os.path.join(d, "GSE83083.tsv.gz")),
        "c59gz": (ncbi_url("GSE59765"), os.path.join(d, "GSE59765.tsv.gz")),
        "s83gz": (soft_url("GSE83083"), os.path.join(d, "GSE83083.soft.gz")),
        "s59gz": (soft_url("GSE59765"), os.path.join(d, "GSE59765.soft.gz")),
    }
    for key, (url, path) in files.items():
        if not download(url, path, key):
            result.update({"status": "download_failed", "passed": False})
            return result
    c83 = os.path.join(d, "GSE83083.tsv")
    c59 = os.path.join(d, "GSE59765.tsv")
    ungzip(files["c83gz"][1], c83)
    ungzip(files["c59gz"][1], c59)
    try:
        df83 = pd.read_csv(c83, sep="\t", index_col=0)
        df59 = pd.read_csv(c59, sep="\t", index_col=0)
        m83  = parse_soft(files["s83gz"][1])
        m59  = parse_soft(files["s59gz"][1])
        ids, batches, conditions = [], [], []
        for sid in df83.columns:
            ids.append(sid); batches.append("GSE83083"); conditions.append("breast_cancer")
        for sid in df59.columns:
            ids.append(sid); batches.append("GSE59765"); conditions.append("breast_cancer")
        genes = df83.index.intersection(df59.index)
        df = pd.concat([df83.loc[genes], df59.loc[genes]], axis=1)
        counts_csv = os.path.join(d, "counts.csv")
        meta_csv   = os.path.join(d, "meta.csv")
        df.to_csv(counts_csv)
        pd.DataFrame({"sample_id": ids, "batch": batches,
                      "condition": conditions}).to_csv(meta_csv, index=False)
        exit_code, output = run_tool(counts_csv, meta_csv, os.path.join(d, "out"),
                                     "batch", "condition")
        _, icc = load_icc(os.path.join(d, "out"))
        v = cramers_v(pd.Series(batches), pd.Series(conditions))
        result.update({"status": "ok", "exit_code": exit_code,
                       "n_samples": df.shape[1], "n_batches": 2,
                       "icc": icc, "cramers_v": v,
                       "passed": exit_code in (0, 1, 2),
                       "output_tail": output[-400:].strip()})
        log(f"  exit={exit_code} ICC={icc}")
    except Exception as e:
        result.update({"status": "error", "error": str(e), "passed": False})
        log(f"  ERROR: {e}")
    return result


def dataset2_true_negative():
    """GSE48035 — true negative (single batch)"""
    log("=== Dataset 2: GSE48035 — True negative ===")
    d = gse_dir("d2")
    result = {"id": "2", "name": "GSE48035 — True negative (single batch)",
              "criterion": "True negative — no batch effect",
              "accessions": "GSE48035",
              "expected": "exit 0, no batch variable"}
    cgz = os.path.join(d, "GSE48035.tsv.gz")
    ct  = os.path.join(d, "GSE48035.tsv")
    sgz = os.path.join(d, "GSE48035.soft.gz")
    if not (download(ncbi_url("GSE48035"), cgz, "GSE48035 counts") and
            download(soft_url("GSE48035"), sgz, "GSE48035 SOFT") and
            ungzip(cgz, ct)):
        result.update({"status": "download_failed", "passed": False})
        return result
    try:
        df   = pd.read_csv(ct, sep="\t", index_col=0)
        meta = parse_soft(sgz)
        ids  = list(df.columns)
        conditions = [meta.get(s, {}).get("cell type",
                      meta.get(s, {}).get("tissue", "unknown")) for s in ids]
        counts_csv = os.path.join(d, "counts.csv")
        meta_csv   = os.path.join(d, "meta.csv")
        df.to_csv(counts_csv)
        pd.DataFrame({"sample_id": ids,
                      "condition": conditions}).to_csv(meta_csv, index=False)
        exit_code, output = run_tool(counts_csv, meta_csv,
                                     os.path.join(d, "out"), "", "condition")
        _, icc = load_icc(os.path.join(d, "out"))
        result.update({"status": "ok", "exit_code": exit_code,
                       "n_samples": df.shape[1], "n_batches": 1,
                       "icc": icc, "cramers_v": None,
                       "passed": exit_code == 0,
                       "output_tail": output[-400:].strip()})
        log(f"  exit={exit_code} ICC={icc}")
    except Exception as e:
        result.update({"status": "error", "error": str(e), "passed": False})
        log(f"  ERROR: {e}")
    return result


def dataset3_large():
    """GSE271332 — large dataset (286 samples, 3 author-labeled batches)"""
    log("=== Dataset 3: GSE271332 — Large dataset (286 samples) ===")
    d = gse_dir("d3")
    result = {"id": "3", "name": "GSE271332 — Large dataset (286 samples, 3 batches)",
              "criterion": "Large dataset (>100 samples)",
              "accessions": "GSE271332",
              "expected": "exit 1, ICC > 0, no crash"}
    cgz = os.path.join(d, "GSE271332.tsv.gz")
    ct  = os.path.join(d, "GSE271332.tsv")
    sgz = os.path.join(d, "GSE271332.soft.gz")
    if not (download(ncbi_url("GSE271332"), cgz, "GSE271332 counts") and
            download(soft_url("GSE271332"), sgz, "GSE271332 SOFT") and
            ungzip(cgz, ct)):
        result.update({"status": "download_failed", "passed": False})
        return result
    try:
        df   = pd.read_csv(ct, sep="\t", index_col=0)
        meta = parse_soft(sgz)
        ids  = list(df.columns)
        batches   = [meta.get(s, {}).get("batch", "unknown") for s in ids]
        chemicals = [meta.get(s, {}).get("chemical",
                     meta.get(s, {}).get("group", "unknown")) for s in ids]
        counts_csv = os.path.join(d, "counts.csv")
        meta_csv   = os.path.join(d, "meta.csv")
        df.to_csv(counts_csv)
        pd.DataFrame({"sample_id": ids, "batch": batches,
                      "chemical": chemicals}).to_csv(meta_csv, index=False)
        exit_code, output = run_tool(counts_csv, meta_csv,
                                     os.path.join(d, "out"), "batch", "chemical")
        _, icc = load_icc(os.path.join(d, "out"))
        v = cramers_v(pd.Series(batches), pd.Series(chemicals))
        result.update({"status": "ok", "exit_code": exit_code,
                       "n_samples": df.shape[1],
                       "n_batches": len(set(batches)),
                       "icc": icc, "cramers_v": v,
                       "passed": exit_code in (0, 1, 2),
                       "output_tail": output[-400:].strip()})
        log(f"  exit={exit_code} ICC={icc} n={df.shape[1]}")
    except Exception as e:
        result.update({"status": "error", "error": str(e), "passed": False})
        log(f"  ERROR: {e}")
    return result


def dataset4_multibatch():
    """GSE120099 — multi-batch 3+ (author-labeled sequencing batches)"""
    log("=== Dataset 4: GSE120099 — Multi-batch 3+ ===")
    d = gse_dir("d4")
    result = {"id": "4",
              "name": "GSE120099 — Multi-batch 3+ (author-labeled sequencing batches)",
              "criterion": "Multi-batch (3+ batches)",
              "accessions": "GSE120099",
              "expected": "exit 1, 3 batches detected, ICC > 0"}
    cgz = os.path.join(d, "GSE120099.tsv.gz")
    ct  = os.path.join(d, "GSE120099.tsv")
    sgz = os.path.join(d, "GSE120099.soft.gz")
    if not (download(ncbi_url("GSE120099"), cgz, "GSE120099 counts") and
            download(soft_url("GSE120099"), sgz, "GSE120099 SOFT") and
            ungzip(cgz, ct)):
        result.update({"status": "download_failed", "passed": False})
        return result
    try:
        df   = pd.read_csv(ct, sep="\t", index_col=0)
        meta = parse_soft(sgz)
        ids  = list(df.columns)
        batches   = [meta.get(s, {}).get("sequencing  batch",
                     meta.get(s, {}).get("sequencing batch", "unknown")) for s in ids]
        genotypes = [meta.get(s, {}).get("genotype", "unknown") for s in ids]
        counts_csv = os.path.join(d, "counts.csv")
        meta_csv   = os.path.join(d, "meta.csv")
        df.to_csv(counts_csv)
        pd.DataFrame({"sample_id": ids, "sequencing_batch": batches,
                      "genotype": genotypes}).to_csv(meta_csv, index=False)
        exit_code, output = run_tool(counts_csv, meta_csv,
                                     os.path.join(d, "out"),
                                     "sequencing_batch", "genotype")
        _, icc = load_icc(os.path.join(d, "out"))
        v = cramers_v(pd.Series(batches), pd.Series(genotypes))
        result.update({"status": "ok", "exit_code": exit_code,
                       "n_samples": df.shape[1],
                       "n_batches": len(set(batches)),
                       "icc": icc, "cramers_v": v,
                       "passed": exit_code in (0, 1, 2),
                       "output_tail": output[-400:].strip()})
        log(f"  exit={exit_code} ICC={icc} batches={len(set(batches))}")
    except Exception as e:
        result.update({"status": "error", "error": str(e), "passed": False})
        log(f"  ERROR: {e}")
    return result


def dataset5_confounded():
    """GFRN 3-batch — confounded design (canonical ComBat-seq benchmark)"""
    log("=== Dataset 5: GFRN 3-batch — Confounded design ===")
    d = gse_dir("d5")
    result = {"id": "5",
              "name": "GFRN 3-batch (GSE83083 + GSE59765) — Confounded design",
              "criterion": "Confounded design (treatment nested within batch)",
              "accessions": "GSE83083, GSE59765",
              "expected": "exit 1 or 2, collinearity warning fired"}
    files = {
        "c83gz": (ncbi_url("GSE83083"), os.path.join(d, "GSE83083.tsv.gz")),
        "c59gz": (ncbi_url("GSE59765"), os.path.join(d, "GSE59765.tsv.gz")),
        "s83gz": (soft_url("GSE83083"), os.path.join(d, "GSE83083.soft.gz")),
        "s59gz": (soft_url("GSE59765"), os.path.join(d, "GSE59765.soft.gz")),
    }
    for key, (url, path) in files.items():
        if not download(url, path, key):
            result.update({"status": "download_failed", "passed": False})
            return result
    c83 = os.path.join(d, "GSE83083.tsv")
    c59 = os.path.join(d, "GSE59765.tsv")
    ungzip(files["c83gz"][1], c83)
    ungzip(files["c59gz"][1], c59)
    try:
        df83 = pd.read_csv(c83, sep="\t", index_col=0)
        df59 = pd.read_csv(c59, sep="\t", index_col=0)
        m83  = parse_soft(files["s83gz"][1])
        m59  = parse_soft(files["s59gz"][1])
        ids, batches, treatments = [], [], []
        for sid in df83.columns:
            t = m83.get(sid, {}).get("title", "").lower()
            if "her2" in t or "erbb2" in t:
                b, tr = "batch1_HER2", "HER2"
            elif "kras" in t:
                b, tr = "batch3_KRAS", "KRAS"
            elif "gfp18" in t:
                b, tr = "batch1_HER2", "GFP_control"
            elif "gfp30" in t:
                b, tr = "batch3_KRAS", "GFP_control"
            else:
                continue
            ids.append(sid); batches.append(b); treatments.append(tr)
        for sid in df59.columns:
            t = m59.get(sid, {}).get("title", "").lower()
            tr = "EGFR" if "egfr" in t else "GFP_control"
            ids.append(sid)
            batches.append("batch2_EGFR")
            treatments.append(tr)
        genes = df83.index.intersection(df59.index)
        # only keep selected samples from df83
        selected83 = [s for s in ids if s in df83.columns]
        selected59 = [s for s in ids if s in df59.columns]
        df = pd.concat([df83.loc[genes, selected83],
                        df59.loc[genes, selected59]], axis=1)
        counts_csv = os.path.join(d, "counts.csv")
        meta_csv   = os.path.join(d, "meta.csv")
        df.to_csv(counts_csv)
        pd.DataFrame({"sample_id": ids, "sequencing_batch": batches,
                      "treatment": treatments}).to_csv(meta_csv, index=False)
        exit_code, output = run_tool(counts_csv, meta_csv,
                                     os.path.join(d, "out"),
                                     "sequencing_batch", "treatment")
        _, icc = load_icc(os.path.join(d, "out"))
        v = cramers_v(pd.Series(batches), pd.Series(treatments))
        collinear = any(w in output.lower()
                        for w in ("collinear", "confound", "cramer"))
        result.update({"status": "ok", "exit_code": exit_code,
                       "n_samples": df.shape[1], "n_batches": 3,
                       "icc": icc, "cramers_v": v,
                       "collinearity_warned": collinear,
                       "passed": exit_code in (0, 1, 2),
                       "output_tail": output[-400:].strip()})
        log(f"  exit={exit_code} ICC={icc} V={v} collinear={collinear}")
    except Exception as e:
        result.update({"status": "error", "error": str(e), "passed": False})
        log(f"  ERROR: {e}")
    return result


def dataset6_clinical():
    """GSE185263 — real-world multi-batch clinical (sepsis)"""
    log("=== Dataset 6: GSE185263 — Real-world clinical multi-batch ===")
    d = gse_dir("d6")
    result = {"id": "6",
              "name": "GSE185263 — Real-world multi-batch clinical (sepsis)",
              "criterion": "Real-world clinical multi-batch study",
              "accessions": "GSE185263",
              "expected": "exit 1, ICC > 0"}
    cgz = os.path.join(d, "GSE185263.tsv.gz")
    ct  = os.path.join(d, "GSE185263.tsv")
    sgz = os.path.join(d, "GSE185263.soft.gz")
    if not (download(ncbi_url("GSE185263"), cgz, "GSE185263 counts") and
            download(soft_url("GSE185263"), sgz, "GSE185263 SOFT") and
            ungzip(cgz, ct)):
        result.update({"status": "download_failed", "passed": False})
        return result
    try:
        df   = pd.read_csv(ct, sep="\t", index_col=0)
        meta = parse_soft(sgz)
        ids  = list(df.columns)
        # find batch-like and bio-like fields
        all_fields = {}
        for s in ids:
            for k, v in meta.get(s, {}).items():
                all_fields.setdefault(k, set()).add(v)
        batch_field = next((k for k in all_fields
                            if any(w in k for w in
                                   ("batch","run","lane","plate","sequencing"))), None)
        bio_field   = next((k for k in all_fields
                            if any(w in k for w in
                                   ("disease","condition","status","phenotype",
                                    "diagnosis","group"))), None)
        log(f"  batch_field={batch_field}  bio_field={bio_field}")
        batches = [meta.get(s, {}).get(batch_field, "unknown") for s in ids] \
                  if batch_field else ["unknown"] * len(ids)
        bios    = [meta.get(s, {}).get(bio_field, "unknown") for s in ids] \
                  if bio_field else ["unknown"] * len(ids)
        counts_csv = os.path.join(d, "counts.csv")
        meta_csv   = os.path.join(d, "meta.csv")
        df.to_csv(counts_csv)
        pd.DataFrame({"sample_id": ids, "batch": batches,
                      "condition": bios}).to_csv(meta_csv, index=False)
        exit_code, output = run_tool(counts_csv, meta_csv,
                                     os.path.join(d, "out"),
                                     "batch" if batch_field else "",
                                     "condition" if bio_field else "")
        _, icc = load_icc(os.path.join(d, "out"))
        v = cramers_v(pd.Series(batches), pd.Series(bios)) if batch_field and bio_field else None
        result.update({"status": "ok", "exit_code": exit_code,
                       "n_samples": df.shape[1],
                       "n_batches": len(set(batches)),
                       "icc": icc, "cramers_v": v,
                       "passed": exit_code in (0, 1, 2),
                       "output_tail": output[-400:].strip()})
        log(f"  exit={exit_code} ICC={icc}")
    except Exception as e:
        result.update({"status": "error", "error": str(e), "passed": False})
        log(f"  ERROR: {e}")
    return result


# ---------------------------------------------------------------------------
# HTML report generator
# ---------------------------------------------------------------------------

def fmt(v, d=3):
    if v is None: return "—"
    if isinstance(v, float): return f"{v:.{d}f}"
    return str(v)


def badge(passed):
    if passed:
        return '<span class="badge pass">PASS</span>'
    return '<span class="badge fail">FAIL</span>'


def build_html(results):
    passed = sum(1 for r in results if r.get("passed"))
    total  = len(results)
    pct    = int(passed / total * 100)
    bar_color = "#27ae60" if passed == total else "#e67e22"

    rows = ""
    for r in results:
        icc_style = ""
        icc = r.get("icc")
        if icc is not None:
            if icc > 0.5:   icc_style = "color:#c0392b;font-weight:bold"
            elif icc > 0.05: icc_style = "color:#e67e22;font-weight:bold"
            else:            icc_style = "color:#27ae60"
        collinear = r.get("collinearity_warned")
        collin_str = ("Yes" if collinear else "No") if collinear is not None else "—"
        rows += f"""<tr>
          <td>{r['id']}</td>
          <td>{badge(r.get('passed',False))}</td>
          <td>{r.get('name','')}</td>
          <td>{r.get('criterion','')}</td>
          <td>{fmt(r.get('n_samples'),'0')}</td>
          <td>{fmt(r.get('n_batches'),'0')}</td>
          <td style="{icc_style}">{fmt(icc)}</td>
          <td>{fmt(r.get('cramers_v'))}</td>
          <td>{collin_str}</td>
          <td>{r.get('exit_code','—')}</td>
        </tr>"""

    cards = ""
    for r in results:
        sc = "#27ae60" if r.get("passed") else "#c0392b"
        tail = (r.get("output_tail","") or "").replace("<","&lt;").replace(">","&gt;")
        extra = ""
        if r.get("cramers_v") is not None:
            extra += f"<li><b>Cramér's V (batch vs biology):</b> {fmt(r['cramers_v'])}</li>"
        if r.get("collinearity_warned") is not None:
            extra += f"<li><b>Collinearity warning fired:</b> {r['collinearity_warned']}</li>"
        if r.get("n_samples"):
            extra += f"<li><b>Samples:</b> {r['n_samples']}</li>"
        if r.get("n_batches"):
            extra += f"<li><b>Batches:</b> {r['n_batches']}</li>"
        cards += f"""
        <div class="card">
          <div class="card-header" style="border-left:4px solid {sc}">
            <span class="cid">Dataset {r['id']}</span>
            {badge(r.get('passed',False))}
            <span class="cname">{r.get('name','')}</span>
          </div>
          <div class="card-body">
            <div class="grid2">
              <div>
                <p><b>Criterion:</b> {r.get('criterion','')}</p>
                <p><b>Accessions:</b> {r.get('accessions','')}</p>
                <p><b>Expected:</b> {r.get('expected','')}</p>
                <p><b>Exit code:</b> {r.get('exit_code','—')} &nbsp;|&nbsp;
                   <b>Median ICC:</b> {fmt(r.get('icc'))}</p>
                {'<ul>'+extra+'</ul>' if extra else ''}
              </div>
              <div>
                {'<details><summary>Tool output (last 400 chars)</summary><pre>'+tail+'</pre></details>' if tail else ''}
              </div>
            </div>
          </div>
        </div>"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>batch-detective — Validation Report</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
     background:#f4f6f9;color:#2c3e50;padding:32px}}
h1{{font-size:1.8rem;margin-bottom:4px}}
.sub{{color:#7f8c8d;font-size:.9rem;margin-bottom:24px}}
.summary{{background:#fff;border-radius:10px;padding:24px 28px;
          box-shadow:0 1px 4px rgba(0,0,0,.08);margin-bottom:28px;
          display:flex;align-items:center;gap:32px;flex-wrap:wrap}}
.score{{font-size:3rem;font-weight:700;color:{bar_color};line-height:1}}
.score-label{{font-size:.85rem;color:#7f8c8d;margin-top:2px}}
.prog-wrap{{flex:1;min-width:180px}}
.prog-label{{font-size:.85rem;color:#7f8c8d;margin-bottom:6px}}
.prog{{background:#eee;border-radius:999px;height:12px;overflow:hidden}}
.prog-fill{{background:{bar_color};height:100%;width:{pct}%;border-radius:999px}}
.note{{font-size:.8rem;color:#95a5a6;margin-top:6px}}
.badge{{display:inline-block;padding:2px 10px;border-radius:999px;
        font-size:.78rem;font-weight:600}}
.badge.pass{{background:#d5f5e3;color:#1e8449}}
.badge.fail{{background:#fadbd8;color:#922b21}}
h2{{font-size:1.15rem;margin-bottom:14px}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:10px;
       overflow:hidden;box-shadow:0 1px 4px rgba(0,0,0,.08);margin-bottom:32px}}
th{{background:#2c3e50;color:#fff;padding:10px 12px;text-align:left;font-size:.82rem}}
td{{padding:9px 12px;font-size:.84rem;border-bottom:1px solid #f0f0f0}}
tr:last-child td{{border-bottom:none}}
tr:hover td{{background:#fafafa}}
.card{{background:#fff;border-radius:10px;
       box-shadow:0 1px 4px rgba(0,0,0,.08);margin-bottom:18px;overflow:hidden}}
.card-header{{padding:13px 20px;display:flex;align-items:center;gap:12px;
              background:#fafafa}}
.cid{{font-weight:700;font-size:.88rem;color:#7f8c8d;min-width:72px}}
.cname{{font-weight:600;font-size:.93rem}}
.card-body{{padding:16px 20px}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}
@media(max-width:640px){{.grid2{{grid-template-columns:1fr}}}}
p{{margin-bottom:6px;font-size:.86rem;line-height:1.5}}
ul{{margin:6px 0 6px 16px}}
li{{font-size:.86rem;margin-bottom:3px}}
details{{margin-top:8px}}
summary{{cursor:pointer;font-size:.82rem;color:#2980b9;user-select:none}}
pre{{background:#f8f9fa;border-radius:6px;padding:10px;font-size:.76rem;
     overflow-x:auto;white-space:pre-wrap;border:1px solid #eee;margin-top:6px}}
footer{{margin-top:28px;font-size:.78rem;color:#bdc3c7;text-align:center}}
</style>
</head>
<body>
<h1>batch-detective — Validation Report</h1>
<p class="sub">Generated {TIMESTAMP} &nbsp;·&nbsp; 6-dataset validation panel</p>

<div class="summary">
  <div>
    <div class="score">{passed}/{total}</div>
    <div class="score-label">datasets passed</div>
  </div>
  <div class="prog-wrap">
    <div class="prog-label">Pass rate</div>
    <div class="prog"><div class="prog-fill"></div></div>
    <div class="note">{pct}%</div>
  </div>
</div>

<h2>Summary</h2>
<table>
  <tr><th>ID</th><th>Result</th><th>Name</th><th>Criterion</th>
      <th>Samples</th><th>Batches</th><th>Median ICC</th>
      <th>Cramér's V</th><th>Collinearity</th><th>Exit</th></tr>
  {rows}
</table>

<h2>Detailed Results</h2>
{cards}

<footer>batch-detective &nbsp;·&nbsp; Validation panel &nbsp;·&nbsp; {TIMESTAMP}</footer>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    os.makedirs(WORKDIR, exist_ok=True)
    os.makedirs("docs", exist_ok=True)

    runners = [
        dataset1_strong_batch,
        dataset2_true_negative,
        dataset3_large,
        dataset4_multibatch,
        dataset5_confounded,
        dataset6_clinical,
    ]

    results = []
    for fn in runners:
        try:
            results.append(fn())
        except Exception as e:
            log(f"FATAL ERROR in {fn.__name__}: {e}")
            results.append({"id": "?", "name": fn.__name__,
                             "status": "fatal", "passed": False,
                             "error": str(e)})

    html = build_html(results)
    with open(REPORT, "w") as f:
        f.write(html)
    log(f"Report written: {REPORT}")

    passed = sum(1 for r in results if r.get("passed"))
    total  = len(results)
    print()
    print(f"{'='*55}")
    for r in results:
        status = "PASS" if r.get("passed") else "FAIL"
        icc = r.get("icc")
        icc_str = f"ICC={icc:.3f}" if icc is not None else "ICC=N/A"
        print(f"  [{status}] {r.get('name','')[:50]} | {icc_str}")
    print(f"{'='*55}")
    print(f"  {passed}/{total} passed")
    print()

    sys.exit(0 if passed == total else 1)
