"""
Generate a clean HTML validation report for batch-detective v0.1.1
from validation_panel/validation_panel_results.json
Usage:
    cd ~/Downloads/batch_detective
    python generate_validation_report.py
"""
import json, datetime, os, sys

RESULTS_JSON = "validation_panel/validation_panel_results.json"
OUT_HTML = "validation_panel/validation_report.html"

with open(RESULTS_JSON) as f:
    data = json.load(f)

datasets = data["datasets"]
summary = data.get("summary", {})
generated = data.get("generated", "unknown")

def badge(passed):
    if passed:
        return '<span class="badge pass">PASS</span>'
    return '<span class="badge fail">FAIL</span>'

def fmt(val, decimals=3):
    if val is None: return "—"
    if isinstance(val, float): return f"{val:.{decimals}f}"
    return str(val)

rows = ""
for d in datasets:
    icc = d.get("batch_icc")
    cv  = (d.get("cramers_v_batch_vs_chemical")
           or d.get("cramers_v_batch_vs_genotype")
           or d.get("cramers_v_batch_vs_treatment")
           or d.get("cramers_v_processing_vs_tissue"))
    collin = d.get("collinearity_warning_fired")
    collin_str = ("Yes" if collin else "No") if collin is not None else "—"
    icc_color = ""
    if icc is not None:
        if icc > 0.5: icc_color = "style='color:#c0392b;font-weight:bold'"
        elif icc > 0.05: icc_color = "style='color:#e67e22;font-weight:bold'"
        else: icc_color = "style='color:#27ae60'"

    rows += f"""
    <tr>
      <td>{d['test_id']}</td>
      <td>{badge(d.get('passed', False))}</td>
      <td class="name">{d.get('name','')}</td>
      <td>{d.get('criterion', d.get('type',''))}</td>
      <td>{fmt(d.get('n_samples'),'0')}</td>
      <td>{fmt(d.get('n_batches') or d.get('n_processing_groups'), '0')}</td>
      <td {icc_color}>{fmt(icc)}</td>
      <td>{fmt(cv)}</td>
      <td>{collin_str}</td>
      <td>{d.get('exit_code','—')}</td>
    </tr>"""

detail_sections = ""
for d in datasets:
    status_color = "#27ae60" if d.get("passed") else "#c0392b"
    notes = d.get("notes", "").replace("<","&lt;").replace(">","&gt;")
    icc_rows = d.get("icc_rows", [])
    icc_table = ""
    if icc_rows:
        icc_table = """<table class='icc-table'>
        <tr><th>Covariate</th><th>N batches</th><th>Median ICC</th><th>Mean ICC</th></tr>"""
        for row in icc_rows:
            icc_table += f"""<tr>
              <td>{row.get('covariate','')}</td>
              <td>{row.get('n_batches','')}</td>
              <td>{fmt(row.get('median_icc'))}</td>
              <td>{fmt(row.get('mean_icc'))}</td>
            </tr>"""
        icc_table += "</table>"

    extra = ""
    for k in ["cramers_v_batch_vs_chemical","cramers_v_batch_vs_genotype",
              "cramers_v_batch_vs_treatment","cramers_v_processing_vs_tissue"]:
        if k in d:
            label = k.replace("_"," ").replace("cramers v","Cramér's V")
            extra += f"<li><b>{label}:</b> {fmt(d[k])}</li>"
    if "collinearity_warning_fired" in d:
        extra += f"<li><b>Collinearity warning fired:</b> {d['collinearity_warning_fired']}</li>"
    if "n_genes" in d:
        extra += f"<li><b>Genes in matrix:</b> {d['n_genes']:,}</li>"

    detail_sections += f"""
    <div class="detail-card">
      <div class="detail-header" style="border-left:4px solid {status_color}">
        <span class="detail-id">Dataset {d['test_id']}</span>
        {badge(d.get('passed',False))}
        <span class="detail-name">{d.get('name','')}</span>
      </div>
      <div class="detail-body">
        <div class="meta-grid">
          <div>
            <p><b>Criterion:</b> {d.get('criterion', d.get('type',''))}</p>
            <p><b>Source:</b> {d.get('dataset','')}</p>
            <p><b>Expected behaviour:</b> {d.get('expected','')}</p>
            <p><b>Samples:</b> {d.get('n_samples','—')} &nbsp;|&nbsp;
               <b>Batches:</b> {d.get('n_batches') or d.get('n_processing_groups','—')} &nbsp;|&nbsp;
               <b>Exit code:</b> {d.get('exit_code','—')}</p>
            {'<ul>' + extra + '</ul>' if extra else ''}
          </div>
          <div>
            <b>ICC table:</b>
            {icc_table if icc_table else '<p>No ICC table generated.</p>'}
          </div>
        </div>
        {'<details><summary>Last 500 chars of tool output</summary><pre>' + notes + '</pre></details>' if notes else ''}
      </div>
    </div>"""

passed = summary.get("passed", sum(1 for d in datasets if d.get("passed")))
total  = summary.get("total", 6)
bar_pct = int(passed / total * 100)
bar_color = "#27ae60" if passed == total else "#e67e22"

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>batch-detective v0.1.1 — Validation Report</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          background: #f4f6f9; color: #2c3e50; padding: 32px; }}
  h1 {{ font-size: 1.8rem; margin-bottom: 4px; }}
  .subtitle {{ color: #7f8c8d; margin-bottom: 24px; font-size: 0.95rem; }}
  .summary-bar {{ background:#fff; border-radius:10px; padding:24px 28px;
                  box-shadow:0 1px 4px rgba(0,0,0,0.08); margin-bottom:28px;
                  display:flex; align-items:center; gap:32px; flex-wrap:wrap; }}
  .big-score {{ font-size:3rem; font-weight:700; color:{bar_color}; line-height:1; }}
  .big-label {{ font-size:0.85rem; color:#7f8c8d; margin-top:2px; }}
  .progress-wrap {{ flex:1; min-width:200px; }}
  .progress-label {{ font-size:0.85rem; color:#7f8c8d; margin-bottom:6px; }}
  .progress {{ background:#eee; border-radius:999px; height:12px; overflow:hidden; }}
  .progress-fill {{ background:{bar_color}; height:100%; width:{bar_pct}%;
                    border-radius:999px; transition:width .5s; }}
  .meta-note {{ font-size:0.82rem; color:#95a5a6; }}
  .badge {{ display:inline-block; padding:2px 10px; border-radius:999px;
            font-size:0.78rem; font-weight:600; letter-spacing:.03em; }}
  .badge.pass {{ background:#d5f5e3; color:#1e8449; }}
  .badge.fail {{ background:#fadbd8; color:#922b21; }}
  table {{ width:100%; border-collapse:collapse; background:#fff;
           border-radius:10px; overflow:hidden;
           box-shadow:0 1px 4px rgba(0,0,0,0.08); margin-bottom:32px; }}
  th {{ background:#2c3e50; color:#fff; padding:10px 12px;
        text-align:left; font-size:0.82rem; font-weight:600; }}
  td {{ padding:9px 12px; font-size:0.85rem; border-bottom:1px solid #f0f0f0; }}
  tr:last-child td {{ border-bottom:none; }}
  tr:hover td {{ background:#fafafa; }}
  td.name {{ max-width:280px; }}
  h2 {{ font-size:1.2rem; margin-bottom:16px; color:#2c3e50; }}
  .detail-card {{ background:#fff; border-radius:10px;
                  box-shadow:0 1px 4px rgba(0,0,0,0.08);
                  margin-bottom:20px; overflow:hidden; }}
  .detail-header {{ padding:14px 20px; display:flex; align-items:center;
                    gap:12px; background:#fafafa; }}
  .detail-id {{ font-weight:700; font-size:0.9rem; color:#7f8c8d; min-width:70px; }}
  .detail-name {{ font-weight:600; font-size:0.95rem; }}
  .detail-body {{ padding:18px 20px; }}
  .meta-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:24px; }}
  @media(max-width:700px){{ .meta-grid{{ grid-template-columns:1fr; }} }}
  p {{ margin-bottom:6px; font-size:0.87rem; line-height:1.5; }}
  ul {{ margin:8px 0 8px 18px; }}
  li {{ font-size:0.87rem; margin-bottom:3px; }}
  .icc-table {{ width:100%; border-collapse:collapse; margin-top:6px; }}
  .icc-table th {{ background:#ecf0f1; color:#2c3e50; padding:6px 10px;
                   font-size:0.8rem; text-align:left; }}
  .icc-table td {{ padding:5px 10px; font-size:0.83rem;
                   border-bottom:1px solid #f0f0f0; }}
  details {{ margin-top:12px; }}
  summary {{ cursor:pointer; font-size:0.83rem; color:#2980b9;
             user-select:none; margin-bottom:6px; }}
  pre {{ background:#f8f9fa; border-radius:6px; padding:10px 14px;
         font-size:0.78rem; overflow-x:auto; white-space:pre-wrap;
         color:#2c3e50; border:1px solid #eee; }}
  footer {{ margin-top:32px; font-size:0.8rem; color:#bdc3c7; text-align:center; }}
</style>
</head>
<body>

<h1>batch-detective v0.1.1 — Validation Report</h1>
<p class="subtitle">Generated {generated} &nbsp;·&nbsp; 6-dataset validation panel</p>

<div class="summary-bar">
  <div>
    <div class="big-score">{passed}/{total}</div>
    <div class="big-label">datasets passed</div>
  </div>
  <div class="progress-wrap">
    <div class="progress-label">Overall pass rate</div>
    <div class="progress"><div class="progress-fill"></div></div>
    <div class="meta-note" style="margin-top:6px">{bar_pct}% &nbsp;·&nbsp; Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
  </div>
</div>

<h2>Summary Table</h2>
<table>
  <tr>
    <th>#</th><th>Result</th><th>Name</th><th>Criterion</th>
    <th>Samples</th><th>Batches</th><th>Median ICC</th>
    <th>Cramér's V</th><th>Collinearity</th><th>Exit</th>
  </tr>
  {rows}
</table>

<h2>Detailed Results</h2>
{detail_sections}

<footer>batch-detective v0.1.1 &nbsp;·&nbsp; Validation panel &nbsp;·&nbsp; {datetime.datetime.now().strftime('%Y-%m-%d')}</footer>
</body>
</html>"""

with open(OUT_HTML, "w") as f:
    f.write(html)

print(f"Report written to: {OUT_HTML}")
