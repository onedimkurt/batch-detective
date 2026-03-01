#!/usr/bin/env python3
"""
tests/validation/run_validation_panel.py
========================================
Validation panel CI entry point for batch-detective.

Runs all 6 validation datasets and exits non-zero if any test fails,
so it can be used as a release gate in CI (GitHub Actions, etc.).

Usage
-----
  # from project root:
  python tests/validation/run_validation_panel.py

  # verbose with full output:
  python tests/validation/run_validation_panel.py --verbose

  # skip downloads (use cached data):
  python tests/validation/run_validation_panel.py --no-download

Exit codes
----------
  0  All tests passed
  1  One or more tests failed
  2  Script error (missing dependency, bad environment, etc.)

CI example (.github/workflows/validate.yml)
-------------------------------------------
  - name: Run validation panel
    run: python tests/validation/run_validation_panel.py
"""

import sys
import os
import argparse
import subprocess
import datetime
import json
import shutil

# ---------------------------------------------------------------------------
# Resolve paths relative to this script so it works from any working directory
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
MAIN_SCRIPT = os.path.join(SCRIPT_DIR, "run_validation_final.py")
REPORT_SCRIPT = os.path.join(SCRIPT_DIR, "generate_validation_report.py")
RESULTS_JSON = os.path.join(RESULTS_DIR, "validation_panel_results.json")


def log(msg, verbose=True):
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    if verbose:
        print(f"[{ts}] {msg}", flush=True)


def check_environment():
    """Verify batch-detective is installed and importable."""
    result = subprocess.run(
        ["batch-detective", "--version"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("ERROR: batch-detective not found or not installed.")
        print("Install with: pip install -e .")
        sys.exit(2)
    return result.stdout.strip()


def run_panel(verbose=False, no_download=False):
    """Execute run_validation_final.py as a subprocess."""
    env = os.environ.copy()
    env["VALIDATION_RESULTS_DIR"] = RESULTS_DIR
    if no_download:
        env["VALIDATION_NO_DOWNLOAD"] = "1"

    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Copy results dir reference into the main script's expected location
    # (the main script writes to validation_panel/ relative to cwd)
    os.chdir(PROJECT_ROOT)

    cmd = [sys.executable, MAIN_SCRIPT]
    if verbose:
        result = subprocess.run(cmd, env=env)
    else:
        result = subprocess.run(cmd, env=env, capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stdout[-2000:])
            print(result.stderr[-500:])

    return result.returncode


def load_results():
    """Load and return the results JSON, or None if missing."""
    json_path = os.path.join(PROJECT_ROOT, "validation_panel",
                             "validation_panel_results.json")
    if not os.path.exists(json_path):
        return None
    with open(json_path) as f:
        return json.load(f)


def print_summary(results, verbose=False):
    """Print a compact CI-friendly summary table."""
    datasets = results.get("datasets", [])
    summary = results.get("summary", {})

    print()
    print("=" * 65)
    print("  batch-detective Validation Panel — Results")
    print("=" * 65)
    print(f"  {'#':<4} {'Result':<8} {'ICC':>7}  {'V':>6}  Name")
    print(f"  {'-'*4} {'-'*8} {'-'*7}  {'-'*6}  {'-'*35}")

    all_passed = True
    for d in datasets:
        passed = d.get("passed", False)
        if not passed:
            all_passed = False
        status = "PASS ✓" if passed else "FAIL ✗"
        icc = d.get("batch_icc")
        icc_str = f"{icc:.3f}" if icc is not None else "  N/A"
        cv = (d.get("cramers_v_batch_vs_chemical")
              or d.get("cramers_v_batch_vs_genotype")
              or d.get("cramers_v_batch_vs_treatment")
              or d.get("cramers_v_processing_vs_tissue"))
        cv_str = f"{cv:.3f}" if cv is not None else "  N/A"
        name = d.get("name", "")[:45]
        print(f"  {d['test_id']:<4} {status:<8} {icc_str:>7}  {cv_str:>6}  {name}")

    passed_count = summary.get("passed", sum(1 for d in datasets if d.get("passed")))
    total = summary.get("total", len(datasets))
    print()
    print(f"  Result: {passed_count}/{total} passed")
    print("=" * 65)
    print()
    return all_passed


def generate_report():
    """Generate the HTML report if the report script exists."""
    if os.path.exists(REPORT_SCRIPT):
        subprocess.run([sys.executable, REPORT_SCRIPT], cwd=PROJECT_ROOT)


def main():
    parser = argparse.ArgumentParser(
        description="Run batch-detective validation panel"
    )
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Stream full output from validation run")
    parser.add_argument("--no-download", action="store_true",
                        help="Skip downloads, use cached data only")
    parser.add_argument("--report-only", action="store_true",
                        help="Skip running tests, just regenerate HTML report")
    args = parser.parse_args()

    version = check_environment()
    log(f"batch-detective version: {version}", verbose=True)
    log(f"Project root: {PROJECT_ROOT}", verbose=True)

    if not args.report_only:
        log("Starting validation panel...", verbose=True)
        rc = run_panel(verbose=args.verbose, no_download=args.no_download)
        if rc not in (0, 1, 2):
            print(f"ERROR: Validation script exited with unexpected code {rc}")
            sys.exit(2)

    results = load_results()
    if results is None:
        print("ERROR: Results JSON not found after run.")
        sys.exit(2)

    all_passed = print_summary(results, verbose=args.verbose)
    generate_report()

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
