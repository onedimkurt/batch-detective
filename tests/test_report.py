"""Tests for report module."""

import pandas as pd
import pytest


def test_traffic_light_red():
    from batch_detective.report import _determine_traffic_light
    icc = pd.DataFrame({
        "covariate": ["batch"],
        "median_icc": [0.45],
        "ci_lower_95": [0.35],
        "ci_upper_95": [0.55],
        "icc_tier": ["moderate"],
        "label": ["technical"],
    })
    assoc = pd.DataFrame({
        "covariate": ["batch"],
        "significant_q05": [True],
    })
    result = _determine_traffic_light(icc, assoc, False, ["batch"])
    assert result["emoji"] == "🔴"


def test_traffic_light_green():
    from batch_detective.report import _determine_traffic_light
    icc = pd.DataFrame({
        "covariate": ["batch"],
        "median_icc": [0.05],
        "ci_lower_95": [0.02],
        "ci_upper_95": [0.08],
        "icc_tier": ["negligible"],
        "label": ["technical"],
    })
    assoc = pd.DataFrame({
        "covariate": ["batch"],
        "significant_q05": [False],
    })
    result = _determine_traffic_light(icc, assoc, False, ["batch"])
    assert result["emoji"] == "🟢"


def test_traffic_light_no_tech_covariates():
    from batch_detective.report import _determine_traffic_light
    icc = pd.DataFrame({"covariate": ["batch"], "median_icc": [0.4], "icc_tier": ["moderate"], "label": ["unknown"]})
    assoc = pd.DataFrame({"covariate": ["batch"], "significant_q05": [True]})
    result = _determine_traffic_light(icc, assoc, False, [])
    assert result["emoji"] == "🔵"
