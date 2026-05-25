"""Tests for nwt_substrate.cosmology.eta_B."""

from nwt_substrate.cosmology.eta_B import (
    ETA_B_PRED,
    ETA_B_DENOM,
    ETA_B_PLANCK,
    ETA_B_BBN,
    eta_B,
    summary,
)
from nwt_substrate.isa.constants import ALPHA_SUBSTRATE, RANK_SO7, N_VERTICES_K7


def test_formula_is_three_alpha4_over_14():
    assert ETA_B_DENOM == 14 == 2 * N_VERTICES_K7
    assert eta_B() == ETA_B_PRED == RANK_SO7 * ALPHA_SUBSTRATE**4 / 14


def test_value_near_6e_minus_10():
    # 3α⁴/14 ≈ 6.08e-10
    assert 6.0e-10 < ETA_B_PRED < 6.2e-10


def test_matches_observation():
    # headline: −0.4% vs Planck; −1.0% vs BBN (light-element abundances)
    assert abs(ETA_B_PRED / ETA_B_PLANCK - 1.0) < 0.005
    assert abs(ETA_B_PRED / ETA_B_BBN - 1.0) < 0.011


def test_summary_keys():
    s = summary()
    assert s["numerator"] == 3
    assert s["denominator"] == 14
    assert 0.99 < s["pred_over_planck"] < 1.0
