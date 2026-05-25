"""Tests for nwt_substrate.cosmology.lambda_cc."""

from nwt_substrate.cosmology.lambda_cc import (
    M_E_OVER_M_PL,
    LAMBDA_EXPONENT,
    RHO_LAMBDA_PRED,
    RHO_LAMBDA_OBS,
    lambda_cc,
    summary,
)
from nwt_substrate.isa.constants import ALPHA_SUBSTRATE, H_COXETER_SO7


def test_exponent_and_me_ratio():
    assert LAMBDA_EXPONENT == 16
    # m_e/M_Pl ≈ 4.19e-23 (Paper 17 Wilson amplitude)
    assert 4.0e-23 < M_E_OVER_M_PL < 4.3e-23


def test_formula():
    expected = M_E_OVER_M_PL**4 * ALPHA_SUBSTRATE**16 * H_COXETER_SO7
    assert lambda_cc() == RHO_LAMBDA_PRED == expected


def test_value_matches_observed_vacuum_energy():
    # ρ_Λ/M_Pl⁴ ≈ 1.2e-123; the famous scale QFT misses by ~10^120
    assert 1.1e-123 < RHO_LAMBDA_PRED < 1.3e-123
    assert abs(RHO_LAMBDA_PRED / RHO_LAMBDA_OBS - 1.0) < 0.01  # within ~1%


def test_summary_keys():
    s = summary()
    assert s["alpha_exponent"] == 16
    assert s["h_cox"] == 6
    assert 0.99 < s["pred_over_obs"] < 1.0
