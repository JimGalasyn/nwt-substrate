"""Tests for nwt_substrate.cosmology.omega_b_c."""

from nwt_substrate.cosmology.omega_b_c import (
    Q_CINQ_SQ,
    OMEGA_B_C_PRED,
    OMEGA_B_C_PLANCK,
    omega_b_c,
    summary,
)
from nwt_substrate.isa.constants import ALPHA_SUBSTRATE, H_V_SO7, RANK_SO7


def test_lo_coeff_is_h_v_squared():
    assert Q_CINQ_SQ == 25 == H_V_SO7**2


def test_formula_is_25alpha_plus_75alpha2():
    expected = 25 * ALPHA_SUBSTRATE + 75 * ALPHA_SUBSTRATE**2
    assert omega_b_c() == OMEGA_B_C_PRED == expected
    # NLO coefficient = RANK_SO7 * q^2 = 3*25 = 75
    assert RANK_SO7 * Q_CINQ_SQ == 75


def test_value_near_planck():
    assert 0.185 < OMEGA_B_C_PRED < 0.187
    assert 0.185 < OMEGA_B_C_PLANCK < 0.187


def test_matches_planck_within_100ppm():
    # actual offset ~67 ppm with the substrate α
    assert abs(OMEGA_B_C_PRED / OMEGA_B_C_PLANCK - 1.0) < 1e-4


def test_summary_keys():
    s = summary()
    assert s["lo_coeff"] == 25
    assert s["nlo_coeff"] == 75
    assert abs(s["ppm_offset"]) < 100
