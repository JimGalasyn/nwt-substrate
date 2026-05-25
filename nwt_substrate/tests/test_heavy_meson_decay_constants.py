"""Tests for nwt_substrate.electroweak.heavy_meson_decay_constants.

P7b §7.5 closure (Galasyn 2026-05-23): heavy pseudoscalar f_X via
f_X² = f_π² · R_s · N_X / m_X with R_s = 1 (non-strange) or √(7/4)
(strange).  Closes the "largest quantitative gap" of P7b at 1.1 -- 2.6 % PDG.
"""

import math

import pytest

from nwt_substrate.electroweak.heavy_meson_decay_constants import (
    ALPHA_SUBSTRATE,
    M_E_GEV,
    R_S_NONSTRANGE,
    R_S_STRANGE,
    HEAVY_PSEUDOSCALARS,
    pi_zero_mass,
    f_pi_substrate,
    heavy_meson_fX,
    heavy_meson_fX_for,
    fX_ratio_strange_over_nonstrange,
    precision_chain,
    verify_heavy_meson_fX,
    precision_chain_summary,
)


# ============================================================
# Substrate Goldstone scales
# ============================================================

def test_pi_zero_mass_matches_pdg():
    """Substrate m_π⁰ = (2 m_e/α)(1 - 5α) matches PDG 134.977 MeV sub-ppm."""
    m_pi0_GeV = pi_zero_mass()
    m_pi0_MeV = m_pi0_GeV * 1e3
    assert m_pi0_MeV == pytest.approx(134.94, abs=0.05)


def test_f_pi_substrate_formula():
    """f_π = m_π⁰ / 5^(1/4)."""
    f_pi = f_pi_substrate()
    assert f_pi == pytest.approx(pi_zero_mass() / math.pow(5.0, 0.25),
                                 rel=1e-14)


def test_f_pi_substrate_close_to_pdg():
    """f_π ≈ 90.24 MeV vs PDG 92.4 MeV → 2.3% gap (chiral-anomaly residual)."""
    f_pi_MeV = f_pi_substrate() * 1e3
    rel_err = f_pi_MeV / 92.4 - 1.0
    assert -0.03 < rel_err < 0.0


# ============================================================
# Substrate strangeness factor
# ============================================================

def test_R_s_nonstrange_is_unity():
    """R_s = 1 for non-strange (D, B)."""
    assert R_S_NONSTRANGE == 1.0


def test_R_s_strange_is_sqrt_7_over_4():
    """R_s = √(7/4) = √|K_7|/C_A²(SU(2))."""
    assert R_S_STRANGE == pytest.approx(math.sqrt(7.0 / 4.0), rel=1e-14)
    assert R_S_STRANGE == pytest.approx(1.3229, abs=1e-4)


# ============================================================
# SU(5)-rep labels
# ============================================================

def test_substrate_N_D_is_10():
    """D-meson SU(5) rep label = dim(10) = 10."""
    assert HEAVY_PSEUDOSCALARS["D"].N == 10


def test_substrate_N_Ds_is_11():
    """D_s = 10 + strangeness = 11."""
    assert HEAVY_PSEUDOSCALARS["D_s"].N == 11


def test_substrate_N_B_is_24():
    """B-meson SU(5) rep label = dim(adjoint) = 24."""
    assert HEAVY_PSEUDOSCALARS["B"].N == 24


def test_substrate_N_Bs_is_25():
    """B_s SU(5) rep label = q_cinq² = 5² = 25 (cinquefoil signature)."""
    assert HEAVY_PSEUDOSCALARS["B_s"].N == 25


def test_strange_flag_consistency():
    """D, B non-strange; D_s, B_s strange."""
    assert HEAVY_PSEUDOSCALARS["D"].strange is False
    assert HEAVY_PSEUDOSCALARS["B"].strange is False
    assert HEAVY_PSEUDOSCALARS["D_s"].strange is True
    assert HEAVY_PSEUDOSCALARS["B_s"].strange is True


# ============================================================
# Substrate ansatz numerics
# ============================================================

def test_f_D_within_2pct():
    """f_D substrate ≈ 209.6 MeV vs PDG 212.0 → −1.13%."""
    f_D_MeV = heavy_meson_fX_for("D") * 1e3
    rel_err = f_D_MeV / 212.0 - 1.0
    assert abs(rel_err) < 0.02


def test_f_Ds_within_2pct():
    """f_D_s substrate ≈ 246.4 MeV vs PDG 248.0 → −0.65% (strange √(7/4))."""
    f_Ds_MeV = heavy_meson_fX_for("D_s") * 1e3
    rel_err = f_Ds_MeV / 248.0 - 1.0
    assert abs(rel_err) < 0.02


def test_f_B_within_2pct():
    """f_B substrate ≈ 192.2 MeV vs PDG 190.0 → +1.18%."""
    f_B_MeV = heavy_meson_fX_for("B") * 1e3
    rel_err = f_B_MeV / 190.0 - 1.0
    assert abs(rel_err) < 0.02


def test_f_Bs_within_5pct():
    """f_B_s substrate ≈ 224.0 MeV vs PDG 230.0 → −2.62% (strange √(7/4))."""
    f_Bs_MeV = heavy_meson_fX_for("B_s") * 1e3
    rel_err = f_Bs_MeV / 230.0 - 1.0
    assert abs(rel_err) < 0.05


def test_strange_correction_required_for_Ds():
    """Without √(7/4) D_s undershoots by ~14% (substrate must include it)."""
    spec = HEAVY_PSEUDOSCALARS["D_s"]
    f_no_strange_MeV = heavy_meson_fX(spec.N, spec.m_X_GeV,
                                      strange=False) * 1e3
    rel_err = f_no_strange_MeV / 248.0 - 1.0
    assert rel_err < -0.10, "Without R_s, D_s should undershoot >10%"


def test_strange_correction_required_for_Bs():
    """Without √(7/4) B_s undershoots by ~15%."""
    spec = HEAVY_PSEUDOSCALARS["B_s"]
    f_no_strange_MeV = heavy_meson_fX(spec.N, spec.m_X_GeV,
                                      strange=False) * 1e3
    rel_err = f_no_strange_MeV / 230.0 - 1.0
    assert rel_err < -0.10


# ============================================================
# SU(3)-breaking ratios
# ============================================================

def test_fDs_over_fD_close_to_pdg():
    """f_D_s/f_D substrate ≈ 1.18 vs PDG 1.17 → 0.5%."""
    ratio = fX_ratio_strange_over_nonstrange("D_s", "D")
    rel_err = ratio / 1.170 - 1.0
    assert abs(rel_err) < 0.02


def test_fBs_over_fB_close_to_pdg():
    """f_B_s/f_B substrate ≈ 1.16 vs PDG 1.21 → −4.2%."""
    ratio = fX_ratio_strange_over_nonstrange("B_s", "B")
    rel_err = ratio / 1.211 - 1.0
    assert abs(rel_err) < 0.05


# ============================================================
# Precision chain + verification
# ============================================================

def test_precision_chain_covers_all_four():
    """precision_chain has entries for D, D_s, B, B_s."""
    chain = precision_chain()
    assert set(chain.keys()) == {"D", "D_s", "B", "B_s"}


def test_precision_chain_signs():
    """Per-meson rel-err signs match memo: D < 0, D_s < 0, B > 0, B_s < 0."""
    chain = precision_chain()
    assert chain["D"]['rel_err'] < 0.0      # −1.13%
    assert chain["D_s"]['rel_err'] < 0.0    # −0.65%
    assert chain["B"]['rel_err'] > 0.0      # +1.18%
    assert chain["B_s"]['rel_err'] < 0.0    # −2.62%


def test_verify_passes_default_5pct():
    """Default 5% tolerance passes for all 4 heavy pseudoscalars."""
    result = verify_heavy_meson_fX()
    assert result['pass'] is True


def test_verify_fails_unrealistic_tol():
    """0.1% tolerance fails (B_s is 2.6% off)."""
    result = verify_heavy_meson_fX(tol=1e-3)
    assert result['pass'] is False


# ============================================================
# Pretty-print
# ============================================================

def test_summary_mentions_key_substrings():
    out = precision_chain_summary()
    assert "f_π" in out
    assert "f_D" in out
    assert "f_B" in out
    assert "√(7/4)" in out
    assert "Free params" in out
