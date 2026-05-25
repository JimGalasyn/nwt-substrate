"""Tests for nwt_substrate.electroweak.vector_meson_decay_constants.

P7b §7.6 closure (Galasyn 2026-05-23): all 11 vector mesons + B_c
substrate-closed via f_X* · m_X* = 7α · m_τ² · C(X*) with substrate
Casimir multipliers C(X*).
"""

from fractions import Fraction

import pytest

from nwt_substrate.electroweak.vector_meson_decay_constants import (
    ALPHA_SUBSTRATE,
    M_E_GEV,
    VECTOR_MESONS,
    tau_mass_substrate,
    vector_binding_scale,
    vector_meson_fX,
    vector_meson_fX_for,
    precision_chain,
    c_ratio_precision_chain,
    verify_vector_meson_fX,
    precision_chain_summary,
)


# ============================================================
# Substrate τ-mass + binding scale
# ============================================================

def test_tau_mass_matches_pdg():
    """Substrate m_τ = 25 m_e/[α(1-α)²] ≈ 1.776 GeV vs PDG 1.77686."""
    m_tau = tau_mass_substrate()
    assert m_tau == pytest.approx(1.77686, abs=2e-3)


def test_vector_binding_scale_close_to_f_rho_m_rho():
    """S_V = 7α m_τ² ≈ 0.161 GeV² vs PDG f_ρ·m_ρ ≈ 0.167 GeV² → ~-3.3%."""
    S_V = vector_binding_scale()
    pdg_S_V = 0.215 * 0.77526
    rel_err = S_V / pdg_S_V - 1.0
    assert -0.05 < rel_err < 0.0


# ============================================================
# Catalog completeness
# ============================================================

def test_catalog_has_11_states():
    """6 light vectors + 4 heavy vectors + B_c = 11 mesons."""
    assert len(VECTOR_MESONS) == 11


def test_catalog_light_vs_heavy_split():
    """6 light, 5 heavy (D*, D_s*, B*, B_s*, B_c)."""
    light = [s for s in VECTOR_MESONS.values() if not s.is_heavy]
    heavy = [s for s in VECTOR_MESONS.values() if s.is_heavy]
    assert len(light) == 6
    assert len(heavy) == 5


# ============================================================
# Substrate Casimir multipliers C(X*) — substrate vocabulary
# ============================================================

def test_C_rho_is_one():
    """ρ is the reference: C = 1."""
    assert VECTOR_MESONS["rho"].C == Fraction(1, 1)


def test_C_omega_is_7_over_8():
    """ω: C = 7/8 = |K_7|/dim(SU(3))."""
    assert VECTOR_MESONS["omega"].C == Fraction(7, 8)


def test_C_K_star_is_7_over_6():
    """K*: C = 7/6 = |K_7|/(2 C_A(SU(3)))."""
    assert VECTOR_MESONS["K*"].C == Fraction(7, 6)


def test_C_phi_is_10_over_7():
    """φ: C = 10/7 = dim(10_SU(5))/|K_7|."""
    assert VECTOR_MESONS["phi"].C == Fraction(10, 7)


def test_C_J_psi_is_15_over_2():
    """J/ψ: C = 15/2."""
    assert VECTOR_MESONS["J/psi"].C == Fraction(15, 2)


def test_C_Upsilon_is_40():
    """Υ: C = 40 = q_cinq · dim(SU(3)) = 5·8."""
    assert VECTOR_MESONS["Upsilon"].C == Fraction(40, 1)


def test_C_D_star_is_3():
    """D*: C = 3 = C_A(SU(3))."""
    assert VECTOR_MESONS["D*"].C == Fraction(3, 1)


def test_C_Ds_star_is_7_over_2():
    """D_s*: C = 7/2 = |K_7|/C_A(SU(2))."""
    assert VECTOR_MESONS["D_s*"].C == Fraction(7, 2)


def test_C_B_star_is_25_over_4():
    """B*: C = 25/4 = q_cinq²/C_A²(SU(2))."""
    assert VECTOR_MESONS["B*"].C == Fraction(25, 4)


def test_C_Bs_star_is_8():
    """B_s*: C = 8 = dim(SU(3))."""
    assert VECTOR_MESONS["B_s*"].C == Fraction(8, 1)


def test_C_B_c_is_16():
    """B_c: C = 16 = 2 dim(SU(3)) (doubly-heavy gluon Casimir)."""
    assert VECTOR_MESONS["B_c"].C == Fraction(16, 1)


# ============================================================
# Substrate formula numerics
# ============================================================

def test_f_X_formula_factored():
    """f_X* = S_V · C(X*) / m_X*."""
    S_V = vector_binding_scale()
    for name, spec in VECTOR_MESONS.items():
        f_direct = vector_meson_fX_for(name)
        f_expected = S_V * float(spec.C) / spec.m_X_GeV
        assert f_direct == pytest.approx(f_expected, rel=1e-14)


def test_f_rho_substrate_matches_pdg_within_5pct():
    """f_ρ substrate ≈ 208 MeV vs PDG 215 → −3.3 %."""
    f_rho_MeV = vector_meson_fX_for("rho") * 1e3
    rel_err = f_rho_MeV / 215.0 - 1.0
    assert abs(rel_err) < 0.05


def test_f_D_star_within_3pct():
    """f_D* substrate ≈ 240.9 MeV vs PDG 245 → −1.7 %."""
    f_MeV = vector_meson_fX_for("D*") * 1e3
    rel_err = f_MeV / 245.0 - 1.0
    assert abs(rel_err) < 0.03


def test_f_Bc_within_5pct():
    """f_{B_c} substrate ≈ 412 MeV vs lattice 427 → −3.6 %."""
    f_MeV = vector_meson_fX_for("B_c") * 1e3
    rel_err = f_MeV / 427.0 - 1.0
    assert abs(rel_err) < 0.05


# ============================================================
# Precision chain - absolute f_X
# ============================================================

def test_precision_chain_covers_all_11():
    """precision_chain has entries for all 11 mesons."""
    chain = precision_chain()
    chain.pop('pass', None)  # not present in precision_chain itself
    assert set(chain.keys()) == set(VECTOR_MESONS.keys())


def test_all_substrate_predictions_below_pdg():
    """The S_V = 7α m_τ² substrate is below PDG f_ρ·m_ρ by ~3 %, so
    every substrate-predicted f_X* is below PDG."""
    chain = precision_chain()
    for name in VECTOR_MESONS:
        assert chain[name]['rel_err'] < 0.0, (
            f"{name} substrate prediction unexpectedly above PDG"
        )


def test_verify_passes_default_tol():
    """Default 7 % tolerance passes for all 11 states."""
    result = verify_vector_meson_fX()
    assert result['pass'] is True


def test_verify_fails_unrealistic_tol():
    """1 % tolerance fails (S_V is already -3.3 % off PDG)."""
    result = verify_vector_meson_fX(tol=0.01)
    assert result['pass'] is False


# ============================================================
# C-ratio precision chain (cleaner substrate test)
# ============================================================

def test_c_ratio_chain_covers_all_11():
    chain = c_ratio_precision_chain()
    assert set(chain.keys()) == set(VECTOR_MESONS.keys())


def test_c_ratio_rho_is_unity():
    """ρ is the reference: C_substrate = C_empirical = 1 exactly."""
    chain = c_ratio_precision_chain()
    assert chain["rho"]['C_substrate'] == 1.0
    assert chain["rho"]['C_empirical'] == pytest.approx(1.0, rel=1e-14)
    assert chain["rho"]['rel_err'] == pytest.approx(0.0, abs=1e-14)


def test_c_ratio_omega_match_within_1pct():
    """ω substrate C = 7/8 = 0.875 vs empirical ≈ 0.878 → ~-0.35 %."""
    chain = c_ratio_precision_chain()
    assert abs(chain["omega"]['rel_err']) < 0.01


def test_c_ratio_all_below_5pct():
    """All 11 states C-ratio match within 5 % (memo claims 0.2 -- 3.6 %)."""
    chain = c_ratio_precision_chain()
    for name in VECTOR_MESONS:
        assert abs(chain[name]['rel_err']) < 0.05, (
            f"{name} C-ratio gap {chain[name]['rel_err']*100:.2f}% "
            f"exceeds 5 %"
        )


# ============================================================
# Pretty-print
# ============================================================

def test_summary_mentions_key_substrings():
    out = precision_chain_summary()
    assert "S_V" in out
    assert "f_ρ·m_ρ" in out
    assert "f_B_c" in out
    assert "Free params" in out
