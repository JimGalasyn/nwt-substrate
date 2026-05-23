"""Tests for nwt_substrate.electroweak.substrate_gf -- substrate G_F closure.

P7b §7.4 closure (Galasyn 2026-05-23): substrate G_F derived in closed form
from substrate α (Paper 17 K_7 Wilson) and m_e (scale anchor) alone.

  G_F = α⁴ / (625 √2 m_e² (1 + 25α/(4√3))²)

Matches PDG G_F = 1.1663787e-5 GeV⁻² at +55 ppm tree-level.
"""

import math

import pytest

from nwt_substrate.electroweak.substrate_gf import (
    ALPHA_SUBSTRATE,
    M_E_GEV,
    V_EW_PDG_GEV,
    v_ew_substrate,
    fermi_constant_substrate,
    fermi_constant_from_vev,
    precision_chain,
    verify_substrate_gf,
    precision_chain_summary,
)
from nwt_substrate.electroweak.constants import G_F_GEV


# ============================================================
# Substrate α  (Paper 17 K_7 Wilson)
# ============================================================

def test_alpha_substrate_paper17_formula():
    """Substrate α = 1/(25π√3 + 1) -- Paper 17 K_7 Wilson amplitude."""
    expected = 1.0 / (25.0 * math.pi * math.sqrt(3.0) + 1.0)
    assert ALPHA_SUBSTRATE == expected


def test_alpha_substrate_matches_pdg_within_8ppm():
    """Substrate α matches PDG α = 1/137.035999084 at +7.6 ppm."""
    pdg = 1.0 / 137.035999084
    ppm = (ALPHA_SUBSTRATE / pdg - 1.0) * 1e6
    assert 5.0 < ppm < 10.0, f"α ppm gap {ppm:+.2f} outside [5, 10]"


# ============================================================
# Substrate v_EW
# ============================================================

def test_v_ew_substrate_factored_formula():
    """v_EW = (25 m_e/α²) · (1 + 25α/(4√3))."""
    alpha = ALPHA_SUBSTRATE
    m_e = M_E_GEV
    sqrt3 = math.sqrt(3.0)
    expected = (25.0 * m_e / alpha ** 2) * (1.0 + 25.0 * alpha / (4.0 * sqrt3))
    assert v_ew_substrate(alpha, m_e) == pytest.approx(expected, rel=1e-14)


def test_v_ew_substrate_matches_pdg_within_30ppm():
    """Substrate v_EW matches PDG v_EW = 246.21965 GeV at ~28 ppm."""
    ppm = (v_ew_substrate() / V_EW_PDG_GEV - 1.0) * 1e6
    assert abs(ppm) < 30.0, f"v_EW ppm gap {ppm:+.2f} exceeds 30"


# ============================================================
# Substrate G_F (closed form)
# ============================================================

def test_g_f_substrate_matches_pdg_within_100ppm():
    """Substrate G_F matches PDG at +55 ppm tree-level."""
    ppm = (fermi_constant_substrate() / G_F_GEV - 1.0) * 1e6
    assert abs(ppm) < 100.0, f"G_F ppm gap {ppm:+.2f} exceeds 100"


def test_g_f_substrate_specific_55ppm():
    """The +55 ppm gap is the substrate-α gap (+7.6 ppm in α) propagated
    through v_EW (∝ 1/α²) → G_F (∝ 1/v² = α⁴) at the closed form."""
    ppm = (fermi_constant_substrate() / G_F_GEV - 1.0) * 1e6
    assert 50.0 < ppm < 60.0, f"G_F ppm gap {ppm:+.2f} outside [50, 60]"


def test_closed_form_equals_tree_level_identity():
    """G_F = α⁴/(625√2 m_e² NLO²) equals 1/(√2 v_EW²) identically."""
    v_sub = v_ew_substrate()
    g_f_closed = fermi_constant_substrate()
    g_f_tree = fermi_constant_from_vev(v_sub)
    assert g_f_closed == pytest.approx(g_f_tree, rel=1e-14)


def test_g_f_inverse_v_squared_scaling():
    """G_F ∝ v_EW⁻² (factor 2 in α perturbation propagates ppm-wise)."""
    # Perturb α by epsilon, check G_F changes by factor 2 in opposite sign
    alpha = ALPHA_SUBSTRATE
    eps = 1e-5
    v0 = v_ew_substrate(alpha)
    v1 = v_ew_substrate(alpha * (1 + eps))
    gf0 = fermi_constant_substrate(alpha)
    gf1 = fermi_constant_substrate(alpha * (1 + eps))
    # Relative changes
    dv_over_v = (v1 - v0) / v0
    dgf_over_gf = (gf1 - gf0) / gf0
    # G_F = 1/(√2 v²) so d ln G_F = -2 d ln v
    assert dgf_over_gf == pytest.approx(-2.0 * dv_over_v, rel=1e-4)


# ============================================================
# Precision chain
# ============================================================

def test_precision_chain_keys():
    """precision_chain returns dict with alpha, v_ew, g_f sub-dicts."""
    chain = precision_chain()
    assert set(chain.keys()) == {'alpha', 'v_ew', 'g_f'}
    for key in chain:
        assert set(chain[key].keys()) == {'substrate', 'pdg', 'ppm'}


def test_precision_chain_ppm_signs():
    """Specific ppm signs: α positive (+7.6), v_EW negative (-28), G_F positive (+55)."""
    chain = precision_chain()
    assert chain['alpha']['ppm'] > 0      # substrate α > PDG α
    assert chain['v_ew']['ppm'] < 0       # substrate v_EW < PDG v_EW
    assert chain['g_f']['ppm'] > 0        # substrate G_F > PDG G_F


def test_precision_chain_ppm_magnitudes():
    """Specific ppm magnitudes match P7b §7.4 documentation."""
    chain = precision_chain()
    assert 5.0 < chain['alpha']['ppm'] < 10.0       # +7.6 ppm
    assert -35.0 < chain['v_ew']['ppm'] < -20.0     # -27.7 ppm
    assert 50.0 < chain['g_f']['ppm'] < 60.0        # +55 ppm


# ============================================================
# Verification
# ============================================================

def test_verify_substrate_gf_passes_default_tol():
    """Default 100 ppm tolerance passes for substrate G_F (~55 ppm gap)."""
    result = verify_substrate_gf()
    assert result['pass'] is True


def test_verify_substrate_gf_fails_unrealistic_tol():
    """A 10 ppm tolerance fails (substrate G_F is +55 ppm)."""
    result = verify_substrate_gf(ppm_tol=10.0)
    assert result['pass'] is False


# ============================================================
# Pretty-print
# ============================================================

def test_precision_chain_summary_contains_key_substrings():
    """Summary string mentions substrate G_F, PDG comparison, and closed form."""
    out = precision_chain_summary()
    assert "G_F" in out
    assert "v_EW" in out
    assert "α" in out
    assert "Closed form" in out
    assert "0" in out  # zero free params line


# ============================================================
# Integration with substrate namespace
# ============================================================

def test_substrate_namespace_exposes_g_f():
    """The substrate namespace at ew.substrate.fermi_constant returns G_F."""
    import nwt_substrate.electroweak as ew
    g_f = ew.substrate.fermi_constant()
    assert g_f == pytest.approx(fermi_constant_substrate(), rel=1e-14)


def test_substrate_namespace_verify_g_f():
    """ew.substrate.verify_g_f() returns dict with pass=True."""
    import nwt_substrate.electroweak as ew
    result = ew.substrate.verify_g_f()
    assert result['pass'] is True
