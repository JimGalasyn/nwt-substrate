"""Tests for nwt_substrate.atomic.hydrogen -- substrate hydrogen spectroscopy.

P8 substrate Coulomb closure (Galasyn 2026-05-23): substrate V(r) = -α ℏc/r
gives the full hydrogen Bohr-Rydberg spectrum at ppm-level PDG match via
substrate α (Paper 17 K_7 Wilson) alone.
"""

import math

import pytest

from nwt_substrate.atomic.hydrogen import (
    M_E_EV, M_P_EV, HBAR_C_EV_PM, G_P,
    rydberg,
    bohr_radius,
    hydrogen_R_H,
    lyman_alpha,
    hyperfine_21cm,
    fine_structure_2P,
    electron_a_e_one_loop,
    lamb_shift_scale,
    precision_chain,
    verify_substrate_hydrogen,
    precision_chain_summary,
)
from nwt_substrate.electroweak.substrate_gf import ALPHA_SUBSTRATE


# ============================================================
# Rydberg (m_e c² α² / 2)
# ============================================================

def test_rydberg_formula():
    """R_y = m_e c² α² / 2."""
    alpha = ALPHA_SUBSTRATE
    expected = 0.5 * M_E_EV * alpha ** 2
    assert rydberg() == pytest.approx(expected, rel=1e-14)


def test_rydberg_matches_pdg_within_20ppm():
    """Substrate Rydberg matches PDG 13.6056931 eV at +15 ppm."""
    R_y = rydberg()
    ppm = (R_y / 13.6056931 - 1.0) * 1e6
    assert 10.0 < ppm < 20.0, f"Rydberg ppm gap {ppm:+.2f} outside [10, 20]"


# ============================================================
# Bohr radius (ℏc / (m_e c² α))
# ============================================================

def test_bohr_radius_formula():
    """a_0 = ℏc / (m_e c² α)."""
    alpha = ALPHA_SUBSTRATE
    expected = HBAR_C_EV_PM / (M_E_EV * alpha)
    assert bohr_radius() == pytest.approx(expected, rel=1e-14)


def test_bohr_radius_matches_pdg_within_20ppm():
    """Substrate Bohr radius matches PDG 52.91772 pm at -8 ppm."""
    a_0 = bohr_radius()
    ppm = (a_0 / 52.91772 - 1.0) * 1e6
    assert -15.0 < ppm < -5.0, f"Bohr radius ppm gap {ppm:+.2f} outside [-15, -5]"


def test_bohr_radius_product_with_rydberg():
    """R_y · a_0 = ℏc · α / 2  (dimensionless ratio R_y a_0 / ℏc = α/2)."""
    R_y = rydberg()
    a_0 = bohr_radius()
    ratio = R_y * a_0 / HBAR_C_EV_PM
    assert ratio == pytest.approx(ALPHA_SUBSTRATE / 2.0, rel=1e-14)


# ============================================================
# Hydrogen R_H (reduced-mass)
# ============================================================

def test_hydrogen_R_H_reduced_mass():
    """R_H = R_y · M_p / (M_p + m_e)."""
    R_y = rydberg()
    expected = R_y * M_P_EV / (M_P_EV + M_E_EV)
    assert hydrogen_R_H() == pytest.approx(expected, rel=1e-14)


def test_hydrogen_R_H_matches_pdg_within_20ppm():
    """Substrate R_H matches PDG 13.5982886 eV at +15 ppm."""
    R_H = hydrogen_R_H()
    ppm = (R_H / 13.5982886 - 1.0) * 1e6
    assert 10.0 < ppm < 20.0


def test_hydrogen_R_H_less_than_rydberg():
    """R_H < R_y by the reduced-mass factor ~1 - m_e/M_p."""
    assert hydrogen_R_H() < rydberg()
    ratio = hydrogen_R_H() / rydberg()
    assert abs(ratio - M_P_EV / (M_P_EV + M_E_EV)) < 1e-14


# ============================================================
# Lyman-α (3/4 R_H)
# ============================================================

def test_lyman_alpha_formula():
    """Lyman-α = (3/4) R_H = R_H · (1 - 1/4)."""
    assert lyman_alpha() == pytest.approx(0.75 * hydrogen_R_H(), rel=1e-14)


def test_lyman_alpha_matches_pdg_within_1pct():
    """Substrate Lyman-α matches PDG 10.20486 eV at -587 ppm."""
    L_a = lyman_alpha()
    ppm = (L_a / 10.20486 - 1.0) * 1e6
    assert -700.0 < ppm < -400.0


# ============================================================
# 21-cm hyperfine
# ============================================================

def test_hyperfine_21cm_formula():
    """ΔE_HF(1S) = (4/3) g_p (m_e/M_p) α⁴ m_e c²."""
    alpha = ALPHA_SUBSTRATE
    expected = (4.0/3.0) * G_P * (M_E_EV / M_P_EV) * alpha**4 * M_E_EV
    assert hyperfine_21cm() == pytest.approx(expected, rel=1e-14)


def test_hyperfine_21cm_matches_pdg_within_1pct():
    """Substrate 21-cm matches PDG 5.87433 μeV at ~0.1%."""
    hf = hyperfine_21cm()
    pct = (hf / 5.87433e-6 - 1.0) * 100
    assert -1.0 < pct < 1.0, f"21-cm pct gap {pct:+.3f} outside ±1%"


def test_hyperfine_21cm_microeV_scale():
    """Hyperfine is at the μeV scale (∼6 × 10⁻⁶ eV)."""
    hf = hyperfine_21cm()
    assert 5e-6 < hf < 7e-6


# ============================================================
# 2P fine structure (m_e c² α⁴ / 32)
# ============================================================

def test_fine_structure_2P_formula():
    """ΔE_FS(2P) = m_e c² α⁴ / 32."""
    alpha = ALPHA_SUBSTRATE
    expected = M_E_EV * alpha ** 4 / 32.0
    assert fine_structure_2P() == pytest.approx(expected, rel=1e-14)


def test_fine_structure_2P_matches_pdg_within_1pct():
    """Substrate 2P FS matches PDG 45.354 μeV at -0.15%."""
    fs = fine_structure_2P()
    pct = (fs / 45.354e-6 - 1.0) * 100
    assert -0.5 < pct < 0.5, f"FS pct gap {pct:+.3f} outside ±0.5%"


# ============================================================
# Electron a_e (Schwinger 1-loop)
# ============================================================

def test_electron_a_e_one_loop_formula():
    """a_e^{(1)} = α / (2π) (Schwinger 1948)."""
    alpha = ALPHA_SUBSTRATE
    expected = alpha / (2.0 * math.pi)
    assert electron_a_e_one_loop() == pytest.approx(expected, rel=1e-14)


def test_electron_a_e_one_loop_matches_pdg_within_15ppm():
    """Substrate a_e^(1) matches PDG 1-loop part 0.00116141 at +8 ppm."""
    a_e = electron_a_e_one_loop()
    ppm = (a_e / 0.00116141 - 1.0) * 1e6
    assert 5.0 < ppm < 15.0, f"a_e ppm gap {ppm:+.2f} outside [5, 15]"


# ============================================================
# Lamb shift scale
# ============================================================

def test_lamb_shift_scale_formula():
    """α⁵ m_e c² scale (without Bethe-log factor)."""
    alpha = ALPHA_SUBSTRATE
    expected = alpha ** 5 * M_E_EV
    assert lamb_shift_scale() == pytest.approx(expected, rel=1e-14)


def test_lamb_shift_scale_microeV():
    """α⁵ m_e c² ≈ 10 μeV (substrate prediction; PDG Lamb shift 4.37 μeV)."""
    L_scale = lamb_shift_scale()
    assert 5e-6 < L_scale < 15e-6


# ============================================================
# Precision chain
# ============================================================

def test_precision_chain_keys():
    """precision_chain returns dict with all 7 expected observables."""
    chain = precision_chain()
    expected_keys = {
        'rydberg', 'bohr_radius', 'hydrogen_R_H', 'hyperfine_21cm',
        'fine_structure_2P', 'a_e_one_loop', 'lyman_alpha',
    }
    assert set(chain.keys()) == expected_keys


def test_precision_chain_substrate_pdg_ppm_structure():
    """Each observable maps to {substrate, pdg, ppm}."""
    chain = precision_chain()
    for name, data in chain.items():
        assert set(data.keys()) == {'substrate', 'pdg', 'ppm'}


def test_precision_chain_ppm_signs():
    """Verify the documented ppm-sign pattern."""
    chain = precision_chain()
    assert chain['rydberg']['ppm'] > 0           # +15 ppm
    assert chain['bohr_radius']['ppm'] < 0       # -8 ppm
    assert chain['hydrogen_R_H']['ppm'] > 0      # +15 ppm
    assert chain['a_e_one_loop']['ppm'] > 0      # +8 ppm
    # 21-cm and FS are at %-level (much larger than ppm)
    assert chain['hyperfine_21cm']['ppm'] > 0    # +0.1% = +1000 ppm
    assert chain['fine_structure_2P']['ppm'] < 0 # -0.15%


# ============================================================
# Verification
# ============================================================

def test_verify_substrate_hydrogen_passes_default():
    """All hydrogen observables pass at default tolerances (100 ppm / 0.5%)."""
    result = verify_substrate_hydrogen()
    assert result['pass'] is True
    for name, passed in result['per_observable_pass'].items():
        assert passed, f"Observable {name} failed default tolerance check"


def test_verify_substrate_hydrogen_fails_unrealistic_ppm():
    """A 1 ppm tolerance fails (Rydberg at +15 ppm doesn't fit)."""
    result = verify_substrate_hydrogen(ppm_tol_ppm=1.0)
    assert result['pass'] is False
    assert result['per_observable_pass']['rydberg'] is False


# ============================================================
# Pretty-print
# ============================================================

def test_precision_chain_summary_contains_key_substrings():
    out = precision_chain_summary()
    assert "Rydberg" in out
    assert "Bohr" in out
    assert "21-cm" in out
    assert "fine structure" in out.lower()
    assert "Coulomb" in out
    assert "0" in out  # zero free params


# ============================================================
# Module + substrate α integration
# ============================================================

def test_uses_same_alpha_as_ew_substrate_gf():
    """The atomic shim uses the same ALPHA_SUBSTRATE as electroweak.substrate_gf."""
    # rydberg with default should equal rydberg with ALPHA_SUBSTRATE passed in
    assert rydberg() == rydberg(alpha=ALPHA_SUBSTRATE)


def test_atomic_module_exports_all_functions():
    """The atomic submodule re-exports all hydrogen formulas."""
    import nwt_substrate.atomic as at
    assert hasattr(at, 'rydberg')
    assert hasattr(at, 'bohr_radius')
    assert hasattr(at, 'hyperfine_21cm')
    assert hasattr(at, 'verify_substrate_hydrogen')
    assert hasattr(at, 'precision_chain')
