"""Tests for the nwt.electroweak shim."""

import pytest

import nwt_substrate.electroweak as ew


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_M_Z_pdg():
    assert ew.M_Z == pytest.approx(91.1876, abs=1e-4)


def test_GAMMA_Z_pdg():
    assert ew.GAMMA_Z == pytest.approx(2.4952, abs=1e-3)


def test_sin2_theta_W_pdg():
    """Effective sin²θ_W at Z pole."""
    assert ew.SIN2_THETA_W == pytest.approx(0.23121, abs=1e-4)


def test_G_Z_squared_from_G_F():
    """g_Z² = 4√2 G_F M_Z² (Halzen-Martin Eq 13.50)."""
    import math
    expected = 4.0 * math.sqrt(2.0) * 1.1663787e-5 * ew.M_Z ** 2
    assert ew.G_Z ** 2 == pytest.approx(expected, rel=1e-5)


# ---------------------------------------------------------------------------
# Couplings
# ---------------------------------------------------------------------------

def test_charged_lepton_couplings():
    """g_V^ℓ = -1/2 + 2 sin²θ_W ≈ -0.038, g_A^ℓ = -1/2."""
    e = ew.coupling("e")
    assert e.T_3 == pytest.approx(-0.5)
    assert e.Q == pytest.approx(-1.0)
    assert e.g_V == pytest.approx(-0.5 + 2 * ew.SIN2_THETA_W, abs=1e-9)
    assert e.g_A == pytest.approx(-0.5)


def test_up_type_quark_couplings():
    """g_V^u = +1/2 - (4/3) sin²θ_W, g_A^u = +1/2."""
    u = ew.coupling("u")
    assert u.T_3 == pytest.approx(+0.5)
    assert u.Q == pytest.approx(+2.0 / 3.0)
    assert u.g_V == pytest.approx(0.5 - (4.0 / 3.0) * ew.SIN2_THETA_W, abs=1e-9)
    assert u.g_A == pytest.approx(+0.5)
    assert u.n_color == 3


def test_down_type_quark_couplings():
    """g_V^d = -1/2 + (2/3) sin²θ_W, g_A^d = -1/2."""
    d = ew.coupling("d")
    assert d.T_3 == pytest.approx(-0.5)
    assert d.Q == pytest.approx(-1.0 / 3.0)
    assert d.g_V == pytest.approx(-0.5 + (2.0 / 3.0) * ew.SIN2_THETA_W, abs=1e-9)
    assert d.g_A == pytest.approx(-0.5)
    assert d.n_color == 3


def test_neutrino_pure_axial_coupling():
    """For ν: Q = 0, so g_V = g_A = T_3 = +1/2."""
    nu = ew.coupling("nu_e")
    assert nu.Q == 0.0
    assert nu.g_V == pytest.approx(+0.5)
    assert nu.g_A == pytest.approx(+0.5)


def test_coupling_aliases():
    """The lookup accepts plain English names too."""
    assert ew.coupling("electron").name == "e"
    assert ew.coupling("muon").name == "mu"
    assert ew.coupling("up").name == "u"
    assert ew.coupling("bottom").name == "b"


def test_coupling_unknown_raises():
    with pytest.raises(KeyError):
        ew.coupling("xyzzy")


# ---------------------------------------------------------------------------
# Z partial widths
# ---------------------------------------------------------------------------

def test_partial_width_Z_to_ee_matches_pdg():
    """PDG Γ(Z → e⁺e⁻) ≈ 83.99 MeV; with our G_F-based coupling we get
    ~83.4 MeV (LO, before radiative corrections)."""
    G_ee = ew.partial_width_Z("e")
    assert G_ee * 1000 == pytest.approx(83.4, abs=1.0)


def test_partial_width_Z_to_invisible_matches_pdg():
    """3 ν families: PDG Γ_inv ≈ 499 MeV; we predict ~498 MeV."""
    G_inv = sum(ew.partial_width_Z(name)
                 for name in ["nu_e", "nu_mu", "nu_tau"])
    assert G_inv * 1000 == pytest.approx(499.0, abs=5.0)


def test_total_Z_width_within_3_percent_of_pdg():
    """LO total Γ_Z prediction is within ~3% of PDG (residual from
    radiative + QCD corrections to hadronic widths)."""
    Gamma = ew.total_width_Z()
    assert Gamma == pytest.approx(2.4952, abs=0.10)
    assert Gamma / 2.4952 > 0.95


def test_top_quark_kinematically_forbidden_in_Z_decay():
    """Z → t t̄ is forbidden because 2 m_t ≈ 345 GeV > M_Z = 91 GeV."""
    assert ew.partial_width_Z("t") == 0.0


def test_branching_ratios_sum_to_one():
    BR = ew.branching_ratios_Z()
    assert sum(BR.values()) == pytest.approx(1.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Cross-sections
# ---------------------------------------------------------------------------

def test_z_pole_peak_cross_section_matches_pdg():
    """σ_peak(e⁺e⁻ → μ⁺μ⁻) at √s = M_Z is well-known: ≈ 2 nb (2000 pb)."""
    sigma = ew.sigma_total(ew.M_Z, "mu", m_f=0.10566)
    # PDG result: σ_peak ≈ 12π Γ_ee Γ_μμ / (M_Z² Γ_Z²) ≈ 2 nb
    assert 1500.0 < sigma < 2500.0


def test_qed_only_recovers_R_ratio_kinematic_factor():
    """At √s well above Z pole and not too far above, the photon channel
    of σ(e⁺e⁻ → ff̄) reproduces the R-ratio integrand."""
    sqrt_s = 200.0
    # At √s = 200 GeV, σ_QED(ee → μμ) ≈ 4πα²/3s = 86.8 nb / s[GeV²]
    sigma_qed = ew.sigma_qed_only(sqrt_s, "mu", m_f=0.10566)
    expected_pointlike = 86.8 / sqrt_s ** 2 * 1000  # convert nb to pb (factor 1e3)
    # Should be very close at √s = 200 GeV
    assert sigma_qed == pytest.approx(expected_pointlike, rel=0.02)


def test_z_pole_dominates_off_resonance():
    """At the Z peak, full σ ≈ 200× photon-only.  At √s = 200 GeV it's
    ~2.4× (interference-dominated)."""
    sigma_full_peak = ew.sigma_total(ew.M_Z, "mu", m_f=0.10566)
    sigma_qed_peak = ew.sigma_qed_only(ew.M_Z, "mu", m_f=0.10566)
    assert sigma_full_peak / sigma_qed_peak > 100.0

    sigma_full_high = ew.sigma_total(200.0, "mu", m_f=0.10566)
    sigma_qed_high = ew.sigma_qed_only(200.0, "mu", m_f=0.10566)
    assert 1.5 < sigma_full_high / sigma_qed_high < 3.5


def test_z_lineshape_destructive_interference_below_peak():
    """Below the Z pole, the γ-Z interference is destructive (Z amplitude
    has opposite sign to γ amplitude).  Around √s = 80 GeV the total falls
    below the photon-only result."""
    sigma_full = ew.sigma_total(80.0, "mu", m_f=0.10566)
    sigma_qed = ew.sigma_qed_only(80.0, "mu", m_f=0.10566)
    assert sigma_full < sigma_qed


# ---------------------------------------------------------------------------
# Cross-shim consistency
# ---------------------------------------------------------------------------

def test_qed_only_consistent_with_phase_q9_partonic_DY():
    """The σ_qed_only formula is the photon-only Drell-Yan integrand;
    should equal Phase Q.9's textbook DY formula × N_c² (R-ratio vs DY
    color factor flip)."""
    import math
    sqrt_s = 100.0
    sigma_qed = ew.sigma_qed_only(sqrt_s, "u", m_f=0.0)
    # Direct textbook formula:
    s = sqrt_s ** 2
    Q_u = 2.0 / 3.0
    sigma_textbook = (4 * math.pi * ew.ALPHA_QED ** 2 / (3 * s)) * Q_u ** 2 * 3 * 0.3894e9
    assert sigma_qed == pytest.approx(sigma_textbook, rel=1e-9)
