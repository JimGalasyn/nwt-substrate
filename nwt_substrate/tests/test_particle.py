"""Tests for nwt_substrate.particles.particle (Particle class)."""

import pytest

from nwt_substrate.particles.particle import Particle


def test_construct_minimal_particle():
    p = Particle(name="X", p=2, q=3, m=5, n_q=3)
    assert p.p == 2 and p.q == 3 and p.m == 5 and p.n_q == 3
    assert p.f == 0
    assert p.sector == "spinor"
    assert p.J == 0.5


def test_sector_to_spin_mapping():
    """scalar -> 0, vector -> 1, spinor -> 1/2."""
    s = Particle(name="s", p=1, q=1, m=2, n_q=0, sector="scalar")
    v = Particle(name="v", p=1, q=1, m=2, n_q=0, sector="vector")
    f = Particle(name="f", p=1, q=1, m=2, n_q=0, sector="spinor")
    assert s.J == 0.0
    assert v.J == 1.0
    assert f.J == 0.5


def test_isospin_from_framing():
    """I_3 = f / 2."""
    p_plus = Particle(name="+", p=1, q=1, m=2, n_q=0, f=2)
    p_zero = Particle(name="0", p=1, q=1, m=2, n_q=0, f=0)
    p_minus = Particle(name="-", p=1, q=1, m=2, n_q=0, f=-2)
    assert p_plus.I3 == 1.0
    assert p_zero.I3 == 0.0
    assert p_minus.I3 == -1.0


def test_hopf_charge_property():
    """Q_H = p * m."""
    p = Particle(name="X", p=3, q=4, m=7, n_q=3)
    assert p.Q_H == 21


def test_baryon_substrate_n_q_mod_2():
    """Substrate baryon number rule: B = n_q mod 2."""
    meson = Particle(name="m", p=2, q=3, m=5, n_q=2)
    baryon = Particle(name="b", p=1, q=3, m=5, n_q=5)
    assert meson.B_substrate == 0
    assert baryon.B_substrate == 1


def test_carrier_property_uses_canonical_map():
    p = Particle(name="X", p=1, q=1, m=2, n_q=3)
    assert p.carrier == "trefoil"


def test_torus_knot_family_classification():
    p = Particle(name="X", p=2, q=3, m=4, n_q=3)
    assert p.torus_knot_family == "trefoil"
    assert p.genus == 1
    assert p.crossings == 3


def test_mass_pred_returns_none_for_unphysical():
    """If m^2 <= p^2, beta is imaginary and mass should be None."""
    p = Particle(name="bad", p=3, q=2, m=2, n_q=0)
    assert p.mass_pred is None


def test_charge_via_extended_gell_mann_nishijima():
    """Q = I_3 + (B + S + C) / 2; verify a few canonical hadrons."""
    # K+: I_3 = 1/2, B = 0, S = 1, C = 0  -> Q = 1/2 + 1/2 = 1
    K_plus = Particle(name="K+", p=2, q=5, m=8, n_q=2, f=1,
                      B_pdg=0, S_pdg=1, C_pdg=0)
    assert K_plus.Q_pred == pytest.approx(1.0)
    # Sigma-: I_3 = -1, B = 1, S = -1, C = 0  -> Q = -1 + 0 = -1
    sig_minus = Particle(name="Sigma-", p=1, q=3, m=6, n_q=5, f=-2,
                         B_pdg=1, S_pdg=-1, C_pdg=0)
    assert sig_minus.Q_pred == pytest.approx(-1.0)
    # Lambda: I_3 = 0, B = 1, S = -1, C = 0  -> Q = -0.5 + 0.5 = 0
    lam = Particle(name="Lambda", p=3, q=4, m=12, n_q=3, f=0,
                   B_pdg=1, S_pdg=-1, C_pdg=0)
    assert lam.Q_pred == pytest.approx(0.0)


def test_summary_runs_without_error():
    """summary() should produce a multi-line string."""
    p = Particle(name="p", p=1, q=3, m=5, n_q=5, f=1, sector="spinor",
                 m_obs=938.27, Q_pdg=1, B_pdg=1, S_pdg=0, C_pdg=0)
    s = p.summary()
    assert "Particle: p" in s
    assert "Mass predicted" in s
    assert "Mass observed" in s


def test_isospin_ceiling():
    """I <= n_q/2 topological ceiling (WRT D17); saturated by pi/rho/Delta."""
    # ceiling = n_q/2
    assert Particle(name="pi+", p=3, q=5, m=5, n_q=2, f=2).isospin_ceiling == 1.0
    assert Particle(name="Delta", p=5, q=4, m=15, n_q=3).isospin_ceiling == 1.5
    # saturating, maximal-isospin states: I = n_q/2
    pip = Particle(name="pi+", p=3, q=5, m=5, n_q=2, f=2, I_pdg=1.0)
    delta = Particle(name="Delta", p=5, q=4, m=15, n_q=3, I_pdg=1.5)
    assert pip.isospin_within_ceiling and pip.I_pdg == pip.isospin_ceiling
    assert delta.isospin_within_ceiling and delta.I_pdg == delta.isospin_ceiling
    # sub-maximal states satisfy the bound strictly
    K = Particle(name="K+", p=2, q=5, m=8, n_q=2, f=1, I_pdg=0.5)
    assert K.isospin_within_ceiling
    # unset isospin -> vacuously within ceiling
    assert Particle(name="X", p=2, q=3, m=5, n_q=3).isospin_within_ceiling
    # a hypothetical violator (I > n_q/2) is flagged but NOT raised
    bad = Particle(name="bad", p=3, q=5, m=5, n_q=2, f=0, I_pdg=2.0)
    assert not bad.isospin_within_ceiling  # 2.0 > 2/2 = 1.0


def test_isospin_ceiling_compendium_consistency():
    """The whole compendium respects I <= n_q/2 (zero violations)."""
    from nwt_substrate.particles.compendium import COMPENDIUM
    from nwt_substrate.particles.particle import Particle as P
    for e in COMPENDIUM:
        if e["n_q"] < 2:
            continue
        part = P(name=e["name"], p=e["p"], q=e["q"], m=e["m"], n_q=e["n_q"],
                 f=e["f"], I_pdg=e["I_pdg"])
        assert part.isospin_within_ceiling, f"{e['name']} violates I <= n_q/2"
