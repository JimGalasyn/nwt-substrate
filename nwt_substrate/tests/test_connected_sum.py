"""Tests for nwt_substrate.compositions.connected_sum."""

import pytest

import nwt_substrate as nwt
from nwt_substrate import CompositeParticle, compose, decompose, LAMBDA_QCD_MEV
from nwt_substrate.particles.particle import Particle


def test_compose_requires_two_factors():
    p = Particle(name="p", p=1, q=3, m=5, n_q=5)
    with pytest.raises(ValueError):
        compose(p)


def test_compose_rejects_unknown_op():
    p = Particle(name="p", p=1, q=3, m=5, n_q=5)
    n = Particle(name="n", p=1, q=3, m=5, n_q=5)
    with pytest.raises(ValueError):
        compose(p, n, op="kruskal")


def test_lambda_qcd_constant():
    """LAMBDA_QCD = 313 MeV per Hopf crossing (Paper 15)."""
    assert LAMBDA_QCD_MEV == 313.0


# ---------------------------------------------------------------------------
# Deuteron: textbook molecular ground state, # rule predicts mass to 0.06%
# ---------------------------------------------------------------------------

def test_deuteron_via_connected_sum():
    """p # n -> d, mass-additive at LO matches AME 2020 to <0.1%."""
    p = Particle(name="p", p=1, q=3, m=5, n_q=5, m_obs=938.272)
    n = Particle(name="n", p=1, q=3, m=5, n_q=5, m_obs=939.565)
    d = compose(p, n, op="#", name="d", m_obs=1875.61)
    assert d.op == "#"
    assert d.n_factors == 2
    assert d.mass_pred is not None
    # 2 * paper6_mass(1,3,5,5) = 2 * 937.24 = 1874.48 MeV vs PDG 1875.61
    assert abs(d.mass_residual) < 0.1


# ---------------------------------------------------------------------------
# X(3872): D0 # D*0 molecular tetraquark
# ---------------------------------------------------------------------------

def test_X3872_via_connected_sum_observed():
    """X(3872) ~ D0 + D*0 to <0.01% using observed constituent masses."""
    D0   = Particle(name="D0",  p=3, q=7, m=7, n_q=2, m_obs=1864.84)
    Dstar0 = Particle(name="D*0", p=4, q=6, m=17, n_q=2, m_obs=2006.85)
    X = compose(D0, Dstar0, op="#", name="X(3872)", m_obs=3871.65)
    # Observed constituent sum is essentially exact for X(3872) (it sits
    # within 40 keV of the D0+D*0 threshold)
    sum_obs = X.mass_constituent_sum_obs
    assert sum_obs is not None
    assert abs(sum_obs - 3871.65) < 0.1


# ---------------------------------------------------------------------------
# Genus and crossing additivity (Schubert)
# ---------------------------------------------------------------------------

def test_genus_additive_under_connected_sum():
    """Schubert: genus(K_A # K_B) = genus(K_A) + genus(K_B)."""
    a = Particle(name="A", p=2, q=3, m=4, n_q=3)  # trefoil
    b = Particle(name="B", p=2, q=3, m=4, n_q=3)  # another trefoil
    c = compose(a, b, op="#")
    assert c.total_genus == a.genus + b.genus


def test_crossings_additive_under_connected_sum():
    a = Particle(name="A", p=2, q=3, m=4, n_q=3)
    b = Particle(name="B", p=2, q=5, m=6, n_q=5)
    c = compose(a, b, op="#")
    assert c.total_crossings == a.crossings + b.crossings


# ---------------------------------------------------------------------------
# Schubert decomposition recovers prime factors
# ---------------------------------------------------------------------------

def test_decompose_recovers_factors():
    a = Particle(name="A", p=2, q=3, m=4, n_q=3)
    b = Particle(name="B", p=1, q=4, m=5, n_q=4)
    c = compose(a, b, op="#")
    primes = decompose(c)
    assert len(primes) == 2
    assert primes[0].name == "A" and primes[1].name == "B"


# ---------------------------------------------------------------------------
# Hopf-linked composition: takes a binding correction
# ---------------------------------------------------------------------------

def test_hopf_link_carries_binding():
    """A Hopf-linked composite uses op='hopf' with explicit binding."""
    delta = Particle(name="Delta", p=5, q=4, m=15, n_q=3, m_obs=1232.0)
    # d*(2380) = 2 Delta, ~ 0.27 Hopf crossings worth of binding
    d_star = compose(delta, delta, op="hopf", name="d*(2380)",
                       binding=-0.27 * LAMBDA_QCD_MEV, m_obs=2380.0)
    assert d_star.op == "hopf"
    assert d_star.binding < 0  # bound state
    assert d_star.mass_pred is not None
    # Pred should be in the right ballpark of observed (within a few %)
    assert abs(d_star.mass_residual) < 5.0


# ---------------------------------------------------------------------------
# Top-level nwt.compose API
# ---------------------------------------------------------------------------

def test_top_level_compose_via_factory():
    """nwt.compose works via the factory built-in particles."""
    p = nwt.particle("p")
    n = nwt.particle("n")
    d = nwt.compose(p, n, op="#", name="d", m_obs=1875.61)
    assert isinstance(d, CompositeParticle)
    assert d.mass_pred is not None
    # Same precision check as the explicit-Particle test above
    assert abs(d.mass_residual) < 0.1
