"""Tests for the 24-particle compendium and Paper 6 mass formula benchmark."""

import numpy as np
import pytest

from nwt_substrate import particle, list_particles
from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.particles.factory import all_particles


def test_compendium_size():
    """Compendium currently holds 25 PDG-anchored particles."""
    assert len(COMPENDIUM) == 25
    assert len(list_particles()) == 25


def test_every_particle_buildable():
    """All 24 compendium entries can be built into a Particle without raising."""
    for name in list_particles():
        p = particle(name)
        assert p.name == name


def test_particle_unknown_raises():
    with pytest.raises(KeyError):
        particle("not-a-real-particle")


def test_proton_substrate_quantum_numbers():
    """Proton uses corrected (1, 3, 5, 5) tuple — Q_H = 5 preserved."""
    p = particle("p")
    assert (p.p, p.q, p.m, p.n_q) == (1, 3, 5, 5)
    assert p.f == 1
    assert p.Q_H == 5
    assert p.B_pdg == 1


def test_electron_anchor_values():
    """Electron is the anchor: (p, q, m) = (2, 1, 3), n_q = 0, observed = m_e."""
    e = particle("e-")
    assert (e.p, e.q, e.m, e.n_q) == (2, 1, 3, 0)
    assert e.m_obs == pytest.approx(0.511, abs=1e-3)


def test_paper6_median_residual_under_2_percent():
    """
    Across the full compendium, the Paper 6 mass formula should give
    a median absolute residual below 2% (currently ~0.76% after the
    nucleon-tuple correction).
    """
    residuals = []
    for p in all_particles():
        if p.mass_residual is None:
            continue
        residuals.append(abs(p.mass_residual))
    median = float(np.median(residuals))
    assert median < 2.0, f"median residual {median:.2f}% > 2%"


def test_paper6_proton_residual():
    """Proton predicted mass should be within 1% of 938.27 MeV."""
    p = particle("p")
    assert p.mass_pred is not None
    assert abs(p.mass_residual) < 1.0


def test_paper6_electron_predicts_itself():
    """Electron is the anchor — its predicted ratio is exactly 1."""
    e = particle("e-")
    assert e.mass_pred == pytest.approx(0.510998928, rel=1e-9)


def test_extended_gell_mann_nishijima_for_hadrons():
    """
    For every hadron whose framing matches a single physical state, the
    predicted charge from extended GM-N matches PDG.

    Excluded:
      - leptons (charge from Paper 7 contact-geometry, not GM-N)
      - J=3/2 multiplet placeholders (Delta, Sigma*) that use f=0 to
        represent the average mass; I=3/2 has no I_3=0 member
    """
    leptons = {"e-", "mu-", "tau-"}
    multiplet_placeholders = {"Delta", "Sigma*"}
    skip = leptons | multiplet_placeholders
    for p in all_particles():
        if p.name in skip:
            continue
        if p.B_pdg is None or p.Q_pdg is None:
            continue
        assert p.Q_pred == pytest.approx(p.Q_pdg, abs=1e-9), (
            f"{p.name}: Q_pred={p.Q_pred} but Q_pdg={p.Q_pdg}"
        )


def test_carrier_assignment_consistent_with_n_q():
    """Every compendium entry's carrier name matches its n_q via the canonical map."""
    expected = {
        0: "unknot", 1: "unknot", 2: "Hopf", 3: "trefoil",
        4: "figure-eight", 5: "cinquefoil", 6: "6_1 knot",
    }
    for p in all_particles():
        assert p.carrier == expected[p.n_q], (
            f"{p.name}: n_q={p.n_q} but carrier={p.carrier}"
        )


def test_isospin_projection_from_framing():
    """I_3 = f / 2 for every compendium entry."""
    for p in all_particles():
        assert p.I3 == p.f / 2.0
