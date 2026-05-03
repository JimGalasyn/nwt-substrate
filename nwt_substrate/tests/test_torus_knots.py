"""Tests for nwt_substrate.topology.torus_knots."""

import pytest

from nwt_substrate.topology.torus_knots import (
    is_torus_knot,
    seifert_genus,
    crossing_number,
    hopf_charge,
    knot_family,
    carrier_for_n_q,
    beta_phase,
)


def test_unknot_when_p_or_q_is_one():
    """T(1, q) and T(p, 1) are the unknot."""
    assert knot_family(1, 5) == "unknot"
    assert knot_family(7, 1) == "unknot"
    assert crossing_number(1, 5) == 0


def test_trefoil_classification():
    """T(2, 3) and T(3, 2) are the trefoil."""
    assert knot_family(2, 3) == "trefoil"
    assert knot_family(3, 2) == "trefoil"
    assert seifert_genus(2, 3) == 1
    assert crossing_number(2, 3) == 3


def test_cinquefoil_classification():
    """T(2, 5) is the cinquefoil knot, genus 2."""
    assert knot_family(2, 5) == "cinquefoil"
    assert seifert_genus(2, 5) == 2
    assert crossing_number(2, 5) == 5


def test_heptafoil_classification():
    """T(2, 7) is the heptafoil (electron carrier in Paper 6)."""
    assert knot_family(2, 7) == "heptafoil"
    assert seifert_genus(2, 7) == 3


def test_torus_link_when_gcd_gt_1():
    """gcd(p, q) > 1 means multi-component link, not knot."""
    assert not is_torus_knot(2, 4)
    assert is_torus_knot(2, 3)
    assert is_torus_knot(3, 5)


def test_hopf_charge():
    """Q_H = p * m."""
    assert hopf_charge(2, 1) == 2
    assert hopf_charge(1, 5) == 5
    assert hopf_charge(3, 7) == 21


def test_beta_phase_unphysical_when_m_le_p():
    """beta = sqrt(m^2/p^2 - 1) returns None when m^2 <= p^2."""
    assert beta_phase(2, 1) is None
    assert beta_phase(2, 2) is None


def test_beta_phase_electron():
    """Electron (p=2, m=3) gives beta = sqrt(5)/2."""
    import numpy as np
    b = beta_phase(2, 3)
    assert b is not None
    assert abs(b - np.sqrt(5) / 2) < 1e-12


def test_carrier_for_n_q_table():
    """Verify the canonical carrier mapping."""
    assert carrier_for_n_q(0) == "unknot"
    assert carrier_for_n_q(1) == "unknot"
    assert carrier_for_n_q(2) == "Hopf"
    assert carrier_for_n_q(3) == "trefoil"
    assert carrier_for_n_q(4) == "figure-eight"
    assert carrier_for_n_q(5) == "cinquefoil"
    assert carrier_for_n_q(6) == "6_1 knot"
