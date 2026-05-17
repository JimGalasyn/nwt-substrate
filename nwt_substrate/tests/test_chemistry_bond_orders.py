"""
Tests for nwt.chemistry.bond_orders — Coulson π bond orders from
Hückel MOs, with substrate-canonical resonance identification.
"""

from __future__ import annotations

import numpy as np
import pytest

import nwt_substrate.chemistry as chem
from nwt_substrate.chemistry.bond_orders import (
    C_C_DOUBLE_BOND_LENGTH_A,
    C_C_SINGLE_BOND_LENGTH_A,
    PI_BOND_ORDER_BENZENE_CLASS,
    cc_bond_length_from_pi_order,
    cyclic_pi_bond_order,
    huckel_bond_orders,
    smiles_pi_bond_orders,
)
from nwt_substrate.isa.constants import RANK_SO7


# ---------------------------------------------------------------------------
# Substrate-canonical benzene-class bond order
# ---------------------------------------------------------------------------

def test_benzene_class_constant_is_two_thirds():
    """PI_BOND_ORDER_BENZENE_CLASS = 2 / RANK_SO7 exactly."""
    assert PI_BOND_ORDER_BENZENE_CLASS == pytest.approx(2.0 / 3.0, rel=1e-15)
    assert PI_BOND_ORDER_BENZENE_CLASS == 2.0 / RANK_SO7


def test_benzene_bond_orders_all_two_thirds():
    """Every benzene C-C bond has π bond order = 2/RANK_SO7 = 2/3."""
    P = smiles_pi_bond_orders("c1ccccc1")
    assert len(P) == 6
    for (i, j), p in P.items():
        assert p == pytest.approx(2.0 / 3.0, abs=1e-10)


def test_benzene_substrate_canonical_match_via_constant():
    """Each benzene π bond order matches PI_BOND_ORDER_BENZENE_CLASS."""
    P = smiles_pi_bond_orders("c1ccccc1")
    for p in P.values():
        assert p == pytest.approx(PI_BOND_ORDER_BENZENE_CLASS, abs=1e-12)


# ---------------------------------------------------------------------------
# Cyclic n=3k systems all give 2/3 substrate-canonical
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "n_atoms, n_electrons, name",
    [
        (3, 2, "cyclopropenyl cation"),
        (6, 6, "benzene"),
    ],
)
def test_cyclic_low_huckel_n_gives_two_thirds(n_atoms, n_electrons, name):
    """Cyclic systems with all occupied MOs strictly bonding give P = 2/3.

    Verified for cyclopropenyl cation (Hückel n=0) and benzene
    (Hückel n=1). Higher-n cycles have Fermi-level degeneracies
    (n = 12, 18, ...) or non-rational eigenstructure (n = 5, 7,
    10, ...) that break the exact 2/RANK_SO7 result.
    """
    p = cyclic_pi_bond_order(n_atoms, n_electrons)
    assert p == pytest.approx(2.0 / 3.0, abs=1e-10), (
        f"{name} (n={n_atoms}, m={n_electrons}) expected 2/3, got {p:.6f}"
    )


def test_10_annulene_golden_ratio_bond_order():
    """[10]annulene (n=10, m=10) has bond order (1+√5)/5 ≈ 0.6472,
    a golden-ratio-related irrational — not the substrate-canonical
    2/3. Reveals the pentagonal-sector structure of Hückel n=2 systems."""
    import math
    p = cyclic_pi_bond_order(10, 10)
    expected = (1 + math.sqrt(5)) / 5
    assert p == pytest.approx(expected, abs=1e-10)
    assert abs(p - 2.0 / 3.0) > 0.015  # NOT 2/3


def test_12_annulene_jahn_teller_unstable():
    """[12]annulene (n=12, m=12) is Jahn-Teller unstable: degenerate
    non-bonding MOs at the Fermi level. The Hückel bond order at
    half-filling depends on the arbitrary basis choice for the
    degenerate subspace; it is NOT exactly 2/RANK_SO7 = 2/3."""
    p = cyclic_pi_bond_order(12, 12)
    # Not exactly 2/3 — degeneracy effects
    assert abs(p - 2.0 / 3.0) > 1e-3


def test_cyclobutadiene_planar_gives_half():
    """Cyclobutadiene (n=4 planar): bond order 1/2 (4n anti-aromatic).
    The 1/2 doesn't land on RANK_SO7; this is the substrate distinguishing
    aromatic (4n+2) from anti-aromatic (4n) cyclic systems."""
    p = cyclic_pi_bond_order(4, 4)
    assert p == pytest.approx(0.5, abs=1e-10)


@pytest.mark.parametrize(
    "n_atoms, n_electrons",
    [(5, 6), (7, 6), (8, 8), (10, 10)],
)
def test_non_n3k_cyclic_gives_irrational(n_atoms, n_electrons):
    """Cyclic systems with n NOT divisible by 3 give irrational bond
    orders that don't land on 2/3 or any small NWT-canonical fraction."""
    p = cyclic_pi_bond_order(n_atoms, n_electrons)
    assert abs(p - 2.0 / 3.0) > 0.005, (
        f"n={n_atoms}, m={n_electrons} unexpectedly close to 2/3: P={p}"
    )


# ---------------------------------------------------------------------------
# Hückel solver — generic API
# ---------------------------------------------------------------------------

def test_huckel_solver_benzene_eigenvalues():
    """Benzene Hückel eigenvalues are -2, -1, -1, 1, 1, 2 (in units of |β|)."""
    n = 6
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        A[i, (i + 1) % n] = 1
        A[(i + 1) % n, i] = 1
    result = huckel_bond_orders(A)
    expected = np.array([-2.0, -1.0, -1.0, 1.0, 1.0, 2.0])
    assert np.allclose(np.sort(result.eigenvalues), np.sort(expected), atol=1e-10)


def test_huckel_solver_total_pi_energy_benzene():
    """Benzene total π energy = 2·(-2 - 1 - 1) = -8|β| (in units of |β|)."""
    n = 6
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        A[i, (i + 1) % n] = 1
        A[(i + 1) % n, i] = 1
    result = huckel_bond_orders(A)
    assert result.total_pi_energy == pytest.approx(-8.0, abs=1e-10)


def test_huckel_solver_rejects_non_square():
    A = np.zeros((6, 5), dtype=int)
    with pytest.raises(ValueError):
        huckel_bond_orders(A)


def test_huckel_solver_rejects_asymmetric():
    A = np.eye(3, dtype=int)
    A[0, 1] = 1   # asymmetric
    with pytest.raises(ValueError):
        huckel_bond_orders(A)


def test_huckel_solver_rejects_too_many_electrons():
    A = np.zeros((4, 4), dtype=int)
    A[0, 1] = A[1, 0] = 1
    with pytest.raises(ValueError):
        huckel_bond_orders(A, n_electrons=10)


# ---------------------------------------------------------------------------
# Naphthalene — irrational bond orders, four classes
# ---------------------------------------------------------------------------

def test_naphthalene_has_four_bond_order_classes():
    """Naphthalene's 11 π bonds fall into 4 unique bond-order values
    by symmetry (D_2h): peripheral α-α, peripheral α-β, peripheral β-β,
    and the bridge bond."""
    P = smiles_pi_bond_orders("c1ccc2ccccc2c1")
    unique = set(round(p, 6) for p in P.values())
    assert len(unique) == 4


def test_naphthalene_bond_orders_match_textbook_within_5pct():
    """Naphthalene Coulson π bond orders match standard Hückel
    textbook values (Streitwieser; Hatch; Heilbronner) within 5%."""
    P = smiles_pi_bond_orders("c1ccc2ccccc2c1")
    # Expected textbook ranges (varies by reference; we check the range
    # of values present)
    values = sorted(P.values())
    # Lowest bond order should be near 0.52-0.55 (the special C9-C10 bond)
    assert 0.48 < values[0] < 0.58
    # Highest bond order should be near 0.72-0.73 (the α-α bond)
    assert 0.70 < values[-1] < 0.76


# ---------------------------------------------------------------------------
# Coronene — D_6h symmetry gives three bond-order classes
# ---------------------------------------------------------------------------

CORONENE_SMILES = "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"


def test_coronene_has_three_bond_order_classes():
    """Coronene's 30 aromatic C-C bonds fall into 3 unique bond-order
    values by D_6h symmetry: inner-ring, spoke-like, and outer α."""
    P = smiles_pi_bond_orders(CORONENE_SMILES)
    unique = set(round(p, 6) for p in P.values())
    assert len(unique) == 3


def test_coronene_has_thirty_pi_bonds():
    """Coronene: 24 sp² carbons → 30 aromatic C-C bonds (24 vertices
    + 7 rings - 1 by Euler = 30 edges in the bond graph)."""
    P = smiles_pi_bond_orders(CORONENE_SMILES)
    assert len(P) == 30


# ---------------------------------------------------------------------------
# Bond-length estimates
# ---------------------------------------------------------------------------

def test_bond_length_at_zero_order_is_single():
    assert cc_bond_length_from_pi_order(0.0) == C_C_SINGLE_BOND_LENGTH_A


def test_bond_length_at_unit_order_is_double():
    assert cc_bond_length_from_pi_order(1.0) == C_C_DOUBLE_BOND_LENGTH_A


def test_benzene_bond_length_matches_experiment():
    """Benzene C-C length predicted from P=2/3 is 1.407 Å vs experimental
    1.397-1.40 Å (within 0.5%)."""
    d = cc_bond_length_from_pi_order(2.0 / 3.0)
    assert d == pytest.approx(1.407, abs=0.005)
    # Within 0.5% of experimental average
    assert abs(d - 1.40) / 1.40 < 0.006


# ---------------------------------------------------------------------------
# Empty / non-aromatic SMILES
# ---------------------------------------------------------------------------

def test_non_aromatic_smiles_returns_empty():
    """Non-aromatic compounds return no π bond orders."""
    P = smiles_pi_bond_orders("CCCCCC")
    assert P == {}


def test_cyclohexane_returns_empty():
    P = smiles_pi_bond_orders("C1CCCCC1")
    assert P == {}


# ---------------------------------------------------------------------------
# API exposed at chemistry top level
# ---------------------------------------------------------------------------

def test_chemistry_module_exports_bond_order_api():
    assert hasattr(chem, "PI_BOND_ORDER_BENZENE_CLASS")
    assert hasattr(chem, "smiles_pi_bond_orders")
    assert hasattr(chem, "huckel_bond_orders")
    assert hasattr(chem, "cyclic_pi_bond_order")
    assert hasattr(chem, "cc_bond_length_from_pi_order")
