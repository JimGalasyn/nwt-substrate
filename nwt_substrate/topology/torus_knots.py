"""
Torus knot T(p, q) topology: classification, invariants, and carriers.

The (p, q) torus knot wraps p times around the major axis and q times
around the minor axis of a torus.  Properties:

  - Genus:      g = (p-1)(q-1)/2  for gcd(p, q) = 1 (knot)
                g = 0             for gcd(p, q) > 1 (link, multi-component)
  - Self-linking:  Lk = pq
  - Crossing number:  c = min(p(q-1), q(p-1))  for p, q >= 2

The "carrier knot" of an NWT particle is one of:
  - unknot       (n_q = 0/1)   leptons
  - Hopf link    (n_q = 2)      mesons
  - trefoil 3_1  (n_q = 3)      most baryons
  - figure-8 4_1 (n_q = 4)      tetraquarks
  - cinquefoil 5_1 (n_q = 5)    nucleon family + pentaquarks
  - 6-crossing   (n_q = 6)      hexaquarks / dibaryons
"""

from __future__ import annotations

from math import gcd


def is_torus_knot(p: int, q: int) -> bool:
    """T(p, q) is a knot (single component) iff gcd(p, q) = 1."""
    return gcd(p, q) == 1


def seifert_genus(p: int, q: int) -> int:
    """
    Seifert genus of T(p, q).

    For gcd(p, q) = 1: g = (p-1)(q-1)/2.
    For multi-component links (gcd > 1): genus is not simply defined; return 0.
    """
    if gcd(p, q) > 1:
        return 0
    return (p - 1) * (q - 1) // 2


def crossing_number(p: int, q: int) -> int:
    """
    Standard crossing number of T(p, q).

    For p = 1 or q = 1 (unknot when gcd = 1): 0.
    For p, q >= 2: min(p(q-1), q(p-1)).
    """
    if p <= 1 or q <= 1:
        return 0
    return min(p * (q - 1), q * (p - 1))


def hopf_charge(p: int, m: int) -> int:
    """L3 Hopf charge Q_H = p * m (per Paper 6 / Paper 7)."""
    return p * m


def knot_family(p: int, q: int) -> str:
    """
    Classify (p, q) into a named carrier-knot family.

    Returns one of: "unknot", "Hopf", "trefoil", "cinquefoil",
    "heptafoil", "septafoil", or generic "T(p,q)".
    """
    if p == 1 or q == 1:
        return "unknot"
    if (p, q) in {(2, 2), (1, 2), (2, 1)}:
        return "Hopf"
    if {p, q} == {2, 3}:
        return "trefoil"
    if {p, q} == {2, 5}:
        return "cinquefoil"
    if {p, q} == {2, 7}:
        return "heptafoil"
    if {p, q} == {2, 9}:
        return "septafoil"
    return f"T({p},{q})"


def carrier_for_n_q(n_q: int) -> str:
    """
    Map n_q to its canonical carrier knot.

    n_q = crossing number of carrier ∈ [0, MAX_CROSSING_NUMBER] = [0, 6];
    this is the inverse map. Names sourced from `isa.CARRIER_NAMES`.

    Substrate identity: there are N_CARRIER_TYPES = 7 distinct carrier
    types = N_VERTICES_K7, one per vertex of the K_7 substrate graph.
    """
    from ..isa.constants import CARRIER_NAMES, MAX_CROSSING_NUMBER
    if 0 <= n_q <= MAX_CROSSING_NUMBER:
        return CARRIER_NAMES[n_q]
    return f"({n_q}-crossing carrier)"


def beta_phase(p: int, m_int: int) -> float | None:
    """
    Aspect-ratio parameter beta = sqrt(m^2 / p^2 - 1).

    From L3 phase closure on the (p, q) torus.  Used in Paper 6 mass formula.
    Returns None if m_int <= p (would give imaginary beta).
    """
    import numpy as np
    if m_int * m_int <= p * p:
        return None
    return float(np.sqrt(m_int * m_int / (p * p) - 1.0))
