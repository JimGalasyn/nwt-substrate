"""Antimatter walk operations on K_7 (Paper 23 §5).

Per Phase J (2026-05-20), antimatter walks are REVERSE-DIRECTION walks:

  matter walk     := closed walk on K_7 with QR-signed steps (d ∈ {1, 2, 4})
  antimatter walk := closed walk on K_7 with NR-signed steps (d ∈ {3, 5, 6})
  matter-antimatter pair := walk and its reverse, k → walk[L - k]

This definition satisfies all four structural tests:

  1. CPT mass equality: |p|, |q|, L all preserved under reversal
  2. LAW 1 (Fano vertex coverage): preserved
  3. LAW 2 (QR/NR sector swap): exact swap of step fractions under reversal
  4. σ-orbit composition: preserved (same edges, opposite directions)

The cleanest case is proton/antiproton: matter proton walk has 100%
QR-signed steps (all d=+1 Hamilton); antiproton walk has 100% NR-signed
steps (all d=+6 = −1 mod 7 reverse-Hamilton).

References:
  - [[bogoliubov-phase-j-antimatter]] — verification
  - Paper 23 §5 (Cosmic Baryon Asymmetry from Octonion Substrate)
  - [[eta-B-substrate-derivation]] — η_B = (3/14)·α⁴ asymmetry amplitude
"""
from __future__ import annotations

from nwt_substrate.condensate.walks import (
    walk_signed_step_distribution,
    walk_QR_NR_fractions,
    walk_sigma_composition,
    closed_walk_winding,
    QR_SIGNED, NR_SIGNED,
)


# ---------------------------------------------------------------------------
# Antimatter constructions
# ---------------------------------------------------------------------------

def antimatter_reverse(walk: list[int]) -> list[int]:
    """The CANONICAL antimatter walk = reverse-direction walk.

    Per Phase J verification: this satisfies CPT mass equality + LAW 1/2
    preservation + σ-orbit composition invariance.  Other candidates
    (negative-winding, edge-orientation flip) are either equivalent or
    rejected.
    """
    return list(reversed(walk))


def antimatter_negate_winding(p: int, q: int) -> tuple[int, int]:
    """Equivalent labeling of the antimatter walk via (-p, -q) homology class.

    Mathematically identical to antimatter_reverse() for closed walks at v_0.
    Provided for clarity in places where (p, q) labels are the natural
    coordinate.
    """
    return -p, -q


# ---------------------------------------------------------------------------
# Matter / antimatter classification
# ---------------------------------------------------------------------------

def matter_direction_fraction(walk: list[int]) -> float:
    """Fraction of steps in QR-signed direction.

    1.0 = pure matter walk (proton-like, all-QR Hamilton).
    0.0 = pure antimatter walk (antiproton-like, all-NR Hamilton).
    Values between indicate mixed walks (mesons, leptons).
    """
    qr, _ = walk_QR_NR_fractions(walk)
    return qr


def is_pure_matter_walk(walk: list[int]) -> bool:
    """True iff every step is QR-signed (d ∈ {1, 2, 4})."""
    steps = walk_signed_step_distribution(walk)
    return all(d in QR_SIGNED for d in steps.keys())


def is_pure_antimatter_walk(walk: list[int]) -> bool:
    """True iff every step is NR-signed (d ∈ {3, 5, 6})."""
    steps = walk_signed_step_distribution(walk)
    return all(d in NR_SIGNED for d in steps.keys())


def cpt_verification(walk: list[int]) -> dict:
    """Verify CPT invariants between walk and its antimatter partner.

    Returns dict with keys:
      walk_pq, antiwalk_pq: (p, q) winding (signs flipped)
      same_modulus: |p_m|+|q_m| == |p_a|+|q_a|
      sigma_composition_same: σ-orbit counts identical
      QR_NR_swap: matter QR fraction = antimatter NR fraction
    """
    anti = antimatter_reverse(walk)
    p_m, q_m = closed_walk_winding(walk)
    p_a, q_a = closed_walk_winding(anti)
    sig_m = walk_sigma_composition(walk)
    sig_a = walk_sigma_composition(anti)
    qr_m, nr_m = walk_QR_NR_fractions(walk)
    qr_a, nr_a = walk_QR_NR_fractions(anti)
    return {
        "walk_pq": (p_m, q_m),
        "antiwalk_pq": (p_a, q_a),
        "same_modulus": (abs(p_m or 0) + abs(q_m or 0)
                          == abs(p_a or 0) + abs(q_a or 0)),
        "sigma_composition_same": sig_m == sig_a,
        "QR_matter": qr_m, "NR_matter": nr_m,
        "QR_antimatter": qr_a, "NR_antimatter": nr_a,
        "QR_NR_swap": abs(qr_m - nr_a) < 1e-9 and abs(nr_m - qr_a) < 1e-9,
    }
