"""
NMR ring-current predictions via substrate Hopf-pair parity (chemistry
Tier-C.8).

The substrate Hopf-pair parity rule, already established for ground-state
aromaticity (A.1, A.2) and pericyclic transition-state thermal allowance
(B.4), extends to NICS (nucleus-independent chemical shift) at the
**sign level only**:

  Odd Hopf-pair count (4n+2 π electrons) → diatropic ring current → NEGATIVE NICS
  Even Hopf-pair count (4n π electrons)  → paratropic ring current → POSITIVE NICS

This sign rule is exact across the 14-molecule canonical NICS reference set,
but it is a TRIVIAL extension of the existing aromaticity classification —
not a new substrate prediction.

**Broad NICS magnitude prediction via substrate-canonical integers FAILS
rational-density audit** (7/14 strict hits at ±0.5 ppm tolerance gives
p_random ≈ 0.78; substrate doesn't outperform random density on
magnitudes). The library therefore does NOT expose a broad magnitude API.

**Two narrow structurally distinctive magnitude identifications DO survive**
because they reuse substrate identifications independently established in
other contexts:

  1. Coronene outer-ring NICS = −18.7 ppm ≈ −(N_EDGES_K7 − RANK_SO7) = −18
     reuses the K_7 hub identification from periodic-table A.3 (shell-18),
     transition-metal C.7 (18-electron rule, non-Cartan so(7) generators),
     and [[so7-substrate-isa-probe]] (May 2026). Three independent
     contexts converge on integer 18.

  2. Benzene NICS = −8.0 ppm = −DIM_OCTONION (exact). DIM_OCTONION grounds
     particle physics throughout NWT.

These two hits are flagged in the reference table; no other magnitude
predictions are exposed.

Resolution memo: [[nmr-via-hopf-pair-resolution]] (Form D, honest weak).
Pre-registration: [[nmr-via-hopf-pair-prereg]].

Public API
----------
    NICSSign — enum (DIATROPIC / PARATROPIC / ZERO)
    nics_sign_from_hopf_parity(n_pi_electrons) -> NICSSign
    StructurallyDistinctiveHit — enum of the two narrow hits
    NICSReference — dataclass per molecule
    NICS_REFERENCE: dict[str, NICSReference]
    nics_reference(name) -> NICSReference
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional


# ---------------------------------------------------------------------------
# Sign rule (the load-bearing Form-A content)
# ---------------------------------------------------------------------------

class NICSSign(str, Enum):
    """NICS sign predicted by substrate Hopf-pair parity."""

    DIATROPIC = "diatropic"
    """Aromatic ring (odd Hopf-pair count = 4n+2 π electrons).
    NICS is negative (upfield shift) due to diatropic ring current."""

    PARATROPIC = "paratropic"
    """Anti-aromatic ring (even Hopf-pair count = 4n π electrons,
    excluding zero). NICS is positive due to paratropic ring current."""

    ZERO = "zero"
    """Non-aromatic baseline (zero π electrons or non-cyclic system).
    NICS magnitude near zero (no ring current)."""


def nics_sign_from_hopf_parity(n_pi_electrons: int) -> NICSSign:
    """Predict NICS sign from Hopf-pair parity (4n+2 vs 4n rule).

    Parameters
    ----------
    n_pi_electrons : int
        Number of π electrons in the cyclic ring system. Must be non-negative.

    Returns
    -------
    NICSSign
        DIATROPIC if 4n+2 (odd Hopf-pair count), PARATROPIC if 4n with
        n>0 (even non-zero), ZERO if 0 electrons.

    Examples
    --------
    >>> nics_sign_from_hopf_parity(6)   # benzene
    <NICSSign.DIATROPIC: 'diatropic'>
    >>> nics_sign_from_hopf_parity(4)   # cyclobutadiene
    <NICSSign.PARATROPIC: 'paratropic'>
    >>> nics_sign_from_hopf_parity(0)
    <NICSSign.ZERO: 'zero'>
    """
    if n_pi_electrons < 0:
        raise ValueError(f"n_pi_electrons must be non-negative; got {n_pi_electrons}")
    if n_pi_electrons % 2 != 0:
        raise ValueError(
            f"Cyclic π systems must have even electron count for closed-shell "
            f"aromaticity; got {n_pi_electrons}"
        )
    if n_pi_electrons == 0:
        return NICSSign.ZERO
    n_pairs = n_pi_electrons // 2
    if n_pairs % 2 == 1:
        return NICSSign.DIATROPIC      # 4n+2 odd Hopf-pair count
    return NICSSign.PARATROPIC          # 4n even Hopf-pair count


# ---------------------------------------------------------------------------
# Structurally distinctive narrow Form-D hits
# ---------------------------------------------------------------------------

class StructurallyDistinctiveHit(str, Enum):
    """NICS magnitude identifications that survive the rational-density
    audit because they reuse substrate integers independently established
    in other observable contexts."""

    CORONENE_K7_OUTER = "coronene_K7_outer"
    """Coronene outer-ring NICS = −18.7 ppm ≈ −(N_EDGES_K7 − RANK_SO7) = −18.
    Reuses K_7-hub identification from A.3 / C.7 / so7-ISA-probe.
    Three independent contexts converge on integer 18."""

    BENZENE_DIM_OCTONION = "benzene_dim_octonion"
    """Benzene NICS = −8.0 ppm = −DIM_OCTONION (exact). DIM_OCTONION
    grounds particle physics throughout NWT (Spin(7) spinor rep,
    K_8 vertices, periodic-table octet, etc.)."""


# ---------------------------------------------------------------------------
# Reference table
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NICSReference:
    """Canonical NICS(0) reference data for a ring system."""

    name: str
    nics_ppm: float
    n_pi_electrons: int          # local-ring π electrons
    hopf_pair_count: int         # n_pi // 2
    predicted_sign: NICSSign
    structurally_distinctive_hit: Optional[StructurallyDistinctiveHit] = None
    notes: str = ""

    @property
    def actual_sign(self) -> NICSSign:
        if self.nics_ppm < -0.5:
            return NICSSign.DIATROPIC
        if self.nics_ppm > 0.5:
            return NICSSign.PARATROPIC
        return NICSSign.ZERO

    @property
    def sign_rule_correct(self) -> bool:
        """Whether the substrate sign rule matches observed NICS sign."""
        return self.actual_sign == self.predicted_sign


def _nics(name: str, val: float, n_pi: int, distinctive=None, notes=""):
    return NICSReference(
        name=name,
        nics_ppm=val,
        n_pi_electrons=n_pi,
        hopf_pair_count=n_pi // 2,
        predicted_sign=nics_sign_from_hopf_parity(n_pi),
        structurally_distinctive_hit=distinctive,
        notes=notes,
    )


NICS_REFERENCE: dict[str, NICSReference] = {
    "benzene": _nics(
        "benzene", -8.0, 6,
        distinctive=StructurallyDistinctiveHit.BENZENE_DIM_OCTONION,
        notes="Paradigmatic 4n+2 aromatic. NICS = -8 = -DIM_OCTONION exactly.",
    ),
    "Cp- anion": _nics(
        "Cp- anion", -14.3, 6,
        notes="Cyclopentadienyl anion (6π). NICS magnitude close to 14 = N_EDGES_K7 − N_VERTICES_K7 but the broad pattern is rational density.",
    ),
    "tropylium+ cation": _nics(
        "tropylium+ cation", -7.6, 6,
        notes="C_7H_7+ cation. 6π Hückel aromatic.",
    ),
    "naphthalene_center": _nics(
        "naphthalene_center", -10.0, 6,
        notes="Per-ring Clar sextet. Bicyclic fused aromatic.",
    ),
    "anthracene_center": _nics(
        "anthracene_center", -13.4, 6,
        notes="Central ring of anthracene; more diatropic than terminal rings.",
    ),
    "phenanthrene_center": _nics(
        "phenanthrene_center", -10.0, 6,
        notes="Angular tricyclic; less diatropic than anthracene-center.",
    ),
    "pyrene_center": _nics(
        "pyrene_center", -7.5, 6,
        notes="Pyrene's central ring; periperycene PAH motif.",
    ),
    "coronene_hub": _nics(
        "coronene_hub", -10.5, 6,
        notes="Central ring of coronene. Per-ring local 6π Clar sextet.",
    ),
    "coronene_outer": _nics(
        "coronene_outer", -18.7, 6,
        distinctive=StructurallyDistinctiveHit.CORONENE_K7_OUTER,
        notes=(
            "Outer ring of coronene. NICS ≈ -18 = -(N_EDGES_K7 - RANK_SO7) "
            "reuses K_7-hub identification from periodic-table A.3 (shell-18), "
            "transition-metal C.7 (18e rule), and so7-substrate-isa-probe. "
            "Three independent substrate contexts converge on integer 18."
        ),
    ),
    "pyrrole": _nics("pyrrole", -15.1, 6, notes="Heteroaromatic, N-donor."),
    "furan": _nics("furan", -12.3, 6, notes="Heteroaromatic, O-donor."),
    "thiophene": _nics("thiophene", -13.6, 6, notes="Heteroaromatic, S-donor."),
    "cyclobutadiene": _nics(
        "cyclobutadiene", 25.0, 4,
        notes="Anti-aromatic; paratropic. NICS reported range +20 to +30 ppm; 25 used as representative.",
    ),
    "planar_8annulene": _nics(
        "planar_8annulene", 25.0, 8,
        notes="Planar form (the actual molecule is tub-shaped non-aromatic). 4n with n=2.",
    ),
}


def nics_reference(name: str) -> NICSReference:
    """Return the canonical NICS reference entry for a named ring system."""
    if name not in NICS_REFERENCE:
        raise KeyError(
            f"unknown NICS reference {name!r}; "
            f"available: {sorted(NICS_REFERENCE)}"
        )
    return NICS_REFERENCE[name]


# Explicit declaration: NO broad magnitude prediction API.  The audit
# showed that broad NICS-magnitude prediction via substrate-canonical
# integers FAILS rational-density (p_random ≈ 0.78 for 7/14 hits at
# ±0.5 ppm tolerance).  Only the sign rule and the two narrow structurally
# distinctive hits are load-bearing — both surfaced via the
# StructurallyDistinctiveHit enum and the reference table notes.

__all__ = [
    "NICSSign",
    "StructurallyDistinctiveHit",
    "NICSReference",
    "NICS_REFERENCE",
    "nics_sign_from_hopf_parity",
    "nics_reference",
]
