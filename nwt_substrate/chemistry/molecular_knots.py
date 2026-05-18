"""
Molecular-knot accessibility predictions (chemistry Tier-B.6).

The NWT carrier-knot table n_q ∈ {0, 1, 2, 3, 4, 5, 6} from 2I (binary
icosahedral) irrep dimensions predicts synthesis accessibility for prime
molecular knots, with two strong substrate signals and one qualitatively
preserved signal:

  P1 — Sweet spot at n_q=5: pentafoil 5_1 is the empirical small-molecule
       synthesis + host-functionality maximum (K(Cl⁻) ≈ 3.6 × 10¹⁰ M⁻¹,
       3 orders of magnitude tunable via metal exchange).
  P2 — Gap at n_q=6: stevedore 6_1 is empirically absent from small-
       molecule synthesis literature as of 2026-05-18.
  P3 — Crossing-number inversion: figure-8 4_1 (n_q=4) is functionally
       below pentafoil 5_1 (n_q=5), despite the naive crossing-number
       expectation that more crossings = harder synthesis.

The carrier-knot table is shared with particle physics — the SAME n_q
table predicts (a) which particle carrier-knots are stable (proton at
n_q=5 is the only stable hadron) and (b) molecular-knot synthesis
accessibility (pentafoil sweet spot, stevedore gap). This cross-scale
isomorphism is the strongest empirical signal in the chemistry sector
of the NWT substrate program.

Resolution memo: [[molecular-knot-accessibility-resolution]] (2026-05-18).
Pre-registration: [[molecular-knot-accessibility-prereg]] (2026-05-18,
Form B win — P1 + P2 confirmed, P3 weakened but not reversed).

Public API
----------
    accessibility_for_knot(alexander_briggs) -> MolecularKnotEntry
    carrier_class_of(crossing_number) -> int | None
    substrate_predicted_tier(crossing_number) -> AccessibilityTier
    KNOT_REFERENCE: dict[str, MolecularKnotEntry]

Sources
-------
- Ayme et al. 2012 *Nat. Chem.* — pentafoil 5_1 first synthesis
- Leigh et al. 2015 *JACS* (10.1021/jacs.5b06340) — pentafoil Cl⁻ binding
- Leigh et al. 2019 *JACS* (10.1021/jacs.8b12548) — metal-exchange tunability
- Fielden / Leigh / Goldup 2022 *Chem. Soc. Rev.* — "Knotting matters" review
- Zhang / Jin 2019 *Nat. Commun.*, 2020 *JACS*, 2024 *JACS* — 4_1 template
- Danon et al. 2017 *Science* — 8_19 "tightest knot"
- Wales 2025 *JCTC* — model-polymer multifunnel landscapes (3 knots)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from nwt_substrate.isa.constants import MAX_CROSSING_NUMBER
from nwt_substrate.topology.torus_knots import carrier_for_n_q


# ---------------------------------------------------------------------------
# Accessibility tiers (Form B framing)
# ---------------------------------------------------------------------------

class AccessibilityTier(str, Enum):
    """Substrate-predicted synthesis-accessibility tier for a prime molecular
    knot, keyed by crossing number c → n_q carrier class."""

    TRIVIAL = "trivial"               # c ∈ {0, 1} — unknotted; degenerate carrier
    ACCESSIBLE = "accessible"          # c ∈ {2, 3} — Hopf + trefoil baseline
    HARD = "hard"                      # c = 4 — figure-8: template only, below sweet spot
    SWEET_SPOT = "sweet_spot"          # c = 5 — pentafoil: empirical maximum
    GAP = "gap"                        # c = 6 — stevedore: empirically absent
    OUTSIDE_TABLE = "outside_table"    # c > 6 — beyond n_q cap; non-carrier mechanisms


# ---------------------------------------------------------------------------
# Host-guest data record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HostGuestData:
    """Empirical host-guest binding data for a molecular knot."""
    guest: str                        # e.g. "Cl-", "PF6-"
    binding_constant_M_inv: float     # K in M^-1 (1 / molar)
    solvent: str
    selectivity_note: str = ""        # e.g. ">1e10:1 Cl-/PF6-"


# ---------------------------------------------------------------------------
# Molecular-knot entry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MolecularKnotEntry:
    """Substrate-canonical record for a molecular knot topology."""

    alexander_briggs: str
    """Knot identifier in Alexander-Briggs notation (e.g., '3_1', '5_1', '8_19')."""

    crossing_number: int
    """Standard crossing number c of the knot."""

    n_q_carrier: Optional[int]
    """Carrier-knot index n_q ∈ {0..6} if c ≤ MAX_CROSSING_NUMBER, else None."""

    carrier_name: str
    """Canonical NWT carrier name (from isa.CARRIER_NAMES) or '(c-crossing carrier)'
    for outside-table knots."""

    accessibility_tier: AccessibilityTier
    """Substrate-predicted accessibility tier."""

    first_synthesis_year: Optional[int]
    """Year of first successful small-molecule synthesis, or None if never."""

    first_synthesis_group: str
    """Lead group on first small-molecule synthesis, or '' if not yet synthesized."""

    best_yield_pct: Optional[float]
    """Documented best yield in % (None if not synthesized / not reported)."""

    representative_host_guest: Optional[HostGuestData]
    """Strongest documented host-guest binding (None if not characterized)."""

    substrate_signals: tuple[str, ...] = field(default_factory=tuple)
    """Which of the substrate signals P1/P2/P3 this entry serves
    (e.g., ('P1_sweet_spot',), ('P2_gap',), ('P3_inversion',))."""

    notes: str = ""
    """Free-text notes (sources, context, caveats)."""


# ---------------------------------------------------------------------------
# Reference table — canonical empirical data (2026-05-18 audit)
# ---------------------------------------------------------------------------

_HOPF_CATENANE = MolecularKnotEntry(
    alexander_briggs="hopf_catenane",
    crossing_number=2,
    n_q_carrier=2,
    carrier_name="Hopf",
    accessibility_tier=AccessibilityTier.ACCESSIBLE,
    first_synthesis_year=1983,
    first_synthesis_group="Sauvage",
    best_yield_pct=90.0,
    representative_host_guest=None,
    notes="Baseline [2]catenane via Cu(I) template (Sauvage 1983). Robust workhorse, ~90% RCM yields. Two-component; technically a link, not a single-component knot — included for n_q=2 carrier-class continuity.",
)

_TREFOIL = MolecularKnotEntry(
    alexander_briggs="3_1",
    crossing_number=3,
    n_q_carrier=3,
    carrier_name="trefoil",
    accessibility_tier=AccessibilityTier.ACCESSIBLE,
    first_synthesis_year=1989,
    first_synthesis_group="Sauvage",
    best_yield_pct=98.0,
    representative_host_guest=None,
    notes='Dietrich-Buchecker & Sauvage 1989 Angew. Chem. — "a terrible synthesis challenge". Initial yield <10%; 98% with FeII template (2014). 19-year gap before next prime knot 5_1.',
)

_FIGURE_EIGHT = MolecularKnotEntry(
    alexander_briggs="4_1",
    crossing_number=4,
    n_q_carrier=4,
    carrier_name="figure-eight",
    accessibility_tier=AccessibilityTier.HARD,
    first_synthesis_year=2014,
    first_synthesis_group="Sanders (DCL discovery); Zhang/Jin 2019 (template)",
    best_yield_pct=None,
    representative_host_guest=None,
    substrate_signals=("P3_inversion_weakened",),
    notes=(
        "Originally discovered via dynamic combinatorial library (Sanders 2014). "
        "Three deliberate template syntheses since: Zhang/Jin 2019 Nat. Commun. "
        "(coordination self-assembly, Cp*Rh), 2020 JACS (NDI quadruple π-stack), "
        "2024 JACS (4_1 → 4_1² supramolecular fusion). Every 4_1 small-molecule "
        "remains metal-coordination self-assembly, NOT active-metal-template "
        "covalent closure. No 4_1 anion-binding K approaches pentafoil's 1e10 M^-1. "
        "Pentafoil retains functional superiority in primary literature."
    ),
)

_PENTAFOIL = MolecularKnotEntry(
    alexander_briggs="5_1",
    crossing_number=5,
    n_q_carrier=5,
    carrier_name="cinquefoil",
    accessibility_tier=AccessibilityTier.SWEET_SPOT,
    first_synthesis_year=2012,
    first_synthesis_group="Leigh",
    best_yield_pct=98.0,
    representative_host_guest=HostGuestData(
        guest="Cl-",
        binding_constant_M_inv=3.6e10,
        solvent="MeCN",
        selectivity_note=">1e10:1 Cl-/PF6-; 3 OOM tunable via metal exchange (Leigh 2019)",
    ),
    substrate_signals=("P1_sweet_spot",),
    notes=(
        "Ayme et al. 2012 Nat. Chem. — pentanuclear FeII helicate template, "
        "44% one-pot. Leigh 2015 JACS — K(Cl-)=3.6e10 M^-1, ten C-H...Cl- H-bonds, "
        "comparable to AgCl. Leigh 2019 JACS — Zn/Co/Ni/Cu transmetalation "
        "analogs, 3 orders of magnitude binding tunability. No competing "
        "prime-knot small-molecule host with comparable affinity has been "
        "published as of 2026-05-18. Headline carrier-knot table prediction."
    ),
)

_STEVEDORE = MolecularKnotEntry(
    alexander_briggs="6_1",
    crossing_number=6,
    n_q_carrier=6,
    carrier_name="6_1 knot",
    accessibility_tier=AccessibilityTier.GAP,
    first_synthesis_year=None,
    first_synthesis_group="",
    best_yield_pct=None,
    representative_host_guest=None,
    substrate_signals=("P2_gap",),
    notes=(
        "NEVER synthesized as a small molecule (2026-05-18 lit audit; "
        "confirmed by Wikipedia, Fielden/Leigh/Goldup 2022 Chem. Soc. Rev.). "
        "Known only as a protein fold (DehI dehalogenase). Leigh group "
        "post-pentafoil progression jumps 5_1 (2012) → 8_19 (2017) → 7_4 → "
        "5_2 → 7_1 → Vernier (2022), consistently SKIPPING 6_1. Cross-scale "
        "analog: d*(2380) hexaquark (n_q=6) decays in ~1e-23 s — extreme "
        "rarity at n_q=6 in BOTH particle and chemistry sectors."
    ),
)

_5_2 = MolecularKnotEntry(
    alexander_briggs="5_2",
    crossing_number=5,
    n_q_carrier=5,
    carrier_name="cinquefoil",
    accessibility_tier=AccessibilityTier.SWEET_SPOT,
    first_synthesis_year=2020,
    first_synthesis_group="Leigh",
    best_yield_pct=None,
    representative_host_guest=None,
    notes=(
        "Non-prime-knot variant at c=5. Synthesized by Leigh group in the "
        "2020 wave. Shares n_q=5 carrier class with pentafoil 5_1; predicted "
        "accessible by substrate but has not been characterized as a "
        "host-guest matrix to the same depth as 5_1."
    ),
)

_7_1 = MolecularKnotEntry(
    alexander_briggs="7_1",
    crossing_number=7,
    n_q_carrier=None,
    carrier_name="(7-crossing carrier)",
    accessibility_tier=AccessibilityTier.OUTSIDE_TABLE,
    first_synthesis_year=2020,
    first_synthesis_group="Leigh",
    best_yield_pct=None,
    representative_host_guest=None,
    notes=(
        "Outside the n_q ∈ {0..6} carrier table (MAX_CROSSING_NUMBER=6). "
        "Synthesizable via non-carrier (template/braid) mechanisms. The "
        "carrier-knot table does NOT predict it cannot be made — only that "
        "its substrate identity is non-carrier."
    ),
)

_7_4 = MolecularKnotEntry(
    alexander_briggs="7_4",
    crossing_number=7,
    n_q_carrier=None,
    carrier_name="(7-crossing carrier)",
    accessibility_tier=AccessibilityTier.OUTSIDE_TABLE,
    first_synthesis_year=2020,
    first_synthesis_group="Leigh",
    best_yield_pct=None,
    representative_host_guest=None,
    notes="3×3 woven grid template (Leigh 2020). Outside carrier-knot table.",
)

_8_19 = MolecularKnotEntry(
    alexander_briggs="8_19",
    crossing_number=8,
    n_q_carrier=None,
    carrier_name="(8-crossing carrier)",
    accessibility_tier=AccessibilityTier.OUTSIDE_TABLE,
    first_synthesis_year=2017,
    first_synthesis_group="Leigh / Danon",
    best_yield_pct=62.0,
    representative_host_guest=None,
    notes=(
        "Danon et al. 2017 Science — 4-component FeII₄ braid, 62% yield. "
        '"Tightest knot" synthesized. Outside carrier-knot table (c=8 > 6). '
        "The 5_1 → 8_19 jump (skipping 6_1/6_2/6_3/7_1) in Leigh's program "
        "is consistent with 6-prime family being LESS accessible than "
        "8-crossing non-prime topologies."
    ),
)

_STAR_OF_DAVID = MolecularKnotEntry(
    alexander_briggs="6_3_2",   # 6²₃ Star-of-David (2-component link)
    crossing_number=6,
    n_q_carrier=None,
    carrier_name="(non-carrier link)",
    accessibility_tier=AccessibilityTier.OUTSIDE_TABLE,
    first_synthesis_year=2014,
    first_synthesis_group="Leigh",
    best_yield_pct=None,
    representative_host_guest=None,
    notes=(
        "Star-of-David is a 2-component LINK (not a single-component knot) "
        "at c=6. Substrate carrier-knot table is for single-component knots; "
        "links sit outside. The 2024 flow-synthesis paper (Walji et al. Cell "
        "Rep. Phys. Sci.) extended both Star-of-David and pentafoil to flow."
    ),
)

_TRISKELION = MolecularKnotEntry(
    alexander_briggs="3_1#3_1#3_1",
    crossing_number=9,
    n_q_carrier=None,
    carrier_name="(composite, outside carrier)",
    accessibility_tier=AccessibilityTier.OUTSIDE_TABLE,
    first_synthesis_year=2022,
    first_synthesis_group="Leigh",
    best_yield_pct=None,
    representative_host_guest=None,
    notes=(
        "Ashbridge et al. 2022 Science — Vernier 3:4 template, 'trefoil-of-"
        "trefoils'. Composite knot at c=9; outside the prime carrier-knot table."
    ),
)


KNOT_REFERENCE: dict[str, MolecularKnotEntry] = {
    "hopf_catenane":    _HOPF_CATENANE,
    "3_1":              _TREFOIL,
    "4_1":              _FIGURE_EIGHT,
    "5_1":              _PENTAFOIL,
    "5_2":              _5_2,
    "6_1":              _STEVEDORE,
    "6_3_2":            _STAR_OF_DAVID,
    "7_1":              _7_1,
    "7_4":              _7_4,
    "8_19":             _8_19,
    "3_1#3_1#3_1":      _TRISKELION,
}
"""Canonical empirical reference for documented molecular knots / links.

Keys are Alexander-Briggs identifiers (with conventional NWT modifications
for composite knots and links: '#' for connected sum, '_n_m' for n_m link)."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def carrier_class_of(crossing_number: int) -> Optional[int]:
    """Map a knot's crossing number to its NWT carrier-class index n_q.

    Returns the carrier index n_q ∈ {0, 1, ..., MAX_CROSSING_NUMBER} for
    knots whose crossing number lies in the carrier-knot table, or None if
    the crossing number exceeds MAX_CROSSING_NUMBER (the substrate's
    2I-irrep-dimension cap).

    Parameters
    ----------
    crossing_number : int
        Standard crossing number (non-negative integer).

    Returns
    -------
    int | None
        n_q index in [0, MAX_CROSSING_NUMBER], or None if outside the
        carrier-knot table.

    Examples
    --------
    >>> carrier_class_of(5)   # pentafoil → n_q=5
    5
    >>> carrier_class_of(8)   # 8_19 → outside table
    >>> carrier_class_of(0)   # unknot
    0
    """
    if crossing_number < 0:
        raise ValueError(f"crossing_number must be non-negative; got {crossing_number}")
    if crossing_number > MAX_CROSSING_NUMBER:
        return None
    return crossing_number


def substrate_predicted_tier(crossing_number: int) -> AccessibilityTier:
    """Map crossing number to substrate-predicted accessibility tier.

    Substrate predictions (per [[molecular-knot-accessibility-resolution]]
    Form B):

      c = 0, 1     → TRIVIAL (unknotted, degenerate carrier)
      c = 2, 3     → ACCESSIBLE (Hopf, trefoil — synthesis baseline)
      c = 4        → HARD (figure-8: template-only, below sweet spot)
      c = 5        → SWEET_SPOT (pentafoil: empirical maximum)
      c = 6        → GAP (stevedore: empirically absent)
      c > 6        → OUTSIDE_TABLE (beyond carrier-knot cap)

    Parameters
    ----------
    crossing_number : int
        Standard crossing number.

    Returns
    -------
    AccessibilityTier
    """
    if crossing_number < 0:
        raise ValueError(f"crossing_number must be non-negative; got {crossing_number}")
    if crossing_number <= 1:
        return AccessibilityTier.TRIVIAL
    if crossing_number in (2, 3):
        return AccessibilityTier.ACCESSIBLE
    if crossing_number == 4:
        return AccessibilityTier.HARD
    if crossing_number == 5:
        return AccessibilityTier.SWEET_SPOT
    if crossing_number == 6:
        return AccessibilityTier.GAP
    return AccessibilityTier.OUTSIDE_TABLE


def accessibility_for_knot(alexander_briggs: str) -> MolecularKnotEntry:
    """Return the substrate-canonical molecular-knot entry for the given knot.

    Parameters
    ----------
    alexander_briggs : str
        Knot identifier matching a key in KNOT_REFERENCE.

    Returns
    -------
    MolecularKnotEntry

    Raises
    ------
    KeyError
        If the knot identifier is not in the reference table.

    Examples
    --------
    >>> r = accessibility_for_knot("5_1")
    >>> r.accessibility_tier
    <AccessibilityTier.SWEET_SPOT: 'sweet_spot'>
    >>> r.representative_host_guest.binding_constant_M_inv
    36000000000.0
    """
    if alexander_briggs not in KNOT_REFERENCE:
        raise KeyError(
            f"unknown molecular knot {alexander_briggs!r}; "
            f"available: {sorted(KNOT_REFERENCE)}"
        )
    return KNOT_REFERENCE[alexander_briggs]


__all__ = [
    "AccessibilityTier",
    "HostGuestData",
    "MolecularKnotEntry",
    "KNOT_REFERENCE",
    "carrier_class_of",
    "substrate_predicted_tier",
    "accessibility_for_knot",
]
