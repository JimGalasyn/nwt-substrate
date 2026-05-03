"""
(p,q)-string view of substrate particles.

A substrate particle with torus-knot label (p,q) maps onto a (p,q)-string
in IIB string theory: an SL(2,Z)-doublet bound state of p F1-strings and
q D1-strings.  The mass-formula contribution m_geom² ∝ (p² + q²) (Paper 6
"three knots → three forces") is *structurally identical* to the
(p,q)-string tension formula T_{(p,q)} = √(p² + q² g_s²) T_{F1} in the
weak-coupling limit g_s → 1.

This is one of the cleaner substrate ↔ string vocabularies: the same
(p,q) integer doublet labels both, and the squared sum (p² + q²) plays
the same kinematic role.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class PQString:
    """
    A substrate particle viewed as a (p,q)-string.

    Attributes
    ----------
    name : str             -- substrate particle name
    p : int                -- F1-string winding number (substrate: knot p)
    q : int                -- D1-string winding number (substrate: knot q)
    knot_type : str        -- T(p,q) torus knot type (= (p,q) when coprime)
    compactification : str -- target manifold name (e.g. "K_7 toroidal")
    mass_geometric : Optional[float]  -- substrate-predicted mass in MeV
    notes : str
    """
    name: str
    p: int
    q: int
    knot_type: str
    compactification: str
    mass_geometric: Optional[float] = None
    notes: str = ""

    @property
    def winding_squared(self) -> int:
        """p² + q²: the substrate gradient-energy term and the squared
        (p,q)-string tension factor."""
        return self.p * self.p + self.q * self.q

    @property
    def is_trefoil(self) -> bool:
        return (self.p, self.q) in [(2, 1), (1, 2), (2, 3), (3, 2)]

    @property
    def sl2z_doublet(self) -> tuple:
        """(p, q) as SL(2,Z) doublet -- same numerical content."""
        return (self.p, self.q)

    def __str__(self) -> str:
        s = f"({self.p}, {self.q})-string '{self.name}' on {self.compactification}"
        s += f"  [winding² = {self.winding_squared}]"
        if self.mass_geometric:
            s += f"  m ≈ {self.mass_geometric} MeV"
        return s


# ===========================================================================
# Canonical substrate-particle ↔ (p,q)-string identifications
# ===========================================================================
#
# (p,q,m,n_q) tuples are the substrate quantum numbers (Paper 6 mass
# formula).  We expose the (p,q) doublet -- this is what couples to the
# string-theoretic SL(2,Z) action.

CANONICAL_PQ_STRINGS = {
    "electron":  PQString(
        "electron", p=2, q=1,
        knot_type="T(2,1) = unknot mod gauge / canonical electron in NWT",
        compactification="K_7 toroidal embedding (Wilson amplitude n=21)",
        mass_geometric=0.5109989,
        notes="Canonical electron: trefoil in 3D, (2,1) winding labels both "
              "the torus-knot carrier and the (p,q)-string.",
    ),
    "muon": PQString(
        "muon", p=2, q=1,
        knot_type="T(2,1)",
        compactification="K_7 toroidal embedding",
        mass_geometric=105.66,
        notes="Same (p,q) doublet as electron; m index in Paper 6 mass "
              "formula distinguishes generations.",
    ),
    "tau": PQString(
        "tau", p=2, q=1,
        knot_type="T(2,1)",
        compactification="K_7 toroidal embedding",
        mass_geometric=1776.86,
        notes="Same (p,q) doublet as electron / muon; tau is a 'stealth "
              "baryon' under substrate topology (memory: nq-is-2I-irrep-dim).",
    ),
    "up": PQString(
        "up", p=1, q=2,
        knot_type="T(1,2) = T(2,1) up to orientation",
        compactification="K_7 toroidal embedding",
        notes="Up-type quark; n_color=3 is internal (not in (p,q)).",
    ),
    "down": PQString(
        "down", p=1, q=2,
        knot_type="T(1,2)",
        compactification="K_7 toroidal embedding (Wilson amplitude n=20)",
        notes="Down-type quark.",
    ),
    "proton": PQString(
        "proton", p=1, q=3,
        knot_type="T(1,3) torus knot",
        compactification="K_7 toroidal embedding (Wilson amplitude n=18)",
        mass_geometric=937.24,
        notes="Nucleon tuple (1,3,5,5) per Paper 6 mass formula update "
              "(memory: mass-formula-refined-2026-04-30).",
    ),
    "neutron": PQString(
        "neutron", p=1, q=3,
        knot_type="T(1,3) torus knot",
        compactification="K_7 toroidal embedding (Wilson amplitude n=18)",
        mass_geometric=939.57,
        notes="Nucleon tuple (1,3,6,5); (p,q) = (1,3) same as proton.",
    ),
    "Z_boson": PQString(
        "Z_boson", p=0, q=0,
        knot_type="(gauge boson; no carrier knot)",
        compactification="K_7 toroidal embedding (Wilson amplitude n=16)",
        mass_geometric=91187.6,
        notes="Z_0 boson at K_7 Wilson tower n=16 (memory: K7-wilson-mass-tower).",
    ),
}


def pq_string(particle_name: str) -> PQString:
    """Look up the (p,q)-string identification of a substrate particle."""
    n = particle_name.lower()
    # Match case-insensitively against catalog keys (also lower-cased)
    catalog_lower = {k.lower(): v for k, v in CANONICAL_PQ_STRINGS.items()}
    if n not in catalog_lower:
        raise KeyError(
            f"Particle '{particle_name}' not in canonical PQ-string catalog. "
            f"Available: {list(CANONICAL_PQ_STRINGS.keys())}"
        )
    return catalog_lower[n]
