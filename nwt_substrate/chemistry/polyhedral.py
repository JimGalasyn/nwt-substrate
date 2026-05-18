"""
Polyhedral cluster substrate identifications (chemistry Tier-B.5).

The Wade-Mingos PSEPT (Polyhedral Skeletal Electron Pair Theory) rules
classify polyhedral boranes, carboranes, and metal clusters by their
skeletal electron pair count (SEPs):

    closo:    n vertices, n+1 SEPs, closed deltahedron (genus 0, S²)
    nido:     n vertices, n+2 SEPs, (n+1)-vertex parent with 1 vertex removed
    arachno:  n vertices, n+3 SEPs, (n+2)-vertex parent with 2 vertices removed
    hypho:    n vertices, n+4 SEPs, (n+3)-vertex parent with 3 vertices removed

The "+k = 3 − χ" offset comes from generic closed-surface MO topology
(Stone's tensor surface harmonic theory) — Form A weak win in the Tier-B.5
audit ([[wade-rules-3d-aromaticity-resolution]]). The substrate-specific
content of Tier B.5 is that the canonical 3D-aromatic small-n closo
boranes B_5 → B_6 → B_7 → B_8 have (n, n+1) vertex/SEP pairs that
traverse the **Spin(7) representation-class ladder**:

    (n, n+1) = (h_v, h_Coxeter)    = (5, 6)   → B_5  trig bipyramid
             = (h_Coxeter, dim_V)   = (6, 7)   → B_6  octahedron
             = (dim_V, dim_S)       = (7, 8)   → B_7  pentag bipyramid
             = (dim_S, N_POS_ROOTS) = (8, 9)   → B_8  snub disphenoid

4-of-4 double-canonical hits in the substrate integer set, with
random-coincidence probability p ≈ 1.5 × 10⁻⁴ under primitive substrate-
integer density ρ ≈ 0.33 in [3, 35].

B_12 H_12^{2-} (icosahedron, most empirically stable closo borane) has
the (n, n+1) = (12, 13) pair canonical via a *separate* substrate
identification:

    12 = K8_PARTITION[2]     (the 12-edge sector of K_8)
    13 = trefoil(p² + q²) = 2² + 3²   (Paper 13 / Paper 20 PMNS NLO closure)

B_9 H_9^{2-} (tricapped trigonal prism) breaks the small-n ladder at the
SEP side (10 is not substrate-canonical), but has a striking *edge-count*
substrate identification:

    deltahedron edges = 3·9 − 6 = 21 = N_EDGES_K7 = dim(so(7))

This is the bonus finding flagged for follow-up in the resolution memo.

The substrate is silent on the closo/nido/arachno discrete classification
beyond standard Euler-χ topology (Form C fails). The "+k" SEP series
derives from generic MO topology, not NWT-specific algebra (Form A weak).

Public API
----------
    wade_classification(n_vertices, n_seps) -> WadeClass
    closo_borane_sep_count(n_vertices) -> int           # = n + 1
    deltahedron_edge_count(n_vertices) -> int           # = 3n − 6
    deltahedron_face_count(n_vertices) -> int           # = 2n − 4
    deltahedron_euler_chi(n_vertices) -> int            # = 2
    closo_borane_substrate_canonical(n_vertices) -> ClosoCanonicalResult

Sources
-------
- Wade, K. *Adv. Inorg. Chem. Radiochem.* 1976, 18, 1–66.
- Mingos, D. M. P. *Acc. Chem. Res.* 1984, 17, 311–319.
- Resolution memo: [[wade-rules-3d-aromaticity-resolution]] (2026-05-18)
- Pre-registration: [[wade-rules-3d-aromaticity-prereg]] (2026-05-18)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Literal

from nwt_substrate.isa.constants import (
    DIM_OCTONION,
    DIM_S_SPIN7,
    DIM_V_SPIN7,
    H_COXETER_SO7,
    H_V_SO7,
    K8_PARTITION,
    N_EDGES_K7,
    N_POS_ROOTS_SO7,
    N_VERTICES_K7,
    RANK_SO7,
)


# ---------------------------------------------------------------------------
# Wade-Mingos classification
# ---------------------------------------------------------------------------

class WadeClass(str, Enum):
    """Wade-Mingos PSEPT cluster classes by SEP offset k = n_seps − n_vertices."""
    CLOSO = "closo"      # k = 1, χ_parent = 2 (closed S²)
    NIDO = "nido"        # k = 2, parent χ = 1 (one cap removed)
    ARACHNO = "arachno"  # k = 3, parent χ = 0 (two caps removed)
    HYPHO = "hypho"      # k = 4, parent χ = -1 (three caps removed)


# Offset → class lookup
_OFFSET_TO_CLASS: dict[int, WadeClass] = {
    1: WadeClass.CLOSO,
    2: WadeClass.NIDO,
    3: WadeClass.ARACHNO,
    4: WadeClass.HYPHO,
}


def wade_classification(n_vertices: int, n_seps: int) -> WadeClass:
    """Classify a polyhedral cluster by Wade-Mingos PSEPT.

    Parameters
    ----------
    n_vertices : int
        Number of cluster vertices (skeletal atoms).
    n_seps : int
        Number of skeletal electron pairs.

    Returns
    -------
    WadeClass
        One of CLOSO (k=1), NIDO (k=2), ARACHNO (k=3), HYPHO (k=4),
        where k = n_seps − n_vertices.

    Raises
    ------
    ValueError
        If the offset k is not in {1, 2, 3, 4} or vertex count is too
        small to support a deltahedron (n_vertices < 4).

    Examples
    --------
    >>> wade_classification(12, 13)  # B_12 H_12^{2-} icosahedron
    <WadeClass.CLOSO: 'closo'>
    >>> wade_classification(5, 7)    # B_5 H_9 pentaborane (nido)
    <WadeClass.NIDO: 'nido'>
    """
    if n_vertices < 4:
        raise ValueError(
            f"n_vertices must be ≥ 4 to form a deltahedron-style cluster; "
            f"got {n_vertices}"
        )
    if n_seps < n_vertices + 1:
        raise ValueError(
            f"n_seps={n_seps} is below the closo minimum (n_vertices+1) "
            f"for n_vertices={n_vertices}"
        )
    offset = n_seps - n_vertices
    if offset not in _OFFSET_TO_CLASS:
        raise ValueError(
            f"SEP offset k = n_seps − n_vertices = {offset} is outside "
            f"the PSEPT closo/nido/arachno/hypho classes {{1, 2, 3, 4}}"
        )
    return _OFFSET_TO_CLASS[offset]


# ---------------------------------------------------------------------------
# Closo deltahedron counts (generic χ topology — Form A weak win)
# ---------------------------------------------------------------------------

def closo_borane_sep_count(n_vertices: int) -> int:
    """SEP count for a closo cluster with `n_vertices`.

    Closo deltahedron has χ = 2 (closed S² surface), giving k = 3 − χ = 1
    via Stone's tensor surface harmonic decomposition:

        n SEPs = n_radial_AOs (n vertices, 1 bonding combination) +
                 n_tangential_pairs (n surface-tangential MOs)
               = 1 + n

    This is generic closed-surface MO topology, NOT NWT-substrate-specific.
    The substrate-specific content is identified via
    :func:`closo_borane_substrate_canonical`.
    """
    if n_vertices < 4:
        raise ValueError(f"n_vertices must be ≥ 4 for a deltahedron; got {n_vertices}")
    return n_vertices + 1


def deltahedron_edge_count(n_vertices: int) -> int:
    """Edge count for a triangulated deltahedron with `n_vertices`.

    For a closed deltahedron (genus 0, all triangular faces), Euler's
    relation V − E + F = 2 combined with 3F = 2E gives:

        E = 3·V − 6
    """
    if n_vertices < 4:
        raise ValueError(f"n_vertices must be ≥ 4 for a deltahedron; got {n_vertices}")
    return 3 * n_vertices - 6


def deltahedron_face_count(n_vertices: int) -> int:
    """Face count for a triangulated deltahedron with `n_vertices`.

        F = 2·V − 4
    """
    if n_vertices < 4:
        raise ValueError(f"n_vertices must be ≥ 4 for a deltahedron; got {n_vertices}")
    return 2 * n_vertices - 4


def deltahedron_euler_chi(n_vertices: int) -> int:
    """Euler characteristic χ = V − E + F = 2 for any closo deltahedron."""
    return 2


# ---------------------------------------------------------------------------
# Substrate-canonical identifications (Form D content)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ClosoCanonicalResult:
    """Substrate-canonical identifications for a closo borane B_n H_n^{2-}.

    Attributes
    ----------
    n_vertices : int
        Cluster vertex count.
    n_seps : int
        SEP count = n_vertices + 1.
    n_canonical_label : str | None
        Substrate-canonical label for `n_vertices`, or None.
    sep_canonical_label : str | None
        Substrate-canonical label for `n_seps`, or None.
    spin7_ladder_position : str | None
        Position in the Spin(7) rep-class ladder for the (n, n+1) pair,
        or None if not on the ladder. The ladder runs B_5 → B_6 → B_7 → B_8
        through {(h_v, h), (h, dim_V), (dim_V, dim_S), (dim_S, N_POS_ROOTS)}.
    edge_canonical_label : str | None
        Substrate-canonical label for the deltahedron's edge count
        (= 3n − 6), or None. B_9 hits N_EDGES_K7 = 21 exactly.
    notes : tuple[str, ...]
        Additional substrate-resonance notes for this cluster.
    """
    n_vertices: int
    n_seps: int
    n_canonical_label: str | None
    sep_canonical_label: str | None
    spin7_ladder_position: str | None
    edge_canonical_label: str | None
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def double_canonical(self) -> bool:
        """Whether both n_vertices and n_seps are substrate-canonical."""
        return self.n_canonical_label is not None and self.sep_canonical_label is not None

    @property
    def on_spin7_ladder(self) -> bool:
        """Whether this cluster is on the Spin(7) rep-class ladder (B_5–B_8)."""
        return self.spin7_ladder_position is not None


# Primitive substrate-canonical integers (locked from the pre-registration).
# Maps integer → display label.
_PRIMITIVE_CANONICAL: dict[int, str] = {
    RANK_SO7:        "RANK_SO7",                     # 3
    H_V_SO7:         "H_V_SO7",                      # 5
    H_COXETER_SO7:   "H_COXETER_SO7",                # 6
    N_VERTICES_K7:   "N_VERTICES_K7",                # 7
    DIM_OCTONION:    "DIM_OCTONION",                 # 8
    N_POS_ROOTS_SO7: "N_POS_ROOTS_SO7",              # 9
    K8_PARTITION[2]: "K8_PARTITION[2]",              # 12
    N_EDGES_K7:      "N_EDGES_K7",                   # 21
}

# Derived substrate-canonical integers (rep-class arithmetic):
#   13 = trefoil p² + q² = 2² + 3²  (Paper 13 / Paper 20 PMNS NLO closure)
_TREFOIL_P2_Q2: int = 2 * 2 + 3 * 3   # = 13


def _canonical_label(k: int) -> str | None:
    """Return the substrate-canonical label for integer `k`, or None.

    Includes the trefoil p²+q² = 13 identification.
    """
    if k in _PRIMITIVE_CANONICAL:
        return _PRIMITIVE_CANONICAL[k]
    if k == _TREFOIL_P2_Q2:
        return "trefoil(p²+q²) = 2²+3² = 13"
    return None


# Spin(7) rep-class ladder positions for B_5 → B_8
_SPIN7_LADDER: dict[int, str] = {
    5: "(h_v, h_Coxeter) = (5, 6)",
    6: "(h_Coxeter, dim_V) = (6, 7)",
    7: "(dim_V, dim_S) = (7, 8)",
    8: "(dim_S, N_POS_ROOTS) = (8, 9)",
}

# Empirically realized canonical closo set polyhedra (for documentation)
CLOSO_POLYHEDRA: dict[int, tuple[str, str]] = {
    5:  ("trigonal bipyramid",        "D_3h"),
    6:  ("octahedron",                "O_h"),
    7:  ("pentagonal bipyramid",      "D_5h"),
    8:  ("snub disphenoid",           "D_2d"),
    9:  ("tricapped trigonal prism",  "D_3h"),
    10: ("bicapped square antiprism", "D_4d"),
    11: ("octadecahedron",            "C_2v"),
    12: ("icosahedron",               "I_h"),
}


def closo_borane_substrate_canonical(n_vertices: int) -> ClosoCanonicalResult:
    """Return the substrate-canonical identifications for closo B_n H_n^{2-}.

    Parameters
    ----------
    n_vertices : int
        Cluster vertex count. Empirically realized canonical closo set is
        n ∈ {5, 6, 7, 8, 9, 10, 11, 12}; the function accepts any n ≥ 4
        but only the canonical set has documented substrate identifications.

    Returns
    -------
    ClosoCanonicalResult
        With substrate-canonical labels for n_vertices, n_seps, and the
        deltahedron edge count, plus Spin(7) rep-class ladder position
        (if any) and additional substrate-resonance notes.

    Examples
    --------
    >>> r = closo_borane_substrate_canonical(7)   # B_7 H_7^{2-}
    >>> r.n_canonical_label
    'N_VERTICES_K7'
    >>> r.sep_canonical_label
    'DIM_OCTONION'
    >>> r.spin7_ladder_position
    '(dim_V, dim_S) = (7, 8)'
    >>> r.double_canonical
    True

    >>> r = closo_borane_substrate_canonical(9)   # B_9 H_9^{2-}: K_7-edge hit
    >>> r.edge_canonical_label
    'N_EDGES_K7'
    """
    if n_vertices < 4:
        raise ValueError(f"n_vertices must be ≥ 4 for a deltahedron; got {n_vertices}")

    n = n_vertices
    sep = n + 1
    edges = 3 * n - 6

    n_lab = _canonical_label(n)
    sep_lab = _canonical_label(sep)
    edge_lab = _canonical_label(edges)
    ladder_pos = _SPIN7_LADDER.get(n)

    notes: list[str] = []
    if ladder_pos is not None:
        notes.append(
            "On the Spin(7) representation-class ladder for B_5–B_8: "
            "(h_v=5, h_Coxeter=6, dim_V=7, dim_S=8, N_POS_ROOTS=9)."
        )
    if n == 9 and edges == N_EDGES_K7:
        notes.append(
            "Deltahedron edge count = 21 = N_EDGES_K7 = dim(so(7)). "
            "The tricapped trigonal prism's 21 edges match the so(7) "
            "generator count — flagged for follow-up as a possible K_7 "
            "edge-graph realization."
        )
    if n == 12 and sep == 13:
        notes.append(
            "(n, n+1) = (12, 13) substrate-canonical via *separate* "
            "identifications: 12 = K_8 partition[2]; 13 = trefoil p²+q² "
            "= 2² + 3² (Paper 13 / Paper 20 PMNS NLO closure). B_12 "
            "is the empirically most stable closo borane (NICS(0) ≈ "
            "−27 to −33 ppm; spherically aromatic icosahedral I_h)."
        )

    return ClosoCanonicalResult(
        n_vertices=n,
        n_seps=sep,
        n_canonical_label=n_lab,
        sep_canonical_label=sep_lab,
        spin7_ladder_position=ladder_pos,
        edge_canonical_label=edge_lab,
        notes=tuple(notes),
    )


__all__ = [
    "WadeClass",
    "wade_classification",
    "closo_borane_sep_count",
    "deltahedron_edge_count",
    "deltahedron_face_count",
    "deltahedron_euler_chi",
    "ClosoCanonicalResult",
    "closo_borane_substrate_canonical",
    "CLOSO_POLYHEDRA",
]
