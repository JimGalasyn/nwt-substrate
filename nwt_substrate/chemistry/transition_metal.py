"""
Transition-metal stable electron-count rules from the substrate Spin(7)
rep-class ladder (chemistry Tier-C.7).

The empirical d-block stable electron-count series {18, 16, 14, 12} is
reproduced by a single substrate rule:

    N_e = N_EDGES_K7 − k

with k traversing the first four Spin(7) rep-class invariants:

    k = RANK_SO7         = 3   →  N_e = 18  (octahedral / td / sandwich)
    k = H_V_SO7          = 5   →  N_e = 16  (square-planar d⁸)
    k = N_VERTICES_K7    = 7   →  N_e = 14  (rare; linear d¹⁰, etc.)
    k = N_POS_ROOTS_SO7  = 9   →  N_e = 12  (bent lanthanide metallocenes)

The {3, 5, 7, 9} ladder is the SAME Spin(7) rep-class triple
(rank, dual Coxeter, vector dim, positive roots) that surfaced in
[[wade-rules-3d-aromaticity-resolution]] for the canonical closo-borane
(n, n+1) pairs at B_5–B_8.

The f-block 32-electron rule (lanthanide / actinide sandwich complexes
like cerocene Ce(COT)₂ and uranocene U(COT)₂) reuses the periodic-table
A.3 identification:

    N_e = K_7_TRIANGLES − RANK_SO7 = 35 − 3 = 32

A cross-arc bonus: Fe(CO)₅ (canonical 18e organometallic, trigonal
bipyramid) and B_5 H_5²⁻ (canonical n=5 closo borane, trigonal bipyramid
deltahedron) share IDENTICAL substrate V/E/F identifications:
V=5=H_V_SO7, E=9=N_POS_ROOTS_SO7, F=6=H_COXETER_SO7.

Form D wins with strong Form-A flavor: the ladder reproduces the
{18, 16, 14, 12, 32} series, but 16 and 14 have alternative α−β forms
(28−12 and 35−21 respectively), preventing strict Form-A uniqueness.
The 18, 32, and 12 identifications use uniquely-determined substrate
algebra. Geometry → count mapping (octahedral → 18, square-planar → 16)
remains aufbau-driven; substrate identifies the count *set*, not the
geometry-to-count *mapping*.

Public API
----------
    electron_count_class(n_electrons) -> ElectronCountClass
    substrate_canonical_form(n_electrons) -> str | None
    ladder_k_for_count(n_electrons) -> int | None
    is_substrate_predicted_stable(n_electrons) -> bool
    OrganometallicEntry, TRANSITION_METAL_REFERENCE

Sources
-------
- Resolution memo: [[transition-metal-18e-rule-resolution]]
- Pre-registration: [[transition-metal-18e-rule-prereg]]
- Cross-arc: [[wade-rules-3d-aromaticity-resolution]] (Spin(7) rep-class ladder)
- Foundation: [[periodic-table-shell-structure-resolution]] (18 = 21−3, 32 = 35−3)
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from nwt_substrate.isa.constants import (
    H_V_SO7,
    K8_PARTITION,
    N_EDGES_K7,
    N_POS_ROOTS_SO7,
    N_TRIANGLES_K7,
    N_VERTICES_K7,
    RANK_SO7,
)


# ---------------------------------------------------------------------------
# Electron-count classification
# ---------------------------------------------------------------------------

class ElectronCountClass(str, Enum):
    """Substrate-predicted stable electron-count class for a transition-metal
    or lanthanide/actinide complex.

    The d-block ladder uses k ∈ {RANK_SO7, H_V_SO7, N_VERTICES_K7,
    N_POS_ROOTS_SO7} = {3, 5, 7, 9}, generating N_e = N_EDGES_K7 − k
    ∈ {18, 16, 14, 12}.  The f-block uses N_e = K_7_TRIANGLES − RANK_SO7
    = 32.
    """
    TWELVE_E = "12e"
    """12-electron count.  k = N_POS_ROOTS_SO7 = 9 (ladder) AND
    K8_PARTITION[2] = 12 (direct).  TWO independent substrate routes.
    Documented: bent lanthanide metallocenes Cp₂Sm, Cp₂Eu; some W/Re clusters."""

    FOURTEEN_E = "14e"
    """14-electron count.  k = N_VERTICES_K7 = 7 (= dim_V Spin(7)).
    Rare empirically: W(CO)₃(PCy₃)₂, Pt(PCy₃)₂, TiCp₂Cl₂.  Alternative
    substrate form: 35−21 = K_7_TRIANGLES − N_EDGES_K7."""

    SIXTEEN_E = "16e"
    """16-electron count.  k = H_V_SO7 = 5 (dual Coxeter).  Square-planar d⁸
    (Pt(II), Pd(II), Wilkinson's, Vaska's).  Alternative substrate form:
    28−12 = N_EDGES_K8 − K8_PARTITION[2]."""

    EIGHTEEN_E = "18e"
    """18-electron count.  k = RANK_SO7 = 3.  UNIQUE substrate identification
    (no other primitive α−β form gives 18 in the canonical integer set).
    Octahedral, tetrahedral, sandwich, trigonal-bipyramid d-block."""

    THIRTYTWO_E = "32e"
    """32-electron count.  K_7_TRIANGLES − RANK_SO7 = 35 − 3.  UNIQUE
    substrate identification.  Lanthanide/actinide bis-cyclooctatetraene
    sandwich (cerocene Ce(COT)₂, uranocene U(COT)₂)."""

    OTHER = "other"
    """Electron count not on the substrate ladder."""


# Mapping from k → (electron count, ElectronCountClass, k-label)
_D_BLOCK_LADDER: dict[int, tuple[int, ElectronCountClass, str]] = {
    RANK_SO7:        (N_EDGES_K7 - RANK_SO7,        ElectronCountClass.EIGHTEEN_E,  "RANK_SO7"),
    H_V_SO7:         (N_EDGES_K7 - H_V_SO7,         ElectronCountClass.SIXTEEN_E,   "H_V_SO7"),
    N_VERTICES_K7:   (N_EDGES_K7 - N_VERTICES_K7,   ElectronCountClass.FOURTEEN_E,  "N_VERTICES_K7"),
    N_POS_ROOTS_SO7: (N_EDGES_K7 - N_POS_ROOTS_SO7, ElectronCountClass.TWELVE_E,    "N_POS_ROOTS_SO7"),
}

# Inverse: electron count → k (for ladder lookup)
_COUNT_TO_K: dict[int, int] = {ne: k for k, (ne, _, _) in _D_BLOCK_LADDER.items()}

# Substrate canonical form for each predicted count
_CANONICAL_FORM: dict[int, str] = {
    18: "N_EDGES_K7 − RANK_SO7 = 21 − 3",
    16: "N_EDGES_K7 − H_V_SO7 = 21 − 5",
    14: "N_EDGES_K7 − N_VERTICES_K7 = 21 − 7",
    12: "N_EDGES_K7 − N_POS_ROOTS_SO7 = 21 − 9",
    32: "K_7_TRIANGLES − RANK_SO7 = 35 − 3",
}


def electron_count_class(n_electrons: int) -> ElectronCountClass:
    """Classify a complex by its valence electron count.

    Parameters
    ----------
    n_electrons : int
        Total valence electron count (d-electrons + ligand donor pairs × 2,
        or analogous f-block count).

    Returns
    -------
    ElectronCountClass

    Examples
    --------
    >>> electron_count_class(18)
    <ElectronCountClass.EIGHTEEN_E: '18e'>
    >>> electron_count_class(16)
    <ElectronCountClass.SIXTEEN_E: '16e'>
    >>> electron_count_class(20)
    <ElectronCountClass.OTHER: 'other'>
    """
    if n_electrons < 0:
        raise ValueError(f"n_electrons must be non-negative; got {n_electrons}")
    for k, (ne, cls, _) in _D_BLOCK_LADDER.items():
        if n_electrons == ne:
            return cls
    if n_electrons == N_TRIANGLES_K7 - RANK_SO7:   # 35 - 3 = 32 (f-block)
        return ElectronCountClass.THIRTYTWO_E
    return ElectronCountClass.OTHER


def ladder_k_for_count(n_electrons: int) -> Optional[int]:
    """Return the d-block ladder index k such that n_electrons = N_EDGES_K7 − k.

    Returns None if the count is not on the d-block ladder
    (e.g., 32-electron f-block, or counts outside the substrate prediction).

    Examples
    --------
    >>> ladder_k_for_count(18)   # k = RANK_SO7
    3
    >>> ladder_k_for_count(16)   # k = H_V_SO7
    5
    >>> ladder_k_for_count(12)   # k = N_POS_ROOTS_SO7
    9
    >>> ladder_k_for_count(32)   # f-block, not on d-block ladder
    """
    return _COUNT_TO_K.get(n_electrons)


def substrate_canonical_form(n_electrons: int) -> Optional[str]:
    """Return the canonical substrate-algebraic form for a stable electron count.

    Returns None if the count is not substrate-predicted stable.

    Examples
    --------
    >>> substrate_canonical_form(18)
    'N_EDGES_K7 − RANK_SO7 = 21 − 3'
    >>> substrate_canonical_form(32)
    'K_7_TRIANGLES − RANK_SO7 = 35 − 3'
    >>> substrate_canonical_form(20)
    """
    return _CANONICAL_FORM.get(n_electrons)


def is_substrate_predicted_stable(n_electrons: int) -> bool:
    """Whether the electron count is on the substrate-predicted stable set
    {12, 14, 16, 18, 32}."""
    if n_electrons < 0:
        raise ValueError(f"n_electrons must be non-negative; got {n_electrons}")
    return n_electrons in _CANONICAL_FORM


# ---------------------------------------------------------------------------
# Canonical organometallic reference set
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OrganometallicEntry:
    """Canonical empirical record for a transition-metal / f-block complex."""

    formula: str
    geometry: str
    d_or_f_config: str
    electron_count: int
    rule_class: ElectronCountClass
    substrate_form: Optional[str]
    notes: str = ""

    @property
    def on_substrate_ladder(self) -> bool:
        return self.substrate_form is not None


def _make(formula, geometry, dconf, ne, notes=""):
    cls = electron_count_class(ne)
    form = substrate_canonical_form(ne)
    return OrganometallicEntry(
        formula=formula,
        geometry=geometry,
        d_or_f_config=dconf,
        electron_count=ne,
        rule_class=cls,
        substrate_form=form,
        notes=notes,
    )


TRANSITION_METAL_REFERENCE: dict[str, OrganometallicEntry] = {
    # 18-electron canonical set
    "Cr(CO)6":           _make("Cr(CO)6",           "octahedral",      "d6",  18, "Group-6 hexacarbonyl, paradigmatic 18e octahedral"),
    "Fe(CO)5":           _make("Fe(CO)5",           "trig bipyramid",  "d8",  18, "Group-8 pentacarbonyl. Trigonal bipyramid (V=5, E=9, F=6) shares substrate V/E/F with B_5 H_5^{2-} from Wade's rules B.5."),
    "Ni(CO)4":           _make("Ni(CO)4",           "tetrahedral",     "d10", 18, "Group-10 tetracarbonyl"),
    "Mn2(CO)10":         _make("Mn2(CO)10",         "tbp dimer",       "d7",  18, "Dimanganese decacarbonyl; M-M bond counted"),
    "Fe(Cp)2":           _make("Fe(Cp)2",           "sandwich",        "d6",  18, "Ferrocene; Cp = η^5-C_5H_5 donating 6e each"),
    "HMn(CO)5":          _make("HMn(CO)5",          "C_4v",            "d6",  18, "Hydridopentacarbonylmanganate"),
    "Co(NH3)6 3+":       _make("Co(NH3)6 3+",       "octahedral",      "d6",  18, "Werner's classical complex"),
    "V(CO)6 -":          _make("V(CO)6 -",          "octahedral",      "d6",  18, "Hexacarbonylvanadate anion"),
    # 16-electron square-planar d⁸ set
    "Pt(NH3)2Cl2":       _make("Pt(NH3)2Cl2",       "square planar",   "d8",  16, "Cisplatin geometry; Pt(II) d⁸"),
    "RhCl(PPh3)3":       _make("RhCl(PPh3)3",       "square planar",   "d8",  16, "Wilkinson's catalyst"),
    "IrCl(CO)(PPh3)2":   _make("IrCl(CO)(PPh3)2",   "square planar",   "d8",  16, "Vaska's complex"),
    "Ni(CN)4 2-":        _make("Ni(CN)4 2-",        "square planar",   "d8",  16, "Tetracyanonickelate"),
    "Pd(PPh3)4":         _make("Pd(PPh3)4",         "tetrahedral",     "d10", 16, "Pd(0) tetrakis(triphenylphosphine)"),
    # 14-electron rare set
    "W(CO)3(PCy3)2":     _make("W(CO)3(PCy3)2",     "trig pyramidal",  "d6",  14, "Coordinatively unsaturated, trans-effect stabilized"),
    "Pt(PCy3)2":         _make("Pt(PCy3)2",         "linear",          "d10", 14, "Coordinatively unsaturated Pt(0)"),
    "TiCp2Cl2":          _make("TiCp2Cl2",          "bent",            "d0",  14, "Titanocene dichloride; bent sandwich"),
    # 12-electron probe set
    "Cp2Sm":             _make("Cp2Sm",             "bent",            "f6",  12, "Bent samarocene; low-VE lanthanide metallocene"),
    "Cp2Eu":             _make("Cp2Eu",             "bent",            "f7",  12, "Bent europocene; low-VE lanthanide metallocene"),
    # 32-electron f-block sandwich
    "Ce(COT)2":          _make("Ce(COT)2",          "sandwich",        "f1",  32, "Cerocene; bis(η^8-cyclooctatetraene)"),
    "U(COT)2":           _make("U(COT)2",           "sandwich",        "f2",  32, "Uranocene; bis(η^8-cyclooctatetraene)"),
}


def transition_metal_entry(formula: str) -> OrganometallicEntry:
    """Return the canonical reference entry for the given organometallic formula.

    Raises
    ------
    KeyError
        If the formula is not in TRANSITION_METAL_REFERENCE.
    """
    if formula not in TRANSITION_METAL_REFERENCE:
        raise KeyError(
            f"unknown transition-metal complex {formula!r}; "
            f"available: {sorted(TRANSITION_METAL_REFERENCE)}"
        )
    return TRANSITION_METAL_REFERENCE[formula]


__all__ = [
    "ElectronCountClass",
    "OrganometallicEntry",
    "TRANSITION_METAL_REFERENCE",
    "electron_count_class",
    "ladder_k_for_count",
    "substrate_canonical_form",
    "is_substrate_predicted_stable",
    "transition_metal_entry",
]
