"""
String-theoretic compactification targets, as substrate objects.

Each target is a substrate-algebraic object (graph, manifold, surgery
sphere, etc.) re-presented in the vocabulary of string compactification.
The substrate provides the construction; the string-theoretic label
(genus, holonomy, ADE type, ...) names the same object in another
community's language.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Optional


@dataclass(frozen=True)
class CompactificationTarget:
    """
    A substrate object viewed as a string-compactification target.

    Attributes
    ----------
    name : str                    -- substrate name (e.g. "K_7", "S³/2I")
    string_label : str            -- string-theoretic name (e.g. "T² embedding")
    n_dim : int                   -- (real) dimension
    holonomy : str                -- "Spin(7)", "G_2", "SU(N)", "trivial", ...
    genus : Optional[int]         -- topological genus if relevant
    n_cycles : Optional[int]      -- # independent cycles (for KK / Wilson lines)
    ade_type : Optional[str]      -- ADE classification if relevant
    substrate_module : str        -- where the substrate construction lives
    notes : str
    """
    name: str
    string_label: str
    n_dim: int
    holonomy: str
    genus: Optional[int] = None
    n_cycles: Optional[int] = None
    ade_type: Optional[str] = None
    substrate_module: str = ""
    notes: str = ""

    def __str__(self) -> str:
        bits = [f"{self.name} ({self.string_label})"]
        bits.append(f"  dim = {self.n_dim}, holonomy = {self.holonomy}")
        if self.genus is not None:
            bits.append(f"  genus = {self.genus}")
        if self.n_cycles is not None:
            bits.append(f"  cycles = {self.n_cycles}")
        if self.ade_type:
            bits.append(f"  ADE type = {self.ade_type}")
        if self.substrate_module:
            bits.append(f"  substrate construction: {self.substrate_module}")
        return "\n".join(bits)


# ===========================================================================
# Canonical compactification targets in the substrate
# ===========================================================================

K7_torus = CompactificationTarget(
    name="K_7",
    string_label="T² embedding (Heffter triangular, V=7 E=21 F=14)",
    n_dim=2,                    # toroidal embedding is 2-dimensional
    holonomy="Spin(7)",         # K_7 sits inside Spin(7) via 2T → Spin(7)
    genus=1,                    # toroidal
    n_cycles=21,                # 21 K_7 edges = 21 Wilson-line cycles
    ade_type="E_8 (via 2I → 2T → Spin(7) → K_7 stabilizer chain)",
    substrate_module=(
        "K_7 graph: 21 edges, 7 vertices; toroidal embedding "
        "(memory: K7-on-torus-matter-gravity-unification)"
    ),
    notes=(
        "Same object as the K_7 graph state used in the heron shim and the "
        "Wilson amplitude tower α^(n/2)·M_Pl for n=1..21.  In string-theoretic "
        "view, the 21 edges are 1-cycles supporting Wilson lines; matter "
        "couples to local (p,q) windings, gravity to global 21-edge amplitude."
    ),
)


poincare_sphere = CompactificationTarget(
    name="S³/2I = Σ(2,3,5)",
    string_label="Brieskorn sphere with E_8 ADE structure",
    n_dim=3,
    holonomy="2I ⊂ Spin(7) via Cl(0,7) octonion-bivector embedding",
    genus=None,
    n_cycles=None,
    ade_type="E_8 (icosahedral case of finite ADE classification)",
    substrate_module=(
        "Brieskorn surgery sphere from trefoil; spectral data "
        "λ₁(S³/2I) = 168 = 7·24 (memory: poincare-eigenvalue-168, "
        "trefoil-uniqueness-lemma)"
    ),
    notes=(
        "UNIQUE spherical Brieskorn integer homology 3-sphere among the "
        "carrier-knot family.  In string theory, this is the ADE singularity "
        "of type E_8 -- matches the McKay correspondence 2I ↔ E_8 used in "
        "F-theory and heterotic compactifications."
    ),
)


T2_torus = CompactificationTarget(
    name="T² (BPS vortex carrier)",
    string_label="Standard T² compactification torus",
    n_dim=2,
    holonomy="trivial (flat T²)",
    genus=1,
    n_cycles=2,
    ade_type=None,
    substrate_module=(
        "Heegaard genus-1 splitting of S³/2I; trefoil lives on the splitting "
        "torus (memory: paper15-b2-session)"
    ),
    notes=(
        "The (p,q) torus-knot windings (p loops around one cycle, q around "
        "the other) become (p,q) string winding numbers under the T² → "
        "compactification dictionary."
    ),
)


# ===========================================================================
# ADE / McKay correspondence catalog
# ===========================================================================

@dataclass(frozen=True)
class ADEEntry:
    """
    Entry in the ADE / McKay-correspondence catalog: a finite subgroup of
    SU(2), the corresponding simply-laced Lie algebra, and the substrate
    object that realises it.
    """
    finite_subgroup: str
    lie_algebra: str
    coxeter_number: int
    substrate_realization: str

    def __str__(self) -> str:
        return (
            f"{self.finite_subgroup} ↔ {self.lie_algebra} "
            f"(h = {self.coxeter_number}); substrate: {self.substrate_realization}"
        )


ADE_CATALOG = [
    ADEEntry("Z_n (cyclic)", "A_{n-1}", -1,
             "(not used in canonical NWT)"),
    ADEEntry("D_n (binary dihedral)", "D_n", -1,
             "(not used in canonical NWT)"),
    ADEEntry("2T (binary tetrahedral, |2T|=24)", "E_6", 12,
             "trefoil 2I-orbit + Spin(7) embedding "
             "(memory: paper15-b2-session, b2.10/b2.12)"),
    ADEEntry("2O (binary octahedral, |2O|=48)", "E_7", 18,
             "(not directly canonical)"),
    ADEEntry("2I (binary icosahedral, |2I|=120)", "E_8", 30,
             "Poincaré sphere S³/2I = Σ(2,3,5); n_q ∈ {1..6} = dim(2I irreps) "
             "(memory: nq-is-2I-irrep-dim, all-6-2I-irrep-dims-populated)"),
]


def ade_lookup(name: str) -> ADEEntry:
    """Look up an ADE entry by finite-subgroup name."""
    name_l = name.lower()
    for e in ADE_CATALOG:
        if name_l in e.finite_subgroup.lower():
            return e
    raise KeyError(f"ADE entry '{name}' not in catalog")
