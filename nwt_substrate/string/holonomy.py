"""
Special-holonomy classifications used in M-theory and string compactifications,
as substrate-algebraic objects.

Substrate construction       String / M-theory label
-----------------------      -----------------------
so(7) on 7-dim V             Spin(7) holonomy on 8-manifolds (M-theory)
g₂ on 7-dim Im(O)            G₂ holonomy on 7-manifolds (M-theory N=1 D=4)
SU(N) on N-dim fundamental   SU(N) holonomy ↔ Calabi-Yau N-fold
2I ⊂ Spin(7)                 finite-subgroup-of-Spin(7) compactification
                              (memory: paper15-b2-session, b2.10/b2.12)

The substrate's K_7 graph + Cl(0,7) octonion structure carries Spin(7)
holonomy through the explicit 2T → Spin(7) embedding (Phase b2.12); the
trefoil on the Heegaard torus carries 2I orbit structure (memory:
trefoil-uniqueness-lemma, b2.13: K_7 edges ↔ so(7) generators).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HolonomyEntry:
    """
    A special-holonomy classification entry: the holonomy group, the
    string/M-theory compactification context where it appears, and the
    substrate construction that realises it.
    """
    group: str
    target_dim: int
    susy_preserved: str           # "N=2 D=4", "N=1 D=4", etc.
    substrate_realization: str
    string_context: str


HOLONOMY_CATALOG = [
    HolonomyEntry(
        group="trivial",
        target_dim=2,
        susy_preserved="(any -- T² compactification)",
        substrate_realization="Heegaard genus-1 splitting of S³/2I; trefoil on T²",
        string_context="T² compactification of any string theory",
    ),
    HolonomyEntry(
        group="SU(2)",
        target_dim=4,
        susy_preserved="N=4 D=4 (heterotic on K3)",
        substrate_realization="su(2) sector of so(7) (3 internal substrate generators)",
        string_context="K3 compactification, Calabi-Yau 2-fold",
    ),
    HolonomyEntry(
        group="SU(3)",
        target_dim=6,
        susy_preserved="N=2 D=4 (heterotic / type II on CY₃)",
        substrate_realization="SU(3) Gell-Mann generators (substrate QCD)",
        string_context="Calabi-Yau 3-fold; mirror symmetry",
    ),
    HolonomyEntry(
        group="G₂",
        target_dim=7,
        susy_preserved="N=1 D=4 (M-theory on G₂ manifold)",
        substrate_realization=(
            "g₂ ⊂ so(7) automorphisms of octonion product (Cl(0,7) "
            "alternative algebra)"
        ),
        string_context="M-theory on G₂ holonomy 7-manifold → N=1 D=4 SUGRA",
    ),
    HolonomyEntry(
        group="Spin(7)",
        target_dim=8,
        susy_preserved="N=1 D=3 (M-theory on Spin(7) manifold)",
        substrate_realization=(
            "Spin(7) on Cl(0,7) octonion left-mult; 21 generators = "
            "6 Lorentz + 3 internal SU(2) + 12 mixed (memory: "
            "substrate-boson-count-so7); 2T → Spin(7) explicit embedding "
            "(b2.12 octonion-Clifford bivectors)"
        ),
        string_context=(
            "M-theory on Spin(7) holonomy 8-manifold → N=1 D=3 SUGRA; "
            "F-theory uplift candidates"
        ),
    ),
    HolonomyEntry(
        group="2I (finite Spin(7))",
        target_dim=4,
        susy_preserved="(finite-group orbifold)",
        substrate_realization=(
            "Binary icosahedral 2I ⊂ Spin(7); |2I|=120, |Irr(2I)|=6 → n_q ∈ {1..6} "
            "fully populated (memory: nq-is-2I-irrep-dim)"
        ),
        string_context=(
            "ADE = E_8 orbifold singularity (McKay correspondence); F-theory "
            "compactification; heterotic-IIA duality on Σ(2,3,5)"
        ),
    ),
]


def holonomy_lookup(group_name: str) -> HolonomyEntry:
    """Look up a holonomy entry by group name (case-insensitive substring;
    underscores ignored to allow 'G_2' / 'G2' / 'g₂' style queries)."""
    def normalize(s: str) -> str:
        return (s.lower()
                .replace("(", "").replace(")", "")
                .replace("_", "")
                .replace(" ", "")
                .replace("₂", "2").replace("₇", "7"))
    n = normalize(group_name)
    for h in HOLONOMY_CATALOG:
        if normalize(h.group) == n or n in normalize(h.group) or normalize(h.group) in n:
            return h
    raise KeyError(f"Holonomy '{group_name}' not in catalog")


def holonomy_summary() -> str:
    """Pretty-print the substrate ↔ holonomy correspondence."""
    lines = ["Substrate ↔ special-holonomy compactification dictionary:"]
    lines.append("")
    for h in HOLONOMY_CATALOG:
        lines.append(f"  {h.group} (dim {h.target_dim}, {h.susy_preserved})")
        lines.append(f"    substrate: {h.substrate_realization}")
        lines.append(f"    string:    {h.string_context}")
        lines.append("")
    return "\n".join(lines)
