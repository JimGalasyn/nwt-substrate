"""
Pericyclic reaction selection rules from the substrate Hopf-pair
parity, applied to transition states.

The Woodward-Hoffmann rules for pericyclic reactions (electrocyclic,
cycloaddition, sigmatropic, cheletropic) are unified by the 4n+2 / 4n
parity of electrons in the cyclic transition state:

  4n+2 electrons in TS → Hückel topology → thermal allowed
                          (suprafacial / disrotatory / supra-supra)
  4n electrons in TS   → Möbius topology → thermal needs one antara
                          (antarafacial / conrotatory / Möbius twist)
                          OR photochemical allowed

This is **exactly the same 4n+2 / 4n Hopf-pair parity rule** that
`nwt_substrate.chemistry.aromaticity_class` already uses for ground-
state aromatic classification — applied to transition-state rings
instead of ground-state rings.

Substrate framing
-----------------
NWT's substrate algebra provides the **parity-level prediction**:
which reactions are thermally allowed (4n+2 / Hopf-pair-odd) vs
thermally forbidden (4n / Hopf-pair-even). The specific geometric
details — which face is suprafacial vs antarafacial, which terminus
disrotates vs conrotates — come from frontier molecular orbital
(FMO) theory and orbital-symmetry conservation, not from substrate
algebra alone.

Pericyclic reaction selectivity is therefore the **transition-state
analog of substrate aromaticity classification**. The cleanest
substrate observation in this module: the same Hopf-pair parity
discriminates ground-state aromaticity AND transition-state thermal
allowance. One substrate rule, two domains.

Verified predictions (against classical Woodward-Hoffmann)
-----------------------------------------------------------
- [4+2] Diels-Alder cycloaddition (6 TS electrons): thermal allowed,
  Hückel suprafacial-suprafacial ✓
- [2+2] cycloaddition (4 TS electrons): thermal forbidden suprafacial,
  needs antara (geometrically difficult); photochemical allowed ✓
- Electrocyclic 4n (e.g., butadiene → cyclobutene, 4 electrons):
  thermal conrotatory (Möbius topology) ✓
- Electrocyclic 4n+2 (e.g., hexatriene → cyclohexadiene, 6 electrons):
  thermal disrotatory (Hückel topology) ✓
- [1,3]-H sigmatropic shift (4 TS electrons): thermal antara required,
  geometrically forbidden ✓
- [1,5]-H sigmatropic shift (6 TS electrons): thermal suprafacial ✓
- [3,3]-Cope/Claisen rearrangement (6 TS electrons): thermal allowed ✓
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


# ---------------------------------------------------------------------------
# Selection-rule core
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PericyclicSelectionRule:
    """Selection rule for a pericyclic reaction at given TS electron count."""

    n_electrons: int
    """Number of electrons in the transition-state cyclic array."""

    hopf_pair_count: int
    """Hopf-pair count = n_electrons / 2 (electrons paired into Hopf-link units)."""

    parity: Literal["odd", "even"]
    """Parity of the Hopf-pair count. odd = 4n+2 Hückel; even = 4n Möbius."""

    thermal_allowed_suprafacial: bool
    """Whether the reaction is thermally allowed in the all-suprafacial mode."""

    thermal_required_topology: Literal["Huckel", "Mobius"]
    """Topology required for thermal allowance: Hückel (all supra) or
    Möbius (one antara, conrotatory, etc.)."""

    photochemical_required_topology: Literal["Huckel", "Mobius"]
    """Complementary topology required for photochemical allowance.
    (Excited-state selectivity is flipped from thermal.)"""


def selection_rule(n_electrons: int) -> PericyclicSelectionRule:
    """Compute the Woodward-Hoffmann selection rule for a pericyclic TS
    with `n_electrons` in the cyclic array.

    The substrate algebra's Hopf-pair parity drives the entire rule:

      n_electrons = 4n + 2  ↔  n_pairs odd   ↔  Hückel TS thermal-allowed
      n_electrons = 4n      ↔  n_pairs even  ↔  Möbius TS thermal-allowed

    Parameters
    ----------
    n_electrons : int
        Number of electrons in the cyclic transition-state array.
        Must be even and non-negative.

    Returns
    -------
    PericyclicSelectionRule

    Raises
    ------
    ValueError
        If n_electrons is negative or odd.
    """
    if n_electrons < 0:
        raise ValueError(f"n_electrons must be non-negative, got {n_electrons}")
    if n_electrons % 2 != 0:
        raise ValueError(
            f"n_electrons must be even (paired); got {n_electrons}. "
            f"Pericyclic TS electrons are always paired."
        )

    n_pairs = n_electrons // 2
    parity: Literal["odd", "even"] = "odd" if n_pairs % 2 == 1 else "even"

    # The 4n+2 / 4n rule:
    if parity == "odd":
        # 4n+2: Hückel suprafacial thermal allowed
        return PericyclicSelectionRule(
            n_electrons=n_electrons,
            hopf_pair_count=n_pairs,
            parity="odd",
            thermal_allowed_suprafacial=True,
            thermal_required_topology="Huckel",
            photochemical_required_topology="Mobius",
        )
    else:
        # 4n: Möbius (antarafacial) thermal allowed; suprafacial requires light
        return PericyclicSelectionRule(
            n_electrons=n_electrons,
            hopf_pair_count=n_pairs,
            parity="even",
            thermal_allowed_suprafacial=False,
            thermal_required_topology="Mobius",
            photochemical_required_topology="Huckel",
        )


# ---------------------------------------------------------------------------
# Canonical pericyclic-reaction electron counts
# ---------------------------------------------------------------------------

# Each entry maps a reaction-class label to its transition-state electron
# count. Used to validate the substrate selection rule against the
# textbook Woodward-Hoffmann predictions.
PERICYCLIC_TS_ELECTRON_COUNT: dict[str, int] = {
    # Cycloadditions: m+n electrons in the new cyclic TS
    "[2+2]_cycloaddition":          4,   # 4n      — thermal forbidden supra
    "[4+2]_cycloaddition_Diels_Alder": 6,   # 4n+2  — thermal allowed supra (DA)
    "[4+4]_cycloaddition":          8,   # 4n      — thermal forbidden supra
    "[6+2]_cycloaddition":          8,
    "[6+4]_cycloaddition":          10,  # 4n+2  — thermal allowed
    # Electrocyclic: n π electrons in the open chain become n in TS
    "electrocyclic_4_pi":           4,   # butadiene → cyclobutene
    "electrocyclic_6_pi":           6,   # hexatriene → cyclohexadiene
    "electrocyclic_8_pi":           8,
    # Sigmatropic [1,k]-H shifts: 2 σ + (k-1) π = k+1 electrons
    "sigmatropic_1_3_H":            4,   # [1,3]-H, 4 electrons → thermal forbidden
    "sigmatropic_1_5_H":            6,   # [1,5]-H, 6 electrons → thermal allowed
    "sigmatropic_1_7_H":            8,   # [1,7]-H, 8 electrons → thermal forbidden
    # Sigmatropic [j,k] C-C shifts: j+k electrons
    "sigmatropic_3_3_Cope_Claisen": 6,   # Cope, Claisen — 6 electrons
    "sigmatropic_5_5":              10,  # 10 electrons → allowed
    # Cheletropic reactions (single atom adds): treat like cycloaddition
    "cheletropic_4_electron":       4,
    "cheletropic_6_electron":       6,
    "cheletropic_8_electron":       8,
}


def reaction_selection_rule(reaction_label: str) -> PericyclicSelectionRule:
    """Return the substrate selection rule for a named pericyclic reaction.

    Parameters
    ----------
    reaction_label : str
        One of the keys in PERICYCLIC_TS_ELECTRON_COUNT.

    Returns
    -------
    PericyclicSelectionRule

    Examples
    --------
    >>> r = reaction_selection_rule("[4+2]_cycloaddition_Diels_Alder")
    >>> r.thermal_allowed_suprafacial
    True
    >>> r.thermal_required_topology
    'Huckel'
    """
    if reaction_label not in PERICYCLIC_TS_ELECTRON_COUNT:
        raise KeyError(
            f"unknown reaction {reaction_label!r}; available: "
            f"{sorted(PERICYCLIC_TS_ELECTRON_COUNT)}"
        )
    n_e = PERICYCLIC_TS_ELECTRON_COUNT[reaction_label]
    return selection_rule(n_e)


# ---------------------------------------------------------------------------
# Electrocyclic rotation mode (con/dis) from substrate parity
# ---------------------------------------------------------------------------

def electrocyclic_rotation_mode(
    n_electrons: int,
    thermal: bool = True,
) -> Literal["conrotatory", "disrotatory"]:
    """Predict the rotation mode for an electrocyclic ring closure.

    Maps the substrate Hopf-pair parity to the textbook con/dis rule:

      Thermal:
        4n   → conrotatory (Möbius TS allowed)
        4n+2 → disrotatory (Hückel TS allowed)
      Photochemical: complementary.

    Parameters
    ----------
    n_electrons : int
        Number of π electrons in the open-chain reactant.
    thermal : bool
        True for thermal conditions; False for photochemical.

    Returns
    -------
    "conrotatory" or "disrotatory"
    """
    rule = selection_rule(n_electrons)
    if thermal:
        return "disrotatory" if rule.parity == "odd" else "conrotatory"
    else:
        return "conrotatory" if rule.parity == "odd" else "disrotatory"


__all__ = [
    "PericyclicSelectionRule",
    "selection_rule",
    "PERICYCLIC_TS_ELECTRON_COUNT",
    "reaction_selection_rule",
    "electrocyclic_rotation_mode",
]
