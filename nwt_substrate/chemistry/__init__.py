"""
nwt_substrate.chemistry — substrate-algebraic chemistry predictions.

Provides O(1) and O(N) algebraic predictions for chemistry observables
that DFT/CCSD compute at O(N^3) to O(N^7), with substrate-vs-DFT
benchmark machinery for efficiency comparisons.

Two tiers of substrate predictions:

  TIER A — Closed form (substrate O(1), DFT O(N^k)):
    The substrate algebra gives the answer outright via group theory;
    no electronic-structure calculation needed.  Speedup factor:
    10^9 - 10^13 vs DFT depending on system size.  Examples: magic
    numbers, irrep labels, selection rules, orbit sizes, allowed
    coordinations.

  TIER B — Substrate skeleton + 1 calibration (substrate O(N), DFT O(N^k)):
    The substrate provides the structural shape; one experimental
    calibration sets magnitude; subsequent predictions are algebraic.
    Examples: aromatic resonance energies via Hopf-pair rule (5/5
    linear acenes within 6% from 1 benzene calibration); per-AA
    hydrophobicity from K_4 vertex labelling.

Public API:

    >>> import nwt_substrate.chemistry as chem
    >>>
    >>> # Tier A: closed-form predictions
    >>> chem.aromaticity_class("benzene")
    AromaticityResult(classification='aromatic', n_pi_pairs=3, hopf=3, parity='odd')
    >>> chem.fullerene_orbit_counts(60)
    FullereneOrbit(name='C_60', vertices=60, pentagons=12, hexagons=20, edges=90)
    >>> chem.c60_vibrational_summary()
    VibrationalSummary(total_modes=174, ir_modes=4, raman_modes=10, silent=160)
    >>>
    >>> # Tier B: skeleton + calibration
    >>> chem.aromatic_resonance_energy("naphthalene", calibration_kcal=36.0)
    ResonanceEnergy(kcal_per_mol=63.7, n_hopf_pairs=5, ...)
    >>>
    >>> # Benchmarks
    >>> chem.benchmark.compare_to_dft("aromaticity", n_molecules=20)
    BenchmarkResult(substrate_us=12.3, dft_minutes_est=180.0, speedup=8.8e8)

Structure:

    nwt_substrate.chemistry.binary_groups   — 2T, 2O, 2I, BD_n, Z_n
    nwt_substrate.chemistry.aromaticity     — Hopf-pair rule, bond order
    nwt_substrate.chemistry.fullerenes      — C_60 family, orbit counts
    nwt_substrate.chemistry.vibrational     — mode decomposition under point groups
    nwt_substrate.chemistry.mckay           — McKay-admissible coordinations
    nwt_substrate.chemistry.benchmark       — substrate-vs-DFT timing comparison
"""
from .aromaticity import (
    aromaticity_class,
    bond_order,
    aromatic_resonance_energy,
    AromaticityResult,
    ResonanceEnergy,
)
from .smiles import (
    parse_smiles,
    smiles_to_aromaticity,
    smiles_resonance_energy,
    smiles_resonance_energy_batch,
    ring_system_to_so7_adjacency,
    Molecule,
    Atom,
    Bond,
    AromaticRingSystem,
    SmilesAromaticityResult,
)
from .fullerenes import (
    fullerene_orbit_counts,
    c60_anion_magic_states,
    FullereneOrbit,
    AnionState,
)
from .vibrational import (
    c60_vibrational_summary,
    point_group_mode_count,
    VibrationalSummary,
)
from .mckay import (
    allowed_coordinations,
    is_admissible_coordination,
    mckay_orbit_size,
    McKayCheck,
)
from .bond_orders import (
    PI_BOND_ORDER_BENZENE_CLASS,
    HuckelBondOrderResult,
    huckel_bond_orders,
    cyclic_pi_bond_order,
    smiles_pi_bond_orders,
    cc_bond_length_from_pi_order,
)
from .pericyclic import (
    PericyclicSelectionRule,
    selection_rule,
    PERICYCLIC_TS_ELECTRON_COUNT,
    reaction_selection_rule,
    electrocyclic_rotation_mode,
)
from .polyhedral import (
    WadeClass,
    wade_classification,
    closo_borane_sep_count,
    deltahedron_edge_count,
    deltahedron_face_count,
    deltahedron_euler_chi,
    ClosoCanonicalResult,
    closo_borane_substrate_canonical,
    CLOSO_POLYHEDRA,
)
from .molecular_knots import (
    AccessibilityTier,
    HostGuestData,
    MolecularKnotEntry,
    KNOT_REFERENCE,
    carrier_class_of,
    substrate_predicted_tier,
    accessibility_for_knot,
)
from .transition_metal import (
    ElectronCountClass,
    OrganometallicEntry,
    TRANSITION_METAL_REFERENCE,
    electron_count_class,
    ladder_k_for_count,
    substrate_canonical_form,
    is_substrate_predicted_stable,
    transition_metal_entry,
)
from .nmr import (
    NICSSign,
    StructurallyDistinctiveHit,
    NICSReference,
    NICS_REFERENCE,
    nics_sign_from_hopf_parity,
    nics_reference,
)
from . import benchmark

__all__ = [
    # Aromaticity (curated table)
    "aromaticity_class",
    "bond_order",
    "aromatic_resonance_energy",
    "AromaticityResult",
    "ResonanceEnergy",
    # SMILES — substrate prediction from arbitrary molecular structure
    "parse_smiles",
    "smiles_to_aromaticity",
    "smiles_resonance_energy",
    "smiles_resonance_energy_batch",
    "ring_system_to_so7_adjacency",
    "Molecule",
    "Atom",
    "Bond",
    "AromaticRingSystem",
    "SmilesAromaticityResult",
    # Fullerenes
    "fullerene_orbit_counts",
    "c60_anion_magic_states",
    "FullereneOrbit",
    "AnionState",
    # Vibrational
    "c60_vibrational_summary",
    "point_group_mode_count",
    "VibrationalSummary",
    # McKay
    "allowed_coordinations",
    "is_admissible_coordination",
    "mckay_orbit_size",
    "McKayCheck",
    # Bond orders (π Hückel)
    "PI_BOND_ORDER_BENZENE_CLASS",
    "HuckelBondOrderResult",
    "huckel_bond_orders",
    "cyclic_pi_bond_order",
    "smiles_pi_bond_orders",
    "cc_bond_length_from_pi_order",
    # Pericyclic Woodward-Hoffmann
    "PericyclicSelectionRule",
    "selection_rule",
    "PERICYCLIC_TS_ELECTRON_COUNT",
    "reaction_selection_rule",
    "electrocyclic_rotation_mode",
    # Polyhedral clusters (Wade-Mingos PSEPT)
    "WadeClass",
    "wade_classification",
    "closo_borane_sep_count",
    "deltahedron_edge_count",
    "deltahedron_face_count",
    "deltahedron_euler_chi",
    "ClosoCanonicalResult",
    "closo_borane_substrate_canonical",
    "CLOSO_POLYHEDRA",
    # Molecular knots (carrier-knot accessibility, Tier B.6)
    "AccessibilityTier",
    "HostGuestData",
    "MolecularKnotEntry",
    "KNOT_REFERENCE",
    "carrier_class_of",
    "substrate_predicted_tier",
    "accessibility_for_knot",
    # Transition-metal electron-count rules (Tier C.7)
    "ElectronCountClass",
    "OrganometallicEntry",
    "TRANSITION_METAL_REFERENCE",
    "electron_count_class",
    "ladder_k_for_count",
    "substrate_canonical_form",
    "is_substrate_predicted_stable",
    "transition_metal_entry",
    # NMR via Hopf-pair parity (Tier C.8)
    "NICSSign",
    "StructurallyDistinctiveHit",
    "NICSReference",
    "NICS_REFERENCE",
    "nics_sign_from_hopf_parity",
    "nics_reference",
    # Benchmarks
    "benchmark",
]
