"""K_7 BEC condensate framework — Bogoliubov spectrum derivation.

Implements the substrate matter-generation mechanism from
`research_plan_bogoliubov_spectrum` (VV memory): K_7 BEC at the
cosmogenic bridge T² undergoes Bogoliubov quasiparticle production +
σ-orbit branching at decoherence.

Phase A: extract abelian Higgs Hamiltonian from Paper 16 NWT Lagrangian
near the bridge T², canonical-quantize, identify condensate background.

Phase B: scalar Bogoliubov diagonalization — recovers ξ_substrate = λ̄_C
healing length and analog-Hawking spectrum near bridge horizon.

Phase C-F: σ-orbit projection, production rates, validation against
Ω_b/Ω_c = 25α + 75α² (Planck 0.007%).

See [[research-plan-bogoliubov-spectrum]] for the full plan.
"""
from nwt_substrate.condensate.abelian_higgs import (
    AbelianHiggsParams,
    bogoliubov_dispersion,
    sound_speed,
    healing_length,
    line_tension_BPS,
)
from nwt_substrate.condensate.vortex_profile import (
    BPSVortexProfile,
    GLVortexProfile,
    bps_rhs,
    gl_vortex_rhs,
    solve_bps_vortex,
    solve_gl_vortex,
)
from nwt_substrate.condensate.walks import (
    edge_to_orbit,
    bfs_shortest_walks,
    closed_walk_winding,
    walk_arc_length_unit_torus,
    walk_signed_step_distribution,
    walk_symmetric_d_distribution,
    walk_QR_NR_fractions,
    walk_sigma_sequence,
    walk_sigma_composition,
    walk_sigma_transitions,
    walk_sigma_complete_coverage,
    walk_invariants,
    walk_to_3d_curve,
    WalkInvariants,
    QR_SIGNED, NR_SIGNED,
)
from nwt_substrate.condensate.antimatter import (
    antimatter_reverse,
    antimatter_negate_winding,
    matter_direction_fraction,
    is_pure_matter_walk,
    is_pure_antimatter_walk,
    cpt_verification,
)
from nwt_substrate.condensate.chiral_magnetic import (
    cme_growth_rate,
    cme_current,
    SOURCE_SIM,
    JOYCE_SHAPOSHNIKOV_NOTE,
)

__all__ = [
    # abelian_higgs
    "AbelianHiggsParams", "bogoliubov_dispersion",
    "sound_speed", "healing_length", "line_tension_BPS",
    # vortex_profile
    "BPSVortexProfile", "bps_rhs", "solve_bps_vortex",
    "GLVortexProfile", "gl_vortex_rhs", "solve_gl_vortex",
    # walks
    "edge_to_orbit", "bfs_shortest_walks", "closed_walk_winding",
    "walk_arc_length_unit_torus", "walk_signed_step_distribution",
    "walk_symmetric_d_distribution", "walk_QR_NR_fractions",
    "walk_sigma_sequence", "walk_sigma_composition",
    "walk_sigma_transitions", "walk_sigma_complete_coverage",
    "walk_invariants", "walk_to_3d_curve", "WalkInvariants",
    "QR_SIGNED", "NR_SIGNED",
    # antimatter
    "antimatter_reverse", "antimatter_negate_winding",
    "matter_direction_fraction",
    "is_pure_matter_walk", "is_pure_antimatter_walk",
    "cpt_verification",
    # chiral_magnetic (dynamical-framing Chern-Simons / CME)
    "cme_growth_rate", "cme_current",
    "SOURCE_SIM", "JOYCE_SHAPOSHNIKOV_NOTE",
]
