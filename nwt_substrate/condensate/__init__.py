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

__all__ = [
    "AbelianHiggsParams",
    "bogoliubov_dispersion",
    "sound_speed",
    "healing_length",
    "line_tension_BPS",
]
