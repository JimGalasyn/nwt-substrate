"""
Baryon-to-CDM ratio Ω_b/Ω_c = 25α + 75α² from the substrate (K₇⊗K₈ partition).

The cosmic baryon-to-cold-dark-matter ratio is fixed by the substrate algebra at
zero free parameters:

    Ω_b/Ω_c = q²·α + RANK_SO7·q²·α²   with   q² = H_V_SO7² = 25
            = 25α + 75α² ≈ 0.18643

against Planck 2018 (Ω_b h² = 0.02237, Ω_c h² = 0.1200 ⟹ ratio 0.18642):
**pred/obs ≈ 1 + 67 ppm** (with the substrate α = 1/(25π√3+1)).

Integer readings (the LO coefficient 25 has two coincident substrate readings):

  q² = 25 = H_V_SO7²                          [VV cross-arc / substrate "DNA"]
         = cinquefoil (5,2)-torus-knot crossings²  [NWT knot reading]
  The cinquefoil is the (5,2) torus knot — 5 = h_v — so q_cinq² = h_v² = 25.

  NLO factor (1 + RANK_SO7·α):
    RANK_SO7 = 3 = N_generations = the Heffter Z₃ σ-orbit count.
    Structurally this is the "at-most-one-loop across the 3 bridging edges"
    truncation (framework_omega_b_c_nlo_truncation): 3 generations × LO × α.

Derivation route: VV's K₇⊗K₈ bridge-partition Wilson amplitudes (Stage-7;
25α LO = q² × per-edge √α, 75α² NLO), consuming the K₈ code
(:mod:`nwt_substrate.algebra.codes.k8`) as the operator handle. NWT confirmed
no independent derivation on their side — this form is canonical.

Provides:
  - OMEGA_B_C_PRED          : substrate prediction (float)
  - OMEGA_B_C_PLANCK        : observed reference ratio
  - Q_CINQ_SQ               : the LO coefficient (= H_V_SO7² = 25)
  - omega_b_c()             : the prediction
  - summary()               : prediction vs observation (dict)

All values verified by tests in tests/test_omega_b_c.py.
"""
from __future__ import annotations

from nwt_substrate.isa.constants import (
    ALPHA_SUBSTRATE,
    H_V_SO7,
    RANK_SO7,
)

__all__ = [
    "Q_CINQ_SQ",
    "OMEGA_B_C_PRED",
    "OMEGA_B_H2_PLANCK",
    "OMEGA_C_H2_PLANCK",
    "OMEGA_B_C_PLANCK",
    "omega_b_c",
    "summary",
]

#: LO coefficient: cinquefoil (5,2)-knot crossings² = h_v² = 25.
Q_CINQ_SQ: int = H_V_SO7**2

#: Substrate prediction Ω_b/Ω_c = q²·α + RANK_SO7·q²·α² = 25α + 75α².
OMEGA_B_C_PRED: float = (
    Q_CINQ_SQ * ALPHA_SUBSTRATE + RANK_SO7 * Q_CINQ_SQ * ALPHA_SUBSTRATE**2
)

#: Planck 2018 best-fit physical densities.
OMEGA_B_H2_PLANCK: float = 0.02237
OMEGA_C_H2_PLANCK: float = 0.1200
#: Observed baryon-to-CDM ratio (h² cancels).
OMEGA_B_C_PLANCK: float = OMEGA_B_H2_PLANCK / OMEGA_C_H2_PLANCK


def omega_b_c() -> float:
    """Substrate baryon-to-CDM ratio Ω_b/Ω_c = 25α + 75α²."""
    return OMEGA_B_C_PRED


def summary() -> dict:
    """Prediction vs observation, as a dict (numbers for tables/papers)."""
    return {
        "omega_b_c_pred": OMEGA_B_C_PRED,
        "formula": "H_V_SO7**2 * alpha + RANK_SO7 * H_V_SO7**2 * alpha**2 = 25a + 75a**2",
        "lo_coeff": Q_CINQ_SQ,
        "nlo_coeff": RANK_SO7 * Q_CINQ_SQ,
        "omega_b_c_planck": OMEGA_B_C_PLANCK,
        "pred_over_planck": OMEGA_B_C_PRED / OMEGA_B_C_PLANCK,
        "ppm_offset": (OMEGA_B_C_PRED / OMEGA_B_C_PLANCK - 1.0) * 1e6,
    }
