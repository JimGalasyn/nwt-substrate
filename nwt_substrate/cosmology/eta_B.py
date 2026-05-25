"""
Baryon asymmetry η_B = 3α⁴/14 from the substrate (Stage-7 K₇⊗K₈ partition).

The baryon-to-photon ratio is fixed by the substrate algebra at zero free
parameters:

    η_B = RANK_SO7 · α⁴ / (2·N_VERTICES_K7) = 3·α⁴/14 ≈ 6.08×10⁻¹⁰

against the observed η_B = 6.10×10⁻¹⁰ (Planck 2018) / 6.14×10⁻¹⁰ (BBN):
**pred/obs ≈ 0.996 (−0.4%)**.

The two integer factors each have *two* substrate readings (they coincide):

  numerator  3  = RANK_SO7 = N_generations = the 3 K₇ X-stabilizers (Cartan
                  generators of so(7)); the Z₃ σ-orbit count.
  denominator 14 = 2·N_VERTICES_K7
                 = dim(G₂) = dim Aut(𝕆)                     [NWT reading]
                 = the 14 single-Fano-point syndromes        [VV reading]
                   (Steane [[7,1,3]] stabilizer cosets; GL(3,2)/12).
  The bridge between the readings is 2·N_VERTICES_K7 = 14 = dim(G₂).

Two independent derivation routes converge on this value:

  * **VV — K₇⊗K₈ bridge partition.** The cosmogenic H⊗7 polarity-selection
    event on the K₇ Steane code commits one logical-coset outcome per event;
    8 √α traversals × 4 qubits per X-stabilizer give the α⁴ amplitude, ×3
    X-stabilizers (Z₃ σ-orbit), /14 single-Fano syndromes (Paper 22).
  * **NWT — Murasugi torus-knot.** An independent (3/14)α⁴ from the carrier
    knot's Murasugi/Jones structure (`analysis/eta_B_substrate_derivation.py`).

The sign sgn(η_B) = sgn(J_parent) (matter over antimatter) is fixed by the
parent black hole's spin via the bridge pointer-basis chain; this module
gives the magnitude.

Provides:
  - ETA_B_PRED          : substrate prediction (float)
  - ETA_B_PLANCK / _BBN : observed reference values
  - eta_B()             : the prediction
  - summary()           : prediction vs observation table (dict)

All values verified by tests in tests/test_eta_B.py.
"""
from __future__ import annotations

from nwt_substrate.isa.constants import (
    ALPHA_SUBSTRATE,
    RANK_SO7,
    N_VERTICES_K7,
)

__all__ = [
    "ETA_B_PRED",
    "ETA_B_PLANCK",
    "ETA_B_BBN",
    "ETA_B_DENOM",
    "eta_B",
    "summary",
]

# denominator: 14 = 2·N_VERTICES_K7 = dim(G₂) = # single-Fano syndromes
ETA_B_DENOM: int = 2 * N_VERTICES_K7

#: Substrate prediction η_B = RANK_SO7·α⁴ / (2·N_VERTICES_K7) = 3α⁴/14.
ETA_B_PRED: float = RANK_SO7 * ALPHA_SUBSTRATE**4 / ETA_B_DENOM

#: Observed baryon-to-photon ratio (Planck 2018 TT,TE,EE+lowE+lensing).
ETA_B_PLANCK: float = 6.10e-10
#: Observed baryon-to-photon ratio (BBN, light-element abundances).
ETA_B_BBN: float = 6.14e-10


def eta_B() -> float:
    """Substrate baryon asymmetry η_B = 3α⁴/14."""
    return ETA_B_PRED


def summary() -> dict:
    """Prediction vs observation, as a dict (numbers for tables/papers)."""
    return {
        "eta_B_pred": ETA_B_PRED,
        "formula": "RANK_SO7 * alpha**4 / (2*N_VERTICES_K7) = 3*alpha**4/14",
        "numerator": RANK_SO7,
        "denominator": ETA_B_DENOM,
        "eta_B_planck": ETA_B_PLANCK,
        "eta_B_bbn": ETA_B_BBN,
        "pred_over_planck": ETA_B_PRED / ETA_B_PLANCK,
        "pred_over_bbn": ETA_B_PRED / ETA_B_BBN,
    }
