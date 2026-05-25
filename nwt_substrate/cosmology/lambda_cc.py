"""
Cosmological constant ρ_Λ/M_Pl⁴ = (m_e/M_Pl)⁴·α¹⁶·h_Cox from the substrate.

The vacuum energy density is fixed by the substrate algebra at zero free
parameters (Stage-7 Coleman-Weinberg vacuum-energy scale):

    ρ_Λ / M_Pl⁴ = (m_e/M_Pl)⁴ · α¹⁶ · h_Cox  ≈ 1.19×10⁻¹²³

against the observed ρ_Λ/M_Pl⁴ ≈ 1.20×10⁻¹²³: **pred/obs ≈ 0.993 (−0.7%)**.
Matching the 10⁻¹²³ scale to ~3 significant figures from substrate primitives is
the headline — cf. the naive QFT vacuum-energy estimate, off by ~10¹²⁰.

Factors:
  (m_e/M_Pl)⁴ — the substrate Compton-mass⁴ (Coleman-Weinberg 1-loop scale).
                m_e/M_Pl is taken from the K₇ Wilson amplitude
                (:func:`nwt_substrate.isa.k7_wilson_amplitude`, Paper 17, NNLO:
                (8/7)·α^(21/2)·(1+α/7+(21/8)α²)) — i.e. m_e is *derived from M_Pl*,
                per the M_Pl-as-cosmogenic-mass-anchor framing.
  α¹⁶        — the Stage-7 substrate vacuum-energy exponent.
  h_Cox = 6  — Coxeter number of so(7) = B₃ (Weyl-group multiplicity;
               |W(B₃)|=48, h_Cox=6, rank=3).

A ~0.5–0.7% slack is inherent (the difference between (m_e/M_Pl)² and the
leading-order α²¹), documented in the Stage-7 Λ_cc derivation.

Provides:
  - RHO_LAMBDA_PRED   : substrate prediction for ρ_Λ/M_Pl⁴ (float)
  - RHO_LAMBDA_OBS    : observed reference value
  - M_E_OVER_M_PL     : substrate m_e/M_Pl (Wilson NNLO)
  - LAMBDA_EXPONENT   : the α-exponent (16)
  - lambda_cc()       : the prediction (ρ_Λ/M_Pl⁴)
  - summary()         : prediction vs observation (dict)

All values verified by tests in tests/test_lambda_cc.py.
"""
from __future__ import annotations

from nwt_substrate.isa.constants import ALPHA_SUBSTRATE, H_COXETER_SO7
from nwt_substrate.isa import k7_wilson_amplitude

__all__ = [
    "M_E_OVER_M_PL",
    "LAMBDA_EXPONENT",
    "RHO_LAMBDA_PRED",
    "RHO_LAMBDA_OBS",
    "lambda_cc",
    "summary",
]

#: α-exponent of the Stage-7 vacuum-energy scale.
LAMBDA_EXPONENT: int = 16

#: m_e/M_Pl from the K₇ Wilson amplitude (Paper 17, NNLO) — m_e derived from M_Pl.
M_E_OVER_M_PL: float = k7_wilson_amplitude(ALPHA_SUBSTRATE, order="NNLO")

#: Substrate prediction ρ_Λ/M_Pl⁴ = (m_e/M_Pl)⁴ · α¹⁶ · h_Cox.
RHO_LAMBDA_PRED: float = (
    M_E_OVER_M_PL**4 * ALPHA_SUBSTRATE**LAMBDA_EXPONENT * H_COXETER_SO7
)

#: Observed vacuum energy density in Planck units, ρ_Λ/M_Pl⁴ (Planck 2018 Λ).
RHO_LAMBDA_OBS: float = 1.20e-123


def lambda_cc() -> float:
    """Substrate vacuum energy density ρ_Λ/M_Pl⁴ = (m_e/M_Pl)⁴·α¹⁶·h_Cox."""
    return RHO_LAMBDA_PRED


def summary() -> dict:
    """Prediction vs observation, as a dict (numbers for tables/papers)."""
    return {
        "rho_lambda_pred": RHO_LAMBDA_PRED,
        "formula": "(m_e/M_Pl)**4 * alpha**16 * H_COXETER_SO7",
        "m_e_over_M_Pl": M_E_OVER_M_PL,
        "alpha_exponent": LAMBDA_EXPONENT,
        "h_cox": H_COXETER_SO7,
        "rho_lambda_obs": RHO_LAMBDA_OBS,
        "pred_over_obs": RHO_LAMBDA_PRED / RHO_LAMBDA_OBS,
    }
