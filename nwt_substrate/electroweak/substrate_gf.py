"""
Substrate Fermi constant G_F from K_7 Wilson amplitude alone.

Closed-form derivation (Galasyn 2026-05-23 P7b §7.4,
analysis/paper21_p7b_GF_substrate.md in null-worldtube-private):

    v_EW^sub = (25 m_e / α²) · (1 + 25 α / (4 √3))      [Paper 13 / 4-formula]
    G_F^sub  = 1 / (√2 v_EW²)                            [SM tree-level identity]
             = α⁴ / (625 √2 m_e² (1 + 25 α / (4 √3))²)  [closed form]

Substrate inputs: α (Paper 17 K_7 Wilson) and the mass scale m_e (itself
derived from the cosmogenic anchor M_Pl; see isa.constants.M_PL_GEV).
Zero free EW parameters.  Combined with P6b (CKM in 4 substrate integers
+ α at 0.03 %), the SM EW + flavor sector reduces to a single mass scale,
which is the cosmogenic M_Pl (m_e = (8/7)·α^(21/2)·(NNLO)·M_Pl).

Precision chain vs PDG (substrate α = 1/(25π√3 + 1)):

    α      substrate 1/137.034952  vs PDG 1/137.035999  → +7.6 ppm
    v_EW   substrate 246.2128 GeV  vs PDG 246.21965 GeV → −27.7 ppm
    G_F    substrate 1.16644e-5    vs PDG 1.16638e-5    → +55 ppm  (tree)
    τ_μ    (G_F + phase-space F + Sirlin Δq)            → −105 ppm  (full chain)

Substrate Wilson reading:

    α⁴      = 4-fermion contact vertex (4 K_7 edges traversed by W)
    625     = q_cinq⁴ = 5⁴ = cinquefoil scale-setting walk factor⁴
    25 / 4  = q_cinq² / C_A²(SU(2))
              (recurs in BOTH v_EW NLO correction AND Sirlin Δq coefficient
              α/(2π) × (25/4 - π²) — structural-consistency check across
              substrate Wilson amplitude levels)
"""

from __future__ import annotations

import math

from .constants import G_F_GEV  # PDG value for comparison
from ..isa.constants import ALPHA_SUBSTRATE, H_V_SO7  # canonical substrate α + h_v(so7)=5=q_cinq (single source of truth)


# ============================================================
# Substrate inputs
# ============================================================

# ALPHA_SUBSTRATE is imported from isa.constants (single source of truth);
# α = 1/(25π√3 + 1), Paper 17 K_7 Wilson amplitude, +7.6 ppm vs PDG.


M_E_GEV: float = 0.5109989461e-3
"""Electron mass in GeV (PDG), the mass scale entering the G_F chain here.
It is NOT itself fundamental: m_e = (8/7)·α^(21/2)·(NNLO bracket)·M_Pl derives
from the cosmogenic mass anchor M_Pl (see isa.constants.M_PL_GEV)."""


V_EW_PDG_GEV: float = 246.21965
"""PDG v_EW = (sqrt(2) G_F)^(-1/2) from muon lifetime, GeV."""


# ============================================================
# Substrate v_EW
# ============================================================

def v_ew_substrate(alpha: float = ALPHA_SUBSTRATE,
                   m_e: float = M_E_GEV) -> float:
    """Substrate Higgs vacuum-expectation value, GeV.

    Formula::

        v_EW = (25 m_e / α²) · (1 + 25 α / (4 √3))

    Two-factor decomposition::

        Leading:  25 m_e / α²       = cinquefoil scale-setting (Paper 13)
        NLO:      1 + 25 α / (4 √3) = cinquefoil / Casimir correction
                                       (25/4 = q_cinq² / C_A²(SU(2)),
                                        √3   = √C_A(SU(3)))

    Matches PDG v_EW = 246.21965 GeV at ~28 ppm with substrate α.
    """
    sqrt3 = math.sqrt(3.0)
    leading = H_V_SO7 ** 2 * m_e / alpha ** 2          # 25 = q_cinq² = H_V_SO7²
    nlo = 1.0 + H_V_SO7 ** 2 * alpha / (4.0 * sqrt3)    # 25/4 = q_cinq²/C_A²(SU2)
    return leading * nlo


# ============================================================
# Substrate G_F (closed form)
# ============================================================

def fermi_constant_substrate(alpha: float = ALPHA_SUBSTRATE,
                             m_e: float = M_E_GEV) -> float:
    """Substrate Fermi constant G_F, GeV⁻², closed form.

    Formula::

        G_F = α⁴ / (625 √2 m_e² (1 + 25 α / (4 √3))²)

    Derivation: SM tree-level identity G_F = 1/(√2 v_EW²) combined with the
    substrate v_EW formula from :func:`v_ew_substrate`.  Pure substrate
    inputs (α, m_e); zero free EW parameters.

    Matches PDG G_F = 1.1663787e-5 GeV⁻² at +55 ppm tree-level (inherited
    from substrate v_EW at -28 ppm via factor 2 since G_F ∝ v⁻²).  Full
    muon-lifetime chain with phase-space F(m_e²/m_μ²) and Sirlin Δq QED
    correction matches PDG τ_μ at -105 ppm.
    """
    sqrt2 = math.sqrt(2.0)
    sqrt3 = math.sqrt(3.0)
    nlo = 1.0 + H_V_SO7 ** 2 * alpha / (4.0 * sqrt3)    # 25/4 = q_cinq²/C_A²(SU2)
    return alpha ** 4 / (H_V_SO7 ** 4 * sqrt2 * m_e ** 2 * nlo ** 2)   # 625 = q_cinq⁴ = H_V_SO7⁴


def fermi_constant_from_vev(v_ew: float) -> float:
    """SM tree-level identity G_F = 1/(√2 v_EW²).

    Provided for direct comparison with the closed-form
    :func:`fermi_constant_substrate`; the two agree identically when
    v_ew == v_ew_substrate(alpha, m_e).
    """
    return 1.0 / (math.sqrt(2.0) * v_ew ** 2)


# ============================================================
# Precision chain + verification
# ============================================================

_ALPHA_PDG: float = 1.0 / 137.035999084


def precision_chain(alpha: float = ALPHA_SUBSTRATE,
                    m_e: float = M_E_GEV) -> dict:
    """Compute substrate-vs-PDG precision chain for the G_F derivation.

    Returns a dict with keys ``alpha``, ``v_ew``, ``g_f``, each mapping to
    a sub-dict ``{'substrate': float, 'pdg': float, 'ppm': float}``.

    The ppm field is signed: positive means substrate > PDG.
    """
    v_sub = v_ew_substrate(alpha, m_e)
    gf_sub = fermi_constant_substrate(alpha, m_e)
    return {
        'alpha': {
            'substrate': alpha,
            'pdg': _ALPHA_PDG,
            'ppm': (alpha / _ALPHA_PDG - 1.0) * 1e6,
        },
        'v_ew': {
            'substrate': v_sub,
            'pdg': V_EW_PDG_GEV,
            'ppm': (v_sub / V_EW_PDG_GEV - 1.0) * 1e6,
        },
        'g_f': {
            'substrate': gf_sub,
            'pdg': G_F_GEV,
            'ppm': (gf_sub / G_F_GEV - 1.0) * 1e6,
        },
    }


def verify_substrate_gf(alpha: float = ALPHA_SUBSTRATE,
                        m_e: float = M_E_GEV,
                        ppm_tol: float = 100.0) -> dict:
    """Verify substrate G_F derivation matches PDG within ppm_tol.

    Default tolerance 100 ppm accommodates the +55 ppm tree-level gap with
    headroom for floating-point noise.  Returns the precision_chain dict
    augmented with a boolean ``pass`` key.
    """
    chain = precision_chain(alpha, m_e)
    chain['pass'] = abs(chain['g_f']['ppm']) <= ppm_tol
    return chain


def precision_chain_summary(alpha: float = ALPHA_SUBSTRATE,
                            m_e: float = M_E_GEV) -> str:
    """Pretty-print the substrate G_F precision chain (substrate vs PDG)."""
    chain = precision_chain(alpha, m_e)
    lines = [
        "Substrate G_F precision chain (substrate vs PDG):",
        "",
        f"  α          = {chain['alpha']['substrate']:.10f}  "
        f"vs PDG {chain['alpha']['pdg']:.10f}  → {chain['alpha']['ppm']:+7.1f} ppm",
        f"  v_EW (GeV) = {chain['v_ew']['substrate']:.5f}     "
        f"vs PDG {chain['v_ew']['pdg']:.5f}    → {chain['v_ew']['ppm']:+7.1f} ppm",
        f"  G_F (GeV⁻²) = {chain['g_f']['substrate']:.5e}    "
        f"vs PDG {chain['g_f']['pdg']:.5e}   → {chain['g_f']['ppm']:+7.1f} ppm",
        "",
        "  Closed form:  G_F = α⁴ / (625 √2 m_e² (1 + 25α/(4√3))²)",
        "  Inputs:       α (K_7 Wilson, Paper 17) + m_e (scale anchor)",
        "  Free params:  0 (EW sector); m_e is the dimensional anchor",
    ]
    return "\n".join(lines)
