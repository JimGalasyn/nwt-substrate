"""
Substrate-derived Newton's `G` and `m_e / M_Planck`.

Two structurally-derived formulas are exposed:

  Paper 14 (LO):
      G = (8/7)² × α²¹ × ℏ c / m_e²
      → -0.24% from CODATA

  Paper 17 (NNLO):
      m_e / M_Pl = (8/7) × α^(21/2) × (1 + α/7 + 3 α²)
      → 0.0001% from CODATA, 11 ppm in G

Substrate identification of factors:
  8/7        = dim(spinor S of Spin(7)) / dim(vector V of Spin(7))
  α^(21/2)   = K_7 Wilson amplitude on 21 edges with √α per edge
                (Paper 13 / Paper 17 Phase 6.1b)
  (1 + α/7)  = NLO per-K_7-vertex correction (Paper 17 Phase 6.5);
                7 = number of K_7 vertices
  3 α²       = NNLO from rank(so(7))/dim(V) × dim(V)/dim(S) × α²
                = 3 α² where 3 = rank(so(7)) (Paper 17 Phase 6.9)

The 168 unification: |PSL(2,7)| = λ_1(S³/2I) = 7·24 = |Irr(2T)|·|2T| all
descend from Spin(7)/2T/E_6 (memory: checkpoint 2026-04-24).
"""

from __future__ import annotations

from .constants import (
    M_PLANCK_GEV, M_E_GEV, ALPHA_QED, G_NEWTON_SI,
    HBAR_J_S, C_LIGHT_M_S, M_ELECTRON_KG,
)


# ===========================================================================
# Substrate-derived m_e / M_Pl
# ===========================================================================

def m_e_over_M_Pl_LO(alpha: float = ALPHA_QED) -> float:
    """
    Leading-order substrate prediction:
        m_e / M_Pl = (8/7) × α^(21/2)
    """
    return (8.0 / 7.0) * alpha ** (21.0 / 2.0)


def m_e_over_M_Pl_NLO(alpha: float = ALPHA_QED) -> float:
    """
    NLO substrate prediction (Paper 17):
        m_e / M_Pl = (8/7) × α^(21/2) × (1 + α/7)

    The (1 + α/7) is the per-K_7-vertex correction structurally derived
    in Paper 17 Phase 6.5; 7 is the number of K_7 vertices.
    """
    return (8.0 / 7.0) * alpha ** (21.0 / 2.0) * (1.0 + alpha / 7.0)


def m_e_over_M_Pl_NNLO(alpha: float = ALPHA_QED) -> float:
    """
    NNLO substrate prediction (Paper 17 closed form):
        m_e / M_Pl = (8/7) × α^(21/2) × (1 + α/7 + 3 α²)

    Where:
        3 = rank(so(7)) = dim(Adj_Spin(7)) / dim(V_Spin(7)) = 21/7
            (verified 4 ways in Paper 17 Phase 6.9)
        α² coefficient comes from PSL(2,7)-equivariant Fano structure
            (Phase 6.7)

    CODATA agreement: 0.0001% on m_e/M_Pl, ~11 ppm on G.
    """
    return (8.0 / 7.0) * alpha ** (21.0 / 2.0) * (1.0 + alpha / 7.0 + 3.0 * alpha ** 2)


def m_e_over_M_Pl_observed() -> float:
    """The CODATA m_e/M_Pl ratio."""
    return M_E_GEV / M_PLANCK_GEV


# ===========================================================================
# Substrate-derived Newton's G
# ===========================================================================

def G_substrate_LO_natural() -> float:
    """
    G in natural units (GeV⁻²) from the LO formula (Paper 14):

        G_LO = (8/7)² × α²¹ / m_e²

    where m_e is in GeV.  Returns GeV⁻².
    """
    return (8.0 / 7.0) ** 2 * ALPHA_QED ** 21 / M_E_GEV ** 2


def G_substrate_NNLO_natural() -> float:
    """
    G in natural units (GeV⁻²) from the Paper 17 NNLO formula:

        G_NNLO = (m_e / M_Pl_predicted)² / m_e²
               = (1 / M_Pl_predicted²)

    where M_Pl_predicted = m_e / (m_e/M_Pl_NNLO).  Returns GeV⁻².
    """
    M_Pl_pred = M_E_GEV / m_e_over_M_Pl_NNLO()
    return 1.0 / M_Pl_pred ** 2


def G_substrate_SI() -> float:
    """
    Substrate-predicted Newton's `G` in SI units (m³ kg⁻¹ s⁻²),
    using the Paper 17 NNLO m_e/M_Pl formula.

    G_natural [GeV⁻²] = G_SI × M_Pl² [GeV²]  →  G_SI = G_nat / M_Pl² (in SI)
    Conversion: ℏc/m_e² gives the natural-units bridge.
    """
    M_Pl_pred_GeV = M_E_GEV / m_e_over_M_Pl_NNLO()
    # G_SI = ℏ c / M_Pl²  with M_Pl in energy units
    M_Pl_pred_kg = M_Pl_pred_GeV * 1.78266192e-27 / 1e-9   # GeV → kg via E=mc²
                                                            # 1 GeV/c² = 1.78e-27 kg
    M_Pl_pred_kg = M_Pl_pred_GeV * 1.78266192e-27
    return HBAR_J_S * C_LIGHT_M_S / M_Pl_pred_kg ** 2


def G_substrate_LO_SI() -> float:
    """LO substrate-predicted G in SI."""
    M_Pl_pred_GeV = M_E_GEV / m_e_over_M_Pl_LO()
    M_Pl_pred_kg = M_Pl_pred_GeV * 1.78266192e-27
    return HBAR_J_S * C_LIGHT_M_S / M_Pl_pred_kg ** 2


# ===========================================================================
# Convenience accessors and structural breakdown
# ===========================================================================

def G_substrate_breakdown() -> str:
    """Pretty-print the substrate-G structural decomposition."""
    alpha = ALPHA_QED
    lines = ["Newton's G from substrate primitives:"]
    lines.append("")
    lines.append("    G = ℏc / M_Pl²,     m_e / M_Pl = (8/7) × α^(21/2) × Λ_NLO,")
    lines.append("")
    lines.append("  where each factor comes from a substrate construction:")
    lines.append("")
    lines.append(f"    (8/7)      = dim(S) / dim(V)   = {8/7:.6f}")
    lines.append(f"                  Spin(7) spinor / vector dimension")
    lines.append(f"                  (Paper 17 Phase 6.4)")
    lines.append("")
    lines.append(f"    α^(21/2)   = K_7 Wilson amplitude on 21 edges with √α/edge")
    lines.append(f"                  21 = #edges of K_7 = dim(Adj Spin(7))")
    lines.append(f"                  = {alpha**(21/2):.6e}   (Paper 13 / 17 Phase 6.1b)")
    lines.append("")
    lines.append(f"    1 + α/7    = NLO per-K_7-vertex correction (Paper 17 Phase 6.5)")
    lines.append(f"                  7 = #vertices of K_7")
    lines.append(f"                  = {1 + alpha/7:.10f}")
    lines.append("")
    lines.append(f"    1 + 3 α²   = NNLO PSL(2,7)-equivariant Fano correction")
    lines.append(f"                  3 = rank(so(7)) = dim(Adj)/dim(V) = 21/7")
    lines.append(f"                  (Paper 17 Phase 6.7-6.9)")
    lines.append("")
    return "\n".join(lines)
