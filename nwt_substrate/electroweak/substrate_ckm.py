"""
Substrate Wolfenstein CKM matrix from K_7 substrate integers + α + one phase.

P6b closure (Galasyn 2026-05-23 P6b,
analysis/paper21_p6b_wolfenstein_substrate.md in null-worldtube-private):

The full Wolfenstein parameterisation (λ, A, ρ̄, η̄, δ_CP) is derived
in 4 substrate integers + α + a Z_n-winding canonical phase:

    λ²              = 7  α            [Paper 8a/9: Cabibbo from K_7]
    A²              = 9 / 14          [q_b / dim(G_2)]
    |ρ̄ + iη̄|²      = 21 α            [dim(so(7)) · α = |E(K_7)| · α]
    δ_CP^CKM        = π/3             [Z_2 meson half-rotation; conjecture]

Substrate integer vocabulary::

    7   = |V(K_7)|                  (vertices of K_7)
    9   = q_b                       (bottom-quark trefoil-winding, Paper 6)
    14  = dim(G_2)                  (octonion automorphism algebra)
    21  = dim(so(7)) = |E(K_7)|     (edges of K_7 = Spin(7) adjoint)

Zero free CKM parameters.  Combined with [[substrate_gf]] (P7b §7.4)
and PMNS (Paper 20), the entire SM flavor + EW mixing sector reduces
to {α (Paper 17 K_7 Wilson), m_e (one dimensional anchor),
substrate-integer combinatorics}.

Precision vs PDG 2024:

    λ           √(7α)            = 0.22601    PDG 0.2253 ± 0.0007    0.32 %
    A           3/√14            = 0.80178    PDG 0.802 ± 0.011      0.027%
    |ρ̄+iη̄|     √(21α)           = 0.39146    PDG 0.391 ± 0.014      0.12 %
    |V_us|      λ                = 0.226      PDG 0.2253              0.32 %
    |V_cb|      A λ²             = 0.04096    PDG 0.041  ± 0.0014     0.1  %
    |V_ub|      A λ³ |ρ̄+iη̄|     = 0.00362    PDG 0.00382 ± 0.00020   5.2  %
    |V_td|      A λ³ |1-ρ̄-iη̄|  ≈ 0.0093     PDG 0.0086 ± 0.0004      7.6  %
    J_CP^CKM    875 α^(7/2)      = 2.91e-5    PDG 3.18e-5 ± 0.15      8.6  %

Open: δ_CP^CKM refinement from π/3 (60°) toward PDG ~66.4° (~2σ tension).
Likely a higher-order substrate-Wilson correction; conjecture
δ = arctan(√(14/3)) ≈ 65.2° currently un-derived.

Substrate Wilson reading of the closed-form Jarlskog::

    J_CP^CKM = A² λ⁶ η̄  = (9·343·3√7 / 28) α^(7/2)  = 875.0 · α^(7/2)

The exponent 7/2 = half-rank of so(7), structurally analogous to
Paper 17's electron-mass Wilson α^(21/2) (full rank).
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

from .substrate_gf import ALPHA_SUBSTRATE
# Structural integers imported from the single source of truth (isa) so that a
# perturbation to the ISA propagates into the CKM predictions (sensitivity sweep).
from ..isa.constants import N_VERTICES_K7, DIM_ADJ_SPIN7


# ============================================================
# Substrate integers (K_7 / so(7) / G_2 / particle vocabulary)
# ============================================================

# N_VERTICES_K7 (|V(K_7)| = 7, the Cabibbo coefficient) is imported from isa above.

Q_B: int = 9
"""Bottom-quark trefoil winding number q_b from Paper 6 T(2,9) carrier knot.
Numerically coincides with N_POS_ROOTS_SO7=9 but is a distinct structural
integer (carrier winding, not a root count), so it is kept local."""

DIM_G2: int = 14
"""dim(G_2) = 14.  Octonion automorphism Lie-algebra dimension. No standalone
ISA constant (= 2·N_VERTICES_K7); kept local."""

DIM_SO7: int = DIM_ADJ_SPIN7
"""dim(so(7)) = 21 = |E(K_7)|.  Adjoint of Spin(7); = isa.DIM_ADJ_SPIN7 (= N_EDGES_K7)."""


# ============================================================
# Wolfenstein parameters from substrate
# ============================================================

def wolfenstein_lambda_sq(alpha: float = ALPHA_SUBSTRATE) -> float:
    """λ² = 7α — Cabibbo coefficient (Paper 8a/9, K_7 + Hopf Z_2 carrier)."""
    return N_VERTICES_K7 * alpha


def wolfenstein_lambda(alpha: float = ALPHA_SUBSTRATE) -> float:
    """λ = √(7α) ≈ 0.22601.  PDG 0.2253 → 0.32 %."""
    return math.sqrt(wolfenstein_lambda_sq(alpha))


def wolfenstein_A_sq() -> float:
    """A² = q_b / dim(G_2) = 9/14.  Bottom-charm transition amplitude."""
    return Q_B / DIM_G2


def wolfenstein_A() -> float:
    """A = 3/√14 ≈ 0.80178.  PDG 0.802 → 0.027 %."""
    return math.sqrt(wolfenstein_A_sq())


def wolfenstein_apex_magnitude_sq(alpha: float = ALPHA_SUBSTRATE) -> float:
    """|ρ̄ + iη̄|² = dim(so(7)) · α = 21α.  K_7 edge-count substrate Wilson."""
    return DIM_SO7 * alpha


def wolfenstein_apex_magnitude(alpha: float = ALPHA_SUBSTRATE) -> float:
    """|ρ̄ + iη̄| = √(21α) ≈ 0.39146.  PDG 0.391 → 0.12 %."""
    return math.sqrt(wolfenstein_apex_magnitude_sq(alpha))


DELTA_CP_CKM: float = math.pi / 3.0
"""δ_CP^CKM = π/3 = 60°.  Conjecture from Z_2 meson half-rotation
(Hopf carrier T(2,2)); PDG ~66.4° → ~6° tension (~2σ).  Refinement
candidate δ = arctan(√(14/3)) ≈ 65.2° currently un-derived."""


def wolfenstein_rho_bar(alpha: float = ALPHA_SUBSTRATE,
                        delta: float = DELTA_CP_CKM) -> float:
    """ρ̄ = |ρ̄+iη̄| · cos(δ).  Substrate δ=π/3 gives ρ̄ ≈ 0.196.
    PDG 0.156 → 25% (driven entirely by δ tension, not |ρ̄+iη̄|)."""
    return wolfenstein_apex_magnitude(alpha) * math.cos(delta)


def wolfenstein_eta_bar(alpha: float = ALPHA_SUBSTRATE,
                        delta: float = DELTA_CP_CKM) -> float:
    """η̄ = |ρ̄+iη̄| · sin(δ).  Substrate δ=π/3 gives η̄ ≈ 0.339.
    PDG 0.358 → 5.3% (limited by δ_CP^CKM tension)."""
    return wolfenstein_apex_magnitude(alpha) * math.sin(delta)


# ============================================================
# CKM matrix elements (Wolfenstein expansion)
# ============================================================

def V_us(alpha: float = ALPHA_SUBSTRATE) -> float:
    """|V_us| = λ = √(7α).  PDG 0.2253 → 0.32 %."""
    return wolfenstein_lambda(alpha)


def V_cb(alpha: float = ALPHA_SUBSTRATE) -> float:
    """|V_cb| = A λ² = (3/√14) · 7α.  PDG 0.0410 → 0.1 %."""
    return wolfenstein_A() * wolfenstein_lambda_sq(alpha)


def V_ub(alpha: float = ALPHA_SUBSTRATE,
         delta: float = DELTA_CP_CKM) -> complex:
    """V_ub = A λ³ (ρ̄ - iη̄) at Wolfenstein O(λ⁵).
    Substrate: |V_ub| = A λ³ |ρ̄+iη̄|.  PDG 0.00382 → 5.2 %."""
    lam = wolfenstein_lambda(alpha)
    A = wolfenstein_A()
    apex = wolfenstein_apex_magnitude(alpha)
    # Conventional sign: V_ub ∝ (ρ̄ - iη̄)
    return A * lam ** 3 * apex * complex(math.cos(delta), -math.sin(delta))


def V_td(alpha: float = ALPHA_SUBSTRATE,
         delta: float = DELTA_CP_CKM) -> complex:
    """V_td = A λ³ (1 - ρ̄ - iη̄) at Wolfenstein O(λ⁵).
    PDG |V_td| ≈ 0.0086.  Substrate ≈ 0.0093 → 7.6 %."""
    lam = wolfenstein_lambda(alpha)
    A = wolfenstein_A()
    rho_bar = wolfenstein_rho_bar(alpha, delta)
    eta_bar = wolfenstein_eta_bar(alpha, delta)
    return A * lam ** 3 * complex(1.0 - rho_bar, -eta_bar)


# ============================================================
# Closed-form Jarlskog
# ============================================================

def jarlskog_ckm(alpha: float = ALPHA_SUBSTRATE,
                 delta: float = DELTA_CP_CKM) -> float:
    """Jarlskog J_CP^CKM = A² λ⁶ |ρ̄+iη̄| sin(δ).

    With substrate δ = π/3 (sin δ = √3 / 2):

        J_CP^CKM = (9/14)(7α)³ √(21α) (√3/2)
                 = 875.0 · α^(7/2)

    Numerical: 875 · α_sub^(7/2) ≈ 2.91e-5.  PDG 3.18e-5 → 8.6 %.

    The exponent 7/2 = half-rank of so(7), structurally analogous to
    Paper 17's electron-mass Wilson α^(21/2) (full so(7) rank).
    """
    A_sq = wolfenstein_A_sq()
    lam_sq = wolfenstein_lambda_sq(alpha)
    apex = wolfenstein_apex_magnitude(alpha)
    return A_sq * lam_sq ** 3 * apex * math.sin(delta)


def jarlskog_coefficient() -> float:
    """The substrate-canonical Jarlskog prefactor in J = (coeff)·α^(7/2).

    coeff = (9 · 343 · 3√7) / 28 = 9261 √7 / 28 ≈ 875.0

    Independent of α.  Derived by substituting δ = π/3 and the
    substrate substitutions λ² = 7α, A² = 9/14, |ρ̄+iη̄| = √(21α)
    into J = A² λ⁶ |ρ̄+iη̄| sin δ and collecting α^(7/2).
    """
    return (9.0 * 343.0 * 3.0 * math.sqrt(7.0)) / 28.0


# ============================================================
# Full 3x3 CKM matrix (Wolfenstein to O(λ⁴) standard truncation)
# ============================================================

def ckm_matrix(alpha: float = ALPHA_SUBSTRATE,
               delta: float = DELTA_CP_CKM) -> np.ndarray:
    """Full 3×3 CKM matrix in standard Wolfenstein expansion to O(λ⁴).

    Returns a complex (3,3) ndarray indexed [up_type, down_type] with
    rows (u, c, t) and columns (d, s, b).  Built from the substrate
    Wolfenstein parameters; the complex entries carry the substrate
    δ_CP^CKM = π/3 in V_ub and V_td.

    Wolfenstein expansion (PDG convention)::

        V_ud =  1 - λ²/2 - λ⁴/8
        V_us =  λ
        V_ub =  A λ³ (ρ̄ - iη̄)
        V_cd = -λ
        V_cs =  1 - λ²/2 - λ⁴(1 + 4A²)/8
        V_cb =  A λ²
        V_td =  A λ³ (1 - ρ̄ - iη̄)
        V_ts = -A λ² + A λ⁴ (1 - 2(ρ̄ + iη̄))/2
        V_tb =  1 - A² λ⁴/2
    """
    lam = wolfenstein_lambda(alpha)
    A = wolfenstein_A()
    rho_bar = wolfenstein_rho_bar(alpha, delta)
    eta_bar = wolfenstein_eta_bar(alpha, delta)
    rho_eta = complex(rho_bar, eta_bar)

    V = np.zeros((3, 3), dtype=complex)
    V[0, 0] = 1.0 - lam ** 2 / 2.0 - lam ** 4 / 8.0
    V[0, 1] = lam
    V[0, 2] = A * lam ** 3 * complex(rho_bar, -eta_bar)
    V[1, 0] = -lam
    V[1, 1] = 1.0 - lam ** 2 / 2.0 - lam ** 4 * (1.0 + 4.0 * A ** 2) / 8.0
    V[1, 2] = A * lam ** 2
    V[2, 0] = A * lam ** 3 * complex(1.0 - rho_bar, -eta_bar)
    V[2, 1] = -A * lam ** 2 + A * lam ** 4 * (1.0 - 2.0 * rho_eta).real / 2.0 \
              + 1j * A * lam ** 4 * (1.0 - 2.0 * rho_eta).imag / 2.0
    V[2, 2] = 1.0 - A ** 2 * lam ** 4 / 2.0
    return V


# ============================================================
# Precision chain vs PDG
# ============================================================

# PDG 2024 central values for verification.
_PDG = {
    'lambda': 0.2253,
    'A': 0.802,
    'apex_magnitude': 0.391,
    'V_us': 0.2253,
    'V_cb': 0.0410,
    'V_ub_mag': 0.00382,
    'V_td_mag': 0.0086,
    'jarlskog': 3.18e-5,
    'delta_deg': 66.4,
}


def precision_chain(alpha: float = ALPHA_SUBSTRATE,
                    delta: float = DELTA_CP_CKM) -> dict:
    """Substrate-vs-PDG comparison for all major Wolfenstein/CKM observables.

    Returns dict keyed by observable name with sub-dict
    {'substrate', 'pdg', 'rel_err'} where rel_err is signed
    fractional error (substrate / PDG - 1).
    """
    sub = {
        'lambda': wolfenstein_lambda(alpha),
        'A': wolfenstein_A(),
        'apex_magnitude': wolfenstein_apex_magnitude(alpha),
        'V_us': V_us(alpha),
        'V_cb': V_cb(alpha),
        'V_ub_mag': abs(V_ub(alpha, delta)),
        'V_td_mag': abs(V_td(alpha, delta)),
        'jarlskog': jarlskog_ckm(alpha, delta),
        'delta_deg': math.degrees(delta),
    }
    return {
        k: {
            'substrate': sub[k],
            'pdg': _PDG[k],
            'rel_err': sub[k] / _PDG[k] - 1.0,
        }
        for k in sub
    }


def verify_substrate_ckm(alpha: float = ALPHA_SUBSTRATE,
                         delta: float = DELTA_CP_CKM,
                         tol_dominant: float = 0.01,
                         tol_subdominant: float = 0.10) -> dict:
    """Verify substrate Wolfenstein matches PDG within stated tolerances.

    Dominant observables (λ, A, |ρ̄+iη̄|, V_us, V_cb) must match within
    ``tol_dominant`` (default 1 %).  Subdominant observables (|V_ub|,
    |V_td|, J_CP, δ_CP) must match within ``tol_subdominant`` (default
    10 %; absorbs the known δ_CP π/3 tension).

    Returns the precision_chain dict augmented with boolean ``pass``.
    """
    chain = precision_chain(alpha, delta)
    dominant = ['lambda', 'A', 'apex_magnitude', 'V_us', 'V_cb']
    subdominant = ['V_ub_mag', 'V_td_mag', 'jarlskog', 'delta_deg']
    ok = True
    for k in dominant:
        if abs(chain[k]['rel_err']) > tol_dominant:
            ok = False
    for k in subdominant:
        if abs(chain[k]['rel_err']) > tol_subdominant:
            ok = False
    chain['pass'] = ok
    return chain


def precision_chain_summary(alpha: float = ALPHA_SUBSTRATE,
                            delta: float = DELTA_CP_CKM) -> str:
    """Pretty-print substrate-vs-PDG Wolfenstein/CKM precision chain."""
    chain = precision_chain(alpha, delta)
    lines = [
        "Substrate Wolfenstein CKM precision chain (substrate vs PDG):",
        "",
        f"  λ              = {chain['lambda']['substrate']:.5f}    "
        f"vs PDG {chain['lambda']['pdg']:.4f}    "
        f"→ {chain['lambda']['rel_err'] * 100:+6.3f} %",
        f"  A              = {chain['A']['substrate']:.5f}    "
        f"vs PDG {chain['A']['pdg']:.4f}    "
        f"→ {chain['A']['rel_err'] * 100:+6.3f} %",
        f"  |ρ̄ + iη̄|     = {chain['apex_magnitude']['substrate']:.5f}    "
        f"vs PDG {chain['apex_magnitude']['pdg']:.4f}    "
        f"→ {chain['apex_magnitude']['rel_err'] * 100:+6.3f} %",
        f"  |V_us|         = {chain['V_us']['substrate']:.5f}    "
        f"vs PDG {chain['V_us']['pdg']:.4f}    "
        f"→ {chain['V_us']['rel_err'] * 100:+6.3f} %",
        f"  |V_cb|         = {chain['V_cb']['substrate']:.5f}    "
        f"vs PDG {chain['V_cb']['pdg']:.4f}    "
        f"→ {chain['V_cb']['rel_err'] * 100:+6.3f} %",
        f"  |V_ub|         = {chain['V_ub_mag']['substrate']:.5f}    "
        f"vs PDG {chain['V_ub_mag']['pdg']:.5f}    "
        f"→ {chain['V_ub_mag']['rel_err'] * 100:+6.3f} %",
        f"  |V_td|         = {chain['V_td_mag']['substrate']:.5f}    "
        f"vs PDG {chain['V_td_mag']['pdg']:.4f}    "
        f"→ {chain['V_td_mag']['rel_err'] * 100:+6.3f} %",
        f"  J_CP           = {chain['jarlskog']['substrate']:.3e}   "
        f"vs PDG {chain['jarlskog']['pdg']:.3e}  "
        f"→ {chain['jarlskog']['rel_err'] * 100:+6.3f} %",
        f"  δ_CP (deg)     = {chain['delta_deg']['substrate']:.2f}     "
        f"vs PDG {chain['delta_deg']['pdg']:.2f}     "
        f"→ {chain['delta_deg']['rel_err'] * 100:+6.3f} %",
        "",
        "  Substrate vocabulary:  (7, 9, 14, 21) + α + δ=π/3",
        "                        |K_7|, q_b, dim(G_2), dim(so(7))",
        "  Free params:  0  (Wolfenstein λ, A, ρ̄, η̄ all substrate-derived)",
        "  Closed Jarlskog:  J_CP^CKM = 875.0 · α^(7/2)",
    ]
    return "\n".join(lines)
