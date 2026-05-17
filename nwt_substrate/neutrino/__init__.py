"""
nwt.neutrino
============

Neutrino-sector view of the substrate algebra (Paper 20: *The Neutrino
Sector from Spin(8) Triality on K_7 / K_8*). Implements the closed-form
predictions of Paper 20 from the substrate Wilson-amplitude machinery
inherited from Papers 17 and 19, extended to K_8 via the
``Spin(7) ⊂ Spin(8)`` stabilized eighth vertex.

The shim is **entirely substrate-driven**: every numerical prediction
follows from K_7 / K_8 edge counts, ``Spin(7)`` representation
dimensions, and ``α_NWT``, with no fitted parameters. The only
non-substrate inputs are the experimental ``Δm²_21`` and ``Δm²_31``
mass-squared differences, used to extrapolate from the lightest
substrate-predicted active mass ``m_1`` to ``m_2`` and ``m_3``.

Predictions exposed:

- :func:`active_masses`: ``(m_1, m_2, m_3) ≈ (14.84, 17.16, 53.00) meV``
- :func:`sterile_masses`: ``(m_N1, m_N2, m_N3) ≈ (61.3, 70.8, 218.8) MeV``
- :func:`sterile_active_mixing_sq`: ``|U_α4|² ≈ α^(9/2) ≈ 2.4×10⁻¹⁰``
- :func:`pmns_angles_leading_order`: tri-bi-maximal +
  ``θ_13 = arcsin(√(3α))``
- :func:`delta_cp_winding`: ``δ_CP = -2π/3`` from ``π_1(PSU(3))`` Z_3
  winding (Baez Fano convention)
- :func:`jarlskog_invariant`: ``|J_CP| ≈ 0.029``

Quick start::

    import nwt_substrate.neutrino as nu

    nu.active_masses()                       # (1.484e-2, 1.716e-2, 5.300e-2) eV
    nu.sterile_masses()                      # ≈ (6.13e7, 7.08e7, 2.19e8) eV
    nu.sterile_active_mixing_sq()            # ≈ 2.4e-10
    nu.pmns_angles_leading_order()           # PMNSAngles(θ_12, θ_23, θ_13)
    nu.delta_cp_winding()                    # -2π/3
    print(nu.substrate_breakdown())          # Substrate-identity table

The :func:`substrate_breakdown` function prints the full Paper 20
identity table relating substrate edge counts to each prediction.
"""

from __future__ import annotations

import math
from typing import NamedTuple

from nwt_substrate.isa.constants import (
    ALPHA_NWT,
    DIM_S_SPIN7,
    K8_ACTIVE_EDGE_COUNT,
    K8_PARTITION,
    K8_SEESAW_EDGE_DIFFERENCE,
    K8_STERILE_EDGE_COUNT,
    N_EDGES_K8,
    N_VERTICES_K7,
    N_VERTICES_K8,
    NLO_VERTEX_COEFFICIENT,
    RANK_SO7,
)
from nwt_substrate.gravity.constants import M_PLANCK_GEV


# ---------------------------------------------------------------------------
# Experimental anchors (used only for m_2, m_3 extrapolation)
# ---------------------------------------------------------------------------

DELTA_M_SQ_21_eV2 = 7.41e-5
"""Solar mass-squared difference (PDG 2024 / NuFIT 6.0), in eV²."""

DELTA_M_SQ_31_eV2 = 2.51e-3
"""Atmospheric mass-squared difference (PDG 2024 / NuFIT 6.0, normal
hierarchy), in eV²."""


# Planck mass in eV (from gravity.constants in GeV)
_M_PLANCK_eV = M_PLANCK_GEV * 1.0e9


# ---------------------------------------------------------------------------
# Substrate Wilson amplitude — the unified Paper-17/19/20 mass formula
# ---------------------------------------------------------------------------

def wilson_mass_eV(
    N_v: int,
    N_e: int,
    alpha: float = ALPHA_NWT,
    order: str = "NNLO",
) -> float:
    """The unified Wilson amplitude mass formula (Paper 17 → 20).

    ``m = (DIM_S_SPIN7 / N_v) · α^(N_e/2) · (1 + α/N_VERTICES_K7) · (1 + 3α²) · m_Pl``

    Parameters
    ----------
    N_v : int
        Vertex count of the graph the phase soliton lives on. 7 for the
        electron (Paper 17), 8 for active neutrinos and sterile partners
        (Paper 20).
    N_e : int
        Edge count of the same graph (or sub-graph for sterile partners).
        21 for electron, 28 for active neutrino, 19 for sterile.
    alpha : float
        Fine structure constant; defaults to ``ALPHA_NWT`` from
        ``isa.constants``.
    order : {"LO", "NLO", "NNLO"}
        Wilson-amplitude expansion order. NLO adds ``(1 + α/7)``;
        NNLO adds ``(1 + 3α²)``.

    Returns
    -------
    float
        Mass in eV.
    """
    if order not in ("LO", "NLO", "NNLO"):
        raise ValueError(f"order must be LO/NLO/NNLO, got {order!r}")

    prefactor = DIM_S_SPIN7 / N_v
    alpha_power = alpha ** (N_e / 2.0)
    mass_LO_natural = prefactor * alpha_power

    nlo = 1.0 if order == "LO" else (1.0 + alpha * NLO_VERTEX_COEFFICIENT)
    nnlo = 1.0 if order != "NNLO" else (1.0 + 3.0 * alpha * alpha)

    return mass_LO_natural * nlo * nnlo * _M_PLANCK_eV


# ---------------------------------------------------------------------------
# Active neutrino masses
# ---------------------------------------------------------------------------

def active_mass_lightest(order: str = "NNLO") -> float:
    """Lightest active neutrino mass ``m_1`` in eV.

    From the substrate K_8 Wilson amplitude (Paper 19 W3.3-J II,
    Paper 20 §3) with ``(N_v, N_e) = (8, K8_ACTIVE_EDGE_COUNT) = (8, 28)``.
    """
    return wilson_mass_eV(
        N_v=N_VERTICES_K8,
        N_e=K8_ACTIVE_EDGE_COUNT,
        order=order,
    )


def active_masses(order: str = "NNLO") -> tuple:
    """Three active neutrino masses (normal hierarchy) in eV.

    ``m_1`` from the substrate K_8 Wilson amplitude (Paper 19 W3.3-J II);
    ``m_2`` and ``m_3`` from the experimental ``Δm²_21`` and ``Δm²_31``
    via ``m_i = √(m_1² + Δm²_i1)``.
    """
    m1 = active_mass_lightest(order=order)
    m2 = math.sqrt(m1 * m1 + DELTA_M_SQ_21_eV2)
    m3 = math.sqrt(m1 * m1 + DELTA_M_SQ_31_eV2)
    return (m1, m2, m3)


def sum_active_masses(order: str = "NNLO") -> float:
    """``Σ m_ν`` in eV. Compared against the cosmological-constant
    joint bound in Paper 20 §11."""
    return sum(active_masses(order=order))


# ---------------------------------------------------------------------------
# Sterile neutrino masses and active-sterile mixing
# ---------------------------------------------------------------------------

def seesaw_ratio() -> float:
    """Active/sterile seesaw ratio (Paper 20 §6).

    ``m_active / m_sterile = α^(K8_SEESAW_EDGE_DIFFERENCE / 2) = α^(9/2)``,
    where the 9 = ``K8_ACTIVE_EDGE_COUNT - K8_STERILE_EDGE_COUNT`` = 28 − 19
    is the difference between active (full K_8) and sterile (K_8 minus
    the 9 cross-orbit edges of the Z_3 ⊂ G_2 decomposition) edge counts.
    """
    return ALPHA_NWT ** (K8_SEESAW_EDGE_DIFFERENCE / 2.0)


def sterile_masses(order: str = "NNLO") -> tuple:
    """Three sterile neutrino masses ``(m_N1, m_N2, m_N3)`` in eV.

    From the unified Wilson recipe with ``(N_v, N_e) = (8, 19)``:
    ``m_N_i = m_ν_i / seesaw_ratio() = m_ν_i × α^(-9/2) ≈ m_ν_i × 4.13×10⁹``.

    Predicted values: ``(61.3, 70.8, 218.8) MeV``, all inside the νMSM
    window targeted by SHiP, PIONEER, NA62, and DUNE-near.
    """
    active = active_masses(order=order)
    ratio = seesaw_ratio()
    return tuple(m / ratio for m in active)


def sterile_active_mixing_sq() -> float:
    """``|U_α4|² ≈ α^(9/2) ≈ 2.4×10⁻¹⁰`` (Paper 20 §6).

    Universal active-sterile mixing-squared from the same edge-count
    difference that fixes the seesaw ratio. At the edge of current
    bounds; testable by SHiP, PIONEER, NA62-HL, DUNE-near in the
    next decade.
    """
    return ALPHA_NWT ** (K8_SEESAW_EDGE_DIFFERENCE / 2.0)


# ---------------------------------------------------------------------------
# PMNS matrix at leading order from Spin(8) triality
# ---------------------------------------------------------------------------

class PMNSAngles(NamedTuple):
    """PMNS mixing angles in radians.

    Tri-bi-maximal-with-θ_13 from Spin(8) triality (Paper 19 W3.3-K
    and Paper 20 §7). ``θ_13`` comes from the ``Spin(7) ⊂ Spin(8)``
    breaking via the ``RANK_SO7 × α`` substrate identity.
    """

    theta_12: float
    """Solar mixing angle ``arctan(1/√2) ≈ 35.26°``."""

    theta_23: float
    """Atmospheric mixing angle ``π/4 = 45°``."""

    theta_13: float
    """Reactor mixing angle ``arcsin(√(rank · α)) ≈ 8.43°`` with
    ``rank = RANK_SO7 = 3``."""


def pmns_angles_leading_order() -> PMNSAngles:
    """PMNS mixing angles at leading order (Paper 19 W3.3-K, Paper 20 §7).

    NLO corrections (``δθ_12 ≈ -1.8°``, ``δθ_23 ≈ +2.7°``) reproduce
    NuFIT 6.0 central values via the ``α``-suppressed ``U_ℓ`` rotation
    ansatz of Paper 20 §7.6; that first-principles derivation from the
    charged-lepton mass matrix remains open and is not exposed here.
    """
    return PMNSAngles(
        theta_12=math.atan(1.0 / math.sqrt(2.0)),
        theta_23=math.pi / 4.0,
        theta_13=math.asin(math.sqrt(RANK_SO7 * ALPHA_NWT)),
    )


# ---------------------------------------------------------------------------
# CP violation from π_1(PSU(3)) winding
# ---------------------------------------------------------------------------

def delta_cp_winding() -> float:
    """Dirac CP phase ``δ_CP = -2π/3`` (Paper 20 §7.5).

    From the ``Z_3`` winding number of a closed loop in
    ``π_1(PSU(3)) = Z(SU(3)) ≅ Z_3``. The rephasing-invariant content
    of the substrate prediction is ``|δ_CP| = 2π/3 = 120°``;
    the negative sign follows the Baez Fano octonion-product convention.
    """
    return -2.0 * math.pi / 3.0


def jarlskog_invariant() -> float:
    """Jarlskog invariant ``|J_CP| ≈ 0.029`` (Paper 20 §7.5).

    Computed from the PMNS leading-order angles + ``δ_CP = -2π/3``:
    ``J_CP = sin(θ_12) cos(θ_12) · sin(θ_23) cos(θ_23) · sin(θ_13) cos²(θ_13) · sin(δ_CP)``.

    The T2K-only ``J_CP ≈ -0.027`` matches at ~0.1σ; the NuFIT 6.0 global
    fit gives ``J_CP ≈ -0.003`` (from ``δ_CP ≈ 177°``) which lies
    ~1.3σ below the substrate magnitude. DUNE and Hyper-K will
    resolve the discrepancy within the next decade.
    """
    a = pmns_angles_leading_order()
    delta = delta_cp_winding()
    s12, c12 = math.sin(a.theta_12), math.cos(a.theta_12)
    s23, c23 = math.sin(a.theta_23), math.cos(a.theta_23)
    s13, c13 = math.sin(a.theta_13), math.cos(a.theta_13)
    return s12 * c12 * s23 * c23 * s13 * c13 * c13 * math.sin(delta)


# ---------------------------------------------------------------------------
# Substrate identity breakdown
# ---------------------------------------------------------------------------

def substrate_breakdown() -> str:
    """Return the Paper 20 substrate-identity table for the neutrino shim.

    Prints the K_8 edge-decomposition, mass predictions, mixing angles,
    and CP-violation phase together with their substrate origins, so a
    reader can trace any number back to the K_7 / K_8 / Spin(7) / so(7)
    structural constants from which it follows.
    """
    am = active_masses()
    sm = sterile_masses()
    pmns = pmns_angles_leading_order()
    delta = delta_cp_winding()
    j = jarlskog_invariant()
    lines = [
        "nwt.neutrino — Paper 20 substrate-identity breakdown",
        "====================================================",
        f"K_7 vertex count             N_VERTICES_K7 = {N_VERTICES_K7}",
        f"K_8 vertex count             N_VERTICES_K8 = {N_VERTICES_K8}",
        f"K_8 edge count               N_EDGES_K8    = {N_EDGES_K8}",
        f"K_8 Z_3 edge partition       {K8_PARTITION}",
        "                               (6 Lorentz + 3 internal-SU(2)"
        f" + 12 SM-flavor + 1 Higgs + 6 Yukawa = {sum(K8_PARTITION)})",
        f"K_8 active edges             K8_ACTIVE_EDGE_COUNT    = {K8_ACTIVE_EDGE_COUNT}",
        f"K_8 sterile edges            K8_STERILE_EDGE_COUNT   = {K8_STERILE_EDGE_COUNT}",
        f"Seesaw edge difference       K8_SEESAW_EDGE_DIFFERENCE = "
        f"{K8_SEESAW_EDGE_DIFFERENCE}  (= 12 SM-flavor − 3 internal)",
        "",
        "Active neutrino masses (eV):",
        f"  m_1 = {am[0]:.4e}  (K_8 Wilson, N_v=8, N_e=28)",
        f"  m_2 = {am[1]:.4e}  (extrapolated via Δm²_21 = {DELTA_M_SQ_21_eV2:.3e} eV²)",
        f"  m_3 = {am[2]:.4e}  (extrapolated via Δm²_31 = {DELTA_M_SQ_31_eV2:.3e} eV²)",
        f"  Σ m_ν = {sum(am):.4e}",
        "",
        "Sterile neutrino masses (eV):",
        f"  m_N1 = {sm[0]:.4e}",
        f"  m_N2 = {sm[1]:.4e}",
        f"  m_N3 = {sm[2]:.4e}",
        f"  Seesaw ratio          α^(9/2) = {seesaw_ratio():.4e}",
        f"  |U_α4|² (mixing²)             = {sterile_active_mixing_sq():.4e}",
        "",
        "PMNS angles (leading order from Spin(8) triality):",
        f"  θ_12 = {math.degrees(pmns.theta_12):.4f}°"
        "  = arctan(1/√2)",
        f"  θ_23 = {math.degrees(pmns.theta_23):.4f}°"
        "  = π/4",
        f"  θ_13 = {math.degrees(pmns.theta_13):.4f}°"
        "  = arcsin(√(3α))   [3 = RANK_SO7]",
        "",
        "Dirac CP phase from π_1(PSU(3)) Z_3 winding:",
        f"  δ_CP   = {math.degrees(delta):.4f}°  = -2π/3 "
        "(Baez Fano octonion-product convention)",
        f"  |J_CP| = {abs(j):.4f}",
    ]
    return "\n".join(lines)


__all__ = [
    # Functions
    "active_mass_lightest",
    "active_masses",
    "delta_cp_winding",
    "jarlskog_invariant",
    "pmns_angles_leading_order",
    "seesaw_ratio",
    "sterile_active_mixing_sq",
    "sterile_masses",
    "substrate_breakdown",
    "sum_active_masses",
    "wilson_mass_eV",
    # Types
    "PMNSAngles",
    # Anchors
    "DELTA_M_SQ_21_eV2",
    "DELTA_M_SQ_31_eV2",
]
