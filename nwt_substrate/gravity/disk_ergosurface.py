"""
Disk-ergosurface matching condition and Spin(7) Coxeter decomposition.

Premise: a Kerr BH's accretion disk is itself a substrate-vortex with
aspect ratio κ_disk. At the inner edge it must match the toroidal
ergosurface (κ_parent). The match isn't 1:1 — the substrate algebra
sets a specific ratio determined by Spin(7)'s octonion components.

The Spin(7) Coxeter decomposition of the 8 octonion components:

    8 = 5 + 3 = H_V_SO7 + RANK_SO7
                ↑           ↑
                steady-     cosmogenic
                state       channel
                (radiation  (throat
                channel)    transfer)

This splits the substrate interface into two channels with distinct
behavior:

  - The 5 dual-Coxeter components form the steady-state radiation
    channel (Bardeen extremal-Kerr 5√α match).
  - The 3 rank components form the cosmogenic transfer channel (f_J
    = (1-√α)³ via 3 K_7 Z_3 bridging edges).

The disk-ergosurface matching condition:

    κ_disk / κ_parent = DIM_OCTONION × (1 - √α) = 8(1 - √α) ≈ 7.32

which matches the observed QPO-implied κ_disk/κ_parent ratio at ~0.5 %.

Compounded with the cosmogenic κ ratio:

    κ_disk / κ_daughter = (κ_disk / κ_parent) × (κ_parent / κ_daughter)
                        = 8(1 - √α) × 3(1 + √α)
                        = 24 (1 - α)

Reference: vortex-vision/analysis/disk_ergosurface_matching.py,
matching_condition_refinement.py (2026-05-15).
"""
from __future__ import annotations

import math

from ..isa import (
    ALPHA_NWT,
    RANK_SO7,
    H_V_SO7,
    DIM_OCTONION,
)


# ---------------------------------------------------------------------------
# Section 1 — Disk / ergosurface aspect-ratio matching
# ---------------------------------------------------------------------------

def kappa_disk_over_parent(alpha: float = ALPHA_NWT) -> float:
    """Inner-disk to outer-ergosurface aspect-ratio match.

    Substrate form:
        κ_disk / κ_parent = DIM_OCTONION × (1 - √α) = 8(1 - √α)

    All 8 octonion components contribute (1 - √α) interfacial
    transmission. Observed ratio from BH-XRB QPO median ≈ 7.28
    (vortex-vision QPO test); substrate prediction 7.32 (deviation 0.5 %).
    """
    return DIM_OCTONION * (1.0 - math.sqrt(alpha))


def kappa_disk_over_daughter(alpha: float = ALPHA_NWT) -> float:
    """Inner-disk to daughter-substrate aspect-ratio match (compounded).

    Substrate form:
        κ_disk / κ_daughter = 8(1 - √α) × 3(1 + √α) = 24(1 - α)

    The (1 - √α)(1 + √α) = 1 - α telescoping is the algebraic
    fingerprint of the matched channel — interferences across the
    cosmogenic / steady-state channels at the substrate level.
    """
    return 8.0 * 3.0 * (1.0 - alpha)


# ---------------------------------------------------------------------------
# Section 2 — Spin(7) Coxeter decomposition of the interface
# ---------------------------------------------------------------------------

def coxeter_decomposition() -> dict:
    """Return the Spin(7) Coxeter decomposition 8 = 5 + 3 with semantics.

    The octonion 8-component substrate interface decomposes as:

        8 = h_v(so(7)) + rank(so(7))
          = 5         + 3
          = (steady-state radiation channel)
          + (cosmogenic transfer channel)

    Under Z_3 ⊂ G_2 acting on octonions:
        1 identity (singlet)
      + 1 Z_3-fixed Im(O) (polar substrate direction)
      + 3 Z_3-orbit A (cosmogenic)
      + 3 Z_3-orbit B (cosmogenic alternate; Z_2-related)
        = 8

    The 5 = 1 + 1 + 3 (one Z_3 orbit + 2 fixed components) gives
    the steady-state channel; the remaining 3 (the other Z_3 orbit)
    gives the cosmogenic channel.
    """
    return {
        "dim_octonion": DIM_OCTONION,
        "h_v_so7_steady_state": H_V_SO7,
        "rank_so7_cosmogenic": RANK_SO7,
        "z3_decomp": {
            "identity": 1,
            "z3_fixed_polar": 1,
            "z3_orbit_A": 3,
            "z3_orbit_B": 3,
            "total": 8,
        },
        "channels": {
            "steady_state": "1 + 1 + 3 (singlet + Z_3-fixed + orbit A)",
            "cosmogenic":   "3 (orbit B = Z_2-flip of orbit A)",
        },
    }


# ---------------------------------------------------------------------------
# Section 3 — Substrate energy budget at the disk-ergosurface interface
# ---------------------------------------------------------------------------

def substrate_energy_budget(alpha: float = ALPHA_NWT) -> dict:
    """Total substrate dissipation at the disk-ergosurface interface
    decomposed by Spin(7) Coxeter channel.

    At extremality, all 8 octonion components dissipate at √α each:
        total = 8√α
        steady_state (5) = 5√α        ↔ Bardeen prograde 1 - 1/√3
        cosmogenic (3)   = 3√α        ↔ Penrose 1 - 1/√2 (leading order)
    """
    sa = math.sqrt(alpha)
    return {
        "total_dissipation":   DIM_OCTONION * sa,
        "steady_state_loss":   H_V_SO7 * sa,
        "cosmogenic_loss":     RANK_SO7 * sa,
        "label": "8√α = 5√α (steady state) + 3√α (cosmogenic)",
    }
