"""
Cosmogenesis: throat selection, energy transfer, and Thorne-equilibrium spin.

The framework hypothesis: a sufficiently spun-up Kerr black hole hosts a
toroidal ergosurface (Beltracchi-class regular rotating BH) that acts as
a substrate-vortex with aspect ratio κ_parent. Under throat formation,
the parent's K_7 modes adiabatically rescale to the daughter universe's
substrate torus at fixed κ_daughter = κ_Macken ≈ 9.844. Mode-overlap
between (κ_parent, κ_daughter) bases determines the energy-transfer
fraction f_J.

Two structurally equivalent forms for f_J:

    EDGE model (rank-3 bridging edges, Paper 17 √α-per-edge rule):
        f_J = (1 - √α) ^ RANK_SO7
            = (1 - √α) ^ 3
            ≈ 0.7437 - small higher-order corrections

    LEADING-ORDER expansion:
        f_J ≈ 1 - 3√α

The "3" is rank(so(7)), the algebraic dimension of the cosmogenic channel
in the Spin(7) Coxeter decomposition 8 = h_v + rank = 5 + 3. The "√α" is
the per-edge Wilson amplitude.

Thorne-equilibrium prediction: astrophysical thin-disk accretion drives
parent BHs to exactly the spin a*_Thorne ≈ 0.998 where the toroidal
ergosurface aspect ratio satisfies κ_parent / κ_daughter = 3(1 + √α).
This identifies Thorne's thin-disk equilibrium spin as the cosmogenic-
viability threshold.

Reference: vortex-vision/analysis/f_J_derivation_from_mode_overlap.py
(2026-05-15).
"""
from __future__ import annotations

import math

from ..isa import (
    ALPHA_NWT,
    KAPPA_MACKEN,
    RANK_SO7,
)


# Daughter substrate's aspect ratio is exactly κ_Macken (closed form).
KAPPA_DAUGHTER: float = KAPPA_MACKEN
"""Daughter-universe substrate aspect ratio = κ_Macken ≈ 9.844.

Closed form: κ_daughter = √((25π√3 + 1) / √2) = √(1/α / √2).
This is the substrate's intrinsic length-scale ratio and the same ratio
that appears in all K_7 mode-spectrum calculations (Paper 6 mass spectrum,
Paper 17 trefoil Aharonov-Bohm)."""


# ---------------------------------------------------------------------------
# Section 1 — Energy-transfer fraction f_J
# ---------------------------------------------------------------------------

def f_J_cosmogenic(alpha: float = ALPHA_NWT,
                   n_bridging_edges: int = RANK_SO7) -> float:
    """Cosmogenic energy-transfer fraction parent → daughter universe.

    Edge model:
        f_J = (1 - √α) ^ n_bridging_edges

    With n_bridging_edges = RANK_SO7 = 3:
        f_J = (1 - √α)³ ≈ 0.7649  (full closed form)
        1 - 3√α        ≈ 0.7437  (leading-order expansion)

    The 3 bridging edges are the rank-3 algebraic dimension of the
    aspect-ratio bridge between κ_parent and κ_daughter (one edge per
    Z_3 orbit of the K_7 Heffter automorphism). Each contributes a
    per-edge amplitude (1 - √α) from Paper 17's √α-per-edge Wilson rule.

    Parameters
    ----------
    alpha : float, default ALPHA_NWT
        Fine-structure constant. Use ALPHA_NWT (closed form) for the
        framework's structural prediction; pass a different α for
        consistency checks against CODATA.
    n_bridging_edges : int, default RANK_SO7
        Number of K_7 edges in the cosmogenic channel. Default 3 = rank(so(7)).
    """
    return (1.0 - math.sqrt(alpha)) ** n_bridging_edges


def f_J_cosmogenic_leading(alpha: float = ALPHA_NWT) -> float:
    """Leading-order expansion of f_J:  f_J ≈ 1 - 3√α.

    Drops the α and α^(3/2) terms in (1 - √α)³.
    Use f_J_cosmogenic() for the full closed-form result.
    """
    return 1.0 - RANK_SO7 * math.sqrt(alpha)


# ---------------------------------------------------------------------------
# Section 2 — Parent / daughter κ ratio for throat-selection viability
# ---------------------------------------------------------------------------

def kappa_parent_over_daughter(alpha: float = ALPHA_NWT) -> float:
    """Cosmogenic-viability ratio of parent ergosurface aspect ratio to
    daughter substrate aspect ratio.

    Substrate form:
        κ_parent / κ_daughter = RANK_SO7 × (1 + √α) = 3(1 + √α)

    Interpretation: the parent's toroidal ergosurface must host the 3
    Z_3 K_7 orbits at aspect ratio 3× the daughter's, with an O(√α)
    correction per orbit from the per-edge amplitude factor. This is
    the throat-selection viability criterion — only BHs with
    κ_parent matching this ratio can act as cosmogenic parents.
    """
    return RANK_SO7 * (1.0 + math.sqrt(alpha))


def kappa_parent_required(alpha: float = ALPHA_NWT,
                          kappa_daughter: float = KAPPA_DAUGHTER) -> float:
    """Required κ_parent for cosmogenic viability.

    κ_parent = 3(1 + √α) × κ_daughter ≈ 32.06
    """
    return kappa_parent_over_daughter(alpha) * kappa_daughter


# ---------------------------------------------------------------------------
# Section 3 — Thorne thin-disk equilibrium spin from cosmogenic viability
#
# Geometric derivation of κ_parent (substrate-vortex circulation matching):
#     R_major = 2M  (equatorial bulge of outer ergosurface)
#     r_minor = (r_+ - r_-) / 2 = √(M² - a²) = M·√(1 - a*²)
#               (horizon half-gap; radial extent of inter-horizon region)
#     κ_parent = R / r = 2 / √(1 - a*²)        — EXACT Kerr identity
#
# Physical principle: a substrate vortex sits in the inter-horizon region,
# centered at the bifurcation 2-sphere. Its minor radius equals the
# horizon half-gap, and its circulation Γ = π(r_+ - r_-) = 4π·M·r_+·κ_surf
# is set by the causal bound v_inner = Γ/(2π r_min) = c at the vortex's
# inner edge. The circulation is proportional to surface gravity.
#
# Cosmogenic-viability condition:
#     κ_parent / κ_daughter = 3(1 + √α)       (K_7 algebra)
# gives κ_parent_required ≈ 32.054.
#
# Inverting κ_parent = 2/√(1-a*²) under this condition:
#     a*_eq = √(1 - 4 / κ_parent_required²) ≈ 0.998052
#
# This is the unique spin where the parent's horizon half-gap geometrically
# realizes the K_7 viability ratio against the daughter's universal κ_Macken.
# Thorne 1974's thin-disk equilibrium hits this spin via disk-as-substrate-
# vortex (see framework_disk_as_substrate_vortex.md), making the K_7-Thorne
# convergence a single geometric identity rather than a coincidence.
#
# Reference: analysis/beltracchi_ergosurface_kappa.py Section 6 contains
# the numerical verification and full derivation.
# ---------------------------------------------------------------------------

def thorne_a_star_equilibrium(alpha: float = ALPHA_NWT,
                              kappa_daughter: float = KAPPA_DAUGHTER) -> float:
    """Cosmogenic-viability spin a*_eq from the κ ratio condition.

    Inverts κ_parent ≈ 2/√(1 - a*²) under
        κ_parent_required = 3(1 + √α) × κ_daughter

    Result: a*_eq ≈ 0.998052 — within 0.005 of Thorne's 1974 thin-disk
    equilibrium spin a*_Thorne ≈ 0.998. Astrophysical accretion drives
    BHs to exactly the cosmogenic-viability threshold.
    """
    kappa_p = kappa_parent_required(alpha, kappa_daughter)
    return math.sqrt(1.0 - (2.0 / kappa_p) ** 2)


# ---------------------------------------------------------------------------
# Section 4 — Mode-overlap transmission between κ_parent and κ_daughter
#
# Each (p, q) Fourier mode on the parent ergosurface has frequency
#     ω_p(p, q) = √((p/κ_parent)² + q²)
# and on the daughter substrate
#     ω_d(p, q) = √((p/κ_daughter)² + q²)
#
# Adiabatic throat evolution gives Bogoliubov coefficients
#     α(p,q) = ½(√(ω_p/ω_d) + √(ω_d/ω_p))
#     β(p,q) = ½(√(ω_p/ω_d) - √(ω_d/ω_p))
# with |α|² - |β|² = 1. The transmission probability is
#     T(p,q) = 1/|α|² = 4 ω_p ω_d / (ω_p + ω_d)²
#
# T = 1 when frequencies match; T < 1 when they don't, with the lost
# fraction going into Hawking-like mode-pair production.
# ---------------------------------------------------------------------------

def mode_overlap(p: int, q: int,
                 kappa_parent: float, kappa_daughter: float) -> float:
    """Bogoliubov mode-transmission probability T(p, q) for a (p, q) K_7
    mode propagating from parent (κ_parent) to daughter (κ_daughter).

    T = 4 ω_p ω_d / (ω_p + ω_d)²
    """
    omega_p = math.sqrt((p / kappa_parent) ** 2 + q ** 2)
    omega_d = math.sqrt((p / kappa_daughter) ** 2 + q ** 2)
    denom = omega_p + omega_d
    if denom < 1e-30:
        return 1.0
    return 4.0 * omega_p * omega_d / denom ** 2


# ---------------------------------------------------------------------------
# Section 5 — Convenience accessors
# ---------------------------------------------------------------------------

def cosmogenesis_summary(alpha: float = ALPHA_NWT) -> dict:
    """Return a dict of the framework's cosmogenesis predictions for
    quick comparison."""
    return {
        "f_J_full": f_J_cosmogenic(alpha),
        "f_J_leading": f_J_cosmogenic_leading(alpha),
        "kappa_parent_over_daughter": kappa_parent_over_daughter(alpha),
        "kappa_parent_required": kappa_parent_required(alpha),
        "thorne_a_star_eq": thorne_a_star_equilibrium(alpha),
        "thorne_a_star_observed": 0.998,
        "label": (
            "f_J = (1-√α)³, κ_p/κ_d = 3(1+√α), a*_eq from κ_p ≈ 2/√(1-a*²)"
        ),
    }
