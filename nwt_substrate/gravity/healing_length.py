"""
Healing length ξ — the substrate's domain-of-applicability marker.

In every superfluid condensate, the healing length is the characteristic
length over which the order parameter recovers from a disturbance. It
sets the smallest scale at which topological-defect descriptions (vortex
cores, sound modes, Wilson amplitudes) remain valid. Beyond ξ, the
substrate algebra fails to bite — observables decouple from substrate
predictions.

For the NWT substrate (universal scale, electron-mass-set):

    ξ_substrate = ℏ / (m_e c) = reduced electron Compton wavelength
                ≈ 3.86 × 10⁻¹³ m

This is the natural length unit of all framework predictions at
particle / atomic / molecular scales.

At cosmic scales, the framework derives a "second-tier" healing length
via Paper 17's Spin(7)_k Chern-Simons level k_CS ≈ 31.73 (which itself
comes from κ_Macken via sin(π/(k + h_v)) = 1/(2^(1/4) κ)):

    ξ_cosmo = ξ_substrate × (1/√α)^k_CS ≈ 100 kpc

This explains why galaxy-scale structures show only inherited-axis
substrate effects rather than internal K_7 saturation, and why
supercluster-scale tests fail: superclusters are at ~100 × ξ_cosmo,
firmly outside the K_7 domain of applicability.

Reference: vortex-vision/analysis/xi_cosmo_derivation.py (2026-05-15).
"""
from __future__ import annotations

import math

from ..isa import (
    ALPHA_NWT,
    KAPPA_MACKEN,
    H_V_SO7,
)
from .constants import (
    HBAR_J_S,
    C_LIGHT_M_S,
    M_ELECTRON_KG,
)


# ---------------------------------------------------------------------------
# Section 1 — Universal substrate healing length
# ---------------------------------------------------------------------------

XI_SUBSTRATE_M: float = HBAR_J_S / (M_ELECTRON_KG * C_LIGHT_M_S)
"""ξ_substrate = ℏ / (m_e c) ≈ 3.8616 × 10⁻¹³ m.

The reduced electron Compton wavelength. This is the NWT substrate's
healing length: the universal length below which all substrate-algebraic
predictions (vortex cores, K_7 Wilson amplitudes, mode spectra) are
load-bearing, and above which they progressively decouple.

Setting m_e (and therefore ξ_substrate) is what Paper 17 / Paper 19's
m_e / M_Pl = (8/7) × α^(21/2) × Λ_NLO derivation accomplishes — the
universe's substrate length scale is fixed by the electron-Compton bridge."""


def xi_substrate_m() -> float:
    """ξ_substrate in metres (alias for XI_SUBSTRATE_M as a function)."""
    return XI_SUBSTRATE_M


# ---------------------------------------------------------------------------
# Section 2 — Paper 17 Spin(7)_k Chern-Simons level k_CS
#
# Derived in Paper 17 from sin(π/(k + h_v)) = 1/(2^(1/4) κ), where
# κ = κ_Macken is the substrate aspect ratio and h_v = 5 is the dual
# Coxeter number of so(7). Solving:
#
#     k_CS = π / arcsin(1 / (2^(1/4) κ)) - h_v ≈ 31.73
# ---------------------------------------------------------------------------

def k_chern_simons(kappa: float = KAPPA_MACKEN,
                   h_v: int = H_V_SO7) -> float:
    """Paper 17 Spin(7)_k Chern-Simons level k_CS.

    Closed form: k_CS = π / arcsin(1 / (2^(1/4) κ)) - h_v
    Numerical value (κ = κ_Macken, h_v = 5): k_CS ≈ 31.73.

    k_CS quantizes the maximum coherent transition the substrate
    condensate can sustain — above k_CS coupling steps, the
    topological description fails."""
    arg = 1.0 / (2.0 ** 0.25 * kappa)
    return math.pi / math.asin(arg) - h_v


# ---------------------------------------------------------------------------
# Section 3 — Cosmic-scale healing length ξ_cosmo
# ---------------------------------------------------------------------------

def xi_cosmo_m(alpha: float = ALPHA_NWT,
               xi_substrate: float = XI_SUBSTRATE_M,
               kappa: float = KAPPA_MACKEN,
               h_v: int = H_V_SO7) -> float:
    """ξ_cosmo from Paper 17 Chern-Simons step-count amplification.

    Closed form: ξ_cosmo = ξ_substrate × (1/√α) ^ k_CS

    With k_CS ≈ 31.73, α = α_NWT, ξ_substrate = ℏ/(m_e c):

        ξ_cosmo ≈ 3.08 × 10²¹ m ≈ 100 kpc (galactic-halo virial scale)

    Note: the precise value depends sensitively on k_CS (exponential
    amplification). The ~100-kpc band brackets galactic halos and
    falls below cluster/supercluster scales, in agreement with the
    framework's mesoscale-band boundary (substrate K_7 active inside
    galaxies, dormant at cluster scales and above).
    """
    k_cs = k_chern_simons(kappa, h_v)
    return xi_substrate * (1.0 / math.sqrt(alpha)) ** k_cs


# ---------------------------------------------------------------------------
# Section 4 — Domain-of-applicability tables
# ---------------------------------------------------------------------------

# Length scales (in metres) of standard physical regimes, used to map
# substrate prediction domains onto observational scales.
_PARSEC_M = 3.0857e16
_KPC_M = 1e3 * _PARSEC_M
_MPC_M = 1e6 * _PARSEC_M


def scale_regime_table(alpha: float = ALPHA_NWT) -> dict:
    """Return a dict of substrate-relevant length scales for reference.

    Anything < ξ_substrate × (1/√α)^k_CS is K_7-active. The "mesoscale
    band" (galaxy halos, dwarf galaxies) sits near the transition,
    and clusters/superclusters are firmly in the K_7-dormant regime.
    """
    xi_c = xi_cosmo_m(alpha)
    return {
        "xi_substrate_m":       XI_SUBSTRATE_M,
        "electron_compton_m":   XI_SUBSTRATE_M,
        "atomic_bohr_radius_m": 5.29e-11,
        "nuclear_fm":           1e-15,
        "molecular_nm":         1e-9,
        "stellar_AU":           1.496e11,
        "galactic_disk_kpc":    30.0 * _KPC_M,
        "xi_cosmo_m":           xi_c,
        "xi_cosmo_kpc":         xi_c / _KPC_M,
        "cluster_Mpc":          10.0 * _MPC_M,
        "supercluster_Mpc":     100.0 * _MPC_M,
        "label": (
            "K_7-active below ξ_cosmo; inherited-axis only between "
            "ξ_cosmo and ~10× ξ_cosmo; K_7-dormant above"
        ),
    }
