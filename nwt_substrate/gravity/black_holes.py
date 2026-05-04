"""
Schwarzschild geometry and Hawking thermodynamics for SM particles.

Uses the substrate-derived `G` (NNLO formula from Paper 17) consistently
to compute classical-GR observables for any particle in the substrate
compendium.
"""

from __future__ import annotations
import math

from .constants import (
    G_NEWTON_SI, HBAR_J_S, C_LIGHT_M_S,
    M_ELECTRON_KG, M_PROTON_KG,
)
from .coupling import G_substrate_SI


# Boltzmann constant (J/K)
K_B_J_K = 1.380649e-23


def schwarzschild_radius_m(mass_kg: float, G: float = G_NEWTON_SI) -> float:
    """
    Schwarzschild radius r_s = 2 G M / c²  [m].
    """
    return 2.0 * G * mass_kg / C_LIGHT_M_S ** 2


def hawking_temperature_K(mass_kg: float, G: float = G_NEWTON_SI) -> float:
    """
    Hawking temperature T_H = ℏ c³ / (8π G M k_B)  [K].
    """
    return (HBAR_J_S * C_LIGHT_M_S ** 3
            / (8.0 * math.pi * G * mass_kg * K_B_J_K))


def hawking_luminosity_W(mass_kg: float, G: float = G_NEWTON_SI) -> float:
    """
    Hawking luminosity L_H = ℏ c⁶ / (15360 π G² M²)  [W].
    """
    return (HBAR_J_S * C_LIGHT_M_S ** 6
            / (15360.0 * math.pi * G ** 2 * mass_kg ** 2))


def evaporation_time_s(mass_kg: float, G: float = G_NEWTON_SI) -> float:
    """
    Black-hole evaporation time t_evap = 5120 π G² M³ / (ℏ c⁴)  [s].
    (Standard textbook formula for a Schwarzschild black hole.)
    """
    return (5120.0 * math.pi * G ** 2 * mass_kg ** 3
            / (HBAR_J_S * C_LIGHT_M_S ** 4))


# Catalog of SM particles (kg) for convenient lookup
PARTICLE_MASSES_KG = {
    "electron": M_ELECTRON_KG,
    "muon": 1.883531627e-28,
    "tau": 3.16754e-27,
    "proton": M_PROTON_KG,
    "neutron": 1.67492749804e-27,
    "W": 1.43269e-25,
    "Z": 1.62558e-25,
    "Higgs": 2.231e-25,
    "top": 3.078e-25,
    "Earth": 5.972e24,
    "Sun": 1.989e30,
    "M_Planck": 2.176e-8,        # 1 Planck mass in kg
}


def black_hole_summary(particle: str, use_substrate_G: bool = False) -> dict:
    """
    Compute Schwarzschild radius, Hawking temperature, luminosity, and
    evaporation time for a given particle treated as a hypothetical
    Schwarzschild black hole.

    Set use_substrate_G=True to use the Paper 17 NNLO substrate prediction
    (-11 ppm from CODATA, inside the +-22 ppm experimental band).
    """
    if particle not in PARTICLE_MASSES_KG:
        raise KeyError(f"Particle '{particle}' not in PARTICLE_MASSES_KG. "
                       f"Available: {list(PARTICLE_MASSES_KG.keys())}")
    M = PARTICLE_MASSES_KG[particle]
    G = G_substrate_SI() if use_substrate_G else G_NEWTON_SI
    return {
        "particle": particle,
        "mass_kg": M,
        "G_used": G,
        "r_s_m": schwarzschild_radius_m(M, G),
        "T_H_K": hawking_temperature_K(M, G),
        "L_H_W": hawking_luminosity_W(M, G),
        "t_evap_s": evaporation_time_s(M, G),
    }
