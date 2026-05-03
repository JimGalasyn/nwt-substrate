"""
Gravitational and cosmological constants (CODATA 2018 / Planck 2018).

The substrate's *prediction* for G is computed in coupling.py from
substrate primitives; this file provides the OBSERVED CODATA values for
comparison.
"""

from __future__ import annotations
import math


# Fundamental constants (SI)
HBAR_J_S = 1.054571817e-34        # J·s
C_LIGHT_M_S = 2.99792458e8         # m/s, exact
G_NEWTON_SI = 6.67430e-11          # m³ kg⁻¹ s⁻²  (CODATA 2018, σ ≈ 22 ppm)

# Particle masses (kg)
M_ELECTRON_KG = 9.1093837015e-31
M_PROTON_KG = 1.67262192369e-27
M_NEUTRON_KG = 1.67492749804e-27

# Planck scales (natural units, GeV)
M_PLANCK_GEV = 1.220890e19         # M_Pl = √(ℏ c / G) in GeV
M_PLANCK_REDUCED_GEV = M_PLANCK_GEV / math.sqrt(8.0 * math.pi)
                                   # M̃_Pl = √(ℏ c / 8π G) ≈ 2.435 × 10¹⁸ GeV

# Cosmological observables (Planck 2018)
H_0_KM_S_MPC = 67.4                # Hubble parameter, km/s/Mpc
H_0_GEV = 1.439e-42                # converted to GeV (= H_0 ℏ / c)
LAMBDA_CC_OBSERVED_GEV4 = 5.6e-47  # Λ_cc ≈ 3 H_0² Ω_Λ (in GeV⁴)
RHO_CRIT_GEV4 = 8.1e-47            # critical density, GeV⁴

# Atomic / lab constants
ALPHA_QED = 1.0 / 137.035999084
M_E_GEV = 0.510998928e-3
M_P_GEV = 0.93827208816
