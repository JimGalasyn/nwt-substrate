"""
Electroweak constants (PDG 2022 values).

The substrate's internal SU(2) (Cl(0,7) bivector triplet, Walk-phase 4a)
is vector under Lorentz, while the SM's electroweak SU(2)_L is chiral.
The connection between the two is one of the substrate's open threads;
the EW shim works at the SM Lagrangian level (chiral V-A couplings) and
labels them with substrate references where the connection is established.

Standard EW parameters:
"""

from __future__ import annotations

# Gauge boson masses (GeV) and widths (GeV)
M_Z = 91.1876
GAMMA_Z = 2.4952
M_W = 80.379
GAMMA_W = 2.085

# Weinberg angle, effective at the Z pole
SIN2_THETA_W = 0.23121
COS2_THETA_W = 1.0 - SIN2_THETA_W

# Couplings.  We define g_Z from the Fermi constant (G_F) because the
# Z-coupling absorbs the running of α from low energy up to M_Z (~6%
# enhancement).  Using low-energy α here would underestimate Γ_Z and
# the Z-pole cross-section by that amount.
#
#   G_F = g_W² / (4 √2 M_W²) = g_Z² / (4 √2 M_Z²)   →   g_Z² = 4 √2 G_F M_Z²
import math as _math
ALPHA_QED = 1.0 / 137.035999084          # low-energy α (used for QED amplitudes)
G_F_GEV = 1.1663787e-5                    # GeV⁻², from muon lifetime
G_Z_SQ = 4.0 * _math.sqrt(2.0) * G_F_GEV * M_Z ** 2
G_Z = G_Z_SQ ** 0.5
G_W_SQ = G_Z_SQ * COS2_THETA_W            # g_W² = g_Z² cos² θ_W
G_W = G_W_SQ ** 0.5
# Effective EW charge at the Z pole (used in electroweak couplings)
E_CHARGE_Z = (G_Z_SQ * SIN2_THETA_W * COS2_THETA_W) ** 0.5
# Low-energy electric charge (for QED amplitude / γ-channel)
E_CHARGE = (4.0 * 3.141592653589793 * ALPHA_QED) ** 0.5

# Higgs VEV
V_HIGGS_GEV = 246.22
M_HIGGS = 125.10
