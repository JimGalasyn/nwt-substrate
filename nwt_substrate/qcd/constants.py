"""
QCD physical constants — PDG-style.  All energies in GeV.
"""

import numpy as np


# Coupling at M_Z
alpha_s = 0.1179                 # PDG 2022
g_s = float(np.sqrt(4.0 * np.pi * alpha_s))   # ~ 1.218

# Color algebra constants (SU(3) fundamental)
N_c = 3
C_F = 4.0 / 3.0
C_A = 3.0
T_R = 0.5                        # Tr(T^a T^b) = T_R delta^ab

# Lambda_QCD (1-loop, n_f=5)
Lambda_QCD_5flavor = 0.087       # GeV (1-loop estimate; PDG ~0.21 with 2-loop)

# Current quark masses (GeV, MS-bar at 2 GeV for u,d,s; pole or m_q(m_q) for c,b,t)
m_u = 2.16e-3
m_d = 4.67e-3
m_s = 0.0934
m_c = 1.27
m_b = 4.18
m_t = 172.7

# Constituent quark masses (effective, for vacuum-polarization-like estimates)
m_u_constituent = 0.336
m_d_constituent = 0.336
m_s_constituent = 0.486

# Reference scales (GeV)
m_Z = 91.1876
m_proton = 0.93828
Lambda_chiral = 4.0 * np.pi * 0.0925  # 4 pi f_pi ~ 1.2 GeV
