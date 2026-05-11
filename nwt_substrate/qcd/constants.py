"""
QCD physical constants — PDG-style.  All energies in GeV.

Color algebra constants (N_c, C_F, C_A, T_R, N_gluons) are sourced
from `nwt_substrate.isa.constants` — the same source-of-truth that
gives gravity its α^(21/2) Wilson exponent and QED its 8×8 Dirac γ
matrices.  See `nwt_substrate.isa.constants` for the substrate-
identity derivations:

    N_c       = N_C_SU3       = RANK_SO7 = 3
    N_gluons  = N_GLUONS      = DIM_OCTONION = DIM_S_SPIN7 = 8
    C_F       = C_F_SU3       = DIM_OCTONION / (2 × N_C_SU3) = 4/3
    C_A       = C_A_SU3       = N_C_SU3 = 3
    T_R       = T_R_SU3       = 1/2

The 8 gluons of QCD and the 8-dim Dirac spinor space of QED are
both the same DIM_OCTONION — this is one substrate identity, not
two coincidences.
"""

import numpy as np

from ..isa.constants import (
    N_C_SU3, N_GLUONS, C_F_SU3, C_A_SU3, T_R_SU3,
)


# Coupling at M_Z
alpha_s = 0.1179                 # PDG 2022
g_s = float(np.sqrt(4.0 * np.pi * alpha_s))   # ~ 1.218

# Color algebra constants (sourced from substrate ISA)
N_c = N_C_SU3                    # = 3 = RANK_SO7
C_F = C_F_SU3                    # = 4/3 = DIM_OCTONION / (2 × N_C_SU3)
C_A = C_A_SU3                    # = 3 = N_C_SU3
T_R = T_R_SU3                    # = 1/2

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
