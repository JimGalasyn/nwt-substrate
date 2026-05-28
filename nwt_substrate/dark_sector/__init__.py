"""
Dark-sector phenomenology from the NWT substrate (K_8 mass tower).

The substrate predicts a multi-species dark matter tower at mass scales
m_n = c_n · α^(N_e/2) · M_Pl  with discrete N_e set by K_8 closed-walk
edge enumeration.  The 98 GeV WIMP is the N_e=16 rung.

Modules:
    wimp_98gev — L_NWT Higgs-portal calculation for the 98 GeV DM
                 (direct-detection sigma_SI + LHC production).
"""

from .wimp_98gev import (
    WIMP_98GeV,
    higgs_portal_lagrangian_text,
    sigma_si_higgs_portal,
    lhc_production_cross_section,
    compare_with_lz_limit,
    predict_all,
)

__all__ = [
    "WIMP_98GeV",
    "higgs_portal_lagrangian_text",
    "sigma_si_higgs_portal",
    "lhc_production_cross_section",
    "compare_with_lz_limit",
    "predict_all",
]
