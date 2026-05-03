"""
QED physical constants — the things you'd look up in a PDG booklet.
All energies in GeV.
"""

import numpy as np


# Coupling
alpha = 1.0 / 137.035999084          # fine-structure constant (Thomson limit)
e_charge = float(np.sqrt(4.0 * np.pi * alpha))   # = 0.30282..., dimensionless

# Lepton masses (GeV)
m_e   = 0.510998928e-3               # electron
m_mu  = 0.10565837                   # muon
m_tau = 1.77686                      # tau

# Electroweak / hadronic reference scales (GeV)
m_Z   = 91.1876
m_W   = 80.379

# Classical electron radius (GeV^-1, natural units)
r_e = alpha / m_e
