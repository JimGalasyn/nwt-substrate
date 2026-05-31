"""
QED physical constants.  The masses and reference scales are genuine measured
inputs (the things you'd look up in a PDG booklet); the coupling ``alpha`` is the
substrate's *derived* value, not the booklet one — it is the one quantity here
NWT predicts rather than measures.  All energies in GeV.
"""

import numpy as np

from ..isa import ALPHA_SUBSTRATE


# Coupling — the substrate's *derived* fine-structure constant (Thomson limit),
# α = 1/(25π√3 + 1) ≈ 1/137.0350, NOT the measured CODATA booklet value.  The QED
# cross sections (σ ∝ α²) must close on the derived coupling so that perturbing it
# moves them — otherwise α never reaches the prediction and these observables read
# as single-sector "hidden SPOFs" (cf. the route-diversity readout).  This is the
# same derived-not-measured discipline o10/predict enforce; the booklet value here
# was a measured-input leak into a quantity the substrate derives.
alpha = float(ALPHA_SUBSTRATE)       # = 1/137.0350 (substrate), 7.6 ppm from CODATA
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
