"""
Amplitude infrastructure for substrate-algebraic Feynman diagrams.

Provides:
  - spinors:        positive/negative energy spinors u(p,s), v(p,s), adjoints
  - polarizations:  transverse polarization vectors for massless / massive bosons
  - propagators:    fermion S_F, scalar D_F, photon D_mu_nu Feynman propagators
  - vertices:       QED gamma^mu, V-A weak gamma^mu(1 - gamma^5), Yukawa
  - kinematics:     2->2 phase space helpers
  - processes:      pre-built standard amplitudes (Compton, etc.)

All amplitudes use the substrate Cl(1,3) Dirac matrices from
nwt_substrate.algebra.dirac.  Spinors are 8-dimensional (Cl(0,7) embedding):
each (p, m) eigenspace is 4-fold degenerate (2 physical spin x 2 internal SU(2)).

Quick start:

    >>> from nwt_substrate.algebra.octonions import make_octonion_table
    >>> from nwt_substrate.algebra.dirac import make_lorentz_gammas
    >>> from nwt_substrate.amplitudes.processes import compton
    >>>
    >>> T = make_octonion_table()
    >>> gammas = make_lorentz_gammas(T)
    >>> M_sq, omega_prime = compton.M_squared_avg(omega=0.1, theta=1.5,
    ...                                            m=0.511, gammas=gammas)
"""

from .spinors import (
    positive_energy_spinors,
    negative_energy_spinors,
    null_spinors,
    adjoint_spinor,
)
from .polarizations import (
    transverse_polarizations,
    massive_polarizations,
)
from .propagators import (
    fermion_propagator,
    scalar_propagator,
    photon_propagator,
)
from .vertices import (
    qed_vertex,
    va_vertex,
    yukawa_vertex,
)
from .cross_sections import (
    cm_momentum,
    sigma_2to2_cm,
    dsigma_dOmega_cm,
    sigma_compton_lab,
    thomson_cross_section,
    sigma_eemumu_high_energy,
    gev_inv_sq_to_mb,
    gev_inv_sq_to_pb,
    gamma_1to2,
    gamma_muon_decay_textbook,
    G_F,
)

__all__ = [
    "positive_energy_spinors",
    "negative_energy_spinors",
    "null_spinors",
    "adjoint_spinor",
    "transverse_polarizations",
    "massive_polarizations",
    "fermion_propagator",
    "scalar_propagator",
    "photon_propagator",
    "qed_vertex",
    "va_vertex",
    "yukawa_vertex",
]
