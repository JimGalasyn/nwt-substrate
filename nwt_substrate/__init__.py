"""
nwt_substrate — substrate-algebraic computation for Null Worldtube Theory.

Provides Particle, walk-phase scattering/decay, mass formula, Hopf-linked
exotic predictions, K_7 graph state, and (eventually) Qiskit interface.

Quick start:

    >>> import nwt_substrate as nwt
    >>> p = nwt.particle("p")          # build proton from compendium
    >>> p.mass_pred                    # Paper 6 mass prediction (MeV)
    937.24...
    >>> p.J, p.P, p.Q, p.B
    (0.5, 1, 1, 1)

Structure:

    nwt_substrate.algebra      — octonions, Clifford, Dirac, 2I irreps
    nwt_substrate.topology     — torus knots, K_7 toroidal embedding
    nwt_substrate.particles    — Particle class, compendium, mass, charge
    nwt_substrate.walk_phase   — substrate-algebraic Compton, decay
    nwt_substrate.compositions — molecular and Hopf-linked composites
    nwt_substrate.gravity      — K_7 Wilson amplitude, G derivation
    nwt_substrate.qiskit       — quantum-hardware interface (TODO)
"""

from .particles.particle import Particle
from .particles.factory import particle, list_particles
from .compositions import (
    CompositeParticle,
    compose,
    decompose,
    LAMBDA_QCD_MEV,
)
from . import diagrams

__version__ = "0.1.0"

__all__ = [
    "Particle",
    "particle",
    "list_particles",
    "CompositeParticle",
    "compose",
    "decompose",
    "LAMBDA_QCD_MEV",
    "diagrams",
]
