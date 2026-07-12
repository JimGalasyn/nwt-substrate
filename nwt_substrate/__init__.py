"""
nwt_substrate — substrate-algebraic computation for Null Worldtube Theory.

.. deprecated:: July 2026
    The NWT program is retired: its physical claims did not survive
    precision audit (31/36 dimensionless claims excluded at experimental σ;
    the fitting pattern matches ~83% of RANDOM targets inside Newton's-G's
    error bar; provenance of the headline constants adjudicated
    post-selected/fitted/motivated by a memory-blind Auditor).  The library
    is preserved unmodified as the audited case-study specimen.  Read the
    README deprecation notice, ``benchmarks/surface.py``,
    ``benchmarks/neighbouring_value.py`` and
    ``docs/BENCHMARK_TRIAGE_2026-07-12.md`` before citing any value it
    computes.  A ``DeprecationWarning`` is emitted on import.

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
    nwt_substrate.cosmology    — Λ, η_B, Ω_b/Ω_c, anisotropy axes
    nwt_substrate.neutrino     — K_8 mass tower, PMNS, sterile spectrum
    nwt_substrate.electroweak  — couplings, CKM, decays, decay constants
    nwt_substrate.atomic       — hydrogen spectroscopy + a_e
    nwt_substrate.qed/.qcd     — gauge-theory shims with substrate α
    nwt_substrate.chemistry    — aromaticity, NICS, vibrational modes
    nwt_substrate.dark_sector  — L_NWT Higgs portal for K_8 DM tower
    nwt_substrate.benchmarks   — substrate vs traditional speed/accuracy
    nwt_substrate.qiskit       — quantum-hardware interface (TODO)
"""

import warnings as _warnings

_warnings.warn(
    "nwt-substrate is DEPRECATED (July 2026): the NWT program is retired — "
    "its physical claims did not survive precision audit. No value this "
    "library computes should be cited as a confirmed prediction. See the "
    "README deprecation notice and docs/BENCHMARK_TRIAGE_2026-07-12.md.",
    DeprecationWarning,
    stacklevel=2,
)

from .particles.particle import Particle
from .particles.factory import particle, list_particles
from .compositions import (
    CompositeParticle,
    compose,
    decompose,
    LAMBDA_QCD_MEV,
)
from . import diagrams
from . import isa

from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    __version__ = _pkg_version("nwt-substrate")
except PackageNotFoundError:
    # Running from a source tree that hasn't been installed (e.g. tests
    # invoked directly from a git checkout without `pip install -e .`).
    __version__ = "0.0.0+unknown"

__all__ = [
    "Particle",
    "particle",
    "list_particles",
    "CompositeParticle",
    "compose",
    "decompose",
    "LAMBDA_QCD_MEV",
    "diagrams",
    "isa",
]
