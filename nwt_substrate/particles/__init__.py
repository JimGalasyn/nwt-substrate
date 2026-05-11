"""Particle representation: Particle class, compendium, mass, charge, factory.

Substrate identity: every NWT particle has a carrier knot of crossing
number n_q ∈ [0, MAX_CROSSING_NUMBER] = [0, 6], giving N_CARRIER_TYPES
= 7 = N_VERTICES_K7 distinct carrier topologies (one per K_7 vertex).
The Paper 6 mass formula `m/m_e = ((p²+q²)/5) × β-ratio × log-ratio ×
n_q^q` operates on the substrate-K_7 representation theory.

See `particles.substrate_breakdown()` for the structural identities.
"""

from .particle import Particle
from .factory import particle, list_particles
from .mass import paper6_mass_ratio, paper6_mass_mev, ME_MEV
from .charge import gell_mann_nishijima

from ..isa import (
    N_CARRIER_TYPES,
    MAX_CROSSING_NUMBER,
    CARRIER_NAMES,
    N_VERTICES_K7,
)


def substrate_breakdown() -> str:
    """Pretty-print the particle-spectrum substrate identities.

    Maps the Paper 6 (p, q, m, n_q) substrate quantum numbers onto the
    K_7 / Spin(7) structural constants from `nwt_substrate.isa`.
    """
    lines = ["Particle spectrum from K_7 substrate (Paper 6):"]
    lines.append("")
    lines.append(f"    Carrier crossing number n_q ∈ [0, MAX_CROSSING_NUMBER]")
    lines.append(f"      = [0, {MAX_CROSSING_NUMBER}], giving")
    lines.append(f"    N_CARRIER_TYPES = {N_CARRIER_TYPES} distinct carrier topologies")
    lines.append(f"      = N_VERTICES_K7 = {N_VERTICES_K7}")
    lines.append(f"    (one per K_7 vertex / Spin(7) vector direction)")
    lines.append("")
    lines.append(f"    Carrier table (n_q → canonical knot):")
    for n_q, name in enumerate(CARRIER_NAMES):
        physics = {
            0: "leptons (unknotted)",
            1: "leptons (extended unknot)",
            2: "mesons (Hopf link)",
            3: "baryons (trefoil)",
            4: "tetraquarks (figure-8)",
            5: "nucleons + pentaquarks (cinquefoil)",
            6: "hexaquarks / dibaryons",
        }[n_q]
        lines.append(f"      n_q = {n_q}:  {name:<14}  → {physics}")
    lines.append("")
    lines.append(f"    Paper 6 mass formula:")
    lines.append(f"      m/m_e = ((p²+q²)/5) × (β/β_e) × (ln(8β)/ln(8β_e)) × n_q^q")
    lines.append(f"    The '5' in the denominator is the electron's (p²+q²)")
    lines.append(f"      = 2² + 1² = 5 (anchor: the simplest stable knot).")
    lines.append(f"    The n_q^q link enhancement counts contributions over the")
    lines.append(f"      carrier's q-fold link structure.")
    lines.append("")
    lines.append(f"    Cross-shim ID: same K_7 substrate that drives:")
    lines.append(f"      chemistry's K_7-hub stabilization (coronene)")
    lines.append(f"      gravity's α^(21/2) Wilson amplitude")
    lines.append(f"      qed's 8×8 Dirac γ matrices")
    lines.append(f"      qcd's 8 gluons and 3 colors")
    lines.append(f"      now selects which 7 carrier topologies are allowed.")
    return "\n".join(lines)


class _SubstrateNamespace:
    """Substrate identities visible from the particles shim."""
    N_CARRIER_TYPES = N_CARRIER_TYPES
    MAX_CROSSING_NUMBER = MAX_CROSSING_NUMBER
    CARRIER_NAMES = CARRIER_NAMES
    N_VERTICES_K7 = N_VERTICES_K7
    breakdown = staticmethod(substrate_breakdown)


substrate = _SubstrateNamespace()


__all__ = [
    "Particle",
    "particle",
    "list_particles",
    "paper6_mass_ratio",
    "paper6_mass_mev",
    "ME_MEV",
    "gell_mann_nishijima",
    "substrate",
    "substrate_breakdown",
]
