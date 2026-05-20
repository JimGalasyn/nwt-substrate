"""σ-orbit Wilson-loop dynamics — bound-state masses on K_7 Heffter.

Phase D-2 of the Bogoliubov spectrum derivation.

Phase D-1 established that the substrate has TWO regimes:
  - RELATIVISTIC gap modes (Phase B): m_σ, m_A at fractions of m_e
  - NR SUPERFLUID (Paper 5 §II): ξ = λ̄_C enforced by c_s = c

Bound-state masses (electron, muon, proton, ...) live in the NR
superfluid limit, where vortex-ring self-energy on (p, q) torus knots
gives the mass spectrum (Paper 11 §III):

    m(p,q,m,n_q) / m_e  =  (p² + q²)/5  ×  β/β_e  ×  ln(8β)/ln(8β_e)  ×  n_q^q

with
    β  = √(m²/p² - 1)              phase-closure aspect ratio
    β_e = √5/2                      electron anchor

The (p,q,m,n_q) quantum numbers per particle are imported from
Paper 6 / Paper 11; the mass formula is anchored by the electron
(m_e, β_e, (2,1,3,0)) and reproduces 25 particle masses spanning
4 orders of magnitude to 0.4% median error.

This module connects the K_7 σ-orbit topology to the Paper 11 mass
formula via topological constraints:

  σ-orbit topology       →  allowed (n_q, p, q) classes
    polar=3 (bridging)    →  n_q ∈ {0, 1} (single-component knots)
                              leptons (e, μ, τ)
    cross=3 (parallel)    →  n_q ∈ {3, 5} (multi-component links)
                              hyperons, nucleons (proton, neutron)
    equatorial triangle  →  n_q ∈ {1, 2} (mesons/unknot)
    twisted cross        →  n_q ∈ {3} with Z_3 chirality
                              tau-family, exotic hadrons

The detailed σ-orbit ↔ (p,q,m,n_q) map requires Paper 11/13 reading
of the K_7-on-torus topology; this module sets up the framework and
the verifiable predictions.

References:
  - Paper 5 §III: vortex-ring energy + (2,1) electron derivation
  - Paper 6, 11: mass formula + (p,q,m,n_q) quantum numbers
  - Paper 13: SM capstone, σ-orbit ↔ knot map for hadrons
  - [[framework_k7_heffter_embedding]] σ-orbit structure
  - [[bogoliubov_phase_d1_complete]] NR-superfluid reconciliation
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from nwt_substrate.condensate.sigma_orbits import (
    SIGMA_ORBITS, orbit_invariants,
)
from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.particles.mass import paper6_mass_mev


# ---------------------------------------------------------------------------
# σ-orbit topological constraints on particle quantum numbers
# ---------------------------------------------------------------------------

# n_q allowed by σ-orbit topology.  These are conjectures based on
# the 3-edge structure of each σ-orbit:
#   - star pattern (3 edges sharing one vertex):  single component, n_q ∈ {0,1}
#   - triangle pattern (closed 3-loop):           single loop, n_q ∈ {1}
#   - parallel pattern (3 disconnected edges):    multi-component, n_q ∈ {3, 5}
#   - twisted cross (3 edges with Z_3 twist):     linked components, n_q ∈ {3, 5}
SIGMA_ORBIT_ALLOWED_N_Q: dict[int, list[int]] = {
    0: [0],       # matter bridging:    leptons (electron, muon, tau if same orbit)
    1: [0],       # Z_2 partner:        dark-matter-class
    2: [2, 3, 5], # τ-fixed cross:      mesons, hyperons, nucleons
    3: [1, 2],    # equatorial upper:   mesons, light triangles
    4: [1, 2],    # equatorial lower:   Z_2 partner of orbit 3
    5: [3],       # twisted cross +:    Z_3+ hyperons
    6: [3],       # twisted cross -:    Z_3- hyperons
}


def compatible_particles(orbit_id: int) -> list[dict]:
    """Return compendium entries whose n_q is allowed by this σ-orbit's
    topology.  This is a candidate map; detailed verification requires
    Paper 13 / Paper 11 cross-checks.
    """
    allowed = SIGMA_ORBIT_ALLOWED_N_Q.get(orbit_id, [])
    return [e for e in COMPENDIUM if e.get("n_q") in allowed]


# ---------------------------------------------------------------------------
# Mass formula wrapper — Paper 11 §III
# ---------------------------------------------------------------------------

def kelvin_saffman_aspect(p: int, m_pc: int) -> float:
    """Phase-closure aspect β = √(m²/p² - 1), Paper 11 eq. (closure)."""
    arg = (m_pc ** 2 / p ** 2) - 1.0
    if arg <= 0:
        return float("nan")
    return math.sqrt(arg)


def particle_mass_MeV(p: int, q: int, m_pc: int, n_q: int) -> float:
    """Paper 11 mass formula, in MeV.  Anchor is electron m_e = 0.511 MeV.

    Returns NaN if (p, q, m, n_q) is invalid (e.g., m ≤ p).
    """
    res = paper6_mass_mev(p, q, m_pc, n_q)
    return float("nan") if res is None else res


# ---------------------------------------------------------------------------
# σ-orbit → ground-state particle assignment
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OrbitParticleAssignment:
    orbit_id: int
    orbit_name: str
    particle_name: str
    p: int
    q: int
    m_pc: int
    n_q: int
    mass_pred_MeV: float
    mass_obs_MeV: float | None
    error_pct: float | None


def ground_state_per_orbit() -> list[OrbitParticleAssignment]:
    """For each σ-orbit, find the lowest-mass compatible compendium
    particle.  This is the ground-state assignment.
    """
    out = []
    for orbit_id in range(7):
        compat = compatible_particles(orbit_id)
        if not compat:
            continue
        # Sort by observed mass; missing m_obs → push to end
        compat_sorted = sorted(
            compat,
            key=lambda e: e.get("m_obs", float("inf")) or float("inf"),
        )
        ground = compat_sorted[0]
        m_obs = ground.get("m_obs")
        m_pred = particle_mass_MeV(
            ground["p"], ground["q"], ground["m"], ground["n_q"])
        if m_obs and not math.isnan(m_pred):
            err = (m_pred - m_obs) / m_obs * 100.0
        else:
            err = None
        out.append(OrbitParticleAssignment(
            orbit_id=orbit_id,
            orbit_name=SIGMA_ORBITS[orbit_id]["name"],
            particle_name=ground["name"],
            p=ground["p"], q=ground["q"],
            m_pc=ground["m"], n_q=ground["n_q"],
            mass_pred_MeV=m_pred,
            mass_obs_MeV=m_obs,
            error_pct=err,
        ))
    return out
