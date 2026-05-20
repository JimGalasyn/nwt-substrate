"""K_7 σ-orbit enumeration and projection of Bogoliubov spectrum.

The K_7 graph has 7 vertices and 21 edges.  The Heffter embedding
([[framework_k7_heffter_embedding]]) singles out a polar vertex P
and 6 equatorial vertices arranged as two triangles E_1,E_2,E_3 and
F_1,F_2,F_3.  The subgroup Z_3 × Z_2 ⊂ S_7 = Aut(K_7) fixes P and
permutes the equatorial vertices, partitioning the 21 edges into
**7 σ-orbits of 3 edges each**:

  orbit  edges                          interpretation
  -----  ----------------------------   --------------------------------
  0      {P-E_1, P-E_2, P-E_3}          bridging (matter σ-orbit)
  1      {P-F_1, P-F_2, P-F_3}          bridging (Universal Z_2 partner)
  2      {E_1-F_1, E_2-F_2, E_3-F_3}    τ-fixed cross (Z_2-neutral)
  3      {E_1-E_2, E_2-E_3, E_3-E_1}    upper triangle (matter triangle)
  4      {F_1-F_2, F_2-F_3, F_3-F_1}    lower triangle (Z_2 partner)
  5      {E_1-F_2, E_2-F_3, E_3-F_1}    cross-twist + (Z_3 forward)
  6      {E_1-F_3, E_2-F_1, E_3-F_2}    cross-twist − (Z_3 backward)

Total: 3 × 7 = 21 edges ✓

The 3 bridging edges of orbit 0 give f_J = (1-√α)³ via Wilson product
(Layer 3 Phase F).  Orbits 0 and 1 are Z_2-related; orbits 5 and 6 are
Z_3-related.  Orbit 2 is fixed by Z_2 polarity.

Each σ-orbit defines a projection of the abelian-Higgs Bogoliubov
spectrum into a substrate quasiparticle species.  The projection
operator is determined by the edge-count and topological invariants
of the σ-orbit; the resulting effective mass is the species' rest mass.

This module:
  - Enumerates the 7 σ-orbits explicitly
  - Computes substrate-natural σ-orbit invariants
  - Provides σ-orbit projection of the abelian-Higgs spectrum
  - Identifies the σ-orbit that yields ξ = λ̄_C (substrate healing length)

References: [[framework_k7_heffter_embedding]],
[[framework_omega_b_c_nlo_truncation]],
[[research_plan_bogoliubov_spectrum]] Phase C.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from nwt_substrate.condensate.abelian_higgs import (
    AbelianHiggsParams, HBARC_EV_M, M_E_EV, LAMBDA_C_M,
)
from nwt_substrate.isa.constants import ALPHA_NWT

# ---------------------------------------------------------------------------
# K_7 vertex labels and σ-orbit edge sets
# ---------------------------------------------------------------------------

# Vertex labels: 0 = P (polar), 1-3 = E_1..E_3, 4-6 = F_1..F_3
P, E1, E2, E3, F1, F2, F3 = 0, 1, 2, 3, 4, 5, 6
VERTEX_LABELS = ["P", "E_1", "E_2", "E_3", "F_1", "F_2", "F_3"]


def _edge(a: int, b: int) -> tuple[int, int]:
    return (min(a, b), max(a, b))


SIGMA_ORBITS: dict[int, dict] = {
    0: {
        "name": "bridging (matter)",
        "edges": [_edge(P, E1), _edge(P, E2), _edge(P, E3)],
        "z2_partner": 1,
        "physical_role": "matter σ-orbit (3 bridging edges, f_J Wilson)",
    },
    1: {
        "name": "bridging (Universal Z_2 partner)",
        "edges": [_edge(P, F1), _edge(P, F2), _edge(P, F3)],
        "z2_partner": 0,
        "physical_role": "Z_2 partner of matter (dark matter candidate)",
    },
    2: {
        "name": "τ-fixed cross",
        "edges": [_edge(E1, F1), _edge(E2, F2), _edge(E3, F3)],
        "z2_partner": 2,  # self-conjugate under Z_2
        "physical_role": "Z_2-neutral channel (radiation/photon candidate)",
    },
    3: {
        "name": "upper triangle (E)",
        "edges": [_edge(E1, E2), _edge(E2, E3), _edge(E3, E1)],
        "z2_partner": 4,
        "physical_role": "equatorial torus, upper",
    },
    4: {
        "name": "lower triangle (F)",
        "edges": [_edge(F1, F2), _edge(F2, F3), _edge(F3, F1)],
        "z2_partner": 3,
        "physical_role": "equatorial torus, lower",
    },
    5: {
        "name": "cross-twist + (Z_3 forward)",
        "edges": [_edge(E1, F2), _edge(E2, F3), _edge(E3, F1)],
        "z2_partner": 6,
        "physical_role": "twisted cross, Z_3 chirality +",
    },
    6: {
        "name": "cross-twist − (Z_3 backward)",
        "edges": [_edge(E1, F3), _edge(E2, F1), _edge(E3, F2)],
        "z2_partner": 5,
        "physical_role": "twisted cross, Z_3 chirality −",
    },
}


def verify_partition() -> dict:
    """Verify the 7 σ-orbits partition the 21 edges of K_7 exactly."""
    all_edges = set()
    for o in SIGMA_ORBITS.values():
        for e in o["edges"]:
            assert e not in all_edges, f"Duplicate edge {e}"
            all_edges.add(e)
    expected = {(i, j) for i in range(7) for j in range(i + 1, 7)}
    return {
        "n_orbits": len(SIGMA_ORBITS),
        "n_edges": len(all_edges),
        "covers_K7": all_edges == expected,
        "missing": expected - all_edges,
        "extra": all_edges - expected,
    }


# ---------------------------------------------------------------------------
# σ-orbit projection of the abelian-Higgs spectrum
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SigmaOrbitInvariants:
    """Topological + algebraic invariants of a σ-orbit.

    n_edges: number of edges (always 3 in K_7 Heffter).
    polar_edges: how many edges include the polar vertex P (0 or 3).
    cross_edges: how many edges connect E_i to F_j (0, 3).
    self_z2_fixed: whether the σ-orbit is fixed by Z_2 polarity (True iff orbit 2).
    z3_chirality: ±1 for the cross-twist orbits, 0 for others.
    """
    orbit_id: int
    n_edges: int
    polar_edges: int
    cross_edges: int
    self_z2_fixed: bool
    z3_chirality: int


def orbit_invariants(orbit_id: int) -> SigmaOrbitInvariants:
    """Compute topological invariants for σ-orbit `orbit_id`."""
    orbit = SIGMA_ORBITS[orbit_id]
    n_edges = len(orbit["edges"])
    polar_edges = sum(1 for (a, b) in orbit["edges"] if P in (a, b))
    cross_edges = sum(
        1 for (a, b) in orbit["edges"]
        if (a in (E1, E2, E3) and b in (F1, F2, F3))
        or (a in (F1, F2, F3) and b in (E1, E2, E3))
    )
    self_z2 = orbit_id == 2
    z3_chirality = 0
    if orbit_id == 5:
        z3_chirality = +1
    elif orbit_id == 6:
        z3_chirality = -1
    return SigmaOrbitInvariants(
        orbit_id=orbit_id,
        n_edges=n_edges,
        polar_edges=polar_edges,
        cross_edges=cross_edges,
        self_z2_fixed=self_z2,
        z3_chirality=z3_chirality,
    )


def wilson_product_factor(orbit_id: int, alpha: float = ALPHA_NWT) -> float:
    """Wilson product over a σ-orbit's edges, in substrate-natural form.

    Each edge contributes a holonomy factor (1 - √α) (Layer 3 Phase F).
    For 3 bridging edges, the product is (1 - √α)³ = f_J ≈ 0.765.

    For OTHER σ-orbits, the same Wilson product applies — the orbit's
    role is set by topology (polar vs cross vs equatorial), not by the
    Wilson factor.  All 7 σ-orbits give Wilson product = (1 - √α)³ = f_J.
    """
    return (1.0 - math.sqrt(alpha)) ** 3


def effective_mass(orbit_id: int,
                   p: AbelianHiggsParams,
                   scheme: str = "polar_weighted") -> float:
    """σ-orbit-projected effective mass (eV).

    Several candidate schemes:
      - 'higgs': just m_σ (no projection)
      - 'gauge': just m_A
      - 'polar_weighted': mass scales with (polar_edges/3) — projects
        bridging (matter) orbits onto m_σ and non-polar onto m_A or
        a combination.
      - 'topological': mass = m_σ × (1 - √α)^{polar_edges} contribution
      - 'BPS_combined': m_σ + (n_cross - n_polar) × (m_A - m_σ)/2

    The PHYSICALLY MEANINGFUL scheme should produce m_e for the matter
    σ-orbit (orbit_id = 0).  We test multiple candidates and report
    which (if any) yields ξ = λ̄_C.

    For Phase C this is exploratory — Phase D will fix the scheme from
    first principles.
    """
    inv = orbit_invariants(orbit_id)
    m_sigma = p.m_sigma_eV
    m_A = p.m_A_eV
    if scheme == "higgs":
        return m_sigma
    if scheme == "gauge":
        return m_A
    if scheme == "polar_weighted":
        # f = (polar_edges / n_edges)
        f = inv.polar_edges / inv.n_edges
        return f * m_A + (1 - f) * m_sigma
    if scheme == "topological":
        # mass = m_σ / (1 - √α)^(polar_edges)
        # boosts polar (bridging) orbits to higher mass
        return m_sigma / (1 - math.sqrt(ALPHA_NWT)) ** inv.polar_edges
    if scheme == "BPS_inverse":
        # σ-orbit projection inverts the BPS factor:
        # m_eff = m_σ × √(2π/α^{polar}) — boost by 1/√α for each polar edge
        return m_sigma * (1.0 / math.sqrt(ALPHA_NWT)) ** (inv.polar_edges / 3.0)
    if scheme == "wilson_normalized":
        # mass = m_σ / Wilson_factor — restores energy lost in Wilson product
        wf = wilson_product_factor(orbit_id)
        return m_sigma / wf
    raise ValueError(f"Unknown scheme '{scheme}'")


def effective_correlation_length(orbit_id: int, p: AbelianHiggsParams,
                                  scheme: str = "polar_weighted") -> float:
    """Correlation length ξ = ℏ/(m_eff c), in meters."""
    m = effective_mass(orbit_id, p, scheme)
    return HBARC_EV_M / m


def find_matter_orbit_mass_match(p: AbelianHiggsParams,
                                  schemes: list[str] | None = None) -> dict:
    """Scan candidate projection schemes; report which produces mass m_e
    for the matter σ-orbit (orbit_id = 0).
    """
    if schemes is None:
        schemes = ["higgs", "gauge", "polar_weighted", "topological",
                   "BPS_inverse", "wilson_normalized"]
    results = {}
    for s in schemes:
        m = effective_mass(0, p, s)
        ratio = m / M_E_EV
        xi = HBARC_EV_M / m
        results[s] = {
            "m_eff_eV": m,
            "m_eff_over_m_e": ratio,
            "xi_m": xi,
            "xi_over_lambda_C": xi / LAMBDA_C_M,
            "matches_m_e": abs(ratio - 1.0) < 0.05,
        }
    return results
