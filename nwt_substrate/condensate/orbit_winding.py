"""K_7 Heffter edge winding on the universal cover of T².

Phase E of the Bogoliubov spectrum derivation.

The K_7 Heffter embedding places vertex k at position
(u_k, v_k) = (k/7, 3k mod 7 / 7) on the unit torus.  Each edge has a
specific minimum-length path on the universal cover; the (Δu, Δv)
displacement of that path is a topological invariant of the edge in
the Heffter embedding.

For a closed walk on K_7, the total (Δu, Δv) gives the walk's
homology class in H_1(T²) = Z².  This is the (p, q) winding of any
torus knot supported on that walk.

Phase E hypothesis: the (p, q) quantum numbers of a particle are
the homology class of its supporting closed walk on K_7.  σ-orbit
topology selects which (p, q) classes are reachable.

References:
  - Paper 11 §III mass formula (p, q, m, n_q inputs)
  - Paper 16 §L_3 Derrick + (p,q) torus-knot ansatz
  - [[bogoliubov_phase_d2_complete]] σ-orbit ↔ particle compatibility
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS


from nwt_substrate.isa.constants import N_VERTICES_K7

# K_7 Heffter vertex coordinates on the unit torus (u, v) ∈ [0, 1)²
HEFFTER_VERT_UV = np.array(
    [[(k % N_VERTICES_K7) / N_VERTICES_K7, (3 * k % N_VERTICES_K7) / N_VERTICES_K7] 
     for k in range(N_VERTICES_K7)],
    dtype=float,
)


def short_displacement(a: int, b: int) -> tuple[float, float]:
    """Shortest-path (Δu, Δv) for the edge from vertex a to vertex b
    on the universal cover of the unit torus, using K_7 Heffter coords.

    Returns the (Δu, Δv) with each component in (-1/2, 1/2].
    """
    u_a, v_a = HEFFTER_VERT_UV[a]
    u_b, v_b = HEFFTER_VERT_UV[b]
    du = u_b - u_a
    dv = v_b - v_a
    # Bring to fundamental domain (-1/2, 1/2]
    while du > 0.5:
        du -= 1.0
    while du <= -0.5:
        du += 1.0
    while dv > 0.5:
        dv -= 1.0
    while dv <= -0.5:
        dv += 1.0
    return du, dv


def edge_winding_class(a: int, b: int) -> tuple[int, int]:
    """Integer winding class of the directed edge a → b in 7-step units.

    Returns (n_u, n_v) where Δu_short = n_u/7, Δv_short = n_v/7,
    with n_u, n_v ∈ {-3, -2, -1, 0, 1, 2, 3}.
    """
    du, dv = short_displacement(a, b)
    return round(du * 7), round(dv * 7)


@dataclass(frozen=True)
class OrbitWinding:
    """Topological signature of a σ-orbit's edge set on the Heffter torus."""
    orbit_id: int
    edge_windings: list[tuple[int, int]]    # per-edge (n_u, n_v) in 1/7 units
    total_winding: tuple[int, int]            # sum over edges
    abs_total: tuple[int, int]                # (|sum_u|, |sum_v|)
    p_candidate: int                          # |Δu_total| / 1, rounded
    q_candidate: int                          # |Δv_total| / 1, rounded


def orbit_winding(orbit_id: int) -> OrbitWinding:
    """Compute the winding signature of a σ-orbit's 3 edges, oriented
    in a CONSISTENT direction along the σ-orbit's Z_3 action.

    Convention: each edge oriented from lower-vertex-index to higher.
    """
    orbit = SIGMA_ORBITS[orbit_id]
    edge_windings = []
    sum_u = 0
    sum_v = 0
    for (a, b) in orbit["edges"]:
        n_u, n_v = edge_winding_class(a, b)
        edge_windings.append((n_u, n_v))
        sum_u += n_u
        sum_v += n_v
    return OrbitWinding(
        orbit_id=orbit_id,
        edge_windings=edge_windings,
        total_winding=(sum_u, sum_v),
        abs_total=(abs(sum_u), abs(sum_v)),
        p_candidate=abs(sum_u),
        q_candidate=abs(sum_v),
    )


def all_orbit_windings() -> list[OrbitWinding]:
    return [orbit_winding(o) for o in range(N_VERTICES_K7)]


def closed_walk_winding(walk: list[int]) -> tuple[int, int]:
    """Total winding (n_u, n_v) for a closed walk on K_7 (list of
    vertex indices, last must equal first).

    Sums the short-displacement winding of each directed edge.
    """
    if walk[0] != walk[-1]:
        raise ValueError("Walk must be closed (last == first)")
    sum_u = 0
    sum_v = 0
    for i in range(len(walk) - 1):
        n_u, n_v = edge_winding_class(walk[i], walk[i + 1])
        sum_u += n_u
        sum_v += n_v
    # Convert to actual torus windings (divide by 7 — the walk closes
    # iff sum_u and sum_v are multiples of 7)
    p = sum_u // 7
    q = sum_v // 7
    return p, q
