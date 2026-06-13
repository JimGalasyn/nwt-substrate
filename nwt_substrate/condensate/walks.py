"""Closed walks on K_7 — homology (p, q), invariants, and BFS enumeration.

Library consolidation of utilities used across the Bogoliubov spectrum
derivation (Phases E-K).  Previously duplicated across analysis scripts.

Key concepts (per Phase E):

  - K_7 walk = sequence of vertex indices [v_0, v_1, ..., v_L] with v_0 = v_L
  - Each edge (v_i, v_{i+1}) is one of 21 K_7 edges
  - The walk's (p, q) torus-knot winding = total Δu and Δv (in 1/7 units) /7
  - The walk's σ-orbit composition partitions edges into 7 σ-orbits

This module supplies:

  - BFS enumeration of SHORTEST closed walks per (|p|, |q|) homology class
  - Walk-level invariants (arc length, step distributions, σ-orbit composition)
  - 3D lift to Heffter torus for pyknotid identification
  - σ-orbit transition analysis (carrier-knot precursor)

For ANTIMATTER walk operations, see ``nwt_substrate.condensate.antimatter``.

References:
  - [[bogoliubov-phase-e-complete]] — (p, q) walk homology, 100% closed
  - [[bogoliubov-phase-i-pyknotid]] — pyknotid 3D-lift identification
  - [[research-plan-bogoliubov-spectrum]] — multi-phase plan
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass

import numpy as np

from nwt_substrate.isa.constants import N_VERTICES_K7
from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS
from nwt_substrate.condensate.orbit_winding import (
    edge_winding_class, HEFFTER_VERT_UV,
)


# ---------------------------------------------------------------------------
# Edge → σ-orbit lookup (used in 4+ analysis scripts)
# ---------------------------------------------------------------------------

def edge_to_orbit(a: int, b: int) -> int:
    """Return the σ-orbit ID (0-6) containing edge (a, b), or -1 if not found."""
    e_norm = (min(a, b), max(a, b))
    for oid, orbit in SIGMA_ORBITS.items():
        if e_norm in orbit["edges"]:
            return oid
    return -1


# ---------------------------------------------------------------------------
# BFS enumeration of shortest closed walks per (|p|, |q|)
# ---------------------------------------------------------------------------

def bfs_shortest_walks(max_length: int = 25,
                       start: int = 0) -> dict[tuple[int, int], list[int]]:
    """BFS in (vertex, accumulated_winding) state space to find the shortest
    closed walk at vertex `start` realizing each (|p|, |q|) homology class.

    State: (current_vertex, sum_Δu, sum_Δv) where sum_Δu, sum_Δv are
    accumulated in 1/7 units.  Walk closes on the universal cover when
    final state = (start, 7p, 7q) for some integer (p, q).

    Returns: dict mapping (|p|, |q|) → walk (list of vertex indices).
    
    Raises
    ------
    ValueError
        If start is not a valid K_7 vertex index (0-6).
    """
    if start < 0 or start >= N_VERTICES_K7:
        raise ValueError(f"start vertex {start} not in range [0, {N_VERTICES_K7-1}]")
    edge_w = {(a, b): edge_winding_class(a, b)
              for a in range(N_VERTICES_K7) for b in range(N_VERTICES_K7) if a != b}
    initial = (start, 0, 0)
    visited = {initial: (0, None)}
    queue = deque([initial])
    while queue:
        state = queue.popleft()
        depth, _ = visited[state]
        if depth >= max_length:
            continue
        v, m_u, m_v = state
        for nxt in range(N_VERTICES_K7):
            if nxt == v:
                continue
            dnu, dnv = edge_w[(v, nxt)]
            new_state = (nxt, m_u + dnu, m_v + dnv)
            if new_state not in visited:
                visited[new_state] = (depth + 1, state)
                queue.append(new_state)
    walks = {}
    for state, (depth, _) in visited.items():
        v, m_u, m_v = state
        if v != start or (m_u, m_v) == (0, 0):
            continue
        if m_u % 7 != 0 or m_v % 7 != 0:
            continue
        p, q = m_u // 7, m_v // 7
        key = (abs(p), abs(q))
        walk = [state[0]]
        cur = state
        while visited[cur][1] is not None:
            cur = visited[cur][1]
            walk.append(cur[0])
        walk.reverse()
        if key not in walks or len(walk) - 1 < len(walks[key]) - 1:
            walks[key] = walk
    return walks


def closed_walk_winding(walk: list[int]) -> tuple[int | None, int | None]:
    """Total (p, q) winding of a closed walk in integer units.

    Returns (None, None) if walk does not close on universal cover.
    """
    sum_u = sum_v = 0
    for i in range(len(walk) - 1):
        nu, nv = edge_winding_class(walk[i], walk[i + 1])
        sum_u += nu
        sum_v += nv
    if sum_u % 7 != 0 or sum_v % 7 != 0:
        return None, None
    return sum_u // 7, sum_v // 7


# ---------------------------------------------------------------------------
# Walk arc length on unit torus (used in Phase F, H, K)
# ---------------------------------------------------------------------------

def walk_arc_length_unit_torus(walk: list[int]) -> float:
    """Total arc length of walk on the unit torus.

    Each edge has displacement (Δu/7, Δv/7); arc length per edge =
    √((Δu/7)² + (Δv/7)²).
    """
    total = 0.0
    for i in range(len(walk) - 1):
        nu, nv = edge_winding_class(walk[i], walk[i + 1])
        total += math.sqrt(nu ** 2 + nv ** 2) / 7.0
    return total


# ---------------------------------------------------------------------------
# Step distributions (QR/NR signed and symmetric d)
# ---------------------------------------------------------------------------

# Quadratic residues mod 7 (signed step direction):
QR_SIGNED = {1, 2, 4}
# Non-residues mod 7:
NR_SIGNED = {3, 5, 6}


def walk_signed_step_distribution(walk: list[int]) -> dict[int, int]:
    """Counter of signed step sizes (b - a) mod 7 ∈ {1..6} along the walk."""
    counts: Counter[int] = Counter()
    for i in range(len(walk) - 1):
        d_signed = (walk[i + 1] - walk[i]) % 7
        counts[d_signed] += 1
    return dict(counts)


def walk_symmetric_d_distribution(walk: list[int]) -> dict[int, int]:
    """Counter of symmetric edge-difference d ∈ {1, 2, 3} along the walk."""
    counts: Counter[int] = Counter()
    for i in range(len(walk) - 1):
        d = (walk[i + 1] - walk[i]) % 7
        counts[min(d, 7 - d)] += 1
    return dict(counts)


def walk_QR_NR_fractions(walk: list[int]) -> tuple[float, float]:
    """Fractions of QR-signed vs NR-signed steps in the walk.

    Per Phase J: matter walks are QR-signed (fraction → 1.0 = pure matter
    Hamilton), antimatter walks are NR-signed (fraction → 0.0).
    """
    steps = walk_signed_step_distribution(walk)
    L = sum(steps.values())
    if L == 0:
        return 0.0, 0.0
    qr = sum(c for d, c in steps.items() if d in QR_SIGNED)
    nr = sum(c for d, c in steps.items() if d in NR_SIGNED)
    return qr / L, nr / L


# ---------------------------------------------------------------------------
# σ-orbit decomposition
# ---------------------------------------------------------------------------

def walk_sigma_sequence(walk: list[int]) -> list[int]:
    """σ-orbit ID at each step of the walk."""
    return [edge_to_orbit(walk[i], walk[i + 1])
            for i in range(len(walk) - 1)]


def walk_sigma_composition(walk: list[int]) -> dict[int, int]:
    """Total edge count in each σ-orbit (0-6).  Sums to len(walk) - 1."""
    counts: Counter[int] = Counter()
    for oid in walk_sigma_sequence(walk):
        counts[oid] += 1
    return dict(counts)


def walk_sigma_transitions(walk: list[int],
                            exclude_self_loops: bool = True
                            ) -> list[tuple[int, int]]:
    """List of consecutive (σ_i, σ_j) pairs in the walk's σ-orbit sequence,
    including wrap-around for closed walks.

    If exclude_self_loops=True (default), drops pairs where σ_i == σ_j.
    """
    sig_seq = walk_sigma_sequence(walk)
    n = len(sig_seq)
    out = []
    for i in range(n):
        cur = sig_seq[i]
        nxt = sig_seq[(i + 1) % n]
        if exclude_self_loops and cur == nxt:
            continue
        out.append((cur, nxt))
    return out


def walk_sigma_complete_coverage(walk: list[int]) -> dict[int, bool]:
    """For each σ-orbit, whether the walk uses ALL 3 of its edges.

    Returns {orbit_id: True/False}.
    """
    edges_used = set()
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        edges_used.add((min(a, b), max(a, b)))
    out = {}
    for oid, orbit in SIGMA_ORBITS.items():
        orbit_edges = set(orbit["edges"])
        out[oid] = (orbit_edges & edges_used) == orbit_edges
    return out


# ---------------------------------------------------------------------------
# Walk topology summary (Phase F-3 invariants)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WalkInvariants:
    """Compact bundle of walk invariants used in Phase F-3 sector matching."""
    L: int
    p: int
    q: int
    arc_length: float
    sigma_counts: dict[int, int]
    n_distinct_orbits: int
    n_sigma_complete: int
    n_distinct_vertices: int
    n_return_to_start: int
    bridging_count: int       # σ_0 + σ_1
    triangle_count: int       # σ_3 + σ_4
    cross_count: int          # σ_2 + σ_5 + σ_6
    tc_ratio: float
    QR_frac: float
    NR_frac: float


def walk_invariants(walk: list[int]) -> WalkInvariants:
    """Compute the full Phase F-3 invariant bundle for a closed walk."""
    L = len(walk) - 1
    p, q = closed_walk_winding(walk)
    arc = walk_arc_length_unit_torus(walk)
    sig = walk_sigma_composition(walk)
    sig_complete = walk_sigma_complete_coverage(walk)
    n_dist_orbits = len([o for o, c in sig.items() if c > 0])
    n_complete = sum(1 for v in sig_complete.values() if v)
    n_dist_verts = len(set(walk))
    n_ret = sum(1 for v in walk[1:-1] if v == walk[0]) + 1
    bridging = sig.get(0, 0) + sig.get(1, 0)
    triangle = sig.get(3, 0) + sig.get(4, 0)
    cross = sig.get(2, 0) + sig.get(5, 0) + sig.get(6, 0)
    tc = triangle / cross if cross > 0 else float("inf")
    qr_frac, nr_frac = walk_QR_NR_fractions(walk)
    return WalkInvariants(
        L=L, p=p if p is not None else 0,
        q=q if q is not None else 0,
        arc_length=arc, sigma_counts=sig,
        n_distinct_orbits=n_dist_orbits,
        n_sigma_complete=n_complete,
        n_distinct_vertices=n_dist_verts,
        n_return_to_start=n_ret,
        bridging_count=bridging, triangle_count=triangle,
        cross_count=cross, tc_ratio=tc,
        QR_frac=qr_frac, NR_frac=nr_frac,
    )


# ---------------------------------------------------------------------------
# 3D lift to Heffter torus (for pyknotid identification — Phase I)
# ---------------------------------------------------------------------------

def walk_to_3d_curve(walk: list[int], n_per_edge: int = 30,
                     R: float = 2.5, r: float = 1.0) -> np.ndarray:
    """Lift a K_7 walk to a 3D space curve on the Heffter torus in R³.

    Vertex k sits at (u_k, v_k) = (k/7, 3k mod 7 / 7) on the unit torus.
    Each edge is traversed along the shortest geodesic; cumulative (u, v)
    is unwrapped to keep the curve continuous on the universal cover.

    Parameters
    ----------
    walk : closed walk vertex sequence
    n_per_edge : number of interpolation points per edge
    R, r : major and minor radii of the embedding torus (R > r > 0)

    Returns
    -------
    pts : ndarray, shape (N, 3), N = L · n_per_edge

    Used by pyknotid.spacecurves.Knot for carrier-knot identification.
    Phase I result: this identifies the walk's INTRINSIC torus-knot type
    T(p, q), NOT the Paper 11 carrier knot (which is a separate
    topological object emerging from σ-orbit structure).
    """
    pts = []
    u_curr = HEFFTER_VERT_UV[walk[0]][0]
    v_curr = HEFFTER_VERT_UV[walk[0]][1]
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        nu, nv = edge_winding_class(a, b)
        du = nu / 7.0
        dv = nv / 7.0
        ts = np.linspace(0, 1, n_per_edge, endpoint=False)
        for t in ts:
            u = u_curr + t * du
            v = v_curr + t * dv
            x = (R + r * np.cos(2 * np.pi * v)) * np.cos(2 * np.pi * u)
            y = (R + r * np.cos(2 * np.pi * v)) * np.sin(2 * np.pi * u)
            z = r * np.sin(2 * np.pi * v)
            pts.append((x, y, z))
        u_curr += du
        v_curr += dv
    return np.array(pts)
