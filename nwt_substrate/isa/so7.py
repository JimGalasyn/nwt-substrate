"""so(7) primitives — basis, edge-graph embedding, trace invariants.

The Lie algebra so(7) has dim 21 = |E(K_7)|. The 21 generators
T^{ab} (a < b) of so(7) correspond one-to-one with the 21 edges of
the complete graph K_7. A graph G on 7 vertices embeds into so(7)
as the antisymmetric 7×7 matrix

    M_G = Σ_{(i,j) ∈ E(G)} T^{ij}

with M_G[i, j] = +1 if (i, j) ∈ E(G) and i < j, else -1 if j < i.

Trace invariants Tr(M^{2k}) are the SO(7)-invariant polynomial
functions on so(7). Odd-order traces vanish identically for any
antisymmetric matrix. Even-order traces count closed walks of
length 2k in the graph.

This module is numpy-only (no torch). The batched versions live in
`batched.py` and dispatch via the backend layer.
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# so(7) basis
# ---------------------------------------------------------------------------

def so7_basis() -> np.ndarray:
    """Return the 21 antisymmetric 7×7 generators of so(7).

    Convention: indexed by ordered pair (a, b) with 0 ≤ a < b ≤ 6,
    ordered lexicographically. So index k = 0 → (0,1), k = 1 → (0,2),
    ..., k = 20 → (5,6).

    Returns
    -------
    np.ndarray of shape (21, 7, 7)
    """
    gens = np.zeros((21, 7, 7), dtype=np.float64)
    k = 0
    for a in range(7):
        for b in range(a + 1, 7):
            gens[k, a, b] = 1.0
            gens[k, b, a] = -1.0
            k += 1
    assert k == 21
    return gens


# Cached basis (read-only)
SO7_BASIS: np.ndarray = so7_basis()
SO7_BASIS.flags.writeable = False


# Edge-index lookup
def edge_index(a: int, b: int) -> int:
    """Return the so(7) basis index for the edge (a, b)."""
    if a == b or a < 0 or b > 6:
        raise ValueError(f"Invalid edge ({a}, {b}); need 0 ≤ a < b ≤ 6")
    if a > b:
        a, b = b, a
    # Lex index: k = sum_{i<a} (6-i) + (b - a - 1)
    # Simpler: precompute
    return _EDGE_INDEX[(a, b)]


_EDGE_INDEX: dict[tuple[int, int], int] = {}
_k = 0
for _a in range(7):
    for _b in range(_a + 1, 7):
        _EDGE_INDEX[(_a, _b)] = _k
        _k += 1


def k7_edge_list() -> list[tuple[int, int]]:
    """Return the 21 edges of K_7 as ordered (i, j) tuples with i < j.

    Use this anywhere you need to iterate over K_7 edges (heron CZ gates,
    chemistry K_7 sub-graph detection, Wilson loop amplitudes, etc.)
    instead of hardcoding `for i in range(7): for j in range(i+1, 7)`.

    Returns
    -------
    list of (i, j) tuples, length exactly N_EDGES_K7 = 21.
    """
    return [(a, b) for a in range(7) for b in range(a + 1, 7)]


# ---------------------------------------------------------------------------
# Graph embedding
# ---------------------------------------------------------------------------

def adjacency_to_so7(adj: np.ndarray) -> np.ndarray:
    """Convert a 7×7 symmetric {0, 1} adjacency matrix to its
    antisymmetric so(7) representative.

    For undirected edge (i, j) with i < j: result[i, j] = +1, result[j, i] = -1.

    Parameters
    ----------
    adj : np.ndarray of shape (7, 7) or (..., 7, 7)
        Symmetric 0/1 adjacency. The diagonal and lower triangle are
        ignored; only the upper triangle is read.

    Returns
    -------
    np.ndarray of same shape as input, antisymmetric.
    """
    adj = np.asarray(adj, dtype=np.float64)
    if adj.shape[-1] != 7 or adj.shape[-2] != 7:
        raise ValueError(f"Expected (..., 7, 7); got {adj.shape}")
    # Upper triangle: take adj[..., i, j] for i < j
    iu = np.triu_indices(7, k=1)
    M = np.zeros_like(adj)
    M[..., iu[0], iu[1]] = adj[..., iu[0], iu[1]]
    M[..., iu[1], iu[0]] = -adj[..., iu[0], iu[1]]
    return M


def edge_list_to_so7(edges: list[tuple[int, int]],
                      n_vertices: int = 7) -> np.ndarray:
    """Build an so(7) element from an edge list.

    Parameters
    ----------
    edges : list of (a, b) tuples with 0 ≤ a < b < n_vertices
    n_vertices : int, default 7
        Number of vertices (must be ≤ 7).

    Returns
    -------
    np.ndarray of shape (7, 7), antisymmetric.
    """
    if n_vertices > 7:
        raise ValueError(f"n_vertices={n_vertices} > 7; not in so(7)")
    M = np.zeros((7, 7), dtype=np.float64)
    for a, b in edges:
        if a > b:
            a, b = b, a
        if a < 0 or b >= n_vertices:
            raise ValueError(f"Edge ({a}, {b}) out of range [0, {n_vertices})")
        M[a, b] = 1.0
        M[b, a] = -1.0
    return M


# ---------------------------------------------------------------------------
# Trace invariants (single-matrix path; batched version in batched.py)
# ---------------------------------------------------------------------------

def trace_invariants(M: np.ndarray, orders: tuple[int, ...] = (2, 4, 6)) -> dict[str, float]:
    """Compute SO(7)-invariant trace polynomials Tr(M^k).

    For antisymmetric M, odd-k traces vanish identically; only even
    k yields a non-trivial invariant.

    Parameters
    ----------
    M : np.ndarray of shape (7, 7), antisymmetric
    orders : tuple of even ints (default (2, 4, 6))

    Returns
    -------
    dict mapping f"Tr_M{k}" → float
    """
    max_order = max(orders) if orders else 0
    powers: dict[int, np.ndarray] = {1: M}
    # Build powers up to max_order via M^k = M^{k-1} @ M
    for k in range(2, max_order + 1):
        # Prefer squaring when possible
        if k % 2 == 0 and (k // 2) in powers:
            powers[k] = powers[k // 2] @ powers[k // 2]
        else:
            powers[k] = powers[k - 1] @ M
    return {f"Tr_M{k}": float(np.trace(powers[k])) for k in orders}


# ---------------------------------------------------------------------------
# Canonical graphs (reference invariants)
# ---------------------------------------------------------------------------

CANONICAL_GRAPHS: dict[str, list[tuple[int, int]]] = {
    "empty": [],
    "single_edge": [(0, 1)],
    "path_3": [(0, 1), (1, 2)],
    "triangle": [(0, 1), (1, 2), (0, 2)],
    "C_7": [(i, (i + 1) % 7) for i in range(7)],
    "K_{1,6}": [(0, i) for i in range(1, 7)],
    "W_6": (
        [(0, i) for i in range(1, 7)] +
        [(i, i + 1) for i in range(1, 6)] + [(1, 6)]
    ),
    "K_7": [(a, b) for a in range(7) for b in range(a + 1, 7)],
}


def canonical_invariants() -> dict[str, dict[str, float]]:
    """Return Tr(M²), Tr(M⁴), Tr(M⁶) for canonical 7-vertex graphs.

    These are the "amino acid signatures" of the substrate ISA —
    every aromatic ring system's 7-vertex projection reduces to one
    of these (or interpolates between them).
    """
    return {
        name: trace_invariants(edge_list_to_so7(edges))
        for name, edges in CANONICAL_GRAPHS.items()
    }
