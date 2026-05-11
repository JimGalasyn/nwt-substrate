"""Observable assembly — polynomials in so(7) trace invariants.

Each chemistry/physics observable is expressed as a polynomial of
the trace invariants Tr(M²), Tr(M⁴), Tr(M⁶), with the molecule/process-
specific data feeding in via the ring-adjacency (chemistry) or
vertex-process (qed/qcd) input matrix.

This is the layer where the substrate's universal alphabet (the
21 so(7) generators) gets translated into observable-specific
languages.

Currently implemented:
    - aromaticity_score: chemistry RE prediction (Tier B)
    - hopf_pair_count: π-pair count from edge structure
    - k7_indicator: boolean K_7-hub detection via Tr(M²) threshold

Each function accepts either a single 7×7 matrix or a batched
(N, 7, 7) array, with optional backend selection.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from .backends import Backend, BackendName, get_backend
from .batched import graphs_to_invariants, trace_invariants_batched
from .constants import (
    DIM_S_SPIN7,
    DIM_V_SPIN7,
    DIM_ADJ_SPIN7,
    N_VERTICES_K7,
    N_EDGES_K7,
    SPINOR_VECTOR_RATIO,
    WILSON_EXPONENT_K7,
    NLO_VERTEX_COEFFICIENT,
    NNLO_BRACKET_COEFFICIENT,
    TR_M2_W6,
)


# Calibration constants (matching the existing chem.smiles_resonance_energy)
CALIBRATION_KCAL_PER_PI_PAIR: float = 12.0
K7_STABILIZATION_KCAL: float = 56.0

# K_7 indicator threshold: W_6 wheel graph has Tr(M²) = -24 = -2 × 12
# (12 edges in coronene's ring-adjacency).  Sourced from substrate
# constants (TR_M2_W6) — DO NOT redefine here.
K7_TR_M2_THRESHOLD: float = float(TR_M2_W6)


# ---------------------------------------------------------------------------
# K_7 hub indicator
# ---------------------------------------------------------------------------

def k7_indicator(adj_batch: np.ndarray,
                  n_vertices_used: np.ndarray | None = None,
                  backend: BackendName | Backend = "numpy") -> np.ndarray:
    """Indicator for K_7-toroidal hub presence.

    Returns 1 if the 7-vertex projection of the ring-adjacency graph
    saturates the W_6 wheel signature: Tr(M²) ≤ -24 (= 12 or more
    edges in the 7-vertex neighborhood), AND all 7 vertices are
    actually used (no zero-padding).

    Parameters
    ----------
    adj_batch : array-like of shape (N, 7, 7)
        Symmetric 0/1 adjacency. Vertices not used by the ring system
        should be zero-padded (all-zero rows + cols).
    n_vertices_used : array-like of shape (N,) or None
        Number of non-trivial vertices per graph. If None, inferred
        from which rows are all-zero.
    backend : backend name or descriptor

    Returns
    -------
    np.ndarray of shape (N,), dtype=int (values 0 or 1).
    """
    be = backend if isinstance(backend, Backend) else get_backend(backend)
    A = np.asarray(adj_batch)
    if A.ndim == 2:
        A = A[None, ...]

    # Compute Tr(M²) via batched ISA kernel
    inv = graphs_to_invariants(A, orders=(2,), backend=be)
    tr_m2 = inv["Tr_M2"]

    # Determine n_vertices_used from adjacency if not given
    if n_vertices_used is None:
        # Row k is "used" iff it has at least one nonzero entry
        nonzero_rows = (A.sum(axis=-1) > 0)
        n_used = nonzero_rows.sum(axis=-1)
    else:
        n_used = np.asarray(n_vertices_used)

    return ((tr_m2 <= K7_TR_M2_THRESHOLD) & (n_used == 7)).astype(np.int64)


# ---------------------------------------------------------------------------
# Aromaticity-RE
# ---------------------------------------------------------------------------

def aromaticity_score(adj_batch: np.ndarray,
                       pi_pairs: np.ndarray,
                       n_vertices_used: np.ndarray | None = None,
                       calibration_kcal: float = CALIBRATION_KCAL_PER_PI_PAIR,
                       k7_kcal: float = K7_STABILIZATION_KCAL,
                       backend: BackendName | Backend = "numpy") -> np.ndarray:
    """Aromatic resonance energy (kcal/mol) via so(7) ISA.

    Substrate decomposition:

        RE = calibration_kcal × n_pi_pairs  +  k7_kcal × K_7_indicator

    where K_7_indicator(G) = 1 iff the ring-adjacency saturates the
    W_6 wheel signature (Tr(M²) ≤ -24 on 7 vertices).

    Parameters
    ----------
    adj_batch : array-like of shape (N, 7, 7)
        Symmetric 0/1 ring-adjacency, 7-vertex projection, zero-padded
        for systems with < 7 rings.
    pi_pairs : array-like of shape (N,)
        Number of π-electron pairs per aromatic ring system (precomputed
        by the shim, since π-count is a chemistry-specific concept
        not a graph invariant).
    n_vertices_used : array-like of shape (N,) or None
    calibration_kcal, k7_kcal : float
        Calibration constants.
    backend : backend name or descriptor

    Returns
    -------
    np.ndarray of shape (N,), RE in kcal/mol.
    """
    be = backend if isinstance(backend, Backend) else get_backend(backend)
    pi = np.asarray(pi_pairs, dtype=np.float64)
    k7 = k7_indicator(adj_batch, n_vertices_used=n_vertices_used, backend=be)
    return pi * calibration_kcal + k7.astype(np.float64) * k7_kcal


# ---------------------------------------------------------------------------
# Hopf-pair count via Tr(M²)  (second observable)
# ---------------------------------------------------------------------------

def hopf_pair_count(adj_batch: np.ndarray,
                     pi_electrons_total: np.ndarray,
                     backend: BackendName | Backend = "numpy") -> dict[str, np.ndarray]:
    """Hopf-pair structural decomposition from ring-adjacency.

    Returns:
        n_pi_pairs   = π_electrons / 2  (= "Hopf count" in substrate language)
        n_edges      = -Tr(M²)/2        (= ring-adjacency edges; structural)
        parity       = n_pi_pairs % 2    (1 = Hückel-aromatic, 0 = anti-)

    Note: n_pi_pairs and parity depend on chemistry data (which atoms
    contribute π electrons); n_edges is a pure substrate-graph invariant.
    Returning both makes the substrate vs chemistry split visible.

    Parameters
    ----------
    adj_batch : array-like of shape (N, 7, 7)
    pi_electrons_total : array-like of shape (N,)
    backend : backend name or descriptor

    Returns
    -------
    dict with keys: n_pi_pairs, n_edges, parity, tr_m2
    """
    be = backend if isinstance(backend, Backend) else get_backend(backend)
    pi = np.asarray(pi_electrons_total, dtype=np.int64)
    inv = graphs_to_invariants(adj_batch, orders=(2,), backend=be)
    tr_m2 = inv["Tr_M2"]
    n_pairs = pi // 2
    return {
        "n_pi_pairs": n_pairs,
        "n_edges": (-tr_m2 / 2).astype(np.int64),
        "parity": (n_pairs % 2).astype(np.int64),
        "tr_m2": tr_m2,
    }


# ---------------------------------------------------------------------------
# K_7 Wilson amplitude — substrate gravity coupling
# ---------------------------------------------------------------------------

def k7_wilson_amplitude(alpha: float | np.ndarray,
                         order: str = "NNLO") -> float | np.ndarray:
    """K_7 Wilson amplitude on the complete graph K_7 (Paper 17).

    Each edge of K_7 carries √α; the 21 edges give α^(21/2) at LO.
    Spin(7) representation theory contributes structural prefactors:

      LO:    W = (8/7) × α^(21/2)
      NLO:   W = (8/7) × α^(21/2) × (1 + α/7)
      NNLO:  W = (8/7) × α^(21/2) × (1 + α/7 + (21/8) α²)

    Each factor traces to substrate primitives:
      8/7    = DIM_S_SPIN7 / DIM_V_SPIN7
      21/2   = N_EDGES_K7 / 2 (√α per edge)
      1/7    = 1 / N_VERTICES_K7 (NLO per-vertex correction)
      21/8   = DIM_ADJ_SPIN7 / DIM_S_SPIN7 (NNLO bracket coefficient)

    This function gives m_e / M_Pl_predicted directly (the Sakharov
    gravity coupling — Paper 17 §6).

    Parameters
    ----------
    alpha : scalar or array
        Fine-structure-like coupling at the relevant scale.
    order : {"LO", "NLO", "NNLO"}, default "NNLO"

    Returns
    -------
    Wilson amplitude W (same type as alpha).
    """
    base = SPINOR_VECTOR_RATIO * (alpha ** WILSON_EXPONENT_K7)
    if order == "LO":
        return base
    if order == "NLO":
        return base * (1.0 + NLO_VERTEX_COEFFICIENT * alpha)
    if order == "NNLO":
        return base * (
            1.0
            + NLO_VERTEX_COEFFICIENT * alpha
            + NNLO_BRACKET_COEFFICIENT * (alpha ** 2)
        )
    raise ValueError(f"Unknown order: {order!r}; expected LO/NLO/NNLO")


def k7_wilson_breakdown(alpha: float, order: str = "NNLO") -> dict[str, float]:
    """Per-factor structural decomposition of the K_7 Wilson amplitude.

    Useful for cross-shim verification: any other code path computing
    m_e/M_Pl, G_Newton, or related substrate gravity observables should
    be able to reproduce these factors EXACTLY from the same ISA constants.

    Returns dict with keys:
      spinor_vector_ratio : (8/7) prefactor from Spin(7)
      bare_K7_amplitude   : α^(21/2)
      nlo_correction      : (1 + α/7) factor
      nnlo_correction     : (1 + α/7 + (21/8)·α²) factor
      total               : full amplitude at requested order
    """
    bare = float(alpha) ** WILSON_EXPONENT_K7
    nlo = 1.0 + NLO_VERTEX_COEFFICIENT * alpha
    nnlo = nlo + NNLO_BRACKET_COEFFICIENT * (alpha ** 2)
    if order == "LO":
        total = SPINOR_VECTOR_RATIO * bare
    elif order == "NLO":
        total = SPINOR_VECTOR_RATIO * bare * nlo
    elif order == "NNLO":
        total = SPINOR_VECTOR_RATIO * bare * nnlo
    else:
        raise ValueError(f"Unknown order: {order!r}")
    return {
        "spinor_vector_ratio": SPINOR_VECTOR_RATIO,
        "bare_K7_amplitude": bare,
        "nlo_correction": nlo,
        "nnlo_correction": nnlo,
        "total": total,
    }


# ---------------------------------------------------------------------------
# Graph classifier — which "amino acid" signature is this?
# ---------------------------------------------------------------------------

def classify_signature(adj_batch: np.ndarray,
                        backend: BackendName | Backend = "numpy"
                        ) -> list[str]:
    """For each graph in the batch, return the canonical signature
    it matches (by graph-invariant closed-walk counts).

    Uses the SYMMETRIC adjacency A's trace powers (= closed walks of
    length k), not the antisymmetric M's. The reason: Tr(A^k) is
    invariant under vertex relabeling (graph isomorphism), while
    Tr(M^k) for k ≥ 4 depends on the labeling-induced edge orientation.

    Tr(M²) and Tr(A²) carry the same information (= 2|E| up to sign),
    and that suffices for the K_7-hub indicator. For finer
    classification we use higher Tr(A^k).

    Returns "unclassified" if no exact match against the canonical
    7-vertex graph library.
    """
    from .so7 import CANONICAL_GRAPHS

    # Build reference Tr(A^k) for each canonical graph
    refs: dict[str, tuple[float, float, float]] = {}
    for name, edges in CANONICAL_GRAPHS.items():
        A_ref = np.zeros((7, 7), dtype=np.float64)
        for a, b in edges:
            A_ref[a, b] = A_ref[b, a] = 1.0
        A2 = A_ref @ A_ref
        A4 = A2 @ A2
        A6 = A4 @ A2
        refs[name] = (
            float(np.trace(A2)),
            float(np.trace(A4)),
            float(np.trace(A6)),
        )

    A = np.asarray(adj_batch, dtype=np.float64)
    if A.ndim == 2:
        A = A[None, ...]
    # Batched trace powers of the SYMMETRIC adjacency
    A2 = np.einsum("nij,njk->nik", A, A)
    A4 = np.einsum("nij,njk->nik", A2, A2)
    A6 = np.einsum("nij,njk->nik", A4, A2)
    tr2 = np.einsum("nii->n", A2)
    tr4 = np.einsum("nii->n", A4)
    tr6 = np.einsum("nii->n", A6)

    out = []
    for n in range(len(tr2)):
        match = "unclassified"
        for name, (r2, r4, r6) in refs.items():
            if (abs(tr2[n] - r2) < 1e-6
                    and abs(tr4[n] - r4) < 1e-6
                    and abs(tr6[n] - r6) < 1e-6):
                match = name
                break
        out.append(match)
    return out
