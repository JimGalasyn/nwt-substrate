"""Batched ISA kernels — the substrate's compute primitive layer.

The single bottleneck operation is:

    given a batch of N antisymmetric 7×7 matrices M,
    return Tr(M^k) for k ∈ {2, 4, 6}.

Implemented as a chain of einsums that run identically under numpy
or torch (CPU or CUDA). This is the kernel that should be deeply
optimized; everything else in the ISA either prepares input for it
or assembles polynomials of its outputs.

This is the "substrate ribosome" — the universal active-encoding
machine that turns substrate-algebra elements into observables.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from .backends import Backend, BackendName, get_backend


# ---------------------------------------------------------------------------
# Batched adjacency → so(7) embedding
# ---------------------------------------------------------------------------

def adjacency_to_so7_batched(adj_batch: np.ndarray,
                              backend: BackendName | Backend = "numpy") -> Any:
    """Convert a batch of 7×7 symmetric {0,1} adjacency matrices to
    their antisymmetric so(7) representatives.

    Parameters
    ----------
    adj_batch : array-like of shape (N, 7, 7)
        Symmetric 0/1. Only upper triangle is read.
    backend : "numpy" | "torch_cpu" | "torch_cuda" | Backend

    Returns
    -------
    Array of shape (N, 7, 7) on the backend's device, antisymmetric.
    """
    be = backend if isinstance(backend, Backend) else get_backend(backend)
    A = be.asarray(adj_batch, dtype=np.float64)
    xp = be.xp
    # upper triangle indices
    iu_i, iu_j = _triu_indices_7(xp, be.device)
    if be.name == "numpy":
        out = np.zeros_like(A)
        out[..., iu_i, iu_j] = A[..., iu_i, iu_j]
        out[..., iu_j, iu_i] = -A[..., iu_i, iu_j]
        return out
    else:
        out = xp.zeros_like(A)
        out[..., iu_i, iu_j] = A[..., iu_i, iu_j]
        out[..., iu_j, iu_i] = -A[..., iu_i, iu_j]
        return out


def _triu_indices_7(xp, device):
    """Cached upper-triangle indices for 7×7 matrix."""
    iu = np.triu_indices(7, k=1)
    if xp is np:
        return iu[0], iu[1]
    # torch
    import torch
    return (
        torch.as_tensor(iu[0], device=device, dtype=torch.long),
        torch.as_tensor(iu[1], device=device, dtype=torch.long),
    )


# ---------------------------------------------------------------------------
# Batched trace invariants — THE kernel
# ---------------------------------------------------------------------------

def trace_invariants_batched(M_batch: Any,
                               orders: tuple[int, ...] = (2, 4, 6),
                               backend: BackendName | Backend = "numpy") -> dict[str, Any]:
    """Compute Tr(M^k) for a batch of 7×7 matrices.

    This is the SUBSTRATE ISA KERNEL — every observable that factors
    through so(7) trace invariants compiles down to this call.

    Parameters
    ----------
    M_batch : array-like of shape (N, 7, 7)
    orders : tuple of int (default (2, 4, 6))
        Trace orders to compute. Only even orders give non-trivial
        results for antisymmetric input.
    backend : "numpy" | "torch_cpu" | "torch_cuda" | Backend

    Returns
    -------
    dict mapping f"Tr_M{k}" → array of shape (N,)
    """
    be = backend if isinstance(backend, Backend) else get_backend(backend)
    M = be.asarray(M_batch, dtype=np.float64)
    xp = be.xp

    # Compute power chain via repeated squaring where possible
    powers: dict[int, Any] = {1: M}
    if 2 in orders or 4 in orders or 6 in orders or 8 in orders:
        powers[2] = _matmul_batched(M, M, xp)
    if 4 in orders or 6 in orders or 8 in orders:
        powers[4] = _matmul_batched(powers[2], powers[2], xp)
    if 6 in orders:
        powers[6] = _matmul_batched(powers[4], powers[2], xp)
    if 8 in orders:
        powers[8] = _matmul_batched(powers[4], powers[4], xp)

    out: dict[str, Any] = {}
    for k in orders:
        if k not in powers:
            # Fallback: build by repeated multiplication
            powers[k] = M
            for _ in range(k - 1):
                powers[k] = _matmul_batched(powers[k], M, xp)
        # Batched trace = sum of diagonal
        out[f"Tr_M{k}"] = _batched_trace(powers[k], xp)
    return out


def _matmul_batched(A: Any, B: Any, xp) -> Any:
    """Batched matmul over leading batch dimension."""
    if xp is np:
        return np.einsum("nij,njk->nik", A, B)
    else:
        # torch.matmul broadcasts batch dims natively
        return xp.matmul(A, B)


def _batched_trace(M: Any, xp) -> Any:
    """Batched trace: sum of diagonal over leading batch dim."""
    if xp is np:
        return np.einsum("nii->n", M)
    else:
        return xp.einsum("nii->n", M)


# ---------------------------------------------------------------------------
# Convenience: full pipeline adj → invariants
# ---------------------------------------------------------------------------

def graphs_to_invariants(adj_batch: np.ndarray,
                          orders: tuple[int, ...] = (2, 4, 6),
                          backend: BackendName | Backend = "numpy",
                          return_numpy: bool = True) -> dict[str, np.ndarray]:
    """One-shot: batch of {0,1} adjacency matrices → trace invariants.

    Combines `adjacency_to_so7_batched` + `trace_invariants_batched`
    in a single call.

    Parameters
    ----------
    adj_batch : array-like of shape (N, 7, 7)
    orders : tuple of int
    backend : backend name or descriptor
    return_numpy : bool, default True
        If True, results are brought back to host as numpy arrays.
        If False, results stay on the backend device (saves transfer
        for downstream batched operations).

    Returns
    -------
    dict of f"Tr_M{k}" → numpy array of shape (N,) (or backend tensor
    if return_numpy=False).
    """
    be = backend if isinstance(backend, Backend) else get_backend(backend)
    M = adjacency_to_so7_batched(adj_batch, backend=be)
    inv = trace_invariants_batched(M, orders=orders, backend=be)
    if return_numpy:
        return {k: be.to_numpy(v) for k, v in inv.items()}
    return inv
