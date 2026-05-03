"""
Spinor wavefunctions for substrate-algebraic amplitude calculations.

In the Cl(0,7) -> Cl(1,3) embedding, the spinor space is 8-dimensional.
For a 4-momentum p with p^2 = m^2 > 0, the operator slash{p} has:
  - 4 eigenvalues = +m  (the positive-energy / particle subspace),
  - 4 eigenvalues = -m  (the negative-energy / antiparticle subspace).

Each eigenspace is 4-dimensional: 2 physical spin states x 2 internal SU(2).

For massless states (p^2 = 0), slash{p} is nilpotent (slash{p}^2 = 0) and has
8 zero eigenvalues; we use SVD to find a clean orthonormal null basis.
"""

from __future__ import annotations

import numpy as np
from numpy import linalg as LA

from ..algebra.dirac import slash_p


def positive_energy_spinors(p_mu: np.ndarray, gammas, m: float,
                            tol: float = 1e-7) -> np.ndarray:
    """
    Orthonormal basis of the +m eigenspace of slash{p}, shape (4, 8).

    Returns rows u_s with (slash p) u_s = m u_s, normalized to <u_s, u_s'> = delta.
    """
    sp = slash_p(p_mu, gammas)
    eigvals, eigvecs = LA.eig(sp)
    pos_idx = np.where(np.abs(eigvals - m) < tol * max(abs(m), 1.0))[0]
    if len(pos_idx) != 4:
        raise ValueError(
            f"Expected 4 positive-energy spinors, got {len(pos_idx)}"
        )
    Q, _ = LA.qr(eigvecs[:, pos_idx])
    return Q.T


def negative_energy_spinors(p_mu: np.ndarray, gammas, m: float,
                            tol: float = 1e-7) -> np.ndarray:
    """
    Orthonormal basis of the -m eigenspace of slash{p}, shape (4, 8).

    These are the antiparticle spinors v_s with (slash p) v_s = -m v_s.
    """
    sp = slash_p(p_mu, gammas)
    eigvals, eigvecs = LA.eig(sp)
    neg_idx = np.where(np.abs(eigvals + m) < tol * max(abs(m), 1.0))[0]
    if len(neg_idx) != 4:
        raise ValueError(
            f"Expected 4 negative-energy spinors, got {len(neg_idx)}"
        )
    Q, _ = LA.qr(eigvecs[:, neg_idx])
    return Q.T


def null_spinors(p_mu: np.ndarray, gammas, tol: float = 1e-9) -> np.ndarray:
    """
    Orthonormal basis of the null space of slash{p} for null p (p^2 = 0).

    Uses SVD because eigenvalue decomposition of a nilpotent matrix is ill-
    conditioned.  Returns 4 spinors in shape (4, 8) when slash{p} has rank 4.
    """
    sp = slash_p(p_mu, gammas)
    U, sigma, Vh = LA.svd(sp)
    # Null basis = right singular vectors with sigma ~ 0
    null_mask = sigma < tol * max(sigma)
    n_null = int(null_mask.sum())
    if n_null == 0:
        # All sigmas above tol -> map by smallest 4
        idx = np.argsort(sigma)[:4]
        null_basis = Vh[idx].conj()
    else:
        null_basis = Vh[null_mask].conj()
    Q, _ = LA.qr(null_basis.T)
    return Q.T[:max(4, n_null)][:4]   # take first 4


def adjoint_spinor(u: np.ndarray, g0: np.ndarray) -> np.ndarray:
    """
    Dirac adjoint:  u_bar = u^dagger gamma^0.

    Used for building bilinears  u_bar Gamma u  in matrix elements.
    """
    return np.conjugate(u) @ g0
