"""
Clifford algebra Cl(0,7) via octonion left-multiplication matrices.

The 8x8 matrices L_i (i = 0..7) representing left-multiplication by e_i
satisfy the Cl(0,7) anticommutation relation:

    {L_i, L_j} = -2 delta_ij I_8    for i, j in 1..7

despite octonions being non-associative, because the associator
[e_i, e_j, x] is antisymmetric in (i, j) and cancels in the
anticommutator.
"""

from __future__ import annotations

import numpy as np

from .octonions import basis_vector, octo_mul


def left_mult_matrix(i: int, T: np.ndarray) -> np.ndarray:
    """
    8x8 matrix L_i such that L_i x = e_i * x for octonion x.

    Operates on octonions (8-dim real vectors).  Satisfies:
      L_0 = I_8
      L_i^2 = -I_8 for i = 1..7
      {L_i, L_j} = -2 delta_ij I_8 (Cl(0,7) Clifford algebra)
    """
    L = np.zeros((8, 8))
    e_i = basis_vector(i)
    for j in range(8):
        e_j = basis_vector(j)
        L[:, j] = octo_mul(e_i, e_j, T)
    return L


def right_mult_matrix(i: int, T: np.ndarray) -> np.ndarray:
    """8x8 matrix R_i such that R_i x = x * e_i."""
    R = np.zeros((8, 8))
    e_i = basis_vector(i)
    for j in range(8):
        e_j = basis_vector(j)
        R[:, j] = octo_mul(e_j, e_i, T)
    return R


def bivector_matrix(i: int, j: int, T: np.ndarray) -> np.ndarray:
    """
    Matrix L_(e_i * e_j): left-multiplication by the bivector product e_i * e_j.

    For distinct i, j in 1..7, e_i * e_j = +/- e_k where (i, j, k) is a Fano
    line.  So L_(e_i * e_j) = +/- L_k -- a single imaginary L matrix.

    Note: L_(e_i * e_j) != L_i @ L_j in general; they differ by the action
    of the octonion associator [e_i, e_j, .].  See octonion_dirac_step4d for
    discussion.
    """
    e_i = basis_vector(i)
    e_j = basis_vector(j)
    bivector = octo_mul(e_i, e_j, T)
    L = np.zeros((8, 8))
    for k in range(8):
        L[:, k] = octo_mul(bivector, basis_vector(k), T)
    return L
