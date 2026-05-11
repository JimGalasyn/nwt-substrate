"""
SU(3) color algebra for QCD on the substrate.

The 8 Gell-Mann matrices lambda^a (a = 1..8) are the standard generators
of SU(3) in the fundamental (3-dim) representation.  Conventions:

    T^a = lambda^a / 2                 (color generators)
    [T^a, T^b] = i f^abc T^c            (Lie algebra)
    Tr(T^a T^b) = (1/2) delta^ab        (orthonormal)
    Tr(T^a) = 0
    {T^a, T^b} = (1/N_c) delta^ab + d^abc T^c   (anticommutator)

Casimir invariants:
    T^a T^a = C_F * I_3      with C_F = 4/3   (fundamental)
    f^acd f^bcd = C_A * delta^ab  with C_A = N_c = 3   (adjoint)

Substrate context:  the substrate's algebraic content already contains
SU(3) as a subalgebra of so(7) via SU(3) ⊂ G_2 ⊂ Spin(7).  This module
provides the explicit 3x3 fundamental representation needed to build
QCD amplitudes (color factors at vertices, Casimir contractions).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Gell-Mann matrices (standard convention)
# ---------------------------------------------------------------------------

def gell_mann_matrices() -> list[np.ndarray]:
    """
    Return the 8 hermitian Gell-Mann matrices lambda^a, a = 1..8.

    Index 0 corresponds to lambda^1 (we use 0-indexed Python).

    Substrate identity (asserted at runtime via len + shape check):
      len(matrices) == N_GLUONS == DIM_OCTONION == 8
      each matrix has shape (N_C_SU3, N_C_SU3) == (3, 3)
    """
    from ..isa.constants import N_GLUONS, N_C_SU3
    matrices = _gell_mann_matrices_raw()
    if len(matrices) != N_GLUONS:
        raise AssertionError(
            f"Substrate identity violated: number of Gell-Mann matrices "
            f"{len(matrices)} != N_GLUONS={N_GLUONS}"
        )
    for k, m in enumerate(matrices):
        if m.shape != (N_C_SU3, N_C_SU3):
            raise AssertionError(
                f"Substrate identity violated: λ^{k+1} has shape {m.shape}, "
                f"expected ({N_C_SU3}, {N_C_SU3})"
            )
    return matrices


def _gell_mann_matrices_raw() -> list[np.ndarray]:
    """The raw 8 Gell-Mann matrices without substrate validation."""
    L1 = np.array([[0, 1, 0],
                   [1, 0, 0],
                   [0, 0, 0]], dtype=complex)
    L2 = np.array([[0, -1j, 0],
                   [1j,  0, 0],
                   [0,   0, 0]], dtype=complex)
    L3 = np.array([[1, 0, 0],
                   [0, -1, 0],
                   [0, 0, 0]], dtype=complex)
    L4 = np.array([[0, 0, 1],
                   [0, 0, 0],
                   [1, 0, 0]], dtype=complex)
    L5 = np.array([[0, 0, -1j],
                   [0, 0,  0],
                   [1j, 0, 0]], dtype=complex)
    L6 = np.array([[0, 0, 0],
                   [0, 0, 1],
                   [0, 1, 0]], dtype=complex)
    L7 = np.array([[0, 0,  0],
                   [0, 0, -1j],
                   [0, 1j, 0]], dtype=complex)
    L8 = np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, -2]], dtype=complex) / np.sqrt(3.0)
    return [L1, L2, L3, L4, L5, L6, L7, L8]


def su3_generators() -> list[np.ndarray]:
    """T^a = lambda^a / 2 (Hermitian, traceless, dim-3 fundamental rep)."""
    return [L / 2.0 for L in gell_mann_matrices()]


# ---------------------------------------------------------------------------
# Structure constants  [T^a, T^b] = i f^abc T^c  and  d^abc symmetric tensor
# ---------------------------------------------------------------------------

def structure_constants() -> np.ndarray:
    """
    Compute f^abc by f^abc = -2 i Tr([T^a, T^b] T^c).

    Returns a (8, 8, 8) real array, totally antisymmetric.
    """
    T = su3_generators()
    f = np.zeros((8, 8, 8), dtype=float)
    for a in range(8):
        for b in range(8):
            comm = T[a] @ T[b] - T[b] @ T[a]
            for c in range(8):
                # [T^a, T^b] = i f^abc T^c  =>  f^abc = -2 i Tr([T^a, T^b] T^c)
                val = -2j * np.trace(comm @ T[c])
                f[a, b, c] = float(val.real)
    return f


def d_constants() -> np.ndarray:
    """
    Compute d^abc by d^abc = 2 Tr({T^a, T^b} T^c).

    Returns a (8, 8, 8) real array, totally symmetric.
    """
    T = su3_generators()
    d = np.zeros((8, 8, 8), dtype=float)
    for a in range(8):
        for b in range(8):
            anti = T[a] @ T[b] + T[b] @ T[a]
            for c in range(8):
                # {T^a, T^b} = (1/3) delta^ab I + d^abc T^c
                # Tr({T^a, T^b} T^c) = (1/3) delta^ab Tr(T^c) + d^abc Tr(T^c T^c-prime)
                #                    = 0 + d^abc * (1/2)
                # so d^abc = 2 Tr({T^a, T^b} T^c)
                val = 2.0 * np.trace(anti @ T[c])
                d[a, b, c] = float(val.real)
    return d


# ---------------------------------------------------------------------------
# Casimirs and standard color factors
# ---------------------------------------------------------------------------

C_F = 4.0 / 3.0     # T^a T^a = C_F I_3 in fundamental
C_A = 3.0           # f^acd f^bcd = C_A delta^ab in adjoint (= N_c)
N_C = 3             # Number of colors


def fundamental_casimir() -> np.ndarray:
    """
    Compute T^a T^a in fundamental rep; should equal C_F * I_3 = (4/3) * I_3.
    """
    T = su3_generators()
    casimir = np.zeros((3, 3), dtype=complex)
    for a in range(8):
        casimir = casimir + T[a] @ T[a]
    return casimir


def adjoint_casimir() -> np.ndarray:
    """
    f^acd f^bcd should equal C_A * delta^ab = 3 * delta^ab.

    Returns the 8x8 matrix.
    """
    f = structure_constants()
    return np.einsum("acd,bcd->ab", f, f)
