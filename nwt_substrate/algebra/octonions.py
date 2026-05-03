"""
Octonion algebra (Cayley-Graves 1843/45) via Fano-plane multiplication table.

Provides:
  - make_octonion_table()  : 8x8x8 multiplication tensor T[i,j,k]
  - octo_mul(a, b, T)      : a * b (octonion product)
  - octo_assoc(a, b, c, T) : associator [a,b,c] = (ab)c - a(bc)
  - octo_norm_sq(a)        : |a|^2 (Hurwitz norm)
  - basis_vector(i)        : i-th basis octonion e_i (i = 0..7)
  - FANO_LINES             : the 7 Fano-plane lines (Baez 2002 §2.2 convention)

All operations verified by tests in tests/test_octonions.py.
"""

from __future__ import annotations

import numpy as np


# Baez's Fano plane convention (Baez, "The Octonions", 2002 §2.2).
# Each tuple (i, j, k) is a line in cyclic order:
#   e_i * e_j = +e_k,  e_j * e_k = +e_i,  e_k * e_i = +e_j
#   e_j * e_i = -e_k,  e_k * e_j = -e_i,  e_i * e_k = -e_j
FANO_LINES = [
    (1, 2, 3),
    (1, 4, 5),
    (1, 7, 6),
    (2, 4, 6),
    (2, 5, 7),
    (3, 4, 7),
    (3, 6, 5),
]


def make_octonion_table() -> np.ndarray:
    """
    Build the 8x8x8 multiplication table T such that
        e_i * e_j = sum_k T[i, j, k] e_k.

    e_0 is the real unit (identity); e_1..e_7 are the imaginary basis.
    """
    T = np.zeros((8, 8, 8))
    # Identity: e_0 acts as 1
    for i in range(8):
        T[0, i, i] = 1.0
        T[i, 0, i] = 1.0
    # Imaginary units square to -1
    for i in range(1, 8):
        T[i, i, 0] = -1.0
    # Fano-line products
    for (i, j, k) in FANO_LINES:
        T[i, j, k] = 1.0
        T[j, k, i] = 1.0
        T[k, i, j] = 1.0
        T[j, i, k] = -1.0
        T[k, j, i] = -1.0
        T[i, k, j] = -1.0
    return T


def basis_vector(i: int) -> np.ndarray:
    """Return the i-th basis octonion as an 8-vector."""
    e = np.zeros(8)
    e[i] = 1.0
    return e


def octo_mul(a: np.ndarray, b: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Octonion product a * b using table T."""
    return np.einsum("i,j,ijk->k", a, b, T)


def octo_assoc(a: np.ndarray, b: np.ndarray, c: np.ndarray,
               T: np.ndarray) -> np.ndarray:
    """Octonion associator [a, b, c] = (ab)c - a(bc)."""
    return octo_mul(octo_mul(a, b, T), c, T) - octo_mul(a, octo_mul(b, c, T), T)


def octo_norm_sq(a: np.ndarray) -> float:
    """Squared Euclidean norm |a|^2 = sum a_i^2 (= Hurwitz norm)."""
    return float(np.dot(a, a))


def octo_conj(a: np.ndarray) -> np.ndarray:
    """Octonion conjugate: flip sign of imaginary parts e_1..e_7."""
    out = a.copy()
    out[1:] = -out[1:]
    return out


def fano_partner(i: int, j: int) -> tuple[int, int]:
    """
    Return (k, sign) such that e_i * e_j = sign * e_k for distinct i, j > 0.

    Raises ValueError if (i, j) not on any Fano line (shouldn't happen for
    distinct i, j in 1..7).
    """
    if i == j:
        raise ValueError(f"e_{i} * e_{i} = -e_0 (no Fano partner for equal indices)")
    for line in FANO_LINES:
        if i in line and j in line:
            k = next(x for x in line if x != i and x != j)
            a, b, c = line
            if (i, j) in [(a, b), (b, c), (c, a)]:
                return k, +1
            else:
                return k, -1
    raise ValueError(f"No Fano line contains both e_{i} and e_{j}")
