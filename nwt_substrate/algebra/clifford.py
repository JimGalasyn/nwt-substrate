"""
Clifford algebra Cl(0,7) via octonion left-multiplication matrices.

The 8x8 matrices L_i (i = 0..7) representing left-multiplication by e_i
satisfy the Cl(0,7) anticommutation relation:

    {L_i, L_j} = -2 delta_ij I_8    for i, j in 1..7

despite octonions being non-associative, because the associator
[e_i, e_j, x] is antisymmetric in (i, j) and cancels in the
anticommutator.

Substrate Pauli theorem (Galasyn 2026-05-23):
    L^2(8_spinor) = 7_vector + 21_adjoint = 28 = |E(K_8)|

The antisymmetric two-fermion subspace of the Spin(7) spinor square is
the substrate K_8 edge space.  This is the substrate-topological origin
of Pauli exclusion -- not an external postulate but a consequence of
the bilinear-symmetry table of the Cl(0,7) octonion-Clifford algebra.

See `pauli_decomposition_dims()` and `kform_symmetry()` below.
"""

from __future__ import annotations

from itertools import combinations

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


# ============================================================
# Substrate Pauli theorem -- bilinear-symmetry table
# ============================================================
#
# For Cl(0,7), each L_i (i = 1..7) is antisymmetric as a real 8x8 matrix:
#   (L_i)^T = -L_i
# This is automatic from the octonion construction: the standard inner
# product on O satisfies <e_i x, y> = -<x, e_i y>.
#
# Consequently the Majorana C-matrix is simply C = I_8.  With C = I_8,
# the bilinear bar(psi) gamma^{mu_1...mu_k} chi = psi^T M chi is
# exchange-SYM iff the matrix M = gamma^{mu_1...mu_k} is matrix-SYM.
#
# Computing the matrix symmetry of all rank-k Clifford products yields
# the standard Cl(0,7) bilinear-symmetry table:
#
#     k=0 (scalar, 1):    SYM
#     k=1 (vector, 7):    ANTISYM
#     k=2 (2-form, 21):   ANTISYM
#     k=3 (3-form, 35):   SYM
#     k=4 (4-form, 35):   SYM       (Hodge dual of k=3)
#     k=5 (5-form, 21):   ANTISYM   (Hodge dual of k=2)
#     k=6 (6-form, 7):    ANTISYM   (Hodge dual of k=1)
#     k=7 (volume, 1):    SYM       (Hodge dual of k=0)
#
# Independent ranks contributing to 8 (x) 8 = 64 are k = 0, 1, 2, 3
# (k = 4..7 are Hodge duals).  Summing dimensions:
#
#     L^2(8) = 7_vector + 21_adjoint = 28 = |E(K_8)|
#     Sym^2(8) = 1_scalar + 35_4-form = 36
#
# The 28-dim antisymmetric two-fermion subspace coincides with the
# substrate K_8 edge count.  This is the substrate Pauli theorem.


def gamma_is_antisymmetric(i: int, T: np.ndarray, tol: float = 1e-12) -> bool:
    """Return True iff L_i is antisymmetric, i.e. (L_i)^T = -L_i.

    For i = 1..7 this is always True (octonion-construction property).
    For i = 0 this returns False (L_0 = I_8 is symmetric).
    """
    L_i = left_mult_matrix(i, T)
    return np.allclose(L_i.T, -L_i, atol=tol)


def kform_matrix(indices: tuple[int, ...], T: np.ndarray) -> np.ndarray:
    """Ordered product L_{i_1} L_{i_2} ... L_{i_k} of distinct gamma matrices.

    For distinct indices in {1, ..., 7}, this is the rank-k Clifford-product
    matrix gamma^{i_1 ... i_k} up to an index-ordering sign (which does not
    affect the matrix symmetry property).  For k = 0 (empty tuple), returns
    the 8x8 identity.
    """
    M = np.eye(8)
    for i in indices:
        M = M @ left_mult_matrix(i, T)
    return M


def kform_symmetry(k: int, T: np.ndarray, tol: float = 1e-12) -> str:
    """Classify the matrix symmetry of all rank-k Clifford products.

    By Spin(7)-equivariance, every rank-k k-form has the same matrix
    symmetry type.  This function verifies that empirically by enumerating
    all C(7, k) k-tuples and checking each ordered product.

    Returns:
        'SYM'     -- gamma^{mu_1...mu_k} is symmetric for all k-tuples
        'ANTISYM' -- gamma^{mu_1...mu_k} is antisymmetric for all k-tuples
        'MIXED (s/a/o)' -- otherwise (shouldn't happen for Cl(0,7))

    The bilinear symmetry of bar(psi) gamma^{mu_1...mu_k} chi under particle
    exchange (with Majorana C = I_8) matches the matrix symmetry.
    """
    if k == 0:
        # The empty product is the identity matrix, which is symmetric.
        return 'SYM'
    sym_count = 0
    antisym_count = 0
    total = 0
    for tup in combinations(range(1, 8), k):
        M = kform_matrix(tup, T)
        if np.allclose(M, M.T, atol=tol):
            sym_count += 1
        elif np.allclose(M, -M.T, atol=tol):
            antisym_count += 1
        total += 1
    if antisym_count == total:
        return 'ANTISYM'
    if sym_count == total:
        return 'SYM'
    mixed = total - sym_count - antisym_count
    return f'MIXED ({sym_count}/{antisym_count}/{mixed})'


# Independent rank contributions to 8 (x) 8 = 1 + 7 + 21 + 35 = 64.
# k = 4..7 are Hodge duals of k = 3..0 in 7D, so don't contribute
# independent dimensions.
_INDEPENDENT_RANKS: dict[int, tuple[str, int]] = {
    0: ('scalar (1)', 1),
    1: ('vector (7 = |V(K_7)|)', 7),
    2: ('adjoint (21 = |E(K_7)| = dim so(7))', 21),
    3: ('4-form (35 = dim self-dual 4-form)', 35),
}


def pauli_decomposition_dims(T: np.ndarray) -> dict:
    """Compute L^2(8_spinor) and Sym^2(8_spinor) dimensions and decompositions.

    For each independent rank k = 0..3 of the Spin(7) tensor decomposition
    8 (x) 8 = 1 + 7 + 21 + 35, classifies the rank's bilinear symmetry via
    `kform_symmetry()` and sums dimensions accordingly.

    Returns a dict with keys:
        'antisym':         total antisymmetric dimension (expected 28 = |E(K_8)|)
        'sym':             total symmetric dimension (expected 36)
        'antisym_decomp':  {rank: name} for antisymmetric ranks
        'sym_decomp':      {rank: name} for symmetric ranks
        'k8_edge_match':   True iff antisym == |E(K_8)| = 28 (substrate Pauli theorem)
    """
    antisym = 0
    sym = 0
    antisym_decomp: dict[int, str] = {}
    sym_decomp: dict[int, str] = {}
    for k, (name, dim) in _INDEPENDENT_RANKS.items():
        sym_type = kform_symmetry(k, T)
        if sym_type == 'ANTISYM':
            antisym += dim
            antisym_decomp[k] = name
        elif sym_type == 'SYM':
            sym += dim
            sym_decomp[k] = name
        else:  # pragma: no cover -- shouldn't happen for Cl(0,7)
            raise RuntimeError(
                f"Unexpected MIXED symmetry at rank {k}: {sym_type}"
            )
    from ..isa.constants import N_EDGES_K8  # |E(K_8)| = 28 (lazy import; algebra idiom)
    K8_EDGES = N_EDGES_K8
    return {
        'antisym': antisym,
        'sym': sym,
        'antisym_decomp': antisym_decomp,
        'sym_decomp': sym_decomp,
        'k8_edge_match': antisym == K8_EDGES,
    }


def verify_substrate_pauli(T: np.ndarray | None = None) -> dict:
    """Verify the substrate Pauli theorem by direct Cl(0,7) computation.

    Substrate Pauli theorem: L^2(8_spinor) = 7_vector + 21_adjoint = 28 =
    |E(K_8)|, identifying the antisymmetric two-fermion subspace with the
    substrate K_8 edge space.

    This is the standalone executable verification that supports the
    substrate-canonical claim "fermion antisymmetry is K_8 edges" from
    Galasyn 2026-05-23 P9.  Runs the full bilinear-symmetry table check.

    Args:
        T: octonion multiplication table (default: built fresh)

    Returns:
        Same shape as pauli_decomposition_dims, plus:
        'all_gammas_antisymmetric': True iff (L_i)^T = -L_i for i=1..7
    """
    from .octonions import make_octonion_table
    if T is None:
        T = make_octonion_table()
    result = pauli_decomposition_dims(T)
    result['all_gammas_antisymmetric'] = all(
        gamma_is_antisymmetric(i, T) for i in range(1, 8)
    )
    return result
