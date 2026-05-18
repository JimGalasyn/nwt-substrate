"""
Explicit AGL(1, 7) → Z_3 ⟷ G_2 → SU(3) bridge on octonion imaginaries.

This module makes the discrete-substrate-vertex and continuous-Lie-algebra
pictures of the cosmogenic three-generation breaking demonstrably the same
breaking — and corrects the previously naive "the two Z_3 subgroups are
identical in G_2" claim.

Structural finding (verified numerically below)
------------------------------------------------
1. There is a Paley-to-Baez relabeling σ: F_7 → {1, ..., 7} that maps the
   Paley (7, 3, 1) QR-translate blocks bijectively onto Baez's FANO_LINES
   convention.

2. Under σ, the AGL(1, 7) Z_3 generator x ↦ 2x mod 7 acts on the imaginary
   octonion basis as the permutation
       (e_1)(e_2 e_4 e_6)(e_3 e_5 e_7).

3. This bare permutation does **not** preserve octonion multiplication; it
   must be lifted to a *signed* permutation in G_2 ⊂ SO(7). A brute-force
   search (128 sign assignments) returns the unique-up-to-trivial-rescaling
   sign vector ε = (+1, +1, +1, +1, −1, +1, −1). The resulting 7 × 7 matrix
   M satisfies:
       M^T M = I,  M^3 = I,  T(M·, M·, M·) = T  ⇒ M ∈ G_2,
       M · e_1 = e_1                          ⇒ M ∈ SU(3) ⊂ G_2.

4. M commutes with J := mult-by-e_1 (the complex structure on the 6-dim
   imaginary subspace perpendicular to e_1) ⇒ M is C-linear ⇒ M ∈ U(3).

5. In the holomorphic basis (the +i eigenvectors of J) M is a 3 × 3 cyclic
   permutation matrix
            [[0, 0, 1],
             [1, 0, 0],
             [0, 1, 0]]
   with trace 0 on C^3 and det +1. So M ∈ SU(3).

6. **Therefore: the residual Z_3 of AGL(1, 7) → Z_3 (under cosmogenic
   spin-axis + polarity breaking) embeds in G_2 as the cyclic permutation
   of the three complex lines of C^3 = (imaginary octonions ⊥ spin axis,
   with the spin axis itself supplying the complex structure J).**

7. The SU(3) center Z_3 = ⟨ω · I_3⟩ is a **different** Z_3 ⊂ SU(3). It
   has trace 3ω on C^3 and trace −2 on R^7, vs. M's trace 0 and +1.

   Both are order-3 subgroups of SU(3), both arise from cosmogenic
   structure, but they are not the same Z_3 in G_2. The AGL residual is
   the "rotate-the-three-generations" Z_3; the center is the "phase-the-
   whole-C^3" Z_3. The former carries the cycle structure of three
   generations; the latter does not.

This module exposes the bridge as numerical primitives so the claim is
testable and reproducible — and ensures we do not over-state the
relationship between the discrete and continuous breakings.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .octonions import make_octonion_table


# ---------------------------------------------------------------------------
# Paley → Baez Fano-isomorphism σ (the labeling that aligns the two designs)
# ---------------------------------------------------------------------------

# σ: F_7 = {0,...,6} → {1,...,7} mapping Paley vertices to Baez octonion
# basis indices, chosen so that each Paley QR-translate block becomes a
# Baez FANO_LINE (as a set). This σ was found by brute-force search over
# bijections; see ``find_paley_to_baez_labeling`` below.
PALEY_TO_BAEZ_LABELING: dict[int, int] = {0: 1, 1: 2, 2: 4, 3: 3, 4: 6, 5: 7, 6: 5}


def find_paley_to_baez_labeling() -> dict[int, int]:
    """Brute-force the bijection σ: F_7 → {1,...,7} aligning the two designs.

    Returns the same labeling as ``PALEY_TO_BAEZ_LABELING``; provided so
    the value can be re-derived from primitives.
    """
    from itertools import permutations
    from .octonions import FANO_LINES

    paley_qr = frozenset({1, 2, 4})
    paley_blocks = {
        frozenset((r + i) % 7 for r in paley_qr) for i in range(7)
    }
    baez_blocks = {frozenset(line) for line in FANO_LINES}

    for perm in permutations(range(1, 8)):
        sigma = {i: perm[i] for i in range(7)}
        mapped = {frozenset(sigma[v] for v in b) for b in paley_blocks}
        if mapped == baez_blocks:
            return sigma
    raise RuntimeError("no Paley-to-Baez labeling found (bug)")


# ---------------------------------------------------------------------------
# Lift the AGL(1, 7) Z_3 generator to G_2 ⊂ SO(7)
# ---------------------------------------------------------------------------

# The sign function ε required to lift the bare permutation
#   π = (e_1)(e_2 e_4 e_6)(e_3 e_5 e_7)
# from S_7 ⊂ O(7) to G_2 ⊂ SO(7) (preserving octonion multiplication).
# Derived by solving the seven Fano-cyclic-orientation constraints; see
# tests/test_g2_bridge for the constraint-by-constraint verification.
AGL_Z3_SIGN_LIFT: dict[int, int] = {1: +1, 2: +1, 3: +1, 4: +1, 5: -1, 6: +1, 7: -1}


def baez_permutation_of_agl_z3() -> dict[int, int]:
    """Return the bare permutation π of Baez basis indices induced by the
    AGL(1, 7) Z_3 generator x ↦ 2x mod 7, via the σ relabeling.

    Cycle structure: (e_1)(e_2 e_4 e_6)(e_3 e_5 e_7).
    """
    sigma = PALEY_TO_BAEZ_LABELING
    return {sigma[i]: sigma[(2 * i) % 7] for i in range(7)}


def signed_permutation_matrix(
    pi: dict[int, int], signs: dict[int, int]
) -> np.ndarray:
    """Build the 7 × 7 signed-permutation matrix M on the imaginary octonion
    7-space with M·e_i = signs[i] · e_{pi[i]} (1-indexed basis labels).
    """
    M = np.zeros((7, 7))
    for i in range(1, 8):
        M[pi[i] - 1, i - 1] = signs[i]
    return M


def agl_z3_g2_matrix() -> np.ndarray:
    """The AGL(1, 7) Z_3 generator x ↦ 2x lifted to G_2 ⊂ SO(7).

    Acts on the 7-dim imaginary-octonion representation as the signed
    permutation (π, ε) with π = (e_1)(e_2 e_4 e_6)(e_3 e_5 e_7) and
    ε = AGL_Z3_SIGN_LIFT.

    Properties:
        M^T M = I,  M^3 = I,  trace(M) = +1
        M fixes e_1 (the cosmogenic spin axis)
        M preserves the octonion multiplication table T
        M lies in SU(3) ⊂ G_2 (the stabilizer of e_1)
    """
    return signed_permutation_matrix(
        baez_permutation_of_agl_z3(), AGL_Z3_SIGN_LIFT
    )


# ---------------------------------------------------------------------------
# Complex structure J = mult-by-e_1 and the C^3 holomorphic basis
# ---------------------------------------------------------------------------

def complex_structure_J() -> np.ndarray:
    """The complex structure J on the imaginary octonion 7-space given by
    left-multiplication by e_1.

    On the 6-dim subspace ⊥ e_1, J satisfies J^2 = -I_6; J · e_1 = 0
    (since e_1 · e_1 = -e_0 is real).
    """
    T = make_octonion_table()
    J = np.zeros((7, 7))
    for j in range(1, 8):
        for k in range(1, 8):
            J[k - 1, j - 1] = T[1, j, k]
    return J


def holomorphic_basis_C3() -> np.ndarray:
    """Return a 6 × 3 complex matrix whose columns are an orthonormal
    holomorphic basis (i.e., +i eigenvectors of J) for C^3 = imaginary
    octonions ⊥ e_1.

    The columns span the 3-dim C-subspace on which the SU(3) stabilizer of
    e_1 acts in its fundamental representation.
    """
    J6 = complex_structure_J()[1:, 1:]  # 6 × 6 part
    eigvals, eigvecs = np.linalg.eig(J6)
    plus_i_idx = np.where(np.isclose(eigvals, 1j))[0]
    if len(plus_i_idx) != 3:
        raise RuntimeError(
            f"expected 3 +i eigenvectors of J; found {len(plus_i_idx)}"
        )
    W = eigvecs[:, plus_i_idx]
    # Orthonormalize
    Q, _ = np.linalg.qr(W)
    return Q


# ---------------------------------------------------------------------------
# The bridge verification
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BridgeVerification:
    """Result of the explicit AGL(1, 7) ⟷ G_2 ⟷ SU(3) bridge construction."""

    agl_z3_matrix_R7: np.ndarray
    """The 7 × 7 G_2 element lifting AGL(1, 7) Z_3."""

    trace_R7: float
    """Trace of agl_z3_matrix_R7 on R^7. Expect +1."""

    fixes_e1: bool
    """Whether agl_z3_matrix_R7 · e_1 = e_1."""

    is_orthogonal: bool
    """Whether M^T M = I_7."""

    order_3: bool
    """Whether M^3 = I_7."""

    preserves_octonion_product: bool
    """Whether (extended to 8x8 by e_0 ↦ e_0) M preserves the
    octonion multiplication table T (G_2 membership)."""

    commutes_with_J: bool
    """Whether [M, J] = 0 (C-linearity ⇒ M ∈ U(3))."""

    agl_z3_matrix_C3: np.ndarray
    """The 3 × 3 complex matrix representing M in the C^3 holomorphic basis.
    Should equal the cyclic permutation matrix on three complex lines."""

    trace_C3: complex
    """Trace of the C^3 representation. Expect 0 (cyclic permutation has
    trace 0 on the standard 3-dim rep)."""

    det_C3: complex
    """Determinant of M_C^3. Expect +1 (so M ∈ SU(3))."""

    su3_center_trace_R7: float
    """For comparison: trace of SU(3) center generator ω·I_3 on R^7. = -2."""

    agl_z3_equals_su3_center: bool
    """Headline: is the AGL Z_3 element the SAME Z_3 as the SU(3) center?
    Expect FALSE — these are different Z_3 subgroups of SU(3) ⊂ G_2."""


def verify_bridge(tol: float = 1e-10) -> BridgeVerification:
    """Construct and verify the explicit AGL(1, 7) ⟷ G_2 ⟷ SU(3) bridge.

    Returns a ``BridgeVerification`` whose fields document each invariant
    of the construction. The headline numerical fact is:

      trace(M on R^7) = +1   (the AGL Z_3 generator)
      trace(SU(3) center on R^7) = -2   (different Z_3 ⊂ SU(3))
    """
    T = make_octonion_table()
    M = agl_z3_g2_matrix()
    J = complex_structure_J()

    # Extend M to 8 × 8 with e_0 → e_0
    M8 = np.zeros((8, 8))
    M8[0, 0] = 1.0
    M8[1:, 1:] = M

    T_transformed = np.einsum("ai,bj,ck,ijk->abc", M8, M8, M8, T)
    preserves_T = bool(np.allclose(T_transformed, T, atol=tol))

    e1 = np.zeros(7); e1[0] = 1.0
    fixes_e1 = bool(np.allclose(M @ e1, e1, atol=tol))

    is_ortho = bool(np.allclose(M.T @ M, np.eye(7), atol=tol))
    order3 = bool(np.allclose(M @ M @ M, np.eye(7), atol=tol))
    commutes = bool(np.allclose(M @ J - J @ M, 0, atol=tol))

    # Express M in C^3 holomorphic basis
    M6 = M[1:, 1:]
    W = holomorphic_basis_C3()
    M_C = W.conj().T @ M6 @ W

    # SU(3) center on R^7
    omega = np.exp(2j * np.pi / 3)
    su3_center_trace = 1.0 + 2 * np.real(3 * omega)  # = 1 + 2·(-3/2) = -2

    # The literal "is M the SU(3) center?" check:
    # SU(3) center on C^3 is ω·I_3, trace = 3ω. M_C has trace 0.
    is_center = bool(np.allclose(np.trace(M_C), 3 * omega, atol=tol))

    return BridgeVerification(
        agl_z3_matrix_R7=M,
        trace_R7=float(np.trace(M)),
        fixes_e1=fixes_e1,
        is_orthogonal=is_ortho,
        order_3=order3,
        preserves_octonion_product=preserves_T,
        commutes_with_J=commutes,
        agl_z3_matrix_C3=M_C,
        trace_C3=complex(np.trace(M_C)),
        det_C3=complex(np.linalg.det(M_C)),
        su3_center_trace_R7=float(su3_center_trace),
        agl_z3_equals_su3_center=is_center,
    )


__all__ = [
    "PALEY_TO_BAEZ_LABELING",
    "find_paley_to_baez_labeling",
    "AGL_Z3_SIGN_LIFT",
    "baez_permutation_of_agl_z3",
    "signed_permutation_matrix",
    "agl_z3_g2_matrix",
    "complex_structure_J",
    "holomorphic_basis_C3",
    "BridgeVerification",
    "verify_bridge",
]
