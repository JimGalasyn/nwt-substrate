"""Tests for the explicit AGL(1, 7) → Z_3 ⟷ G_2 → SU(3) bridge.

Foundation-strengthening verification: the AGL Z_3 residual and the SU(3)
center are **different** Z_3 subgroups of SU(3) ⊂ G_2. Both arise from
substrate breaking patterns but they have distinct group-theoretic roles:

  AGL Z_3   = cyclic permutation of the three complex lines of C^3
              (= the three SM generations as substrate objects)
  SU(3) ctr = scalar phase rotation on C^3
              (= chromatic-class action, not generation-cycling)
"""
from __future__ import annotations

import numpy as np
import pytest

from nwt_substrate.algebra import (
    AGL_Z3_SIGN_LIFT,
    PALEY_TO_BAEZ_LABELING,
    agl_z3_g2_matrix,
    baez_permutation_of_agl_z3,
    complex_structure_J,
    holomorphic_basis_C3,
    signed_permutation_matrix,
    verify_bridge,
)
from nwt_substrate.algebra.g2_bridge import find_paley_to_baez_labeling
from nwt_substrate.algebra.octonions import FANO_LINES, make_octonion_table


# ---------------------------------------------------------------------------
# Paley → Baez labeling σ
# ---------------------------------------------------------------------------

class TestPaleyBaezLabeling:
    def test_constant_matches_brute_force(self):
        """The hardcoded labeling matches a fresh brute-force search."""
        assert find_paley_to_baez_labeling() == PALEY_TO_BAEZ_LABELING

    def test_labeling_is_bijection(self):
        values = list(PALEY_TO_BAEZ_LABELING.values())
        assert sorted(values) == [1, 2, 3, 4, 5, 6, 7]

    def test_labeling_maps_paley_lines_to_baez_lines(self):
        sigma = PALEY_TO_BAEZ_LABELING
        paley_qr = {1, 2, 4}
        paley_blocks = {
            frozenset((r + i) % 7 for r in paley_qr) for i in range(7)
        }
        baez_blocks = {frozenset(line) for line in FANO_LINES}
        mapped = {frozenset(sigma[v] for v in b) for b in paley_blocks}
        assert mapped == baez_blocks


# ---------------------------------------------------------------------------
# Bare permutation π and signed-permutation lift
# ---------------------------------------------------------------------------

class TestBaezPermutation:
    def test_pi_fixes_e1(self):
        pi = baez_permutation_of_agl_z3()
        assert pi[1] == 1

    def test_pi_cycle_structure(self):
        """π should have cycles (e_1)(e_2 e_4 e_6)(e_3 e_5 e_7)."""
        pi = baez_permutation_of_agl_z3()
        # Cycle through e_2 → e_4 → e_6 → e_2
        assert pi[2] == 4 and pi[4] == 6 and pi[6] == 2
        # Cycle through e_3 → e_5 → e_7 → e_3
        assert pi[3] == 5 and pi[5] == 7 and pi[7] == 3

    def test_pi_has_order_3(self):
        pi = baez_permutation_of_agl_z3()
        # Apply π three times → identity
        for i in range(1, 8):
            j = pi[pi[pi[i]]]
            assert j == i


class TestSignedPermutationLift:
    def test_sign_lift_constants(self):
        """The hardcoded sign vector (+,+,+,+,-,+,-) is one of two
        admissible lifts (differ by a Z_2 sign rescaling)."""
        assert AGL_Z3_SIGN_LIFT == {
            1: +1, 2: +1, 3: +1, 4: +1, 5: -1, 6: +1, 7: -1,
        }

    def test_M_orthogonal(self):
        M = agl_z3_g2_matrix()
        assert np.allclose(M.T @ M, np.eye(7))

    def test_M_order_3(self):
        M = agl_z3_g2_matrix()
        assert np.allclose(M @ M @ M, np.eye(7))

    def test_M_fixes_e1(self):
        M = agl_z3_g2_matrix()
        e1 = np.zeros(7); e1[0] = 1.0
        assert np.allclose(M @ e1, e1)

    def test_M_preserves_octonion_product(self):
        """The headline G_2-membership test: T transforms as T."""
        M = agl_z3_g2_matrix()
        T = make_octonion_table()
        M8 = np.zeros((8, 8))
        M8[0, 0] = 1.0
        M8[1:, 1:] = M
        T_transformed = np.einsum("ai,bj,ck,ijk->abc", M8, M8, M8, T)
        assert np.allclose(T_transformed, T)


# ---------------------------------------------------------------------------
# Complex structure J and C^3
# ---------------------------------------------------------------------------

class TestComplexStructureJ:
    def test_J_squares_to_minus_identity_perp_e1(self):
        J = complex_structure_J()
        J2_perp = (J @ J)[1:, 1:]
        assert np.allclose(J2_perp, -np.eye(6))

    def test_J_kills_e1(self):
        """e_1 · e_1 = -e_0 is real, so J(e_1) = 0 in the imaginary 7-space."""
        J = complex_structure_J()
        e1 = np.zeros(7); e1[0] = 1.0
        assert np.allclose(J @ e1, 0)

    def test_J_antisymmetric_perp_e1(self):
        J = complex_structure_J()
        J6 = J[1:, 1:]
        assert np.allclose(J6.T, -J6)


class TestHolomorphicBasis:
    def test_three_holomorphic_directions(self):
        W = holomorphic_basis_C3()
        assert W.shape == (6, 3)

    def test_W_is_orthonormal(self):
        W = holomorphic_basis_C3()
        assert np.allclose(W.conj().T @ W, np.eye(3))

    def test_W_columns_are_plus_i_eigenvectors_of_J(self):
        W = holomorphic_basis_C3()
        J6 = complex_structure_J()[1:, 1:]
        for k in range(3):
            v = W[:, k]
            assert np.allclose(J6 @ v, 1j * v)


# ---------------------------------------------------------------------------
# The headline bridge verification
# ---------------------------------------------------------------------------

class TestBridgeVerification:
    def setup_method(self):
        self.result = verify_bridge()

    def test_M_is_in_G2(self):
        assert self.result.preserves_octonion_product

    def test_M_orthogonal(self):
        assert self.result.is_orthogonal

    def test_M_has_order_3(self):
        assert self.result.order_3

    def test_M_fixes_spin_axis(self):
        assert self.result.fixes_e1

    def test_M_is_C_linear(self):
        """[M, J] = 0 → M ∈ U(3)."""
        assert self.result.commutes_with_J

    def test_M_trace_on_R7_is_plus_one(self):
        """Trace of AGL Z_3 generator on R^7 = +1."""
        assert abs(self.result.trace_R7 - 1.0) < 1e-10

    def test_M_trace_on_C3_is_zero(self):
        """Trace of cyclic-permutation matrix on C^3 = 0."""
        assert abs(self.result.trace_C3) < 1e-10

    def test_M_determinant_on_C3_is_plus_one(self):
        """det(M_C) = +1 ⇒ M ∈ SU(3) (not just U(3))."""
        assert abs(self.result.det_C3 - 1.0) < 1e-10

    def test_M_is_cyclic_permutation_in_holomorphic_basis(self):
        """The structural payoff: M in C^3 is the cyclic permutation matrix
        of three complex lines (up to a unitary basis choice)."""
        # The eigenvalues of M_C must be {1, ω, ω²} (those of a 3-cycle)
        eigvals = np.linalg.eigvals(self.result.agl_z3_matrix_C3)
        omega = np.exp(2j * np.pi / 3)
        expected = {1.0, omega, omega.conjugate()}
        # Match up to ordering and complex conjugation pairing
        eigs_sorted = sorted(eigvals, key=lambda z: np.angle(z))
        # Should contain one element near each cube root of unity
        cube_roots = sorted([1.0 + 0j, omega, omega.conjugate()], key=lambda z: np.angle(z))
        for e, c in zip(eigs_sorted, cube_roots):
            assert abs(e - c) < 1e-9, f"eigenvalue {e} vs cube root {c}"

    def test_su3_center_trace_on_R7_is_minus_two(self):
        assert abs(self.result.su3_center_trace_R7 - (-2.0)) < 1e-10

    def test_agl_z3_is_NOT_su3_center(self):
        """Foundation-strengthening: the AGL Z_3 residual and the SU(3)
        center are different Z_3 subgroups in SU(3) ⊂ G_2."""
        assert self.result.agl_z3_equals_su3_center is False

    def test_traces_are_distinct(self):
        """Trace = +1 (AGL Z_3 on R^7) ≠ -2 (SU(3) center on R^7)."""
        assert abs(self.result.trace_R7 - self.result.su3_center_trace_R7) > 1.0


# ---------------------------------------------------------------------------
# Cross-checks against NWT substrate constants
# ---------------------------------------------------------------------------

class TestNWTSubstrateConsistency:
    def test_three_complex_lines_match_three_generations(self):
        """The C^3 ⊥ spin axis has exactly 3 complex lines = the structural
        count of SM fermion generations and rank(so(7))."""
        from nwt_substrate.isa.constants import RANK_SO7
        W = holomorphic_basis_C3()
        n_complex_lines = W.shape[1]
        assert n_complex_lines == 3 == RANK_SO7

    def test_R6_real_dim_equals_2_times_three_generations(self):
        """The 6 imaginary octonions ⊥ spin axis = 2 × 3 (real-dim per
        complex line × number of generations)."""
        J = complex_structure_J()
        # 6-dim subspace ⊥ e_1 (where J^2 = -I)
        J6 = J[1:, 1:]
        assert J6.shape == (6, 6)
