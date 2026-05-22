"""Q15 — Spin(7) acting on Cl(0,7) spinor S=8: derivation of n_q sectors.

Q13/Q14 ruled out so(3)_diag ⊂ so(7) as the substrate locus of the
carrier-knot multiplicity n_q ∈ {1, 2, 3, 5} (= consecutive Fibonacci
F_2..F_5), because so(3)_diag has only ODD-dimensional irreps {1, 3, 5}
in the relevant decomposition of so(7) (the (3,3) bi-fundamental of
so(3)_{123} × so(3)_{456}). The "missing 2" for the meson (Hopf) sector
sits structurally OUTSIDE the so(3)_diag rep theory of so(7).

The remaining theoretical lead, per the Q13/Q14 negative memo, is to
test Spin(7) acting on its 8-dim spinor S = Cl(0,7) qubit space. Under
SU(2) ⊂ Spin(7), the 8-dim spinor can decompose with half-integer-spin
content, giving access to dim-2 (spin-1/2) irreps.

We have existing Cl(0,7) infrastructure:
  - octonions.make_octonion_table() → T[i,j,k] 8x8x8
  - clifford.left_mult_matrix(i, T) → L_i 8x8 satisfying {L_i, L_j} = -2δ_ij I
  - g2_bridge.PALEY_TO_BAEZ_LABELING → vertex k (∈ {0..6}) ↔ octonion idx (∈ {1..7})

This script:

  Step 1: Build Cl(0,7) gammas γ_i = L_i (i=1..7), real 8×8 skew-symmetric.
          Verify Clifford anticommutation to machine precision.

  Step 2: Build 21 Spin(7) generators Σ_ij = (1/2) γ_i γ_j for i<j.
          These are real 8×8 skew-symmetric and generate Spin(7) ⊂ SO(8).
          Map K_7 edges to Σ via the Paley→Baez labeling.

  Step 3: Identify several candidate SU(2) ⊂ Spin(7) embeddings:
          (A) so(3)_diag-on-S: octonion analog of Q14's so(3)_{diag}
              with L^A_x = Σ_{23} + Σ_{56}, etc. — test whether S=8
              decomposes differently than V_7.
          (B) Three-generation Z_3-commuting SU(2): SU(2) inside G_2
              centralizer of AGL Z_3 generator.
          (C) Joint-spin SU(2): J^B_α = (1/2) Σ_i γ_i ... candidate
              total-spin operator.
          (D) Cl(0,7)-spinor v-cycle SU(2): align with Heffter v-direction.

          For each, compute the 8-dim spinor decomposition. Look for
          irrep dimensions matching {1, 2, 3, 5}.

  Step 4: For each compendium walk, build A^S = Σ_e Σ_{e} (algebra
          element on S=8) and v-weighted variant A_v^S. Project onto
          each SU(2) candidate's irrep eigenspaces. Test correlation
          of fractional weights with Paper 11 n_q.

If a clean clustering emerges (e.g., mesons in dim-2 spinor, nucleons
in dim-5 quintet), the n_q derivation gap closes structurally.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q15_spin7_s8_nq_derivation.py
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import numpy as np

from nwt_substrate.algebra.clifford import left_mult_matrix
from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.g2_bridge import PALEY_TO_BAEZ_LABELING
from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.particles.compendium import COMPENDIUM


N_VERT = 7              # K_7 vertices {0..6}
SPINOR_DIM = 8          # Cl(0,7) spinor space S = R^8


# ============================================================================
# Step 1: Cl(0,7) gamma matrices on S=8
# ============================================================================

def build_gammas() -> list[np.ndarray]:
    """Return [γ_1, ..., γ_7] as real 8×8 skew-symmetric matrices.

    Convention: γ_i = L_i (octonion left-multiplication).
    Satisfies {γ_i, γ_j} = -2 δ_ij I_8 (Cl(0,7) signature).
    """
    T = make_octonion_table()
    return [left_mult_matrix(i, T) for i in range(1, 8)]


def verify_clifford(gammas: list[np.ndarray], tol: float = 1e-10) -> None:
    """Verify {γ_i, γ_j} = -2 δ_ij I_8 for all i, j in 1..7."""
    I = np.eye(SPINOR_DIM)
    for i in range(7):
        for j in range(7):
            anti = gammas[i] @ gammas[j] + gammas[j] @ gammas[i]
            expected = -2.0 * (1.0 if i == j else 0.0) * I
            if not np.allclose(anti, expected, atol=tol):
                raise AssertionError(
                    f"Clifford failure at (γ_{i+1}, γ_{j+1}): "
                    f"max|diff| = {np.max(np.abs(anti - expected)):.2e}"
                )
    # γ_i skew-symmetric
    for i, g in enumerate(gammas):
        if not np.allclose(g.T, -g, atol=tol):
            raise AssertionError(f"γ_{i+1} not skew-symmetric")
    print("  ✓ {γ_i, γ_j} = -2δ_ij I_8 verified for all 49 pairs")
    print("  ✓ All γ_i real skew-symmetric (anti-Hermitian on R^8)")


# ============================================================================
# Step 2: Spin(7) generators Σ_ij = (1/2) γ_i γ_j on S=8
# ============================================================================

def build_sigma_generators(gammas: list[np.ndarray]) -> dict[tuple[int, int], np.ndarray]:
    """Return {(i, j): Σ_ij} for 1 ≤ i < j ≤ 7.

    Σ_ij = (1/4)[γ_i, γ_j] = (1/2) γ_i γ_j  (since γ_i γ_j = -γ_j γ_i for i≠j).
    These are real skew-symmetric 8×8 matrices generating Spin(7) ⊂ SO(8).
    """
    sigmas = {}
    for i in range(7):
        for j in range(i + 1, 7):
            sigmas[(i + 1, j + 1)] = 0.5 * (gammas[i] @ gammas[j])
    return sigmas


def verify_so7_commutators(sigmas: dict, tol: float = 1e-10) -> None:
    """Spot-check so(7) commutation relations:
       [Σ_ij, Σ_jk] = Σ_ik for distinct i, j, k.

    Specifically: [γ_i γ_j / 2, γ_j γ_k / 2] = γ_i γ_k / 2 · (γ_j² / -1)
        γ_j² = -I, so [Σ_ij, Σ_jk] = (1/4)(γ_i γ_j γ_j γ_k - γ_j γ_k γ_i γ_j)
                                   = (1/4)(- γ_i γ_k - γ_j γ_k γ_i γ_j)
        γ_j γ_k γ_i γ_j: anticommute γ_j past γ_i (i ≠ j): γ_j γ_k γ_i γ_j
            = - γ_j γ_k γ_j γ_i = + γ_j γ_j γ_k γ_i = - γ_k γ_i = γ_i γ_k
        So [Σ_ij, Σ_jk] = (1/4)(- γ_i γ_k - γ_i γ_k) = -(1/2) γ_i γ_k = -Σ_ik
    Convention check.
    """
    g = build_gammas()
    samples = [(1, 2, 3), (3, 4, 5), (2, 5, 7), (1, 6, 7)]
    for (i, j, k) in samples:
        S_ij = sigmas[(min(i, j), max(i, j))]
        if i > j: S_ij = -S_ij
        S_jk = sigmas[(min(j, k), max(j, k))]
        if j > k: S_jk = -S_jk
        comm = S_ij @ S_jk - S_jk @ S_ij
        # Expected -Σ_ik by the derivation above
        S_ik = sigmas[(min(i, k), max(i, k))]
        if i > k: S_ik = -S_ik
        if not np.allclose(comm, -S_ik, atol=tol):
            # Try +Σ_ik
            if not np.allclose(comm, +S_ik, atol=tol):
                raise AssertionError(
                    f"so(7) commutator failure at ({i},{j},{k}): "
                    f"max|[Σ_ij, Σ_jk] ∓ Σ_ik| = "
                    f"{min(np.max(np.abs(comm - S_ik)), np.max(np.abs(comm + S_ik))):.2e}"
                )
    print(f"  ✓ so(7) commutator structure [Σ_ij, Σ_jk] ∝ Σ_ik on {len(samples)} samples")


# ============================================================================
# Step 3: K_7 edge → Σ generator (via Paley→Baez labeling)
# ============================================================================

def k7_edge_to_sigma_label(a: int, b: int) -> tuple[int, int]:
    """Map a directed K_7 edge (a, b) with a, b ∈ {0..6} to a Σ-label
    (i, j) with i, j ∈ {1..7} (octonion-imaginary basis).

    Uses PALEY_TO_BAEZ_LABELING so that the Paley QR-design on K_7
    aligns with the Baez Fano-plane convention used by the octonion
    table.

    Returns (i, j) sorted i < j. Caller must apply a sign for direction.
    """
    i = PALEY_TO_BAEZ_LABELING[a]
    j = PALEY_TO_BAEZ_LABELING[b]
    return (min(i, j), max(i, j))


def walk_algebra_elements_S(walk: list[int],
                              sigmas: dict[tuple[int, int], np.ndarray]
                              ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the spinor-space algebra elements for a closed walk:
      A^S      = Σ_e (sign_e) Σ_e        (unweighted)
      A_u^S    = Σ_e (sign_e Δu_e) Σ_e   (u-direction-weighted)
      A_v^S    = Σ_e (sign_e Δv_e) Σ_e   (v-direction-weighted)
    Output: three 8×8 real matrices (skew-symmetric).
    """
    A_walk = np.zeros((SPINOR_DIM, SPINOR_DIM))
    A_u = np.zeros((SPINOR_DIM, SPINOR_DIM))
    A_v = np.zeros((SPINOR_DIM, SPINOR_DIM))
    for k in range(len(walk) - 1):
        a, b = walk[k], walk[k + 1]
        ij = k7_edge_to_sigma_label(a, b)
        S = sigmas[ij]
        # Sign convention: edge a→b carries +Σ when (Paley→Baez)(a) <
        # (Paley→Baez)(b); -Σ otherwise.
        sign = +1.0 if PALEY_TO_BAEZ_LABELING[a] < PALEY_TO_BAEZ_LABELING[b] else -1.0
        n_u, n_v = edge_winding_class(a, b)
        if a > b:
            n_u, n_v = -n_u, -n_v
        A_walk += sign * S
        A_u += sign * n_u * S
        A_v += sign * n_v * S
    return A_walk, A_u, A_v


# ============================================================================
# Step 4: SU(2) ⊂ Spin(7) candidate embeddings
# ============================================================================

def su2_so3_diag_octonion(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (A): the octonion analog of Q14's so(3)_diag.

    Define generators on S=8:
      J^A_x = Σ_{2,3} + Σ_{5,6}
      J^A_y = Σ_{1,3} + Σ_{4,6}
      J^A_z = Σ_{1,2} + Σ_{4,5}

    These act on the 7-dim imaginary octonion space the same way as Q14's
    so(3)_diag acts on so(7) via the adjoint, BUT here they act on the
    8-dim SPINOR space — and the spinor decomposition can have a different
    Casimir spectrum.
    """
    Jx = sigmas[(2, 3)] + sigmas[(5, 6)]
    Jy = sigmas[(1, 3)] + sigmas[(4, 6)]
    Jz = sigmas[(1, 2)] + sigmas[(4, 5)]
    return Jx, Jy, Jz


def su2_so3_anti_diag_octonion(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (B): anti-diagonal of so(3)_{123} - so(3)_{456}.

      J^B_x = Σ_{2,3} - Σ_{5,6}
      J^B_y = Σ_{1,3} - Σ_{4,6}
      J^B_z = Σ_{1,2} - Σ_{4,5}
    """
    Jx = sigmas[(2, 3)] - sigmas[(5, 6)]
    Jy = sigmas[(1, 3)] - sigmas[(4, 6)]
    Jz = sigmas[(1, 2)] - sigmas[(4, 5)]
    return Jx, Jy, Jz


def su2_so3_123(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (C): so(3) on octonion triple {e_1, e_2, e_3} only.

      J^C_x = Σ_{2,3}, J^C_y = Σ_{1,3}, J^C_z = Σ_{1,2}
    """
    return sigmas[(2, 3)], sigmas[(1, 3)], sigmas[(1, 2)]


def su2_so3_456(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (D): so(3) on octonion triple {e_4, e_5, e_6} only."""
    return sigmas[(5, 6)], sigmas[(4, 6)], sigmas[(4, 5)]


def su2_principal_so7(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (E): 'principal' chain SU(2) using consecutive
    generators along a fixed permutation of the 7 imaginaries.

      J^E_x = Σ_{1,2} + Σ_{3,4} + Σ_{5,6}
      J^E_y = Σ_{2,3} + Σ_{4,5} + Σ_{6,7}
      J^E_z = Σ_{1,3} + Σ_{2,4} + ... (Cartan-like)

    Test whether 8 decomposes non-trivially.
    """
    Jx = sigmas[(1, 2)] + sigmas[(3, 4)] + sigmas[(5, 6)]
    Jy = sigmas[(2, 3)] + sigmas[(4, 5)] + sigmas[(6, 7)]
    Jz = sigmas[(1, 3)] + sigmas[(2, 4)] + sigmas[(3, 5)] + sigmas[(4, 6)] + sigmas[(5, 7)]
    return Jx, Jy, Jz


def su2_diag_L_Sp2(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (G): diagonal SU(2) inside SU(2)_L × SU(2)_R-principal-of-Sp(2).

    Spin(7) ⊃ SU(2)_L × Sp(2) with spinor 8 = (2_L, 4_R). The principal SU(2)
    of Sp(2) makes 4_R = spin-3/2. The diagonal SU(2)_diag = SU(2)_L + SU(2)_R^princ
    decomposes 2 ⊗ 4 = 3 ⊕ 5 — i.e., j=1 ⊕ j=2.

    Empirically construct candidate SU(2)_L acting on a pair of γ's
    (e.g., {γ_1, γ_2}) and SU(2)_R acting on the orthogonal 4-block.
    Take their sum.

    Concretely: SU(2)_L generated by Σ_{12}; we need 3 generators that
    close to SU(2). Use the Cartan-block construction:
      J_x = (1/2)(Σ_{34} + Σ_{56} + Σ_{1,7})  — guess
    This is a SEARCH candidate; will check closure.
    """
    # First guess: 3 mutually-orthogonal Cartan-direction sums.
    Jx = sigmas[(1, 2)] + sigmas[(3, 4)] + sigmas[(5, 6)]
    Jy = sigmas[(1, 3)] + sigmas[(2, 4)] + sigmas[(5, 7)]
    Jz = sigmas[(1, 4)] - sigmas[(2, 3)] + sigmas[(6, 7)]
    return Jx, Jy, Jz


def su2_subregular_so7(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (H): subregular SU(2) in so(7), constructed via
    sl(2)-triple of the subregular nilpotent.

    For B_3 = so(7), the subregular nilpotent orbit has Bala-Carter label
    B_2 (one row of length 5, two rows of length 1 in the partition [5,1,1]).
    Under the corresponding SU(2), the 7-vector decomposes as 5 + 1 + 1,
    and the 8-spinor decomposes as ... (to be measured numerically).

    Construction: nilpotent N = Σ_{12} + Σ_{34} + Σ_{56} (using three
    'consecutive' commuting bivectors) — but those commute, giving abelian
    not nilpotent. Try instead a step-up chain:
      N = Σ_{12} + Σ_{23} + Σ_{34} + Σ_{45} + Σ_{56}
    This nilpotent has the cyclic structure for a principal-in-sl(5)-block
    embedding. Check empirically.
    """
    # Nilpotent step-up chain
    Jplus = (sigmas[(1, 2)] + sigmas[(2, 3)] + sigmas[(3, 4)]
              + sigmas[(4, 5)] + sigmas[(5, 6)])
    Jminus = -Jplus.T  # since Σ_ij are skew-symmetric, J^T = -J
    Jx = 0.5 * (Jplus + Jminus)
    Jy = -0.5 * (Jplus - Jminus)  # skew part
    # Need Jz = [Jx, Jy]/c for closure
    Jz = Jx @ Jy - Jy @ Jx
    return Jx, Jy, Jz


def su2_spin4_L_plus_spin3(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (I): diagonal of SU(2)_L (⊂ Spin(4) on γ_1..γ_4)
    and SU(2)_3 (= Spin(3) on γ_5..γ_7).

    Spin(7) ⊃ Spin(4) × Spin(3). The 8-spinor = S(so(4)) ⊗ S(so(3))
    = (2_L ⊕ 2_R) ⊗ 2_3 = (2_L, 2_3) ⊕ (2_R, 2_3).

    Under SU(2)_L (only): each 2_L is spin-1/2; 2_R untouched. So
    SU(2)_L gives (2_L,1)⊕(2_L,1) ⊕ (1,2_R) ⊕ ... hmm just gets 4·(j=1/2).

    Under DIAGONAL of SU(2)_L and SU(2)_3:
      (2_L, 2_3) → 2_L ⊗ 2_3 = 1 ⊕ 3   (j=0 + j=1, dim 1 + 3 = 4)
      (2_R, 2_3) → 2_3 (since SU(2)_L doesn't act on 2_R; 2_R untouched)
                 So this is a 4-dim block where SU(2)_diag acts only
                 via 2_3, giving 2·(2) (two doublets).

    Total: S=8 = 1 ⊕ 3 ⊕ 2 ⊕ 2 = {1, 2, 2, 3}. ★

    Spin(4) sub-SU(2)_L = J^L_x = Σ_{12} + Σ_{34}, etc. (using indices 1..4)
    Spin(3) on indices {5,6,7}: J^3_x = Σ_{67}, J^3_y = Σ_{57}, J^3_z = Σ_{56}.
    """
    # SU(2)_L: self-dual generators of so(4) on indices {1,2,3,4}
    JLx = sigmas[(1, 2)] + sigmas[(3, 4)]
    JLy = sigmas[(1, 3)] - sigmas[(2, 4)]
    JLz = sigmas[(1, 4)] + sigmas[(2, 3)]
    # SU(2)_3 = Spin(3) on indices {5,6,7}
    J3x = sigmas[(6, 7)]
    J3y = sigmas[(5, 7)]
    J3z = sigmas[(5, 6)]
    # Normalize: SU(2)_L is built from sum of 2 Σ's (each ε(j=1/2)), so
    # commutator [JLx, JLy] = 2·JLz (need to halve for canonical), while
    # [J3x, J3y] = 1·J3z. So J_diag = (1/2) JL + J3 for matched scaling.
    Jx = 0.5 * JLx + J3x
    Jy = 0.5 * JLy + J3y
    Jz = 0.5 * JLz + J3z
    return Jx, Jy, Jz


def su2_spin4_R_plus_spin3(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (J): diagonal of SU(2)_R (anti-self-dual ⊂ Spin(4))
    and SU(2)_3 (= Spin(3) on γ_5..γ_7).

    Same decomposition expected: 8 = 1 ⊕ 2 ⊕ 2 ⊕ 3.
    """
    JRx = sigmas[(1, 2)] - sigmas[(3, 4)]
    JRy = sigmas[(1, 3)] + sigmas[(2, 4)]
    JRz = sigmas[(1, 4)] - sigmas[(2, 3)]
    J3x = sigmas[(6, 7)]
    J3y = sigmas[(5, 7)]
    J3z = sigmas[(5, 6)]
    Jx = 0.5 * JRx + J3x
    Jy = 0.5 * JRy + J3y
    Jz = 0.5 * JRz + J3z
    return Jx, Jy, Jz


def su2_spin4L_only(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2)_L of Spin(4) ⊂ Spin(7) on first 4 indices only.

    The 8 = (2_L, 2_3) ⊕ (2_R, 2_3) restricted to SU(2)_L only:
      (2_L, 2_3) → 2_L · 2_3 = 2 blocks of j=1/2 (i.e., 2·(2))
      (2_R, 2_3) → 1_L · 2_R · 2_3 = 4-dim block fixed under SU(2)_L → 4·(j=0)
    Total: 2·(2) + 4·(1) = 4 + 4 = 8.
    """
    Jx = 0.5 * (sigmas[(1, 2)] + sigmas[(3, 4)])
    Jy = 0.5 * (sigmas[(1, 3)] - sigmas[(2, 4)])
    Jz = 0.5 * (sigmas[(1, 4)] + sigmas[(2, 3)])
    return Jx, Jy, Jz


def su2_g2_compatible(sigmas: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SU(2) embedding (F): G_2-compatible SU(2) commuting with AGL Z_3
    lift permutation (e_1)(e_2 e_4 e_6)(e_3 e_5 e_7).

    Z_3 orbits on the 21 Σ-generators:
      fixed (touching e_1 with the other in the e_1 cycle): only Σ_{i,i}? no
      orbit_e1: Σ pairs (1, 2), (1, 4), (1, 6) form one Z_3 orbit;
                          (1, 3), (1, 5), (1, 7) form another
      ... etc.

    The diagonal sum over a Z_3-orbit is Z_3-invariant.
      J^F_x = Σ_{1,2} + Σ_{1,4} + Σ_{1,6}    (Z_3-symmetrized "polar X")
      J^F_y = Σ_{1,3} + Σ_{1,5} + Σ_{1,7}    (Z_3-symmetrized "polar Y")
      J^F_z = Σ_{2,3} + Σ_{4,5} + Σ_{6,7}    (Z_3-symmetrized "equatorial Z")

    Won't necessarily close as SU(2); test [J_x, J_y] = c J_z.
    """
    Jx = sigmas[(1, 2)] + sigmas[(1, 4)] + sigmas[(1, 6)]
    Jy = sigmas[(1, 3)] + sigmas[(1, 5)] + sigmas[(1, 7)]
    Jz = sigmas[(2, 3)] + sigmas[(4, 5)] + sigmas[(6, 7)]
    return Jx, Jy, Jz


def check_su2_closure(Jx, Jy, Jz, tol: float = 1e-8) -> tuple[bool, float]:
    """Verify [J_x, J_y] = c · J_z for some scalar c. Return (closes, c)."""
    comm = Jx @ Jy - Jy @ Jx
    # Find scaling c
    norm_jz = np.linalg.norm(Jz)
    if norm_jz < tol:
        return False, 0.0
    c = np.trace(comm @ Jz.T) / np.trace(Jz @ Jz.T)
    residual = np.linalg.norm(comm - c * Jz)
    return (residual < tol * max(1.0, np.linalg.norm(comm))), float(c)


def diag_casimir_decomposition(Jx, Jy, Jz, label: str) -> dict:
    """Compute Casimir C = -(J_x² + J_y² + J_z²) and report spectrum.

    Negative sign because J_a are anti-Hermitian (skew-symmetric); their
    squares are negative-semidefinite, so -J_a² is positive.

    For SU(2) irrep of dim 2j+1, C = j(j+1) on that irrep (with proper
    normalization).
    """
    C = -(Jx @ Jx + Jy @ Jy + Jz @ Jz)
    C = 0.5 * (C + C.T)
    eigs = np.linalg.eigvalsh(C)
    # Group by rounded eigenvalue
    groups = defaultdict(int)
    for e in eigs:
        groups[round(float(e), 4)] += 1
    return {
        "label": label,
        "casimir": C,
        "eigenvalues": [round(float(e), 4) for e in eigs],
        "multiplicities": dict(groups),
    }


def fit_j_from_casimir(cas_eigval: float) -> float | None:
    """Given a Casimir eigenvalue c, return j such that c = κ · j(j+1)
    for the normalization κ deduced empirically. Tries a few standard
    normalizations.

    For raw Σ_ij = (1/2) γ_i γ_j, the "fundamental" sigma has eigenvalues
    ±(i/2). Sum of three Σ for a so(3) gives Casimir with specific
    normalization. Just report whether eigval looks like κ · j(j+1) for
    j ∈ {0, 1/2, 1, 3/2, 2}.
    """
    if cas_eigval < 1e-8:
        return 0.0
    return None  # Inspector-only, see numeric output


def build_irrep_projectors(C: np.ndarray, tol: float = 1e-6) -> dict[float, np.ndarray]:
    """Spectral decomposition of Casimir C; return {eigval: projector}."""
    eigvals, eigvecs = np.linalg.eigh(C)
    groups = defaultdict(list)
    for k, ev in enumerate(eigvals):
        groups[round(float(ev), 4)].append(k)
    projectors = {}
    for ev, idx in groups.items():
        P = np.zeros_like(C)
        for k in idx:
            v = eigvecs[:, k:k+1]
            P += v @ v.T
        projectors[ev] = P
    return projectors


def project_matrix_onto_subspace(A: np.ndarray, P: np.ndarray) -> float:
    """Compute fractional weight ||P A P||_F^2 / ||A||_F^2.

    For a Lie-algebra element A acting on S=8 broken into SU(2) irrep
    subspaces by projector P, this fraction measures how much of A's
    'energy' lives in that irrep block.
    """
    blocked = P @ A @ P
    tot = float(np.trace(A @ A.T))
    if tot < 1e-12:
        return 0.0
    return float(np.trace(blocked @ blocked.T)) / tot


# ============================================================================
# Main analysis
# ============================================================================

def get_compendium_walks() -> list[dict]:
    walks_dict = bfs_shortest_walks(max_length=25)
    rows = []
    seen = set()
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key in seen or key not in walks_dict:
            continue
        seen.add(key)
        rows.append({
            "name": entry["name"], "p": key[0], "q": key[1],
            "n_q": entry["n_q"], "sector": entry.get("sector", ""),
            "walk": walks_dict[key],
            "L": len(walks_dict[key]) - 1,
        })
    return rows


def main():
    print("=" * 78)
    print("Q15 — Spin(7) on Cl(0,7) spinor S=8: n_q derivation pilot")
    print("=" * 78)
    print()

    # ---- Step 1: Build and verify Cl(0,7) gammas --------------------------
    print("[Step 1] Building Cl(0,7) gamma matrices γ_1..γ_7 on S=8...")
    gammas = build_gammas()
    verify_clifford(gammas)
    print()

    # ---- Step 2: Build Spin(7) generators Σ_ij ----------------------------
    print("[Step 2] Building 21 Spin(7) generators Σ_ij = (1/2)γ_i γ_j...")
    sigmas = build_sigma_generators(gammas)
    print(f"  ✓ Built {len(sigmas)} generators (expect 21 = C(7,2))")
    verify_so7_commutators(sigmas)
    print()

    # ---- Step 3: SU(2) ⊂ Spin(7) candidates and S=8 decomposition --------
    print("[Step 3] Candidate SU(2) ⊂ Spin(7) embeddings — S=8 decomposition")
    print("-" * 78)
    candidates = [
        ("A: so(3)_diag (octo)",     su2_so3_diag_octonion(sigmas)),
        ("B: so(3)_anti-diag",       su2_so3_anti_diag_octonion(sigmas)),
        ("C: so(3)_{123}",           su2_so3_123(sigmas)),
        ("D: so(3)_{456}",           su2_so3_456(sigmas)),
        ("E: principal-chain",       su2_principal_so7(sigmas)),
        ("F: G_2-Z_3-symmetrized",   su2_g2_compatible(sigmas)),
        ("G: diag SU(2)_L×Sp(2)",    su2_diag_L_Sp2(sigmas)),
        ("H: subregular-chain",      su2_subregular_so7(sigmas)),
        ("I: Spin(4)_L ⊕ Spin(3)",   su2_spin4_L_plus_spin3(sigmas)),
        ("J: Spin(4)_R ⊕ Spin(3)",   su2_spin4_R_plus_spin3(sigmas)),
        ("K: Spin(4)_L only",        su2_spin4L_only(sigmas)),
    ]

    decompositions = {}
    for label, (Jx, Jy, Jz) in candidates:
        closes, c_jz = check_su2_closure(Jx, Jy, Jz)
        decomp = diag_casimir_decomposition(Jx, Jy, Jz, label)
        mults_str = ", ".join(f"{e:.3f}×{m}" for e, m in
                                sorted(decomp['multiplicities'].items()))
        closure_str = f"CLOSED (c={c_jz:+.3f})" if closes else f"OPEN (c={c_jz:+.3f})"
        print(f"  {label:<28} | {closure_str:<24} | Casimir spec: {mults_str}")
        decompositions[label] = (decomp, closes, (Jx, Jy, Jz))
    print()

    # ---- Spinor decomposition under each SU(2) ----------------------------
    print("[Step 3b] Spinor-irrep content per SU(2) candidate")
    print("-" * 78)
    print(f"  {'Candidate':<28} {'Casimir eigenvalues × multiplicity':<48}")
    print("-" * 78)
    for label, (decomp, closes, JJJ) in decompositions.items():
        groups = decomp['multiplicities']
        parts = [f"{ev:.3f}×{mult}" for ev, mult in sorted(groups.items())]
        marker = " ★" if closes else ""
        # Try to interpret as j(j+1) → 2j+1 irreps
        if closes:
            irreps = []
            for ev, mult in sorted(groups.items()):
                if abs(ev) < 1e-4:
                    irreps.append(f"{mult}·(j=0,d=1)")
                else:
                    # Try j(j+1) = ev → j = (-1 + sqrt(1+4ev))/2
                    j = (-1 + np.sqrt(1 + 4 * ev)) / 2
                    d = int(round(2 * j + 1))
                    n_copies = mult // d if d > 0 and mult % d == 0 else None
                    if n_copies:
                        irreps.append(f"{n_copies}·(j={j:.2f},d={d})")
                    else:
                        irreps.append(f"mult={mult}@C={ev:.2f}")
            irrep_str = " ⊕ ".join(irreps)
        else:
            irrep_str = "(not SU(2))"
        print(f"  {label:<28} {', '.join(parts):<28}  {irrep_str}{marker}")
    print()

    # ---- Step 4: project compendium walks onto SU(2) irreps -------------
    print("[Step 4] Project walks onto SU(2) irrep blocks, correlate with n_q")
    print("=" * 78)
    rows = get_compendium_walks()
    print(f"  Loaded {len(rows)} compendium walks "
          f"(distinct (|p|,|q|) classes)\n")

    for label, (decomp, closes, (Jx, Jy, Jz)) in decompositions.items():
        if not closes:
            continue  # only analyze SU(2)-closed candidates
        print(f"  === Candidate {label} ===")
        C = decomp['casimir']
        projectors = build_irrep_projectors(C)
        ev_list = sorted(projectors.keys())
        # Each ev corresponds to one irrep block. Dim of block = mult.
        block_dims = {ev: int(np.trace(P).round()) for ev, P in projectors.items()}
        print(f"    Block structure: " +
              ", ".join(f"C={ev:.2f}(dim {block_dims[ev]})" for ev in ev_list))

        # Collect per-walk projection data
        per_walk = []
        for r in rows:
            A_walk, A_u, A_v = walk_algebra_elements_S(r['walk'], sigmas)
            frac_walk = {ev: project_matrix_onto_subspace(A_walk, P)
                          for ev, P in projectors.items()}
            frac_v = {ev: project_matrix_onto_subspace(A_v, P)
                       for ev, P in projectors.items()}
            per_walk.append({
                **r,
                "frac_walk": frac_walk,
                "frac_v": frac_v,
            })

        # Correlation analysis
        n_q_arr = np.array([r['n_q'] for r in per_walk])
        for ev in ev_list:
            dim_block = block_dims[ev]
            xw = np.array([r['frac_walk'][ev] for r in per_walk])
            xv = np.array([r['frac_v'][ev] for r in per_walk])
            rw = np.corrcoef(xw, n_q_arr)[0, 1] if xw.std() > 0 else 0.0
            rv = np.corrcoef(xv, n_q_arr)[0, 1] if xv.std() > 0 else 0.0
            print(f"    Block C={ev:>+6.2f} (dim {dim_block}): "
                  f"r(frac_walk, n_q) = {rw:+.3f}, "
                  f"r(frac_v, n_q) = {rv:+.3f}")

        # Mean fractional weight per n_q sector
        print(f"    Mean fractional weight (A_v) per n_q sector:")
        by_nq = defaultdict(list)
        for r in per_walk:
            by_nq[r['n_q']].append(r)
        header = "    sector " + " ".join(f"C={ev:>+5.2f}(d{block_dims[ev]})"
                                            for ev in ev_list)
        print(header)
        for nq in sorted(by_nq.keys()):
            vals = []
            for ev in ev_list:
                m = np.mean([r['frac_v'][ev] for r in by_nq[nq]])
                vals.append(f"     {m:.3f}      ")
            print(f"    n_q={nq:<2} ({len(by_nq[nq]):>2})  " + " ".join(vals))
        print()


def direct_spinor_eigenvalue_analysis(sigmas: dict) -> None:
    """Skip SU(2) decomposition: just look at A^S eigenvalues directly.

    A^S = Σ_e (sign) Σ_e is real skew-symmetric 8×8, so its eigenvalues
    come in pairs ±iλ_k with 4 distinct |λ_k| (or fewer with multiplicity).
    The multiset {|λ_1|, ..., |λ_4|} is a 4-tuple invariant per walk.

    Test correlations of various invariants of this multiset with n_q:
      - sum of |λ_k|: trace norm
      - max |λ_k|: dominant eigenvalue
      - count of distinct |λ_k| values: spectral richness
      - product of |λ_k|: pfaffian-related
      - ratios |λ_1| / |λ_2| etc.
    """
    print("=" * 78)
    print("[Direct] Eigenvalue analysis of A^S = Σ_e Σ_e on S=8 (no SU(2) projection)")
    print("=" * 78)
    rows = get_compendium_walks()
    data = []
    for r in rows:
        A_walk, A_u, A_v = walk_algebra_elements_S(r['walk'], sigmas)
        # Eigenvalues of skew-symmetric A: pairs ±iλ_k
        eigs_walk = np.linalg.eigvals(A_walk)
        # Take positive imaginary parts
        lam_walk = sorted([abs(e.imag) for e in eigs_walk if e.imag > 1e-8],
                            reverse=True)
        eigs_v = np.linalg.eigvals(A_v)
        lam_v = sorted([abs(e.imag) for e in eigs_v if e.imag > 1e-8],
                        reverse=True)
        # Pad to length 4 with zeros
        while len(lam_walk) < 4:
            lam_walk.append(0.0)
        while len(lam_v) < 4:
            lam_v.append(0.0)
        # Count distinct |λ| values (rounded)
        n_distinct_walk = len({round(x, 3) for x in lam_walk if x > 1e-6})
        n_distinct_v = len({round(x, 3) for x in lam_v if x > 1e-6})
        # Trace of A^2 = -2 Σ |λ_k|^2 (for skew-symmetric)
        tr_A2_walk = -float(np.trace(A_walk @ A_walk))
        tr_A2_v = -float(np.trace(A_v @ A_v))
        # Pfaffian (4-form invariant) — for 8×8 skew, Pf(A)^2 = det(A)
        det_walk = float(np.linalg.det(A_walk).real)
        det_v = float(np.linalg.det(A_v).real)
        data.append({
            **r,
            "lam_walk": lam_walk, "lam_v": lam_v,
            "max_walk": lam_walk[0], "max_v": lam_v[0],
            "n_distinct_walk": n_distinct_walk,
            "n_distinct_v": n_distinct_v,
            "tr_A2_walk": tr_A2_walk, "tr_A2_v": tr_A2_v,
            "det_walk": det_walk, "det_v": det_v,
        })

    print(f"\n  {'(p,q)':<8} {'name':<10} {'n_q':<4} {'L':<3} "
          f"{'|λ_v| sorted':<35} {'dist':<5} {'|A_v|²':<10}")
    print("-" * 100)
    for d in data:
        lam_str = "[" + ", ".join(f"{x:.2f}" for x in d['lam_v']) + "]"
        print(f"  ({d['p']:>2},{d['q']:>2})  {d['name']:<10} {d['n_q']:<4} "
              f"{d['L']:<3} {lam_str:<35} {d['n_distinct_v']:<5} "
              f"{d['tr_A2_v']:<10.3f}")
    print()

    # Correlations
    n_q = np.array([d['n_q'] for d in data])
    print("  Correlations with n_q:")
    for key in ["max_walk", "max_v", "n_distinct_walk", "n_distinct_v",
                "tr_A2_walk", "tr_A2_v", "det_walk", "det_v"]:
        x = np.array([d[key] for d in data])
        if x.std() > 0:
            r = np.corrcoef(x, n_q)[0, 1]
            print(f"    r({key:<18}, n_q) = {r:+.4f}")

    # Group by n_q
    print()
    print("  Mean |λ_v| spectrum per n_q sector:")
    by_nq = defaultdict(list)
    for d in data:
        by_nq[d['n_q']].append(d)
    print(f"    {'sector':<12} {'count':<5} {'mean |λ_v|_1':<14} "
          f"{'mean |λ_v|_2':<14} {'mean |λ_v|_3':<14} {'mean |λ_v|_4':<14} "
          f"{'mean n_distinct':<15}")
    for nq in sorted(by_nq.keys()):
        entries = by_nq[nq]
        means = [np.mean([d['lam_v'][k] for d in entries]) for k in range(4)]
        nd = np.mean([d['n_distinct_v'] for d in entries])
        print(f"    n_q={nq:<8} {len(entries):<5} "
              f"{means[0]:<14.3f} {means[1]:<14.3f} {means[2]:<14.3f} "
              f"{means[3]:<14.3f} {nd:<15.2f}")


def spinor_wilson_loop_holonomy(sigmas: dict) -> None:
    """Compute G_W = exp(A^S · t) on the 8-dim spinor and report eigenvalue
    spectrum classification by n_q sector.

    Parallel to Q8 (Wilson loop on V_7 = vector rep), but in S=8 = spinor rep.
    The 8 eigenvalues lie on the unit circle (Spin(7) ⊂ SO(8) → unitary
    action on spinor). The multiset of phases and their multiplicities
    give a richer fingerprint than the linear projection used in Step 4.
    """
    from nwt_substrate.isa.constants import ALPHA_NWT
    from scipy.linalg import expm
    print()
    print("=" * 78)
    print("[Spinor Wilson] G_W^S=8 = exp(A^S · t) eigenvalue spectrum, sweep t")
    print("=" * 78)
    rows = get_compendium_walks()

    # Sweep couplings (same as Q8)
    couplings = [("α", ALPHA_NWT), ("√α", np.sqrt(ALPHA_NWT)),
                  ("0.5", 0.5), ("1.0", 1.0), ("π/L", None)]

    for cname, t_fixed in couplings:
        print(f"\n  --- coupling t = {cname} ---")
        per = []
        for r in rows:
            A_walk, A_u, A_v = walk_algebra_elements_S(r['walk'], sigmas)
            L = r['L']
            t = np.pi / L if t_fixed is None else t_fixed
            G = expm(A_walk * t)
            eigs = np.linalg.eigvals(G)
            phases = sorted(abs(np.angle(e)) for e in eigs)
            n_one = sum(1 for e in eigs if abs(np.angle(e)) < 1e-6)
            n_distinct = len({round(p, 4) for p in phases if p > 1e-6})
            tr_G = float(np.trace(G).real)
            casimir_2 = float(sum(np.angle(e) ** 2 for e in eigs))
            casimir_4 = float(sum(np.angle(e) ** 4 for e in eigs))
            max_phase = float(max(abs(np.angle(e)) for e in eigs))
            nontrivial = 8 - n_one
            per.append({
                **r, "t": t, "n_eig_one": n_one,
                "n_distinct": n_distinct, "trace": tr_G,
                "casimir_2": casimir_2, "casimir_4": casimir_4,
                "max_phase": max_phase, "nontrivial_dim": nontrivial,
            })

        # Correlations vs n_q
        nq = np.array([r['n_q'] for r in per])
        for obs in ["nontrivial_dim", "trace", "casimir_2",
                     "casimir_4", "max_phase", "n_distinct"]:
            x = np.array([r[obs] for r in per])
            if x.std() > 0:
                rc = np.corrcoef(x, nq)[0, 1]
                print(f"    r({obs:<20}, n_q) = {rc:+.4f}")

        # Also do v-weighted version
        per_v = []
        for r in rows:
            A_walk, A_u, A_v = walk_algebra_elements_S(r['walk'], sigmas)
            L = r['L']
            t = np.pi / L if t_fixed is None else t_fixed
            G = expm(A_v * t)
            eigs = np.linalg.eigvals(G)
            casimir_2 = float(sum(np.angle(e) ** 2 for e in eigs))
            casimir_4 = float(sum(np.angle(e) ** 4 for e in eigs))
            max_phase = float(max(abs(np.angle(e)) for e in eigs))
            per_v.append({**r, "casimir_2": casimir_2, "casimir_4": casimir_4,
                           "max_phase": max_phase})
        for obs in ["casimir_2", "casimir_4", "max_phase"]:
            x = np.array([r[obs] for r in per_v])
            if x.std() > 0:
                rc = np.corrcoef(x, nq)[0, 1]
                print(f"    [v-weighted] r({obs:<14}, n_q) = {rc:+.4f}")


def brute_force_su2_search(sigmas: dict, target_decomp: list[int]) -> None:
    """Search over linear combinations of 3 Σ-generators for an SU(2)
    sub-algebra whose action on S=8 gives the specified target decomposition
    (multiset of irrep dimensions summing to 8).

    Restricted search: each J_a is a sum of K Σ_ij with coefficients ±1.
    """
    from itertools import combinations as combs
    print()
    print("=" * 78)
    print(f"[Brute force] Searching for SU(2) ⊂ Spin(7) with S=8 decomp = "
          f"{target_decomp}")
    print("=" * 78)
    target_C_values = sorted([(d - 1) / 2 * ((d - 1) / 2 + 1)
                                for d in target_decomp])
    print(f"  Target Casimir eigenvalues (each d-dim irrep → j(j+1)): "
          f"{[round(c, 3) for c in target_C_values]}")

    # Enumerate single-generator SU(2)s (boring), then 2-sum, 3-sum
    sigma_keys = sorted(sigmas.keys())

    matches = []
    # Search 3-sum: each J_a is sum of 3 distinct Σs with ±1 signs
    # That's too combinatorial — restrict to a sample
    # Try sums of 1, 2, or 3 generators per J_a
    import random
    random.seed(42)
    n_tries = 5000
    for _ in range(n_tries):
        # Random 3 Σs for each Ja
        Jx = sum(random.choice([-1, 1]) * sigmas[k]
                  for k in random.sample(sigma_keys, k=random.choice([1, 2, 3])))
        Jy = sum(random.choice([-1, 1]) * sigmas[k]
                  for k in random.sample(sigma_keys, k=random.choice([1, 2, 3])))
        Jz = sum(random.choice([-1, 1]) * sigmas[k]
                  for k in random.sample(sigma_keys, k=random.choice([1, 2, 3])))
        closes, _ = check_su2_closure(Jx, Jy, Jz, tol=1e-6)
        if not closes:
            continue
        decomp = diag_casimir_decomposition(Jx, Jy, Jz, "rand")
        # Compare Casimir eigvalues
        evs = sorted(decomp['multiplicities'].items())
        # Build multiset of irrep dimensions
        irrep_dims = []
        for ev, mult in evs:
            if abs(ev) < 1e-4:
                irrep_dims.extend([1] * mult)
            else:
                j = (-1 + np.sqrt(1 + 4 * ev)) / 2
                d = int(round(2 * j + 1))
                if d > 0 and mult % d == 0:
                    irrep_dims.extend([d] * (mult // d))
                else:
                    irrep_dims.append(None)
        if sorted(filter(None, irrep_dims)) == sorted(target_decomp):
            matches.append((Jx, Jy, Jz, decomp))
            if len(matches) <= 5:
                print(f"  ✓ MATCH found! mult: {decomp['multiplicities']}")
        if len(matches) >= 20:
            break
    print(f"  Total matches found: {len(matches)} (in {n_tries} tries)")
    return matches


if __name__ == "__main__":
    main()

    print()
    print()

    # Direct eigenvalue analysis
    T = make_octonion_table()
    gammas = build_gammas()
    sigmas = build_sigma_generators(gammas)
    direct_spinor_eigenvalue_analysis(sigmas)

    # Wilson loop holonomy on the spinor (parallel to Q8 on V_7)
    spinor_wilson_loop_holonomy(sigmas)

    print()

    # Brute-force search for SU(2) with target 1+2+5 = 8 decomposition
    brute_force_su2_search(sigmas, [1, 2, 5])

    # And for 3+5 = 8 decomposition
    brute_force_su2_search(sigmas, [3, 5])
