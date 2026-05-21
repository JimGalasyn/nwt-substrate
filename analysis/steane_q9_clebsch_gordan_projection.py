"""Q9 — Clebsch-Gordan projection of compendium walks under so(3)_diag.

Following Q8's modest Wilson-loop signal, this script implements the
framework memo's intended Step 4: explicit projection of each walk's
algebra-element A = Σ_e J_e onto the so(3)_diag irrep pieces of so(7)
(the 21-dim adjoint rep).

so(7) under so(3)_diag decomposes as:

  σ_0 (polar matter bridging):    (3)
  σ_1 (polar antimatter bridging): (3)
  σ_2 (so(3)_{123} adjoint):      (3)
  σ_3 (so(3)_{456} adjoint):      (3)
  σ_4 ⊕ σ_5 ⊕ σ_6 (cross-block):  (3) ⊗ (3) = (1) ⊕ (3) ⊕ (5)

Total:  5·(3) ⊕ (1) ⊕ (5) = 15 + 1 + 5 = 21 ✓

Carrier-knot candidate mapping (per framework memo):
  n_q = 0 (lepton):    (1) singlet — Casimir eigenvalue 0
  n_q = 3 (hyperon):   (3) vector  — Casimir eigenvalue 2
  n_q = 5 (nucleon):   (5) sym-traceless — Casimir eigenvalue 6
  n_q = 2 (meson):     does NOT fit so(3)_diag irreps cleanly (so(3) has no 2-dim rep)
                       → empirical test: where do mesons land?

Pipeline:
  1. Build so(3)_diag generators J_x, J_y, J_z on R^7 (7x7 matrices).
  2. Compute the Casimir C = ad_{J_x}² + ad_{J_y}² + ad_{J_z}² as a
     21x21 operator on so(7) in the {J_{ij} : i<j} basis.
  3. Spectrally decompose C — should have eigenvalues {0, 2, 6} with
     multiplicities {1, 15, 5}.
  4. Build spectral projectors P_0, P_1, P_2.
  5. For each compendium walk, compute A = Σ ε_e J_e (signed by
     direction) and project onto each P_j.
  6. Tabulate ||P_j A||² weights and look for clustering by Paper 11 n_q.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q9_clebsch_gordan_projection.py
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.particles.compendium import COMPENDIUM


N = 7


def so7_basis() -> list[np.ndarray]:
    """The 21 generators J_{ij} = E_{ij} - E_{ji} for 0 ≤ i < j ≤ 6,
    each normalized to ⟨J, J⟩ = -Tr(J²)/2 = 1.
    """
    basis = []
    for i, j in combinations(range(N), 2):
        J = np.zeros((N, N), dtype=float)
        J[i, j] = +1.0
        J[j, i] = -1.0
        # ⟨J, J⟩ = -(1/2) Tr(J²) = 1, so J is already orthonormal under
        # this inner product. Verify: J² has -1 at (i,i) and (j,j), so
        # Tr(J²) = -2, and -(1/2)*-2 = 1. ✓
        basis.append(J)
    return basis


SO7_BASIS = so7_basis()
EDGE_INDEX = {}
for k, (i, j) in enumerate(combinations(range(N), 2)):
    EDGE_INDEX[(i, j)] = k
    EDGE_INDEX[(j, i)] = k  # also map reversed


def so3_diag_generators() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """J_x, J_y, J_z for so(3)_diag = diagonal embedding of so(3) into
    so(3)_{123} × so(3)_{456} ⊂ so(7).

    so(3)_{123} acts on V = span{e_1, e_2, e_3}:
      L_x = J_{2,3} = E_{23} - E_{32}
      L_y = J_{3,1} = E_{31} - E_{13}
      L_z = J_{1,2} = E_{12} - E_{21}

    so(3)_{456} acts on V' = span{e_4, e_5, e_6}:
      L_x' = J_{5,6}
      L_y' = J_{6,4}
      L_z' = J_{4,5}

    so(3)_diag generators:
      J_x = L_x + L_x'
      J_y = L_y + L_y'
      J_z = L_z + L_z'
    """
    def gen(i, j):
        M = np.zeros((N, N))
        M[i, j] = 1.0
        M[j, i] = -1.0
        return M

    # so(3) convention: [L_x, L_y] = +L_z. Requires L_y = E_{13} - E_{31}
    # (not E_{31} - E_{13}).
    L_x = gen(2, 3)
    L_y = gen(1, 3)  # = E_{13} - E_{31}, makes [L_x, L_y] = +L_z
    L_z = gen(1, 2)
    L_xp = gen(5, 6)
    L_yp = gen(4, 6)
    L_zp = gen(4, 5)

    J_x = L_x + L_xp
    J_y = L_y + L_yp
    J_z = L_z + L_zp
    return J_x, J_y, J_z


def adjoint_action(J: np.ndarray, basis: list[np.ndarray]) -> np.ndarray:
    """Return the matrix of ad_J = [J, ·] on so(7) in the given basis.
    The matrix is (|basis|, |basis|): adjoint[k, l] = ⟨basis[k], [J, basis[l]]⟩.
    """
    n = len(basis)
    M = np.zeros((n, n))
    for l in range(n):
        commutator = J @ basis[l] - basis[l] @ J
        for k in range(n):
            # Inner product ⟨A, B⟩ = -(1/2) Tr(AB)
            M[k, l] = -0.5 * np.trace(basis[k] @ commutator)
    return M


def walk_to_algebra_element(walk: list[int]) -> np.ndarray:
    """Compute A = Σ_{directed edges} J_{a→b} where J_{a→b} is the
    so(7) generator with E_{ab} = +1, E_{ba} = -1. This is the LEADING
    so(7) algebra element of the walk's Wilson loop at small coupling.

    Returns the 21-dim coordinate vector of A in the SO7_BASIS basis
    (i.e., A = Σ_k coords[k] * SO7_BASIS[k]).
    """
    coords = np.zeros(len(SO7_BASIS))
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        k = EDGE_INDEX[(a, b)]
        i_, j_ = min(a, b), max(a, b)
        # If direction a→b matches basis convention i<j → +1, else -1
        sign = +1 if a < b else -1
        coords[k] += sign
    return coords


def main():
    print("=" * 78)
    print("Q9 — Clebsch-Gordan projection of compendium walks under so(3)_diag")
    print("=" * 78)
    print()

    # ---- Build so(3)_diag generators ------------------------------------
    J_x, J_y, J_z = so3_diag_generators()
    print("so(3)_diag commutation check [J_x, J_y] = J_z?")
    comm_xy = J_x @ J_y - J_y @ J_x
    err = np.linalg.norm(comm_xy - J_z)
    print(f"  ||[J_x, J_y] - J_z|| = {err:.6f}  (should be 0)")
    print()

    # ---- Build Casimir operator on so(7) -------------------------------
    print("Building Casimir C = ad_{J_x}² + ad_{J_y}² + ad_{J_z}² on so(7)...")
    ad_x = adjoint_action(J_x, SO7_BASIS)
    ad_y = adjoint_action(J_y, SO7_BASIS)
    ad_z = adjoint_action(J_z, SO7_BASIS)
    # Casimir = ad_x² + ad_y² + ad_z² but with a SIGN:
    # The Casimir for so(3) is sum L_i² where L_i are the so(3) reps.
    # In adjoint, ad_J : so(7) → so(7), and Casimir = sum (ad_J)² acts
    # on (j)-irreps as -j(j+1) (the minus sign comes from anti-Hermiticity).
    # We'll take the absolute value as eigenvalues.
    C = ad_x @ ad_x + ad_y @ ad_y + ad_z @ ad_z
    # Symmetrize for numerical stability
    C = 0.5 * (C + C.T)

    eigvals, eigvecs = np.linalg.eigh(C)
    print(f"Casimir eigenvalues (sorted):")
    # Group by rounded eigenvalue
    eig_groups = defaultdict(list)
    for k, ev in enumerate(eigvals):
        eig_groups[round(ev, 4)].append(k)
    for ev in sorted(eig_groups.keys()):
        # j(j+1) = -ev (since ad-squared is negative-definite on irreps)
        # eigenvalues should be 0, -2, -6 (j=0, 1, 2)
        abs_ev = abs(ev)
        # j(j+1) = abs_ev → j = (-1 + sqrt(1 + 4 abs_ev)) / 2
        if abs_ev < 1e-6:
            j = 0
        else:
            j = (-1 + np.sqrt(1 + 4 * abs_ev)) / 2
        dim = 2 * round(j) + 1  # SU(2)-style dim
        print(f"  eigenvalue = {ev:+.4f}  →  j ≈ {j:.3f}, dim 2j+1 ≈ {dim}, "
              f"multiplicity = {len(eig_groups[ev])}")
    print()

    # ---- Identify spectral projectors -----------------------------------
    # Build projectors onto each eigenvalue subspace
    projectors = {}
    for ev, indices in eig_groups.items():
        P = np.zeros_like(C)
        for k in indices:
            v = eigvecs[:, k:k + 1]
            P += v @ v.T
        projectors[ev] = P

    # ---- For each compendium walk, project A and tabulate weights -------
    walks = bfs_shortest_walks(max_length=25)

    seen = set()
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key in seen or key not in walks:
            continue
        seen.add(key)
        walk = walks[key]
        A_coords = walk_to_algebra_element(walk)
        total_norm = np.linalg.norm(A_coords) ** 2
        weights = {}
        for ev, P in projectors.items():
            w = float(np.linalg.norm(P @ A_coords) ** 2)
            weights[ev] = w
        rows.append({
            "name": entry["name"], "p": key[0], "q": key[1],
            "n_q": entry["n_q"], "L": len(walk) - 1,
            "total_norm_sq": float(total_norm),
            "weights": weights,
        })

    # ---- Print table -----------------------------------------------------
    sorted_evs = sorted(eig_groups.keys())
    print(f"{'(p,q)':<8} {'name':<10} {'n_q':<4} {'L':<3} {'total':<8}", end='')
    for ev in sorted_evs:
        abs_ev = abs(ev)
        if abs_ev < 1e-6: j_label = "j=0(1)"
        elif abs(abs_ev - 2) < 0.01: j_label = "j=1(3)"
        elif abs(abs_ev - 6) < 0.01: j_label = "j=2(5)"
        else: j_label = f"ev={ev:.2f}"
        print(f"{j_label:<10}", end='')
    print(" dominant")
    print("-" * 110)

    for r in rows:
        # Find dominant projector
        dom_ev = max(r['weights'], key=r['weights'].get)
        abs_ev = abs(dom_ev)
        if abs_ev < 1e-6: dom_label = "j=0"
        elif abs(abs_ev - 2) < 0.01: dom_label = "j=1"
        elif abs(abs_ev - 6) < 0.01: dom_label = "j=2"
        else: dom_label = f"ev={dom_ev:.2f}"
        print(f"({r['p']:>2},{r['q']:>2})  {r['name']:<10} {r['n_q']:<4} "
              f"{r['L']:<3} {r['total_norm_sq']:<8.2f}", end='')
        for ev in sorted_evs:
            w = r['weights'][ev]
            frac = w / r['total_norm_sq'] if r['total_norm_sq'] > 0 else 0
            print(f"{frac:<10.3f}", end='')
        print(f"  {dom_label}")
    print()

    # ---- Cluster by n_q -------------------------------------------------
    print("=" * 78)
    print("Clustering: dominant projector by n_q sector")
    print("=" * 78)
    print()
    by_nq = defaultdict(list)
    for r in rows:
        dom_ev = max(r['weights'], key=r['weights'].get)
        abs_ev = abs(dom_ev)
        if abs_ev < 1e-6: dom_label = "j=0 (1-dim)"
        elif abs(abs_ev - 2) < 0.01: dom_label = "j=1 (3-dim)"
        elif abs(abs_ev - 6) < 0.01: dom_label = "j=2 (5-dim)"
        else: dom_label = f"ev={dom_ev:.2f}"
        by_nq[r['n_q']].append((r['name'], dom_label,
                                  r['weights'], r['total_norm_sq']))

    for nq in sorted(by_nq.keys()):
        print(f"n_q = {nq}:")
        for name, dom_label, weights, total in by_nq[nq]:
            # Show fractional weights
            fracs = {}
            for ev, w in weights.items():
                abs_ev = abs(ev)
                if abs_ev < 1e-6: key = "j=0"
                elif abs(abs_ev - 2) < 0.01: key = "j=1"
                elif abs(abs_ev - 6) < 0.01: key = "j=2"
                else: key = f"ev={ev:.2f}"
                fracs[key] = w / total if total > 0 else 0
            frac_str = ', '.join(f"{k}={v:.2f}" for k, v in
                                  sorted(fracs.items()))
            print(f"  {name:<10} dominant: {dom_label}   ({frac_str})")
        print()

    # ---- Correlation check ----------------------------------------------
    print("=" * 78)
    print("Correlations (fractional weight per j) with n_q")
    print("=" * 78)
    print()
    for ev in sorted_evs:
        abs_ev = abs(ev)
        if abs_ev < 1e-6: j_label = "j=0"
        elif abs(abs_ev - 2) < 0.01: j_label = "j=1"
        elif abs(abs_ev - 6) < 0.01: j_label = "j=2"
        else: j_label = f"ev={ev:.2f}"
        x = np.array([r['weights'][ev] / r['total_norm_sq']
                       if r['total_norm_sq'] > 0 else 0 for r in rows])
        y = np.array([r['n_q'] for r in rows])
        if x.std() > 0 and y.std() > 0:
            corr = np.corrcoef(x, y)[0, 1]
            print(f"  Pearson r(frac_{j_label}, n_q) = {corr:+.4f}")
    print()


if __name__ == "__main__":
    main()
