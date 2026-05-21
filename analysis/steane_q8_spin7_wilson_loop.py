"""Q8 — Spin(7) Wilson-loop pipeline Steps 3+4.

Continuation of [[spin7-rep-theory-carrier-knot-framework]]. Steps 1+2
were prototyped: carrier-knot {unknot, Hopf, trefoil, cinquefoil}
n_q ∈ {0, 2, 3, 5} → so(3)_diag piece dims {1, 2, 3, 5} of so(7)
sub-algebra decomposition; explicit σ-orbit → so(7) generator map
(21 generators / 7 σ-orbits).

This script implements Steps 3 and 4:

  Step 3: For each compendium walk W = v_0 → ... → v_L → v_0, construct
          the Spin(7) Wilson-loop product
              G_W = exp(J_{v_{L-1}, v_L} t) ... exp(J_{v_1, v_2} t) exp(J_{v_0, v_1} t)
          where J_{i,j} = E_{ij} - E_{ji} is the so(7) skew-symmetric
          generator and t is the coupling (sweep t = {α, √α, 1}).

  Step 4: Compute G_W's eigenvalue spectrum. Classify by:
            - rank(G_W - I) = dimension of non-trivial eigen-subspace
            - eigenvalue-phase multiplicity structure
            - trace Tr(G_W) and characteristic-polynomial signature
          Test whether the classification matches Paper 11 n_q sectors.

Honest framing: Step 4 has THREE open theoretical choices (per the
framework memo): coupling t, path-ordering convention, knot-from-G_W
identification. This script tries multiple combinations and reports
what the data say.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q8_spin7_wilson_loop.py
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.linalg import expm

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.isa.constants import ALPHA_NWT


N = 7  # K_7 vertices


def so7_generator(i: int, j: int) -> np.ndarray:
    """J_{ij} = E_{ij} - E_{ji}, skew-symmetric 7×7 matrix.
    For directed edge a→b, use J_{a, b} (sign depends on direction)."""
    J = np.zeros((N, N), dtype=float)
    J[i, j] = 1.0
    J[j, i] = -1.0
    return J


def wilson_loop(walk: list[int], t: float, ordering: str = "left") -> np.ndarray:
    """Build the path-ordered Wilson-loop product.

    Conventions tried:
      'left'    : G = exp(J_L t) ... exp(J_2 t) exp(J_1 t)
                  (path-ordered from start; multiplied left-to-right in time)
      'right'   : G = exp(J_1 t) exp(J_2 t) ... exp(J_L t)
      'sum'     : leading-order term G = exp(t · Σ J_e) — non-path-ordered
                  (the t → 0 limit of the path-ordered product, where
                  commutators don't matter)
    """
    G = np.eye(N)
    edges = [(walk[i], walk[i + 1]) for i in range(len(walk) - 1)]
    if ordering == "sum":
        A = np.zeros((N, N), dtype=float)
        for a, b in edges:
            A[a, b] += 1.0
            A[b, a] += -1.0
        return expm(A * t)
    if ordering == "cross_block_sum":
        # Restrict to the 6-dim cross-block ⊥ polar vertex 0.
        # Drop edges that touch the polar vertex; keep only intra-equatorial.
        A = np.zeros((N, N), dtype=float)
        for a, b in edges:
            if a == 0 or b == 0:
                continue  # skip polar bridging
            A[a, b] += 1.0
            A[b, a] += -1.0
        # Embed in N×N (polar block fixed; cross-block dynamic)
        return expm(A * t)
    for a, b in edges:
        # Directed generator: J with E_{ab} = +1, E_{ba} = -1
        J = np.zeros((N, N), dtype=float)
        J[a, b] = +1.0
        J[b, a] = -1.0
        U = expm(J * t)
        if ordering == "left":
            G = U @ G
        elif ordering == "right":
            G = G @ U
        else:
            raise ValueError(ordering)
    return G


def classify_eigenvalues(G: np.ndarray, tol: float = 1e-6) -> dict:
    """Compute eigenvalue spectrum and structural classifiers."""
    eigs = np.linalg.eigvals(G)
    # Sort by angle
    angles = np.angle(eigs)
    angles_sorted = sorted(angles, key=lambda a: (round(abs(a), 4), a))
    abs_angles = sorted(abs(a) for a in angles)

    # Count eigenvalues near +1 (= angle 0)
    n_one = sum(1 for a in angles if abs(a) < tol)
    n_minus_one = sum(1 for a in angles if abs(abs(a) - np.pi) < tol)
    n_complex = N - n_one - n_minus_one

    # Distinct angle multiset (treating ±θ as one pair)
    pos_angles = sorted(set(round(abs(a), 4) for a in angles
                             if abs(a) > tol and abs(abs(a) - np.pi) > tol))
    # Multiplicities
    angle_mults = {}
    for a_round in pos_angles:
        angle_mults[a_round] = sum(1 for a in angles
                                     if abs(round(abs(a), 4) - a_round) < tol)

    # rank(G - I): dimension of non-trivial subspace = 2 * # nontrivial pairs
    nontrivial_dim = N - n_one

    # Casimir-like invariants on A = -i log(G) (Lie-algebra element)
    # For so(7), C_2 = -(1/2) Tr(A²) and C_4 = -(1/4) Tr(A⁴), etc.
    # Use the algebra-element directly when ordering='sum' (A = Σ J_e):
    A_real = np.zeros((N, N))  # placeholder; caller can pass A if needed
    # We can recover A approximately via:  A = (G - G.T) / (2 * t) at small t
    # but for now compute C_2 from eigenvalue angles (sum of squared angles)
    casimir_2 = float(sum(a * a for a in angles))
    casimir_4 = float(sum(a ** 4 for a in angles))
    max_angle = float(max(abs(a) for a in angles))

    return {
        "n_eig_one": n_one,
        "n_eig_neg_one": n_minus_one,
        "n_eig_complex": n_complex,
        "nontrivial_dim": nontrivial_dim,
        "n_distinct_angle_pairs": len(pos_angles),
        "angle_mults": angle_mults,
        "trace": float(np.trace(G).real),
        "angles_sorted": [round(a, 4) for a in angles_sorted],
        "casimir_2": casimir_2,
        "casimir_4": casimir_4,
        "max_angle": max_angle,
    }


def main():
    print("=" * 78)
    print("Q8 — Spin(7) Wilson-loop pipeline (Steps 3 + 4)")
    print("=" * 78)
    print()
    print(f"ALPHA_NWT = {ALPHA_NWT:.10f}")
    print()

    walks = bfs_shortest_walks(max_length=25)

    # Diagnostic table: one walk per (|p|,|q|) class, n_q from compendium
    seen = set()
    walks_for_test = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key in seen or key not in walks:
            continue
        seen.add(key)
        walks_for_test.append({
            "name": entry["name"], "p": key[0], "q": key[1],
            "n_q": entry["n_q"], "walk": walks[key],
        })

    # Sweep couplings: small (α), intermediate, large
    couplings = [
        ("α",       ALPHA_NWT),
        ("√α",      np.sqrt(ALPHA_NWT)),
        ("α/π",     ALPHA_NWT / np.pi),
        ("0.5",     0.5),
        ("1.0",     1.0),
        ("π/L",     None),  # walk-length-dependent
    ]

    # Also sweep ordering conventions
    orderings = ["left", "sum", "cross_block_sum"]

    for cname, t_val in couplings:
     for ordering in orderings:
        print("=" * 78)
        print(f"Coupling t = {cname}  | ordering = {ordering}")
        print("=" * 78)
        print()
        print(f"{'(p,q)':<8} {'name':<10} {'n_q':<4} {'L':<3} "
              f"{'rank(G-I)':<11} {'#one':<5} {'#cplx pairs':<12} "
              f"{'trace':<8} {'angle mults'}")
        print("-" * 110)
        rows_at_t = []
        for d in walks_for_test:
            t = (np.pi / (len(d['walk']) - 1)) if t_val is None else t_val
            G = wilson_loop(d['walk'], t, ordering=ordering)
            cl = classify_eigenvalues(G)
            mults_str = ', '.join(f"{a:.3f}×{m}" for a, m in cl['angle_mults'].items())
            print(f"({d['p']:>2},{d['q']:>2})  {d['name']:<10} {d['n_q']:<4} "
                  f"{len(d['walk'])-1:<3} {cl['nontrivial_dim']:<11} "
                  f"{cl['n_eig_one']:<5} {cl['n_distinct_angle_pairs']:<12} "
                  f"{cl['trace']:<8.3f} {mults_str}")
            rows_at_t.append({**d, **cl, "t": t})
        print()

        # Test classifications against n_q across multiple observables
        observables = ["nontrivial_dim", "trace", "casimir_2",
                        "casimir_4", "max_angle", "n_distinct_angle_pairs"]
        y = np.array([r['n_q'] for r in rows_at_t])
        for obs in observables:
            x = np.array([r[obs] for r in rows_at_t])
            if x.std() > 0 and y.std() > 0:
                corr = np.corrcoef(x, y)[0, 1]
                print(f"    Pearson r({obs:<24}, n_q) = {corr:+.4f}")
            else:
                print(f"    {obs}: constant (cannot correlate)")
        print()


if __name__ == "__main__":
    main()
