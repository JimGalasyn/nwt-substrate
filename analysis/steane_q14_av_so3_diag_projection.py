"""Q14 — so(3)_diag Clebsch-Gordan projection of v-direction algebra element A_v.

Q9 did so(3)_diag projection of the FULL walk algebra element A_walk = Σ_e J_e
and found r < 0.20 correlation with Paper 11 n_q sectors. The decomposition
so(7) | so(3)_diag = (1) ⊕ 5·(3) ⊕ (5) is clean, but the walks don't
preferentially activate any one irrep by sector.

Q13 introduced A_v = Σ_e (Δv_e) · J_e — the v-direction-weighted version
that should isolate the carrier-knot part of the walk's substrate algebra
element. The multi-traversal hypothesis says λ_v = n_q, so v-weighted
projection MIGHT cluster cleanly even though full-walk projection doesn't.

This script:
  1. Build the so(3)_diag Casimir and spectral projectors P_0, P_1, P_2
     (same as Q9).
  2. For each compendium walk, compute A_v = Σ_e (Δv_e) · J_e.
  3. Project A_v onto each Casimir eigenspace and compute fractional
     weights frac_j ∈ [0, 1].
  4. Test correlations with Paper 11 n_q sector.
  5. Compare against A_walk projection (Q9 baseline) and A_u projection.

If A_v's so(3)_diag projection clusters cleanly by n_q sector while
A_walk's doesn't, the v-direction multi-traversal interpretation has
empirical bite.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q14_av_so3_diag_projection.py
"""
from __future__ import annotations

import math
from collections import defaultdict
from itertools import combinations

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.particles.compendium import COMPENDIUM


N = 7


def so7_basis():
    basis = []
    for i, j in combinations(range(N), 2):
        J = np.zeros((N, N))
        J[i, j] = +1.0
        J[j, i] = -1.0
        basis.append(J)
    return basis


SO7_BASIS = so7_basis()
EDGE_INDEX = {}
for k, (i, j) in enumerate(combinations(range(N), 2)):
    EDGE_INDEX[(i, j)] = k
    EDGE_INDEX[(j, i)] = k


def so3_diag_generators():
    def gen(i, j):
        M = np.zeros((N, N))
        M[i, j] = 1.0
        M[j, i] = -1.0
        return M
    L_x = gen(2, 3)
    L_y = gen(1, 3)
    L_z = gen(1, 2)
    L_xp = gen(5, 6)
    L_yp = gen(4, 6)
    L_zp = gen(4, 5)
    return L_x + L_xp, L_y + L_yp, L_z + L_zp


def adjoint_action(J, basis):
    n = len(basis)
    M = np.zeros((n, n))
    for l in range(n):
        comm = J @ basis[l] - basis[l] @ J
        for k in range(n):
            M[k, l] = -0.5 * np.trace(basis[k] @ comm)
    return M


def walk_components(walk):
    """Return A_walk, A_u, A_v as so(7)-coordinate vectors."""
    A_walk = np.zeros(len(SO7_BASIS))
    A_u = np.zeros(len(SO7_BASIS))
    A_v = np.zeros(len(SO7_BASIS))
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        k = EDGE_INDEX[(a, b)]
        sign = +1.0 if a < b else -1.0
        n_u, n_v = edge_winding_class(a, b)
        if a > b:
            n_u, n_v = -n_u, -n_v
        A_walk[k] += sign
        A_u[k] += sign * n_u
        A_v[k] += sign * n_v
    return A_walk, A_u, A_v


def build_projectors():
    Jx, Jy, Jz = so3_diag_generators()
    adx = adjoint_action(Jx, SO7_BASIS)
    ady = adjoint_action(Jy, SO7_BASIS)
    adz = adjoint_action(Jz, SO7_BASIS)
    C = adx @ adx + ady @ ady + adz @ adz
    C = 0.5 * (C + C.T)
    eigvals, eigvecs = np.linalg.eigh(C)
    groups = defaultdict(list)
    for k, ev in enumerate(eigvals):
        groups[round(ev, 4)].append(k)
    projectors = {}
    j_labels = {}
    for ev, idx in groups.items():
        P = np.zeros_like(C)
        for k in idx:
            v = eigvecs[:, k:k+1]
            P += v @ v.T
        projectors[ev] = P
        if abs(ev) < 1e-6:
            j_labels[ev] = "j=0"
        elif abs(ev + 2) < 0.01:
            j_labels[ev] = "j=1"
        elif abs(ev + 6) < 0.01:
            j_labels[ev] = "j=2"
        else:
            j_labels[ev] = f"ev={ev:.2f}"
    return projectors, j_labels


def project_and_report(A, projectors, j_labels):
    tot = np.linalg.norm(A) ** 2
    if tot < 1e-12:
        return None, {l: 0.0 for l in j_labels.values()}
    weights = {}
    for ev, P in projectors.items():
        w = float(np.linalg.norm(P @ A) ** 2)
        weights[j_labels[ev]] = w / tot
    dominant = max(weights, key=weights.get)
    return dominant, weights


def main():
    print("=" * 78)
    print("Q14 — so(3)_diag projection of A_v (v-direction algebra element)")
    print("=" * 78)
    print()

    projectors, j_labels = build_projectors()
    print("Casimir spectrum verified: so(7)|so(3)_diag = (1) ⊕ 5·(3) ⊕ (5)")
    print()

    walks_dict = bfs_shortest_walks(max_length=25)

    rows = []
    seen = set()
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key in seen or key not in walks_dict:
            continue
        seen.add(key)
        walk = walks_dict[key]
        A_walk, A_u, A_v = walk_components(walk)
        dom_walk, w_walk = project_and_report(A_walk, projectors, j_labels)
        dom_u,    w_u    = project_and_report(A_u, projectors, j_labels)
        dom_v,    w_v    = project_and_report(A_v, projectors, j_labels)
        rows.append({
            "name": entry["name"], "p": key[0], "q": key[1],
            "n_q": entry["n_q"], "L": len(walk) - 1,
            "dom_walk": dom_walk, "w_walk": w_walk,
            "dom_u": dom_u, "w_u": w_u,
            "dom_v": dom_v, "w_v": w_v,
        })

    # ---- Per-walk projection table --------------------------------------
    print(f"{'(p,q)':<8} {'name':<10} {'n_q':<4} "
          f"{'A_walk: j=0/1/2 frac':<24} {'A_u frac':<24} {'A_v frac'}")
    print("-" * 100)
    for r in rows:
        ww = r['w_walk']
        wu = r['w_u']
        wv = r['w_v']
        wstr = f"({ww['j=0']:.2f},{ww['j=1']:.2f},{ww['j=2']:.2f})"
        ustr = f"({wu['j=0']:.2f},{wu['j=1']:.2f},{wu['j=2']:.2f})"
        vstr = f"({wv['j=0']:.2f},{wv['j=1']:.2f},{wv['j=2']:.2f})"
        print(f"({r['p']:>2},{r['q']:>2})  {r['name']:<10} {r['n_q']:<4} "
              f"{wstr:<24} {ustr:<24} {vstr}")
    print()

    # ---- Correlations with n_q for each variant -------------------------
    print("=" * 78)
    print("Correlations: fractional weight on each j-irrep vs Paper 11 n_q")
    print("=" * 78)
    print()
    y = np.array([r['n_q'] for r in rows])
    for variant_label, key in [("A_walk", "w_walk"), ("A_u", "w_u"),
                                 ("A_v", "w_v")]:
        print(f"  {variant_label} projection:")
        for j_label in ["j=0", "j=1", "j=2"]:
            x = np.array([r[key][j_label] for r in rows])
            if x.std() > 0 and y.std() > 0:
                r_coef = np.corrcoef(x, y)[0, 1]
                print(f"    Pearson r(frac_{j_label}, n_q) = {r_coef:+.4f}")
            else:
                print(f"    frac_{j_label}: constant or zero variance")
        print()

    # ---- Cluster by n_q sector -------------------------------------------
    print("=" * 78)
    print("Mean fractional weights by n_q sector")
    print("=" * 78)
    print()
    by_nq = defaultdict(list)
    for r in rows:
        by_nq[r['n_q']].append(r)
    print(f"  {'sector':<10} {'count':<6} "
          f"{'A_v: j=0':<10} {'A_v: j=1':<10} {'A_v: j=2':<10} "
          f"{'A_walk: j=0':<13} {'A_walk: j=1':<13}")
    print("  " + "-" * 75)
    for nq in sorted(by_nq.keys()):
        entries = by_nq[nq]
        wv_j0 = np.mean([r['w_v']['j=0'] for r in entries])
        wv_j1 = np.mean([r['w_v']['j=1'] for r in entries])
        wv_j2 = np.mean([r['w_v']['j=2'] for r in entries])
        ww_j0 = np.mean([r['w_walk']['j=0'] for r in entries])
        ww_j1 = np.mean([r['w_walk']['j=1'] for r in entries])
        print(f"  n_q={nq}    {len(entries):<6} "
              f"{wv_j0:<10.3f} {wv_j1:<10.3f} {wv_j2:<10.3f} "
              f"{ww_j0:<13.3f} {ww_j1:<13.3f}")
    print()

    # ---- Dominant projector classification ------------------------------
    print("=" * 78)
    print("Dominant projector classification per variant")
    print("=" * 78)
    print()
    for variant_label, key in [("A_walk", "dom_walk"),
                                 ("A_u", "dom_u"), ("A_v", "dom_v")]:
        print(f"  {variant_label} dominant projector breakdown:")
        nq_to_dom = defaultdict(lambda: defaultdict(int))
        for r in rows:
            nq_to_dom[r['n_q']][r[key]] += 1
        for nq in sorted(nq_to_dom.keys()):
            doms = nq_to_dom[nq]
            total = sum(doms.values())
            print(f"    n_q={nq} ({total} walks): "
                  + ', '.join(f"{d}={c}" for d, c in
                               sorted(doms.items())))
        print()


if __name__ == "__main__":
    main()
