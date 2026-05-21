"""Q13 — v-cycle Wilson eigenvalues vs n_q sectors.

Multi-traversal hypothesis (Jim's intuition for why exponent = q in
n_q^q): each v-revolution of a (p, q) walk on the Heffter torus
traces the carrier-knot once; q v-revolutions accumulate the
carrier-knot factor q times, giving n_q^q.

Substrate framing:
  - Heffter torus T² has u-cycle (longitudinal, vertex-cycling) and
    v-cycle (latitudinal, twisted by ×3).
  - A walk's Wilson holonomy decomposes as
        W_(p,q) = (u-holonomy)^p · (v-holonomy)^q
  - Mass formula `(p²+q²)/5 · β-factor · n_q^q` puts the carrier-knot
    factor in the v-position. So the hypothesis is:
        λ_v = n_q
  - Each compendium walk's v-direction Wilson loop should have
    eigenvalue magnitude ~n_q (sector-determined).

This script tests:
  1. Compute A_v = Σ_e (Δv_e) · J_e (v-direction-weighted so(7) sum)
     where Δv_e is the v-displacement of edge e in 1/7 units.
  2. Compute W_v = exp(A_v · t) for various coupling t.
  3. Compute observables (trace, log|trace|, max eigenvalue phase,
     determinant subspaces) and test whether they cluster at n_q
     sector values {1, 2, 3, 5}.
  4. Compare against u-direction A_u and total A_walk.

If v-direction observables cluster by n_q sector and u-direction
observables don't, the multi-traversal hypothesis is supported. If
neither separates, the n_q^q form has a more subtle origin.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q13_v_cycle_wilson.py
"""
from __future__ import annotations

import math
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.linalg import expm

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


def walk_components(walk: list[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return three so(7)-coordinate vectors per walk:
      A_walk: Σ_e (sign_e) · J_e (each edge counted with directional sign)
      A_u:    Σ_e (sign_e · Δu_e) · J_e (u-direction-weighted)
      A_v:    Σ_e (sign_e · Δv_e) · J_e (v-direction-weighted)

    Δu_e, Δv_e are the integer winding components (in 1/7 units, i.e.,
    n_u, n_v ∈ {-3..3} for each edge).
    """
    A_walk = np.zeros(len(SO7_BASIS))
    A_u = np.zeros(len(SO7_BASIS))
    A_v = np.zeros(len(SO7_BASIS))
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        k = EDGE_INDEX[(a, b)]
        sign = +1.0 if a < b else -1.0
        n_u, n_v = edge_winding_class(a, b)
        # Note: edge_winding_class returns the directed winding from a→b
        # not from min(a,b)→max(a,b), so we need to flip sign if needed
        if a > b:
            n_u, n_v = -n_u, -n_v
            sign = -sign  # cancels: directed J with directed winding
            # Actually: J_{a,b} with a < b is the basis element; the
            # walk's contribution is sign · n_u · J_{i,j}. Let's just
            # use the directed J + directed n_u/n_v consistently.
        A_walk[k] += sign
        A_u[k] += sign * n_u
        A_v[k] += sign * n_v
    return A_walk, A_u, A_v


def coords_to_matrix(coords: np.ndarray) -> np.ndarray:
    """Build the 7×7 so(7) matrix from basis coordinates."""
    M = np.zeros((N, N))
    for k, c in enumerate(coords):
        if c != 0:
            M += c * SO7_BASIS[k]
    return M


def wilson_observables(A: np.ndarray, t: float) -> dict:
    """Build W = exp(A · t) and extract observables."""
    A_mat = coords_to_matrix(A)
    W = expm(A_mat * t)
    eigs = np.linalg.eigvals(W)
    # SO(N) eigenvalues are on the unit circle, real
    trace = float(np.trace(W).real)
    # Max angle (= max |arg(eig)|)
    angles = np.angle(eigs)
    max_angle = float(max(abs(a) for a in angles))
    # Sum of squared angles (~ trace of A² · t², a Casimir proxy)
    sum_angle_sq = float(sum(a * a for a in angles))
    # Norm of A coordinates (substrate-level "size")
    A_norm = float(np.linalg.norm(A))
    A_sq_norm = float(np.sum(A * A))
    return {
        "trace": trace,
        "max_angle": max_angle,
        "sum_angle_sq": sum_angle_sq,
        "A_norm": A_norm,
        "A_sq_norm": A_sq_norm,
    }


def main():
    print("=" * 78)
    print("Q13 — v-cycle Wilson eigenvalues vs n_q sectors")
    print("=" * 78)
    print()
    print("Hypothesis: λ_v = n_q, so v-direction Wilson loop / A_v Casimir")
    print("            should cluster at {1, 2, 3, 5} by Paper 11 n_q sector.")
    print()

    walks = bfs_shortest_walks(max_length=25)

    rows = []
    seen = set()
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key in seen or key not in walks:
            continue
        seen.add(key)
        walk = walks[key]
        A_walk, A_u, A_v = walk_components(walk)
        rows.append({
            "name": entry["name"], "p": key[0], "q": key[1],
            "n_q": entry["n_q"], "L": len(walk) - 1,
            "A_walk": A_walk, "A_u": A_u, "A_v": A_v,
        })

    # ---- Per-walk A-norm table ------------------------------------------
    print(f"{'(p,q)':<8} {'name':<10} {'n_q':<4} {'L':<3} "
          f"{'|A_walk|':<10} {'|A_u|':<10} {'|A_v|':<10} {'|A_v|/q':<10}")
    print("-" * 80)
    for r in rows:
        nw = np.linalg.norm(r['A_walk'])
        nu = np.linalg.norm(r['A_u'])
        nv = np.linalg.norm(r['A_v'])
        nvq = nv / r['q'] if r['q'] > 0 else 0
        print(f"({r['p']:>2},{r['q']:>2})  {r['name']:<10} {r['n_q']:<4} "
              f"{r['L']:<3} {nw:<10.3f} {nu:<10.3f} {nv:<10.3f} {nvq:<10.3f}")
    print()

    # ---- Cluster |A_v|/q by n_q sector ----------------------------------
    print("=" * 78)
    print("|A_v| / q  clustered by Paper 11 n_q sector")
    print("=" * 78)
    print()
    print("If λ_v = n_q multi-traversal hypothesis holds, |A_v|/q should")
    print("cluster at sector-determined values, NOT necessarily equal to n_q")
    print("itself but distinguishing the four sectors.")
    print()
    by_nq = defaultdict(list)
    for r in rows:
        nv = np.linalg.norm(r['A_v'])
        by_nq[r['n_q']].append((r['name'], nv / r['q'] if r['q'] > 0 else 0))
    for nq in sorted(by_nq.keys()):
        vals = [v for _, v in by_nq[nq]]
        mean_v = np.mean(vals)
        std_v = np.std(vals)
        print(f"  n_q = {nq}: mean(|A_v|/q) = {mean_v:.3f} ± {std_v:.3f}")
        for name, v in by_nq[nq]:
            print(f"    {name:<10} |A_v|/q = {v:.3f}")
        print()

    # ---- Wilson loop observables at various t ---------------------------
    print("=" * 78)
    print("Wilson loop observables W_v = exp(A_v · t) clustered by n_q")
    print("=" * 78)
    print()
    for t_label, t in [("α", 0.0073), ("π/L_avg", math.pi / 13),
                        ("1.0", 1.0)]:
        print(f"--- t = {t_label} ({t:.4f}) ---")
        print(f"  {'sector':<10} "
              f"{'avg trace':<12} {'avg max_angle':<15} "
              f"{'avg sum_angle²':<15}")
        sector_obs = defaultdict(list)
        for r in rows:
            obs = wilson_observables(r['A_v'], t)
            sector_obs[r['n_q']].append(obs)
        for nq in sorted(sector_obs.keys()):
            traces = [o['trace'] for o in sector_obs[nq]]
            angles = [o['max_angle'] for o in sector_obs[nq]]
            sums = [o['sum_angle_sq'] for o in sector_obs[nq]]
            print(f"  n_q={nq}      "
                  f"{np.mean(traces):<12.3f} {np.mean(angles):<15.3f} "
                  f"{np.mean(sums):<15.3f}")
        print()

    # ---- log of substrate observables vs q · log(n_q) -------------------
    print("=" * 78)
    print("Test: which observable f satisfies log(f) ≈ q · log(n_q) ?")
    print("=" * 78)
    print()
    print("Multi-traversal hypothesis predicts an observable with this scaling.")
    print()
    log_target = np.array([r['q'] * math.log(max(r['n_q'], 1)) for r in rows])
    print(f"  Target: q · log(n_q) (= log of Paper 6's R_paper6 = n_q^q)")
    print(f"  Range: [{log_target.min():.3f}, {log_target.max():.3f}]")
    print()

    candidates = [
        ("|A_walk|²",   lambda r: np.sum(r['A_walk']**2)),
        ("|A_u|²",      lambda r: np.sum(r['A_u']**2)),
        ("|A_v|²",      lambda r: np.sum(r['A_v']**2)),
        ("|A_v|²/q²",   lambda r: np.sum(r['A_v']**2) / max(r['q']**2, 1)),
        ("|A_walk|·q",  lambda r: np.linalg.norm(r['A_walk']) * r['q']),
        ("|A_v|·q",     lambda r: np.linalg.norm(r['A_v']) * r['q']),
        ("L",           lambda r: r['L']),
        ("L · q",       lambda r: r['L'] * r['q']),
    ]

    print(f"  {'candidate':<14} {'Pearson r':<11} {'best-fit slope':<16} "
          f"{'RMS residual'}")
    print("  " + "-" * 70)
    for name, fn in candidates:
        vals = np.array([math.log(max(fn(r), 1e-10)) for r in rows])
        if vals.std() < 1e-10 or log_target.std() < 1e-10:
            print(f"  {name:<14} (constant)")
            continue
        r_coef = np.corrcoef(vals, log_target)[0, 1]
        # Linear fit log(f) ≈ a · q·log(n_q) + b
        X = np.column_stack([log_target, np.ones_like(log_target)])
        coefs, *_ = np.linalg.lstsq(X, vals, rcond=None)
        slope, intercept = coefs
        resid = vals - (slope * log_target + intercept)
        rms = float(np.sqrt(np.mean(resid**2)))
        print(f"  {name:<14} {r_coef:<+11.4f} {slope:<+16.4f} {rms:.4f}")
    print()

    # ---- u vs v direction: which one carries the n_q signal? ------------
    print("=" * 78)
    print("u-direction vs v-direction Wilson observables — which carries n_q?")
    print("=" * 78)
    print()
    for direction, key in [("u-direction", "A_u"), ("v-direction", "A_v")]:
        # Compute |A|² normalized by L (per-step Casimir) and check
        # if it clusters by n_q sector.
        print(f"  {direction} (A_norm²/L by sector):")
        sec = defaultdict(list)
        for r in rows:
            val = np.sum(r[key]**2) / max(r['L'], 1)
            sec[r['n_q']].append(val)
        for nq in sorted(sec.keys()):
            vals = sec[nq]
            print(f"    n_q={nq}: mean {np.mean(vals):.3f} ± {np.std(vals):.3f}  "
                  f"(values: {[f'{v:.2f}' for v in vals]})")
        print()


if __name__ == "__main__":
    main()
