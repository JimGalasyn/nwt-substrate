"""Q6 — ordering-dependent walk invariants predicting n_q sector.

Q1 + Q4 + Q5 established that σ-sig (direction-blind edge multiplicities)
cannot determine n_q. Walk-reversal preserves σ-sig but inverts the
d-sequence (QR ↔ NR), so direction-blind invariants lose matter/
antimatter sector and likely the carrier-knot sector.

This script tests ordering-dependent invariants:

  1. d-multiset: count of each d ∈ {1, ..., 6} along the walk
  2. QR-pure / NR-pure / mixed flags
  3. d-sequence run patterns (longest run, run-length stats)
  4. σ-orbit transition matrix (7×7)
  5. d-balance: |#QR − #NR|
  6. d-run stats: longest QR-run, longest NR-run

Target: n_q sector ∈ {0, 2, 3, 5} (= unknot, Hopf, trefoil, cinquefoil).

Per [[sector-hamilton-verification]]:
  - nucleon (n_q=5): EXCLUSIVELY d=1
  - meson (n_q=2):  MIXED d=1 + d=3 (no full d=2 Hamilton)
  - hyperon (n_q=3): EXCLUSIVELY d=3 (with some no-Hamilton seeds)
  - lepton (n_q=0):  no Hamilton cycle structure

This script attempts to derive that classification rigorously from
walk-level d-sequence features.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q6_walk_ordering_to_nq.py
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import groupby

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.walks import edge_to_orbit
from nwt_substrate.particles.compendium import COMPENDIUM


QR = {1, 2, 4}
NR = {3, 5, 6}


def d_sequence(walk: list[int]) -> tuple[int, ...]:
    return tuple((walk[i+1] - walk[i]) % 7 for i in range(len(walk) - 1))


def d_multiset(d_seq: tuple[int, ...]) -> tuple[int, ...]:
    counts = [0] * 7
    for d in d_seq:
        counts[d] += 1
    return tuple(counts)  # index 0 unused (d=0 is self-loop, never)


def qr_nr_balance(d_seq: tuple[int, ...]) -> tuple[int, int, int]:
    """Returns (#QR, #NR, |QR−NR|)."""
    n_qr = sum(1 for d in d_seq if d in QR)
    n_nr = sum(1 for d in d_seq if d in NR)
    return n_qr, n_nr, abs(n_qr - n_nr)


def longest_runs(d_seq: tuple[int, ...]) -> dict:
    """Compute longest runs by:
      (a) same d value
      (b) QR class (d ∈ QR)
      (c) NR class (d ∈ NR)
    """
    same_d = max(sum(1 for _ in g) for k, g in groupby(d_seq))
    qr_runs = [sum(1 for _ in g) for k, g in groupby(d_seq, key=lambda d: d in QR) if k]
    nr_runs = [sum(1 for _ in g) for k, g in groupby(d_seq, key=lambda d: d in NR) if k]
    return {
        "same_d": same_d,
        "qr_run": max(qr_runs) if qr_runs else 0,
        "nr_run": max(nr_runs) if nr_runs else 0,
        "n_qr_runs": len(qr_runs),
        "n_nr_runs": len(nr_runs),
    }


def hamilton_d(walk: list[int]) -> int | None:
    """If the walk traverses K_7 as a Hamilton cycle with uniform d-step,
    return that d ∈ {1..6}; else None.
    """
    if len(walk) - 1 != 7:
        return None
    if set(walk[:-1]) != set(range(7)):
        return None
    d0 = (walk[1] - walk[0]) % 7
    for i in range(1, 7):
        if (walk[i+1] - walk[i]) % 7 != d0:
            return None
    return d0


def hamilton_containment(walk: list[int]) -> dict:
    """Does the walk contain a Hamilton-cycle sub-pattern in d ∈ {1, 2, 3}?
    Returns dict {d: max_consecutive_d_run_length}.
    """
    d_seq = d_sequence(walk)
    return {d: max((sum(1 for _ in g) for k, g in groupby(d_seq) if k == d),
                    default=0)
            for d in (1, 2, 3, 4, 5, 6)}


def sigma_signature(walk: list[int]) -> tuple[int, ...]:
    counts = [0] * 7
    for i in range(len(walk) - 1):
        oid = edge_to_orbit(walk[i], walk[i+1])
        if 0 <= oid < 7:
            counts[oid] += 1
    return tuple(counts)


def sigma_trajectory(walk: list[int]) -> tuple[int, ...]:
    return tuple(edge_to_orbit(walk[i], walk[i+1])
                  for i in range(len(walk) - 1))


def main() -> None:
    print("=" * 78)
    print("Q6 — ordering-dependent walk invariants predicting n_q sector")
    print("=" * 78)
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
        d_seq = d_sequence(walk)
        n_qr, n_nr, bal = qr_nr_balance(d_seq)
        runs = longest_runs(d_seq)
        hc = hamilton_containment(walk)
        rows.append({
            "name": entry["name"], "p": key[0], "q": key[1], "L": len(walk) - 1,
            "n_q": entry["n_q"],
            "walk": walk, "d_seq": d_seq,
            "d_mult": d_multiset(d_seq),
            "n_qr": n_qr, "n_nr": n_nr, "qr_nr_bal": bal,
            "hamilton_d": hamilton_d(walk),
            "longest_same_d": runs["same_d"],
            "longest_qr_run": runs["qr_run"],
            "longest_nr_run": runs["nr_run"],
            "n_qr_runs": runs["n_qr_runs"],
            "n_nr_runs": runs["n_nr_runs"],
            "hc_d1": hc[1], "hc_d2": hc[2], "hc_d3": hc[3],
            "hc_d4": hc[4], "hc_d5": hc[5], "hc_d6": hc[6],
            "sig": sigma_signature(walk),
        })

    # -- Per-walk dump ----------------------------------------------------
    print(f"{'(p,q)':<8} {'n_q':<3} {'L':<3} {'d-mult (d=1..6)':<24} "
          f"{'#QR':<4} {'#NR':<4} {'max same-d':<11} "
          f"{'max d=1':<8} {'max d=3':<8}")
    print("-" * 80)
    for r in rows:
        dm = ' '.join(str(c) for c in r['d_mult'][1:])
        print(f"({r['p']:>2},{r['q']:>2})  {r['n_q']:<3} {r['L']:<3} "
              f"[{dm:<20}] {r['n_qr']:<4} {r['n_nr']:<4} "
              f"{r['longest_same_d']:<11} "
              f"{r['hc_d1']:<8} {r['hc_d3']:<8}")
    print()

    # -- Hypothesis: n_q sector classification by d-pattern ---------------
    print("=" * 78)
    print("HYPOTHESIS: n_q sector ← d-direction Hamilton mixing pattern")
    print("=" * 78)
    print()
    print("Rule from [[sector-hamilton-verification]]:")
    print("  nucleon (n_q=5): exclusive d=1 Hamilton  →  hc_d1 ≥ 7, hc_d3 == 0")
    print("  hyperon (n_q=3): exclusive d=3 Hamilton  →  hc_d3 ≥ 7, hc_d1 == 0")
    print("  meson   (n_q=2): MIXED d=1 + d=3         →  hc_d1, hc_d3 both > 0")
    print("  lepton  (n_q=0): no full Hamilton        →  hc_d1, hc_d3 both < 7")
    print()
    print(f"{'(p,q)':<8} {'name':<10} {'n_q sec':<8} {'hc_d1':<6} {'hc_d3':<6} "
          f"{'predicted sector'}")
    print("-" * 70)
    correct = 0
    for r in rows:
        hd1, hd3 = r['hc_d1'], r['hc_d3']
        if hd1 >= 7 and hd3 == 0:
            pred = "n_q=5 (nucleon)"
            pred_nq = 5
        elif hd3 >= 7 and hd1 == 0:
            pred = "n_q=3 (hyperon)"
            pred_nq = 3
        elif hd1 > 0 and hd3 > 0:
            pred = "n_q=2 (meson)"
            pred_nq = 2
        else:
            pred = "n_q=0 (lepton)"
            pred_nq = 0
        ok = pred_nq == r['n_q']
        if ok:
            correct += 1
        marker = "✓" if ok else "✗"
        print(f"({r['p']:>2},{r['q']:>2})  {r['name']:<10} {r['n_q']:<8} "
              f"{hd1:<6} {hd3:<6} {pred:<20} {marker}")
    print(f"\n  Sector classification: {correct}/{len(rows)} correct")
    print()

    # -- Find walks where simple rule fails -------------------------------
    print("=" * 78)
    print("FAILURE ANALYSIS — which walks defy the simple Hamilton rule?")
    print("=" * 78)
    print()
    for r in rows:
        hd1, hd3 = r['hc_d1'], r['hc_d3']
        if hd1 >= 7 and hd3 == 0:
            pred_nq = 5
        elif hd3 >= 7 and hd1 == 0:
            pred_nq = 3
        elif hd1 > 0 and hd3 > 0:
            pred_nq = 2
        else:
            pred_nq = 0
        if pred_nq != r['n_q']:
            print(f"  ✗ {r['name']:<10} (p,q)=({r['p']},{r['q']}) L={r['L']}")
            print(f"       n_q_actual = {r['n_q']}, predicted = {pred_nq}")
            print(f"       d-multiset (d=1..6): {r['d_mult'][1:]}")
            print(f"       hc_d (d=1..6): "
                  f"{[r[f'hc_d{d}'] for d in range(1,7)]}")
            print(f"       walk: {r['walk']}")
            print(f"       σ-trajectory: {sigma_trajectory(r['walk'])}")
            print()

    # -- d-multiset → n_q linear fit -------------------------------------
    print("=" * 78)
    print("REGRESSION: d-multiset vs n_q")
    print("=" * 78)
    print()
    X = np.array([list(r['d_mult'][1:]) for r in rows])  # 6 features
    y = np.array([r['n_q'] for r in rows])
    A = np.hstack([X, np.ones((X.shape[0], 1))])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    print(f"  d-multiset → n_q   R² = {r2:.4f}")
    print(f"  Coefficients: " + ' '.join(f"d{d+1}={beta[d]:+.3f}"
                                          for d in range(6)) +
          f"  intercept={beta[6]:+.3f}")
    print()

    # -- Refined hypothesis: σ-orbit transition matrix? ------------------
    print("=" * 78)
    print("σ-orbit transition matrix per sector")
    print("=" * 78)
    print()
    by_sector = defaultdict(list)
    for r in rows:
        by_sector[r['n_q']].append(r)
    for nq in sorted(by_sector.keys()):
        print(f"  n_q = {nq} sector  ({len(by_sector[nq])} classes)")
        # Build a transition matrix averaged over walks in this sector
        # Tij[i][j] = # of i→j σ-orbit transitions, normalized by walk
        T = np.zeros((7, 7), dtype=int)
        for r in by_sector[nq]:
            traj = sigma_trajectory(r['walk'])
            for k in range(len(traj) - 1):
                if 0 <= traj[k] < 7 and 0 <= traj[k+1] < 7:
                    T[traj[k]][traj[k+1]] += 1
        print(f"    σ-orbit transition matrix (rows=from, cols=to):")
        print("           " + ' '.join(f"o{j}" for j in range(7)))
        for i in range(7):
            print(f"      o{i}:  " + ' '.join(f"{T[i,j]:>2}" for j in range(7)))
        print()


if __name__ == "__main__":
    main()
