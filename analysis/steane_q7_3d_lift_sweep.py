"""Q7 — 3D-lift parameter exploration for carrier-knot identification.

Q5/Q6 left an open question: pyknotid identifies most compendium walks
as unknot (0_1), in conflict with Paper 11's n_q sector labels (0, 2, 3,
5). Two possible diagnoses:

  (A) The 3D Heffter-torus lift at R=2.5, r=1, n_per_edge=30 is
      insufficient to manifest the carrier knot — try different R/r,
      finer sampling, higher torus genus.
  (B) Paper 11's n_q label is NOT a torus-knot crossing number — the
      mathematical fact is that a (p, q) torus knot with min(|p|,|q|)=1
      is the unknot (e.g. (1, 3) walks are *mathematically* unknots,
      so pyknotid is correct).

This script tests both by:
  1. Computing the mathematical (p, q)-torus knot min_crossings for
     each compendium walk: tk = min(|p|(|q|-1), |q|(|p|-1)) when
     gcd(|p|, |q|) = 1.
  2. Sweeping (R, r, n_per_edge) across a small grid and running
     pyknotid for diagnostic walks one per sector.
  3. Comparing three columns: pyknotid_observed, tk_predicted,
     paper11_n_q.

If pyknotid ≈ tk → 3D lift is correct, Paper 11 uses a different
convention.
If pyknotid ≠ tk → 3D lift has issues that parameter tweaks might fix.
If pyknotid ≈ paper11_n_q at some parameter set → Paper 11 IS the
carrier knot, just need the right parameters.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q7_3d_lift_sweep.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# numpy compat patches for pyknotid
np.float = float
np.int = int
np.complex = complex
np.object = object
np.str = str
np.long = int

from pyknotid.spacecurves import Knot

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.orbit_winding import (
    edge_winding_class, HEFFTER_VERT_UV,
)
from nwt_substrate.particles.compendium import COMPENDIUM


def walk_to_3d_curve(walk: list[int], n_per_edge: int,
                     R: float, r: float) -> np.ndarray:
    pts = []
    u_curr = HEFFTER_VERT_UV[walk[0]][0]
    v_curr = HEFFTER_VERT_UV[walk[0]][1]
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        nu, nv = edge_winding_class(a, b)
        du = nu / 7.0
        dv = nv / 7.0
        ts = np.linspace(0, 1, n_per_edge, endpoint=False)
        for t in ts:
            u = u_curr + t * du
            v = v_curr + t * dv
            x = (R + r * np.cos(2 * np.pi * v)) * np.cos(2 * np.pi * u)
            y = (R + r * np.cos(2 * np.pi * v)) * np.sin(2 * np.pi * u)
            z = r * np.sin(2 * np.pi * v)
            pts.append((x, y, z))
        u_curr += du
        v_curr += dv
    return np.array(pts)


def identify_knot(curve: np.ndarray) -> tuple[str, int | None]:
    try:
        k = Knot(curve, verbose=False)
        cands = k.identify()
        if not cands:
            return "?", None
        best = min(cands, key=lambda c: (c.min_crossings
                    if hasattr(c, 'min_crossings') and c.min_crossings is not None
                    else 999))
        return best.identifier, best.min_crossings
    except Exception as e:
        return f"ERR:{type(e).__name__}", None


def torus_knot_min_crossings(p: int, q: int) -> int | None:
    """For a (p, q) torus knot with gcd(p, q) = 1, the minimum crossing
    number is min(|p|(|q|-1), |q|(|p|-1)).
    Returns None if gcd ≠ 1 (then it's a torus link, not a knot).
    """
    p, q = abs(p), abs(q)
    if p == 0 or q == 0:
        return 0
    if math.gcd(p, q) != 1:
        return None
    return min(p * (q - 1), q * (p - 1))


def main():
    print("=" * 78)
    print("Q7 — 3D-lift parameter exploration")
    print("=" * 78)
    print()

    # ---- Predict torus-knot crossings for ALL compendium classes -------
    walks = bfs_shortest_walks(max_length=25)

    print(f"{'particle':<10} {'(p,q)':<8} {'gcd':<3} {'L':<3} "
          f"{'paper11 n_q':<13} {'(p,q)-torus tk':<14}")
    print("-" * 64)
    seen = set()
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key in seen:
            continue
        seen.add(key)
        tk = torus_knot_min_crossings(key[0], key[1])
        tk_str = "LINK" if tk is None else f"{tk}"
        walk = walks.get(key)
        L = len(walk) - 1 if walk else "?"
        gcd_pq = math.gcd(key[0], key[1])
        print(f"{entry['name']:<10} ({key[0]:>2},{key[1]:>2})  {gcd_pq:<3} "
              f"{L:<3} {entry['n_q']:<13} {tk_str:<14}")
    print()

    # ---- Diagnostic walks: one per sector ----------------------------------
    print("=" * 78)
    print("Parameter sweep on diagnostic walks (one per sector)")
    print("=" * 78)
    print()

    diagnostics = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if entry["name"] in ("e-", "p", "pi+", "tau-", "K0", "Upsilon"):
            if key in walks:
                diagnostics.append({"name": entry["name"], "key": key,
                                     "walk": walks[key], "n_q": entry["n_q"],
                                     "tk": torus_knot_min_crossings(*key)})

    # Sweep parameter grid: (R, r) with various aspect ratios
    param_grid = [
        (2.5, 1.0, 30),     # default
        (2.5, 1.0, 100),    # finer sampling
        (4.0, 1.0, 30),     # less curvy torus
        (10.0, 1.0, 30),    # near-cylindrical
        (2.5, 2.0, 30),     # fatter (still embedded since r < R)
        (3.0, 0.5, 30),     # thin torus
        (3.0, 1.5, 60),     # mid-fat, finer sampling
    ]

    print(f"{'name':<8} {'(p,q)':<8} {'L':<3} {'expected':<22}", end='')
    for R, r, n in param_grid:
        print(f"R={R} r={r} n={n}".ljust(15), end='')
    print()
    print("-" * 130)

    for d in diagnostics:
        expected = f"n_q={d['n_q']}, tk={d['tk']}"
        print(f"{d['name']:<8} ({d['key'][0]:>2},{d['key'][1]:>2})  "
              f"{len(d['walk'])-1:<3} {expected:<22}", end='')
        for R, r, n in param_grid:
            curve = walk_to_3d_curve(d['walk'], n, R, r)
            kid, nc = identify_knot(curve)
            nc_str = "?" if nc is None else str(nc)
            print(f"{kid}({nc_str})".ljust(15), end='')
        print()
    print()

    # ---- Match table: which parameter regime matches Paper 11 vs tk? -------
    print("=" * 78)
    print("Match summary: pyknotid match vs Paper 11 n_q AND vs (p,q)-tk")
    print("=" * 78)
    print()

    print(f"{'parameters':<22} {'matches paper11':<18} "
          f"{'matches (p,q)-tk':<18}")
    print("-" * 60)
    for R, r, n in param_grid:
        n_paper = 0
        n_tk = 0
        for d in diagnostics:
            curve = walk_to_3d_curve(d['walk'], n, R, r)
            kid, nc = identify_knot(curve)
            if nc == d['n_q']:
                n_paper += 1
            if d['tk'] is not None and nc == d['tk']:
                n_tk += 1
        print(f"R={R:<4} r={r:<4} n={n:<4}      "
              f"{n_paper}/{len(diagnostics)}              "
              f"{n_tk}/{len(diagnostics)}")
    print()


if __name__ == "__main__":
    main()
