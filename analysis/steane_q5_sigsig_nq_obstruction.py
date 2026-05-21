"""Q5 — verify σ-sig cannot determine n_q (walk ordering matters).

Q4 showed σ-sig → n_q regression has R² = 0.34, suggesting σ-sig
does NOT determine n_q. This script gives a direct existence proof:
enumerate all closed walks of fixed length L starting at vertex 0
with a given (|p|, |q|) winding; group by σ-sig; if two walks share
σ-sig but lift to different carrier knots on the Heffter torus,
σ-sig cannot determine n_q.

Strategy: for tractable small classes (e.g. (1, 3) at L=7 = nucleons,
(2, 1) at L=5 = electron), enumerate all closed walks via brute-force
DFS; group by σ-sig; run pyknotid on representatives.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q5_sigsig_nq_obstruction.py
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt   # before numpy compat patches

# numpy compat patches for pyknotid (uses deprecated aliases)
np.float = float
np.int = int
np.complex = complex
np.object = object
np.str = str
np.long = int

from pyknotid.spacecurves import Knot

from nwt_substrate.condensate.orbit_winding import (
    edge_winding_class, HEFFTER_VERT_UV,
)
from nwt_substrate.condensate.walks import edge_to_orbit


R_MAJOR = 2.5
R_MINOR = 1.0


def walk_to_3d_curve(walk: list[int], n_per_edge: int = 30,
                     R: float = R_MAJOR, r: float = R_MINOR) -> np.ndarray:
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


def identify_knot(walk: list[int]) -> tuple[str, int | None]:
    try:
        curve = walk_to_3d_curve(walk)
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


def sigma_signature(walk: list[int]) -> tuple[int, ...]:
    counts = [0] * 7
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        oid = edge_to_orbit(a, b)
        if 0 <= oid < 7:
            counts[oid] += 1
    return tuple(counts)


def enumerate_closed_walks(target_pq: tuple[int, int], L: int,
                            start: int = 0) -> list[list[int]]:
    """DFS-enumerate all closed walks of length L starting/ending at `start`
    with winding (p, q) = target_pq (mod 7) on universal cover.

    For tractability, do NOT remove the cyclic/rotation equivalence —
    we want all sequences as visited.
    """
    edge_w = {(a, b): edge_winding_class(a, b)
              for a in range(7) for b in range(7) if a != b}
    p_target, q_target = target_pq
    results: list[list[int]] = []

    def dfs(path: list[int], dpu: int, dpv: int):
        if len(path) - 1 == L:
            if path[-1] == start and dpu == 7 * p_target and dpv == 7 * q_target:
                results.append(list(path))
            return
        v = path[-1]
        for nxt in range(7):
            if nxt == v:
                continue
            dnu, dnv = edge_w[(v, nxt)]
            path.append(nxt)
            dfs(path, dpu + dnu, dpv + dnv)
            path.pop()

    dfs([start], 0, 0)
    return results


def main() -> None:
    print("=" * 78)
    print("Q5 — σ-sig cannot determine n_q: ordering matters")
    print("=" * 78)
    print()
    print("Strategy: for a fixed (|p|,|q|) class at fixed L, enumerate ALL")
    print("closed walks. Group by σ-sig. If two walks share σ-sig but")
    print("differ in carrier knot, the σ-sig → n_q map is multi-valued.")
    print()

    # Both signs of (p, q) are equivalent for σ-sig; pick +p, +q. Take
    # the (p, q) classes from the smallest particle compendium classes.
    targets = [
        ("electron", (2, 1), 5),    # (2,1) L=5  — 6^5 = 7776 worst case
        ("nucleon",  (1, 3), 7),    # (1,3) L=7  — 6^7 = 279k worst case
    ]

    for label, (p, q), L in targets:
        print(f"--- {label}: (|p|,|q|) = ({p},{q}), L = {L} ---")
        all_walks = enumerate_closed_walks((p, q), L)
        print(f"  enumerated {len(all_walks)} closed walks")
        # Group by σ-sig
        by_sig = defaultdict(list)
        for w in all_walks:
            by_sig[sigma_signature(w)].append(w)
        print(f"  distinct σ-sigs: {len(by_sig)}")
        # For each σ-sig group, sample up to 5 walks and identify knot
        for sig, group in sorted(by_sig.items(),
                                  key=lambda kv: -len(kv[1])):
            sig_str = '[' + ' '.join(str(c) for c in sig) + ']'
            print(f"  σ-sig {sig_str}: {len(group)} walks")
            knot_set = set()
            sample = group[: min(len(group), 3)]
            for w in sample:
                kid, nc = identify_knot(w)
                knot_set.add((kid, nc))
                print(f"    walk {w} → carrier knot = {kid} (min_crossings = {nc})")
            if len(knot_set) > 1:
                print(f"  ✓ MULTIPLE KNOT TYPES same σ-sig — σ-sig CANNOT determine n_q")
            else:
                print(f"  → single knot type observed for this σ-sig (within sample)")
        print()


if __name__ == "__main__":
    main()
