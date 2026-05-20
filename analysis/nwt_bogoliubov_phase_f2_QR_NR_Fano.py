"""Bogoliubov Phase F-2 — verify NWT structural laws + use them for m, n_q.

NWT session this morning verified two laws against my Phase E-3 walk data:

  LAW 1 (Hamilton-containment ↔ full Fano-coverage):
    13/16 compendium walks touch all 7 Fano triangles.
    3 exceptions are exactly the non-Hamilton seed walks: e⁻, Λ, π⁺.

  LAW 2 (Three K_7 Hamilton cycles by edge-difference d ∈ {1, 2, 3}):
    d = 1 (QR) — proton walk 0→1→2→3→4→5→6→0 = matter scaffold
    d = 2 (QR) — secondary 0→2→4→6→1→3→5→0
    d = 3 (NR) — 0→3→6→2→5→1→4→0 = σ_6 / CP-violation channel

    NR-Hamilton (d=3) dominates in σ_6 trio: π⁰, Ω⁻, K⁰.

Phase F-2: verify both laws independently from my walk data, then use
the QR/NR decomposition as candidate predictors for n_q and m.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_f2_QR_NR_Fano.py
"""
from __future__ import annotations

import math
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_f2_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


# ---------------------------------------------------------------------------
# Paley QR / NR on Z_7
# ---------------------------------------------------------------------------

# Quadratic residues mod 7: {1, 2, 4}
QR_DIFFS = {1, 2, 4}
NR_DIFFS = {3, 5, 6}
# Symmetric (unsigned) difference classes:
#   |a-b| mod 7 ∈ {1, 2, 3}; QR-symmetric = {1, 2}, NR-symmetric = {3}
# This is the convention NWT uses for "Hamilton by edge-difference"


def edge_signed_diff(a: int, b: int) -> int:
    """Signed step (b - a) mod 7, in {1..6}."""
    return (b - a) % 7


def edge_symmetric_d(a: int, b: int) -> int:
    """Symmetric edge-difference, in {1, 2, 3}."""
    d = (b - a) % 7
    return min(d, 7 - d)


def hamilton_skip_d(d: int) -> list[int]:
    """The K_7 'skip-d' Hamilton cycle starting at 0:
       0 → d → 2d → 3d → 4d → 5d → 6d → 0 (all mod 7).
    Works for d coprime to 7, i.e., d ∈ {1, 2, 3, 4, 5, 6}.
    """
    walk = [0]
    for k in range(1, 7):
        walk.append((k * d) % 7)
    walk.append(0)
    return walk


# Verify the 3 distinct Hamilton-by-d structures
HAMILTON_CYCLES = {
    1: hamilton_skip_d(1),   # QR — proton's scaffold
    2: hamilton_skip_d(2),   # QR — secondary
    3: hamilton_skip_d(3),   # NR — σ_6 / CP channel
}


# ---------------------------------------------------------------------------
# Fano plane on K_7: 7 triangles (Steiner (7, 3, 1) system)
# ---------------------------------------------------------------------------

# Standard Fano-plane line incidence: 7 lines, each containing 3 of the 7 points.
# One standard realization (Heffter rotation-compatible):
FANO_TRIANGLES = [
    {0, 1, 3},
    {1, 2, 4},
    {2, 3, 5},
    {3, 4, 6},
    {4, 5, 0},
    {5, 6, 1},
    {6, 0, 2},
]


def walk_fano_coverage(walk: list[int]) -> tuple[int, list[int]]:
    """Number of distinct Fano triangles touched by the walk's vertex set
    (a triangle is 'touched' iff at least one of its 3 vertices is in the walk).
    Returns (count, list of triangle indices touched).
    """
    verts = set(walk)
    touched = []
    for i, tri in enumerate(FANO_TRIANGLES):
        if tri & verts:
            touched.append(i)
    return len(touched), touched


def walk_fano_full_coverage(walk: list[int]) -> bool:
    """All 7 Fano triangles touched?"""
    return walk_fano_coverage(walk)[0] == 7


# ---------------------------------------------------------------------------
# QR / NR step analysis of a walk
# ---------------------------------------------------------------------------

def walk_step_sequence(walk: list[int]) -> list[int]:
    """List of signed step sizes (b - a) mod 7 ∈ {1..6} along the walk."""
    return [edge_signed_diff(walk[i], walk[i + 1])
            for i in range(len(walk) - 1)]


def walk_max_run_of_step(walk: list[int], target_d: int) -> int:
    """Longest consecutive run of step-d in the walk."""
    steps = walk_step_sequence(walk)
    max_run = 0
    cur = 0
    for s in steps:
        if s == target_d or s == (7 - target_d):
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 0
    return max_run


def walk_contains_hamilton_subcycle(walk: list[int], d: int) -> int:
    """Number of complete K_7-Hamilton-d sub-cycles embedded in the walk.

    A Hamilton-d sub-cycle is 7 consecutive steps all of size ±d (mod 7),
    visiting all 7 vertices and returning to start.
    """
    steps = walk_step_sequence(walk)
    n = len(steps)
    count = 0
    i = 0
    while i + 7 <= n:
        if all(s == d or s == (7 - d) for s in steps[i:i+7]):
            count += 1
            i += 7
        else:
            i += 1
    return count


# ---------------------------------------------------------------------------
# Load shortest walks (re-compute from Phase E-3 BFS)
# ---------------------------------------------------------------------------

def load_shortest_walks(max_length: int = 25) -> dict:
    from collections import deque
    from nwt_substrate.condensate.orbit_winding import edge_winding_class
    edge_w = {(a, b): edge_winding_class(a, b)
              for a in range(7) for b in range(7) if a != b}
    initial = (0, 0, 0)
    visited = {initial: (0, None)}
    queue = deque([initial])
    while queue:
        state = queue.popleft()
        depth, _ = visited[state]
        if depth >= max_length:
            continue
        v, m_u, m_v = state
        for nxt in range(7):
            if nxt == v:
                continue
            dnu, dnv = edge_w[(v, nxt)]
            new_state = (nxt, m_u + dnu, m_v + dnv)
            if new_state not in visited:
                visited[new_state] = (depth + 1, state)
                queue.append(new_state)
    walks = {}
    for state, (depth, _) in visited.items():
        v, m_u, m_v = state
        if v != 0 or (m_u, m_v) == (0, 0):
            continue
        if m_u % 7 != 0 or m_v % 7 != 0:
            continue
        p, q = m_u // 7, m_v // 7
        key = (abs(p), abs(q))
        walk = [state[0]]
        cur = state
        while visited[cur][1] is not None:
            cur = visited[cur][1]
            walk.append(cur[0])
        walk.reverse()
        if key not in walks or len(walk) - 1 < len(walks[key]) - 1:
            walks[key] = walk
    return walks


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE F-2 — verify NWT laws + use for m, n_q")
    print("=" * 78)
    print()

    # ---- Verify the 3 Hamilton cycles by edge-difference ----------------
    print("Step 1 — Verify NWT's three K_7 Hamilton cycles (skip-d cycles):")
    print(f"  d   polarity   cycle")
    print(f"  -   --------   -----")
    for d, cycle in HAMILTON_CYCLES.items():
        polarity = "QR" if d in (1, 2) else "NR"  # d=3 (and 4=7-3 too) is NR
        role = ""
        if d == 1: role = " ← proton walk = matter scaffold"
        if d == 3: role = " ← σ_6 / CP channel"
        cycle_str = '→'.join(str(v) for v in cycle)
        print(f"  {d}   {polarity}         {cycle_str}{role}")
    print()

    # ---- Load shortest walks --------------------------------------------
    print("Step 2 — Load shortest K_7 walks (Phase E-3 BFS, L ≤ 25)…")
    walks = load_shortest_walks(25)
    print(f"  {len(walks)} (|p|, |q|) classes reachable")
    print()

    # ---- Build analysis table -------------------------------------------
    print("Step 3 — Per-compendium-particle: QR/NR/Fano analysis")
    print()
    print(f"  {'particle':<12} {'(p,q,m,n_q)':<14} {'L':<3} "
          f"{'QR1':<4} {'QR2':<4} {'NR3':<4} "
          f"{'maxR_d3':<8} {'Fano':<5} {'Ham⊂?'}")
    print("  " + "-" * 80)
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        steps = walk_step_sequence(walk)
        L = len(steps)
        count_d = Counter()
        for s in steps:
            d = min(s, 7 - s)
            count_d[d] += 1
        max_d3_run = walk_max_run_of_step(walk, 3)
        fano_count, _ = walk_fano_coverage(walk)
        # Hamilton subcycles
        ham_d1 = walk_contains_hamilton_subcycle(walk, 1)
        ham_d2 = walk_contains_hamilton_subcycle(walk, 2)
        ham_d3 = walk_contains_hamilton_subcycle(walk, 3)
        ham_str = ''
        if ham_d1: ham_str += f'd1×{ham_d1} '
        if ham_d2: ham_str += f'd2×{ham_d2} '
        if ham_d3: ham_str += f'd3×{ham_d3} '
        if not ham_str: ham_str = '(none)'
        rows.append({
            **entry,
            "walk": walk,
            "L": L,
            "count_d1": count_d.get(1, 0),
            "count_d2": count_d.get(2, 0),
            "count_d3": count_d.get(3, 0),
            "max_d3_run": max_d3_run,
            "fano_count": fano_count,
            "ham_d1": ham_d1, "ham_d2": ham_d2, "ham_d3": ham_d3,
        })
        print(f"  {entry['name']:<12} "
              f"({entry['p']},{entry['q']},{entry['m']},{entry['n_q']})".ljust(16) +
              f" {L:<3} {count_d.get(1, 0):<4} {count_d.get(2, 0):<4} "
              f"{count_d.get(3, 0):<4} "
              f"{max_d3_run:<8} {fano_count:<5} {ham_str}")
    print()

    # ---- Verify LAW 1: Hamilton-containment ↔ Fano-full ----------------
    print("=" * 78)
    print("LAW 1 VERIFICATION — Hamilton containment ↔ Fano coverage")
    print("=" * 78)
    print()
    ham_full = [r for r in rows if (r['ham_d1'] or r['ham_d2'] or r['ham_d3']) > 0]
    ham_none = [r for r in rows if r['ham_d1'] == 0 and r['ham_d2'] == 0
                and r['ham_d3'] == 0]
    fano_full = [r for r in rows if r['fano_count'] == 7]
    fano_partial = [r for r in rows if r['fano_count'] < 7]
    print(f"  Walks containing ANY skip-d Hamilton sub-cycle: {len(ham_full)}/{len(rows)}")
    print(f"  Walks with full Fano coverage (7/7):            {len(fano_full)}/{len(rows)}")
    print()
    print(f"  Partial-Fano (< 7) walks: {[r['name'] for r in fano_partial]}")
    print(f"  No-Hamilton-subcycle walks: {[r['name'] for r in ham_none]}")
    print()
    overlap = set(r['name'] for r in fano_partial) & set(r['name'] for r in ham_none)
    print(f"  Both partial-Fano AND no-Hamilton: {sorted(overlap)}")
    if set(r['name'] for r in fano_partial) == set(r['name'] for r in ham_none):
        print(f"  ★ LAW 1 CONFIRMED: Hamilton-containment ⟺ full Fano-coverage")
    else:
        print(f"  LAW 1 partial: sets overlap but differ")
        print(f"    Fano-partial only: "
              f"{set(r['name'] for r in fano_partial) - overlap}")
        print(f"    Ham-none only: "
              f"{set(r['name'] for r in ham_none) - overlap}")

    # ---- Verify LAW 2: NR-d3 dominance in σ_6 trio ---------------------
    print()
    print("=" * 78)
    print("LAW 2 VERIFICATION — NR-Hamilton (d=3) dominates σ_6 trio")
    print("=" * 78)
    print()
    rows_by_d3run = sorted(rows, key=lambda r: -r['max_d3_run'])
    print(f"  Walks ranked by max d=3 run length:")
    print(f"  {'particle':<12} {'max d3-run':<11} {'σ_6 trio?'}")
    sigma6_trio = {'pi0', 'Omega-', 'K0'}
    for r in rows_by_d3run[:10]:
        marker = "★" if r['name'] in sigma6_trio else ""
        print(f"  {r['name']:<12} {r['max_d3_run']:<11} {marker}")
    print()

    # ---- n_q candidate: QR/NR Hamilton sub-cycle decomposition ---------
    print("=" * 78)
    print("PHASE F-2 n_q CANDIDATE — QR/NR Hamilton sub-cycle decomposition")
    print("=" * 78)
    print()
    print("Hypothesis: n_q = ham_d1 + ham_d2 + ham_d3 + non-Hamilton-component-count")
    print()
    print(f"  {'particle':<12} {'(p,q,m,n_q)':<14} "
          f"{'ham_d1':<7} {'ham_d2':<7} {'ham_d3':<7} "
          f"{'sum':<5} {'n_q':<5} {'Δ'}")
    n_q_pred_arr = []
    n_q_obs_arr = []
    for r in rows:
        n_q_pred = r['ham_d1'] + r['ham_d2'] + r['ham_d3']
        n_q_pred_arr.append(n_q_pred)
        n_q_obs_arr.append(r['n_q'])
        diff = r['n_q'] - n_q_pred
        print(f"  {r['name']:<12} "
              f"({r['p']},{r['q']},{r['m']},{r['n_q']})".ljust(16) +
              f" {r['ham_d1']:<7} {r['ham_d2']:<7} {r['ham_d3']:<7} "
              f"{n_q_pred:<5} {r['n_q']:<5} {diff:+}")
    r_corr = np.corrcoef(n_q_obs_arr, n_q_pred_arr)[0, 1]
    print()
    print(f"  Pearson r(n_q, ham_d1+ham_d2+ham_d3) = {r_corr:+.3f}")
    print()

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) Fano coverage vs Hamilton containment
    ax = axes[0, 0]
    ham_yes = [int(r['ham_d1'] + r['ham_d2'] + r['ham_d3'] > 0) for r in rows]
    fano_arr = [r['fano_count'] for r in rows]
    colors = ['C2' if h else 'C3' for h in ham_yes]
    ax.scatter(ham_yes, fano_arr, c=colors, s=120, alpha=0.7,
               edgecolor='k')
    for r, hy in zip(rows, ham_yes):
        ax.annotate(r['name'], (hy, r['fano_count']),
                    xytext=(np.random.uniform(-15, 15),
                            np.random.uniform(-15, 15)),
                    textcoords='offset points', fontsize=7, alpha=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['no Ham sub-cycle', 'has Ham sub-cycle'])
    ax.set_ylabel('Fano triangles touched (/7)')
    ax.set_title(f'LAW 1: Hamilton-containment ⟺ Fano-coverage\n'
                 f'(green = both, red = neither expected)')
    ax.set_ylim(4, 7.5)
    ax.grid(alpha=0.3)

    # (b) max d=3 run by particle (σ_6 trio highlighted)
    ax = axes[0, 1]
    names = [r['name'] for r in rows]
    d3_runs = [r['max_d3_run'] for r in rows]
    bar_colors = ['C3' if n in sigma6_trio else 'C0' for n in names]
    ax.bar(range(len(rows)), d3_runs, color=bar_colors, alpha=0.8)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_ylabel('max d=3 (NR) run length')
    ax.set_title(f'LAW 2: NR-Hamilton dominance\n(red = σ_6 trio: π⁰, Ω⁻, K⁰)')
    ax.grid(alpha=0.3, axis='y')

    # (c) n_q candidate vs observed
    ax = axes[1, 0]
    ax.scatter(n_q_pred_arr, n_q_obs_arr, s=120, c='C0', alpha=0.7,
               edgecolor='k')
    for r, p_ in zip(rows, n_q_pred_arr):
        ax.annotate(r['name'], (p_, r['n_q']),
                    xytext=(4, 4), textcoords='offset points',
                    fontsize=7, alpha=0.7)
    ax.plot([0, max(n_q_obs_arr) + 1], [0, max(n_q_obs_arr) + 1],
            'k--', alpha=0.4, label='1:1')
    ax.set_xlabel('n_q predicted = ham_d1 + ham_d2 + ham_d3')
    ax.set_ylabel('n_q observed (Paper 11)')
    ax.set_title(f'n_q candidate from Hamilton sub-cycle count\n'
                 f'Pearson r = {r_corr:.3f}')
    ax.legend()
    ax.grid(alpha=0.3)

    # (d) QR vs NR step balance per particle
    ax = axes[1, 1]
    qr_arr = [r['count_d1'] + r['count_d2'] for r in rows]
    nr_arr = [r['count_d3'] for r in rows]
    width = 0.4
    xs = np.arange(len(rows))
    ax.bar(xs - width/2, qr_arr, width, label='QR (d=1,2)', color='C2', alpha=0.7)
    ax.bar(xs + width/2, nr_arr, width, label='NR (d=3)', color='C3', alpha=0.7)
    ax.set_xticks(xs)
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_ylabel('# of steps')
    ax.set_title('QR vs NR step counts per particle')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    fig.suptitle(
        f"Phase F-2 — NWT structural laws applied to K_7 walks\n"
        f"LAW 1 (Hamilton ⟺ Fano) + LAW 2 (NR dominance in σ_6 trio)",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_f2_QR_NR_Fano.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_f2_QR_NR_Fano.npz",
             names=np.array([r['name'] for r in rows]),
             p=np.array([r['p'] for r in rows]),
             q=np.array([r['q'] for r in rows]),
             m=np.array([r['m'] for r in rows]),
             n_q_obs=np.array([r['n_q'] for r in rows]),
             n_q_pred=np.array(n_q_pred_arr),
             count_d1=np.array([r['count_d1'] for r in rows]),
             count_d2=np.array([r['count_d2'] for r in rows]),
             count_d3=np.array([r['count_d3'] for r in rows]),
             max_d3_run=np.array([r['max_d3_run'] for r in rows]),
             fano_count=np.array([r['fano_count'] for r in rows]),
             ham_d1=np.array([r['ham_d1'] for r in rows]),
             ham_d2=np.array([r['ham_d2'] for r in rows]),
             ham_d3=np.array([r['ham_d3'] for r in rows]))
    print(f"  data saved {OUT_DIR / 'phase_f2_QR_NR_Fano.npz'}")


if __name__ == "__main__":
    main()
