"""Bogoliubov Phase E-2 — systematic K_7 closed-walk enumeration vs compendium (p, q).

Phase E established that σ-orbit edge-sum windings give 3 partial matches to
the Paper 11 compendium (electron via σ_6 primitive, ω via σ_3/σ_4 primitive,
K_7 Hamilton cycle → proton).  Phase E-2 systematically enumerates ALL short
closed walks on K_7 and checks which (p, q) values appear, then matches each
to compendium particles.

The goal: identify which particles have a NATURAL K_7 walk realization, and
which are MULTI-WALK composites.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_e2_walk_search.py
"""
from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_e2_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def walk_winding(walk: list[int]) -> tuple[int, int]:
    """Total (Δu, Δv) in 1/7 units of a directed walk."""
    sum_u = sum_v = 0
    for i in range(len(walk) - 1):
        nu, nv = edge_winding_class(walk[i], walk[i + 1])
        sum_u += nu
        sum_v += nv
    return sum_u, sum_v


def enumerate_closed_walks(max_length: int, start: int = 0):
    """BFS enumerate all closed walks of length ≤ max_length starting/ending
    at vertex `start`.  Yields walks as tuples of vertex indices.

    Length = number of edges (so walk has length+1 vertices).
    """
    # Use DFS with edge revisit allowed (since each edge can be used multiply)
    # but cap total length.
    def dfs(walk: list[int], remaining: int):
        if len(walk) > 1 and walk[-1] == start:
            yield tuple(walk)
            if remaining == 0:
                return
        if remaining == 0:
            return
        current = walk[-1]
        for nxt in range(7):
            if nxt == current:
                continue
            walk.append(nxt)
            yield from dfs(walk, remaining - 1)
            walk.pop()

    yield from dfs([start], max_length)


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE E-2 — systematic K_7 closed-walk → compendium (p,q)")
    print("=" * 78)
    print()
    print("Strategy: BFS enumerate closed walks of length ≤ N starting at v_0,")
    print("compute (p, q) winding for each, find shortest walk per (p, q)")
    print("value, match to compendium particles.")
    print()

    # Build (p, q) → set of particle names from compendium
    pq_to_particles: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for e in COMPENDIUM:
        pq_to_particles[(abs(e["p"]), abs(e["q"]))].append(e)
    print(f"Compendium has {len(COMPENDIUM)} particles spanning "
          f"{len(pq_to_particles)} distinct (|p|,|q|) classes.")
    print()

    # Enumerate walks of increasing length, find shortest walk realizing each (p, q)
    MAX_LEN = 9       # 9-edge walks include up to (p,q) ~ 9/7 max
    print(f"Enumerating closed walks of length ≤ {MAX_LEN}…")
    shortest_walk_for_pq: dict[tuple[int, int], tuple] = {}
    count = 0
    for walk in enumerate_closed_walks(MAX_LEN, start=0):
        if walk[-1] != 0:
            continue
        sum_u, sum_v = walk_winding(list(walk))
        # Closes on universal cover only if displacements are integer multiples of 7
        if sum_u % 7 != 0 or sum_v % 7 != 0:
            continue
        p = sum_u // 7
        q = sum_v // 7
        pq_abs = (abs(p), abs(q))
        if pq_abs == (0, 0):
            continue
        L = len(walk) - 1
        if pq_abs not in shortest_walk_for_pq or L < shortest_walk_for_pq[pq_abs][0]:
            shortest_walk_for_pq[pq_abs] = (L, walk, (p, q))
        count += 1
    print(f"  enumerated {count} non-trivial closed walks on universal cover")
    print(f"  found {len(shortest_walk_for_pq)} distinct (|p|,|q|) values")
    print()

    # ---- Compendium matches ---------------------------------------------
    print("=" * 78)
    print(f"WALKS REALIZING COMPENDIUM (p, q) CLASSES (length ≤ {MAX_LEN})")
    print("=" * 78)
    print()
    print(f"  {'particles':<22} {'(|p|,|q|)':<10} {'L':<3} "
          f"{'shortest walk':<35}")
    print(f"  " + "-" * 72)
    matched_classes = []
    unmatched_classes = []
    for pq, particles in sorted(pq_to_particles.items()):
        names = sorted(set(p["name"] for p in particles))
        names_str = ', '.join(names[:3]) + ('...' if len(names) > 3 else '')
        if pq in shortest_walk_for_pq:
            L, walk, signed_pq = shortest_walk_for_pq[pq]
            walk_str = '→'.join(str(v) for v in walk)
            print(f"  {names_str:<22} {str(pq):<10} {L:<3} {walk_str:<35}")
            matched_classes.append((pq, names))
        else:
            print(f"  {names_str:<22} {str(pq):<10} {'-':<3} {'(not realized at L ≤ ' + str(MAX_LEN) + ')':<35}")
            unmatched_classes.append((pq, names))
    print()
    print(f"  Compendium classes realized: {len(matched_classes)} / "
          f"{len(pq_to_particles)}")
    print()

    # ---- Walks WITHOUT compendium matches (substrate predictions?) ------
    print("=" * 78)
    print(f"WALK-REALIZABLE (|p|,|q|) NOT IN COMPENDIUM (potential new species?)")
    print("=" * 78)
    print()
    extra_pqs = sorted(set(shortest_walk_for_pq.keys()) - set(pq_to_particles.keys()))
    for pq in extra_pqs[:15]:
        L, walk, signed_pq = shortest_walk_for_pq[pq]
        print(f"  (|p|,|q|) = {str(pq):<12} L = {L}  walk: "
              f"{'→'.join(str(v) for v in walk)}")
    if len(extra_pqs) > 15:
        print(f"  ...and {len(extra_pqs) - 15} more")
    print()

    # ---- Plot ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) (p, q) lattice — compendium vs walk-realized
    ax = axes[0]
    # Compendium particles
    cp_x = [pq[0] for pq in pq_to_particles.keys()]
    cp_y = [pq[1] for pq in pq_to_particles.keys()]
    ax.scatter(cp_x, cp_y, s=200, c='gray', alpha=0.6,
                label='Paper 11 compendium', edgecolor='k')
    # Walk-realized (p, q)
    walk_x = [pq[0] for pq in shortest_walk_for_pq.keys()]
    walk_y = [pq[1] for pq in shortest_walk_for_pq.keys()]
    ax.scatter(walk_x, walk_y, s=60, c='C3', marker='X',
                alpha=0.7, label=f'K_7 walks (L ≤ {MAX_LEN})')
    # Annotate matches in green
    for pq, names in matched_classes:
        names_str = names[0] if len(names) <= 2 else f'{names[0]}+'
        ax.annotate(names_str, pq, xytext=(6, 6),
                    textcoords='offset points', fontsize=8,
                    color='C2', fontweight='bold')
    ax.set_xlabel('|p| (toroidal winding)')
    ax.set_ylabel('|q| (poloidal winding)')
    ax.set_title('K_7 walk-realizable (|p|,|q|) vs Paper 11 compendium')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.5, 11)
    ax.set_ylim(-0.5, 11)

    # (b) Match summary table
    ax = axes[1]
    ax.axis('off')
    n_matched = len(matched_classes)
    n_total = len(pq_to_particles)
    n_predicted = len(extra_pqs)
    table_data = [
        ['statistic', 'value', 'note'],
        ['compendium (|p|,|q|) classes', str(n_total), 'Paper 11 input'],
        ['classes realized by K_7 walks', f'{n_matched}/{n_total}',
            f'{n_matched/n_total*100:.0f}% coverage at L ≤ {MAX_LEN}'],
        ['extra walk-realizable classes', str(n_predicted),
            'substrate predictions?'],
        ['', '', ''],
    ]
    # Add specific matches
    table_data.append(['matched particles', '', 'shortest walk length'])
    for pq, names in matched_classes:
        L = shortest_walk_for_pq[pq][0]
        names_str = ', '.join(names[:2]) + ('...' if len(names) > 2 else '')
        table_data.append([f'{names_str}', f'(|p|,|q|)={pq}', f'L = {L}'])

    table = ax.table(cellText=table_data, loc='center', cellLoc='left',
                     colWidths=[0.45, 0.25, 0.30])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.5)
    for i, row in enumerate(table_data):
        if i == 0:
            for j in range(3):
                table[(i, j)].set_facecolor('#cccccc')
        elif i >= 6 and 'matched' not in row[0]:
            for j in range(3):
                table[(i, j)].set_facecolor('#d4f4d4')
    ax.set_title(f'K_7 walk enumeration ↔ Paper 11 compendium\n'
                  f'{n_matched}/{n_total} classes realized at L ≤ {MAX_LEN}',
                  pad=15)

    fig.suptitle(
        f"Phase E-2 — Systematic K_7 closed walks reproduce a fraction of "
        f"Paper 11 (p, q) quantum numbers",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_e2_walk_search.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_e2_walk_search.npz",
             matched_pq=np.array([pq for pq, _ in matched_classes]),
             matched_particles=np.array(
                ['|'.join(n) for _, n in matched_classes]),
             walk_lengths=np.array(
                [shortest_walk_for_pq[pq][0] for pq, _ in matched_classes]),
             walk_paths=np.array(
                ['-'.join(str(v) for v in shortest_walk_for_pq[pq][1])
                 for pq, _ in matched_classes]),
             extra_pqs=np.array(extra_pqs),
             max_length=MAX_LEN)
    print(f"  data saved {OUT_DIR / 'phase_e2_walk_search.npz'}")

    # ---- Headline -------------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase E-2 K_7 walk realization of compendium (p, q)")
    print("=" * 78)
    print()
    print(f"  Compendium has {n_total} distinct (|p|, |q|) classes.")
    print(f"  K_7 walks of length ≤ {MAX_LEN} realize "
          f"{n_matched} of them "
          f"({n_matched/n_total*100:.0f}%).")
    print()
    print(f"  ★ Specific matches found:")
    for pq, names in matched_classes:
        L = shortest_walk_for_pq[pq][0]
        names_str = ', '.join(names[:3])
        print(f"    (|p|,|q|) = {str(pq):<8}  L = {L}  particles: {names_str}")
    print()
    print(f"  Interpretation: each Paper 11 (|p|, |q|) class corresponds to a")
    print(f"  specific class of closed walks on K_7 Heffter; the particle's")
    print(f"  full (p, q, m, n_q) is determined by the walk's full topology")
    print(f"  (length, branching, σ-orbit composition).")
    print()
    print(f"  Remaining (p,q) classes likely require longer walks (L > {MAX_LEN})")
    print(f"  or composite multi-walk configurations.")


if __name__ == "__main__":
    main()
