"""Bogoliubov Phase E-3 — BFS in (vertex, winding) state space to L ≤ 20.

Phase E-2 enumerated walks exhaustively at L ≤ 9, finding 2/16 compendium
matches.  Brute enumeration to L ≤ 20 is intractable (~10^15 walks), but
the actual STATE SPACE is tiny: (current_vertex, accumulated_Δu, accumulated_Δv)
with all components bounded.  Use BFS on that state space to find the
shortest walk realizing each compendium (p, q).

State: (v, m_u, m_v) where v ∈ {0..6}, m_u, m_v ∈ ℤ in 1/7 units.
Initial: (0, 0, 0).
Goal for particle (p, q): state (0, 7p, 7q) reached via closed walk.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_e3_walk_BFS.py
"""
from __future__ import annotations

import math
from collections import deque, defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_e3_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)

MAX_LENGTH = 25


def bfs_shortest_walks(max_length: int = MAX_LENGTH, start: int = 0):
    """BFS in (vertex, m_u, m_v) state space. Returns dict mapping
    (m_u, m_v) (in 1/7 units, integer) reachable at vertex `start` to
    (shortest_length, walk_path) where walk_path is a list of vertex
    indices.

    State space is bounded: |m_u|, |m_v| ≤ max_length * 3.
    """
    # Precompute edge windings
    edge_w = {}
    for a in range(7):
        for b in range(7):
            if a != b:
                edge_w[(a, b)] = edge_winding_class(a, b)

    # BFS state: (v, m_u, m_v), with parents for path reconstruction
    initial = (start, 0, 0)
    visited = {initial: (0, None, None)}  # state → (depth, parent_state, last_v)
    queue = deque([initial])

    while queue:
        state = queue.popleft()
        depth, _, _ = visited[state]
        if depth >= max_length:
            continue
        v, m_u, m_v = state
        for nxt in range(7):
            if nxt == v:
                continue
            dnu, dnv = edge_w[(v, nxt)]
            new_state = (nxt, m_u + dnu, m_v + dnv)
            if new_state not in visited:
                visited[new_state] = (depth + 1, state, v)
                queue.append(new_state)

    # For each (m_u, m_v) reachable at vertex `start` (closed walk on cover),
    # extract shortest walk
    closed_walks: dict[tuple[int, int], tuple[int, list[int]]] = {}
    for state, (depth, parent, _) in visited.items():
        v, m_u, m_v = state
        if v != start:
            continue
        if (m_u, m_v) == (0, 0):
            continue   # trivial
        if m_u % 7 != 0 or m_v % 7 != 0:
            continue   # doesn't close on universal cover (not integer winding)
        p = m_u // 7
        q = m_v // 7
        # Reconstruct walk
        walk = [state[0]]
        cur = state
        while visited[cur][1] is not None:
            cur = visited[cur][1]
            walk.append(cur[0])
        walk.reverse()
        # walk now starts at `start` and ends at `start`
        key = (abs(p), abs(q))
        if key not in closed_walks or depth < closed_walks[key][0]:
            closed_walks[key] = (depth, walk, (p, q))
    return closed_walks


def main() -> None:
    print("=" * 78)
    print(f"BOGOLIUBOV PHASE E-3 — BFS K_7 walk enumeration to L ≤ {MAX_LENGTH}")
    print("=" * 78)
    print()

    # Build compendium (|p|, |q|) lookup
    pq_to_particles: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for e in COMPENDIUM:
        pq_to_particles[(abs(e["p"]), abs(e["q"]))].append(e)
    print(f"Compendium has {len(COMPENDIUM)} particles in "
          f"{len(pq_to_particles)} distinct (|p|, |q|) classes.")
    print()

    print(f"Running BFS in (vertex, winding) state space, max depth {MAX_LENGTH}…")
    closed = bfs_shortest_walks(MAX_LENGTH)
    print(f"  found {len(closed)} reachable (|p|, |q|) values "
          f"(non-trivial closed walks on universal cover)")
    print()

    # ---- Compendium matches ---------------------------------------------
    print("=" * 78)
    print(f"COMPENDIUM (|p|, |q|) REALIZATION at L ≤ {MAX_LENGTH}")
    print("=" * 78)
    print()
    print(f"  {'(|p|,|q|)':<10} {'L':<3} {'particles':<25} "
          f"{'shortest walk'}")
    print(f"  " + "-" * 90)
    matched = []
    unmatched = []
    for pq, particles in sorted(pq_to_particles.items()):
        names = sorted(set(p["name"] for p in particles))
        names_str = ', '.join(names[:3]) + ('...' if len(names) > 3 else '')
        if pq in closed:
            L, walk, signed_pq = closed[pq]
            walk_str = '→'.join(str(v) for v in walk)
            print(f"  {str(pq):<10} {L:<3} {names_str:<25} {walk_str}")
            matched.append((pq, names, L, walk, signed_pq))
        else:
            print(f"  {str(pq):<10} {'-':<3} {names_str:<25} "
                  f"(not reached at L ≤ {MAX_LENGTH})")
            unmatched.append((pq, names))
    print()
    print(f"  Coverage: {len(matched)}/{len(pq_to_particles)} classes "
          f"({len(matched)/len(pq_to_particles)*100:.0f}%)")
    print()

    # ---- σ-orbit composition of each particle's walk --------------------
    from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS

    def edge_to_orbit(a: int, b: int) -> int:
        """Find which σ-orbit contains edge (min(a,b), max(a,b))."""
        e_norm = (min(a, b), max(a, b))
        for oid, orbit in SIGMA_ORBITS.items():
            if e_norm in orbit["edges"]:
                return oid
        return -1

    print("=" * 78)
    print("σ-ORBIT COMPOSITION of shortest walks")
    print("=" * 78)
    print()
    print(f"  {'particle':<22} {'(p,q)':<8} {'L':<3} {'σ-orbit composition'}")
    print(f"  " + "-" * 76)
    for pq, names, L, walk, signed_pq in matched:
        orbits_used = defaultdict(int)
        for i in range(len(walk) - 1):
            o = edge_to_orbit(walk[i], walk[i + 1])
            orbits_used[o] += 1
        orbit_str = ', '.join(f'σ_{o}×{n}' for o, n in sorted(orbits_used.items()))
        names_str = names[0] if len(names) <= 2 else f'{names[0]}+'
        print(f"  {names_str:<22} {str(pq):<8} {L:<3} {orbit_str}")
    print()

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    # (a) (|p|, |q|) lattice — coverage at L ≤ 20
    ax = axes[0]
    # Compendium
    for pq in pq_to_particles.keys():
        ax.scatter(*pq, s=250, c='gray', alpha=0.4,
                    edgecolor='k', zorder=3)
    # Matched
    for pq, names, L, walk, signed_pq in matched:
        names_str = names[0] if len(names) <= 2 else f'{names[0]}+'
        ax.scatter(*pq, s=400, c='C2', alpha=0.85,
                    edgecolor='k', linewidth=2, zorder=10)
        ax.annotate(f'{names_str}\nL={L}', pq, xytext=(6, 6),
                    textcoords='offset points', fontsize=8,
                    color='C2', fontweight='bold')
    # Unmatched
    for pq, names in unmatched:
        names_str = names[0] if len(names) <= 2 else f'{names[0]}+'
        ax.scatter(*pq, s=250, c='C3', alpha=0.6,
                    edgecolor='k', marker='X', zorder=5)
        ax.annotate(f'{names_str}', pq, xytext=(6, 6),
                    textcoords='offset points', fontsize=7,
                    color='C3')
    # All reachable but no compendium match
    for pq in closed.keys():
        if pq in pq_to_particles:
            continue
        ax.scatter(*pq, s=40, c='C0', alpha=0.4,
                    marker='s', zorder=2)
    ax.set_xlabel('|p| (toroidal winding)')
    ax.set_ylabel('|q| (poloidal winding)')
    ax.set_title(f'K_7 walk realizability at L ≤ {MAX_LENGTH}\n'
                 f'green = compendium match, red X = no walk, blue □ = walk but no particle')
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.5, 11)
    ax.set_ylim(-0.5, 11)

    # (b) Walk length distribution
    ax = axes[1]
    lengths = [L for _, _, L, _, _ in matched]
    bins = np.arange(0, MAX_LENGTH + 2) - 0.5
    ax.hist(lengths, bins=bins, color='C2', edgecolor='black', alpha=0.7)
    # Annotate which particles at which L
    L_to_names = defaultdict(list)
    for pq, names, L, walk, signed_pq in matched:
        L_to_names[L].append(names[0])
    for L, namelist in L_to_names.items():
        names_str = ', '.join(namelist[:3])
        ax.text(L, 0.2, names_str, rotation=90, ha='center', va='bottom',
                fontsize=8, color='black')
    ax.set_xlabel('shortest walk length L')
    ax.set_ylabel('# (|p|, |q|) classes')
    ax.set_title(f'Walk-length distribution for {len(matched)} matched compendium classes')
    ax.set_xticks(range(0, MAX_LENGTH + 1, 2))
    ax.grid(alpha=0.3, axis='y')

    fig.suptitle(
        f"Phase E-3 — Extended K_7 walk enumeration to L ≤ {MAX_LENGTH}\n"
        f"Coverage: {len(matched)}/{len(pq_to_particles)} compendium classes "
        f"({len(matched)/len(pq_to_particles)*100:.0f}%)",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_e3_walk_BFS.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    # Save data
    np.savez(OUT_DIR / "phase_e3_walk_BFS.npz",
             max_length=MAX_LENGTH,
             matched_pq=np.array([pq for pq, _, _, _, _ in matched]),
             matched_particles=np.array(
                ['|'.join(n) for _, n, _, _, _ in matched]),
             walk_lengths=np.array(
                [L for _, _, L, _, _ in matched]),
             walk_paths=np.array(
                ['-'.join(str(v) for v in walk)
                 for _, _, _, walk, _ in matched]),
             unmatched_pq=np.array([pq for pq, _ in unmatched]),
             unmatched_particles=np.array(
                ['|'.join(n) for _, n in unmatched]))
    print(f"  data saved {OUT_DIR / 'phase_e3_walk_BFS.npz'}")

    # ---- Headline -------------------------------------------------------
    print()
    print("=" * 78)
    print(f"HEADLINE — Phase E-3 BFS to L ≤ {MAX_LENGTH}")
    print("=" * 78)
    print()
    print(f"  {len(matched)}/{len(pq_to_particles)} compendium (|p|, |q|) classes "
          f"realized at L ≤ {MAX_LENGTH}")
    print(f"  ({len(matched)/len(pq_to_particles)*100:.0f}% coverage)")
    print()
    if unmatched:
        print(f"  Unmatched ({len(unmatched)}):")
        for pq, names in unmatched:
            print(f"    {str(pq):<10} {', '.join(names)}")
    else:
        print(f"  ★ ALL compendium (|p|, |q|) classes realized as K_7 walks!")


if __name__ == "__main__":
    main()
