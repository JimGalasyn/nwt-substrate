"""Bogoliubov Phase F — derive m and n_q quantum numbers from K_7 walk topology.

Phase E established that Paper 11's (p, q) torus-knot quantum numbers
ARE K_7 walk homology classes (100% compendium coverage at L ≤ 25).
Phase F tests whether m (phase-closure integer) and n_q (linked-
component count) are also derivable from the SAME walks' geometry.

m candidate
-----------
Paper 11 phase closure:    β = √(m²/p² − 1),  i.e.,  m² = p²(β² + 1).
For a (p, q) torus knot with major radius R and minor radius ξ:
  - vortex arc length ℓ = 2π · √((pR)² + (qξ)²)
  - phase closure: m · 2π = k · ℓ with k = 1/ξ
  - ⇒ m = √((pR/ξ)² + q²) = √(p² β² + q²)
  - Compare to Paper 11: m² = p²(β² + 1) = p²β² + p²
  - ⇒ Paper 11 implicitly identifies q² ↔ p² for the phase-closure
     condition, OR the geometry differs.

Test: compute K_7 walk arc length, derive β_walk, predict m_walk, compare
to compendium m.

n_q candidate
-------------
Number of LINKED COMPONENTS of the vortex.  Candidates from K_7 walk:
  - Number of distinct edges used (vs total walk length)
  - Edge multiplicity (max repeat count of any edge)
  - Self-crossings of the walk in the Heffter projection
  - Number of vertices visited

Test: for each compendium particle, compute these K_7 invariants and
look for the relation to Paper 11's n_q.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_f_m_nq_derivation.py
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.orbit_winding import (
    edge_winding_class, HEFFTER_VERT_UV,
)
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_f_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def walk_arc_length_on_unit_torus(walk: list[int]) -> float:
    """Arc length of a walk on the unit torus, where each edge's
    displacement is its (Δu, Δv)/7 on the universal cover.

    Length = Σ_i √((Δu_i/7)² + (Δv_i/7)²)
    """
    total = 0.0
    for i in range(len(walk) - 1):
        nu, nv = edge_winding_class(walk[i], walk[i + 1])
        total += math.sqrt(nu ** 2 + nv ** 2) / 7.0
    return total


def walk_edge_multiplicity(walk: list[int]) -> dict:
    """Multiset of undirected edges used by the walk."""
    counter = Counter()
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        e = (min(a, b), max(a, b))
        counter[e] += 1
    return dict(counter)


def walk_vertex_multiplicity(walk: list[int]) -> dict:
    """Vertices visited and their counts."""
    counter = Counter(walk[:-1])  # don't double-count closing vertex
    return dict(counter)


def predict_m_from_geometry(walk: list[int], p_target: int,
                              q_target: int) -> dict:
    """Phase F-1: derive m candidate from walk arc length.

    Theory:
      For a (p, q) torus knot with aspect β = R/ξ:
        m_theoretical² = p² (β² + 1) = (pβ)² + p²       (Paper 11)
        m_theoretical² = (pβ)² + q²                       (geometric)

      For K_7 walk:
        arc_length = 2π · √((pR)² + (qξ)²) / (some unit length)
        ⇒ (arc_length × scale)² = (pR)² + (qξ)²
        ⇒ pR · scale ≡ arc_length component in u direction
        ⇒ qξ · scale ≡ arc_length component in v direction
    """
    arc = walk_arc_length_on_unit_torus(walk)
    # Two candidate m formulas:
    # (1) m = arc * scale_1     (linear in arc length)
    # (2) m² = p² + (arc * scale_2)²    (Paper 11 Pythagorean)
    return {
        "arc_length": arc,
        "m_paper11_target": None,    # filled in by caller
        "p_target": p_target,
        "q_target": q_target,
    }


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE F — derive m, n_q from K_7 walk topology")
    print("=" * 78)
    print()

    # ---- Build BFS for shortest walks (re-use Phase E-3 logic) ---------
    print("Step 1 — load shortest K_7 walks for each compendium particle.")
    from nwt_substrate.condensate.orbit_winding import edge_winding_class as ewc
    from collections import deque

    def bfs_shortest_walks_for(max_length: int = 25):
        edge_w = {(a, b): ewc(a, b) for a in range(7) for b in range(7) if a != b}
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
        # Find shortest walks to each (m_u, m_v) at v=0
        walks: dict[tuple[int, int], list[int]] = {}
        for state, (depth, _) in visited.items():
            v, m_u, m_v = state
            if v != 0 or (m_u, m_v) == (0, 0):
                continue
            if m_u % 7 != 0 or m_v % 7 != 0:
                continue
            p, q = m_u // 7, m_v // 7
            key = (abs(p), abs(q))
            # Reconstruct walk
            walk = [state[0]]
            cur = state
            while visited[cur][1] is not None:
                cur = visited[cur][1]
                walk.append(cur[0])
            walk.reverse()
            if key not in walks or len(walk) - 1 < len(walks[key]) - 1:
                walks[key] = walk
        return walks

    shortest_walks = bfs_shortest_walks_for(25)
    print(f"  Found {len(shortest_walks)} (|p|,|q|) classes reachable at L ≤ 25")
    print()

    # ---- Build per-particle data table ---------------------------------
    print("Step 2 — compute walk invariants for each compendium particle.")
    print()
    rows = []
    for entry in COMPENDIUM:
        p_q = (abs(entry["p"]), abs(entry["q"]))
        if p_q not in shortest_walks:
            continue
        walk = shortest_walks[p_q]
        L = len(walk) - 1
        arc = walk_arc_length_on_unit_torus(walk)
        edge_mult = walk_edge_multiplicity(walk)
        vert_mult = walk_vertex_multiplicity(walk)
        n_distinct_edges = len(edge_mult)
        n_distinct_vertices = len(vert_mult)
        max_edge_repeat = max(edge_mult.values())
        max_vert_repeat = max(vert_mult.values())
        rows.append({
            "name": entry["name"],
            "p": entry["p"],
            "q": entry["q"],
            "m": entry["m"],
            "n_q": entry["n_q"],
            "L": L,
            "arc_length": arc,
            "n_distinct_edges": n_distinct_edges,
            "n_distinct_vertices": n_distinct_vertices,
            "max_edge_repeat": max_edge_repeat,
            "max_vert_repeat": max_vert_repeat,
            "walk": walk,
        })

    # ---- Phase F-1: m derivation candidates ----------------------------
    print("=" * 78)
    print("PHASE F-1: m (phase-closure integer) candidates")
    print("=" * 78)
    print()
    print(f"  {'particle':<12} {'(p,q,m,n_q)':<14} {'L':<4} "
          f"{'arc':<7} "
          f"{'arc·7/p':<9} {'√(arc²·49+p²)':<14} {'L-p-q':<7} {'L_int':<7}")
    print(f"  " + "-" * 95)
    for r in rows:
        ratio = r['arc_length'] * 7 / abs(r['p']) if r['p'] != 0 else float('nan')
        # Paper 11 m: m² = p² + (p·β)² where β = R/r is a "free" geometric param
        # K_7 walk's arc/7 gives an effective β·p ≈ arc·7 in some normalization
        m_geom_candidate = math.sqrt((r['arc_length'] * 7) ** 2 + r['p'] ** 2)
        print(f"  {r['name']:<12} "
              f"({r['p']},{r['q']},{r['m']},{r['n_q']})".ljust(16) +
              f" {r['L']:<4} {r['arc_length']:<7.3f} "
              f"{ratio:<9.3f} {m_geom_candidate:<14.3f} "
              f"{r['L']-abs(r['p'])-abs(r['q']):<7} {int(r['L']/abs(r['p'])) if r['p'] else 0:<7}")
    print()

    # ---- Phase F-2: n_q derivation candidates --------------------------
    print("=" * 78)
    print("PHASE F-2: n_q (linked-component count) candidates")
    print("=" * 78)
    print()
    print(f"  {'particle':<12} {'(p,q,m,n_q)':<14} {'L':<3} "
          f"{'#dist_edges':<13} {'#dist_verts':<13} "
          f"{'max_e_rep':<10} {'max_v_rep':<10}")
    print(f"  " + "-" * 90)
    for r in rows:
        print(f"  {r['name']:<12} "
              f"({r['p']},{r['q']},{r['m']},{r['n_q']})".ljust(16) +
              f" {r['L']:<3} "
              f"{r['n_distinct_edges']:<13} "
              f"{r['n_distinct_vertices']:<13} "
              f"{r['max_edge_repeat']:<10} "
              f"{r['max_vert_repeat']:<10}")
    print()

    # ---- Statistical fits ------------------------------------------------
    print("=" * 78)
    print("CORRELATION ANALYSIS")
    print("=" * 78)
    print()

    # m vs walk arc length / L / etc.
    m_obs = np.array([r['m'] for r in rows], dtype=float)
    L_arr = np.array([r['L'] for r in rows], dtype=float)
    arc_arr = np.array([r['arc_length'] for r in rows], dtype=float)
    p_arr = np.array([abs(r['p']) for r in rows], dtype=float)
    q_arr = np.array([abs(r['q']) for r in rows], dtype=float)
    n_q_obs = np.array([r['n_q'] for r in rows], dtype=float)
    ne_arr = np.array([r['n_distinct_edges'] for r in rows], dtype=float)
    nv_arr = np.array([r['n_distinct_vertices'] for r in rows], dtype=float)
    me_arr = np.array([r['max_edge_repeat'] for r in rows], dtype=float)
    mv_arr = np.array([r['max_vert_repeat'] for r in rows], dtype=float)

    print("Correlation of m with K_7 walk invariants (Pearson r):")
    for name, x in [
        ("L (walk length)", L_arr),
        ("arc_length (unit torus)", arc_arr),
        ("L - p - q", L_arr - p_arr - q_arr),
        ("L / p", L_arr / p_arr),
        ("arc · 7 / p", arc_arr * 7 / p_arr),
        ("√((arc·7)² + p²)", np.sqrt((arc_arr * 7) ** 2 + p_arr ** 2)),
        ("p · √(L²/p² - 1)", p_arr * np.sqrt(np.maximum(L_arr ** 2 / p_arr ** 2 - 1, 0))),
    ]:
        r_corr = np.corrcoef(m_obs, x)[0, 1] if len(set(x)) > 1 else float('nan')
        print(f"    r(m, {name:<30}) = {r_corr:+.4f}")
    print()

    print("Correlation of n_q with K_7 walk invariants (Pearson r):")
    for name, x in [
        ("L (walk length)", L_arr),
        ("# distinct edges", ne_arr),
        ("# distinct vertices", nv_arr),
        ("max edge repeat", me_arr),
        ("max vertex repeat", mv_arr),
        ("L - n_distinct_edges", L_arr - ne_arr),
        ("L / n_distinct_edges", L_arr / ne_arr),
    ]:
        r_corr = np.corrcoef(n_q_obs, x)[0, 1] if len(set(x)) > 1 else float('nan')
        print(f"    r(n_q, {name:<30}) = {r_corr:+.4f}")
    print()

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) m vs candidate predictors
    ax = axes[0, 0]
    candidate_x = np.sqrt((arc_arr * 7) ** 2 + p_arr ** 2)
    ax.scatter(candidate_x, m_obs, s=120, c='C0', alpha=0.7,
                edgecolor='k')
    for r in rows:
        x = math.sqrt((r['arc_length'] * 7) ** 2 + r['p'] ** 2)
        ax.annotate(r['name'], (x, r['m']), xytext=(4, 4),
                    textcoords='offset points', fontsize=7, alpha=0.7)
    ax.set_xlabel('√((arc·7)² + p²)')
    ax.set_ylabel('m (Paper 11)')
    r_val = np.corrcoef(m_obs, candidate_x)[0, 1]
    ax.set_title(f'm vs walk arc-length geometric candidate\n'
                 f'Pearson r = {r_val:.3f}')
    # 1:1 line
    lims = [min(min(candidate_x), min(m_obs)) - 1,
            max(max(candidate_x), max(m_obs)) + 1]
    ax.plot(lims, lims, 'k--', alpha=0.4, label='1:1')
    ax.legend()
    ax.grid(alpha=0.3)

    # (b) n_q vs candidate predictors
    ax = axes[0, 1]
    ax.scatter(me_arr, n_q_obs, s=120, c='C3', alpha=0.7,
                edgecolor='k')
    for r in rows:
        ax.annotate(r['name'], (r['max_edge_repeat'], r['n_q']),
                    xytext=(4, 4), textcoords='offset points',
                    fontsize=7, alpha=0.7)
    ax.set_xlabel('max edge repeat in walk')
    ax.set_ylabel('n_q (Paper 11)')
    r_val = np.corrcoef(n_q_obs, me_arr)[0, 1]
    ax.set_title(f'n_q vs max edge repeat\n'
                 f'Pearson r = {r_val:.3f}')
    ax.grid(alpha=0.3)

    # (c) m residual = m - candidate
    ax = axes[1, 0]
    candidate_m = np.sqrt((arc_arr * 7) ** 2 + p_arr ** 2)
    residual = m_obs - candidate_m
    ys = np.arange(len(rows))
    ax.barh(ys, residual, color='C0', alpha=0.7)
    ax.set_yticks(ys)
    ax.set_yticklabels([r['name'] for r in rows], fontsize=8)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel('m - √((arc·7)² + p²)')
    ax.set_title('m residual against geometric candidate')
    ax.grid(alpha=0.3, axis='x')

    # (d) Walk-structure stats by particle
    ax = axes[1, 1]
    width = 0.2
    xs = np.arange(len(rows))
    ax.bar(xs - 1.5*width, L_arr, width=width, label='L', alpha=0.8)
    ax.bar(xs - 0.5*width, ne_arr, width=width, label='# dist edges')
    ax.bar(xs + 0.5*width, nv_arr, width=width, label='# dist verts')
    ax.bar(xs + 1.5*width, me_arr, width=width, label='max e_rep')
    ax.set_xticks(xs)
    ax.set_xticklabels([r['name'] for r in rows], rotation=90, fontsize=7)
    ax.set_ylabel('count')
    ax.set_title('Walk structural invariants per particle')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis='y')

    fig.suptitle(
        f"Phase F — derive m, n_q from K_7 walk topology\n"
        f"({len(rows)} compendium particles with shortest-walk data)",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_f_m_nq.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_f_m_nq.npz",
             names=np.array([r['name'] for r in rows]),
             p=np.array([r['p'] for r in rows]),
             q=np.array([r['q'] for r in rows]),
             m=m_obs,
             n_q=n_q_obs,
             L=L_arr,
             arc_length=arc_arr,
             n_distinct_edges=ne_arr,
             n_distinct_vertices=nv_arr,
             max_edge_repeat=me_arr,
             max_vert_repeat=mv_arr)
    print(f"  data saved {OUT_DIR / 'phase_f_m_nq.npz'}")


if __name__ == "__main__":
    main()
