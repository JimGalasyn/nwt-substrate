"""Bogoliubov Phase F-3 — σ-orbit linkage for n_q closure.

NWT's Phase F-2 follow-up established four sector-Hamilton categories:
  - Exclusive d=1 → nucleons (n_q=5)
  - Exclusive d=3 → hyperons (n_q=3)
  - Mixed d=1 + d=3 → mesons (n_q=2)
  - No Hamilton → leptons (n_q=0) OR no-Hamilton seeds

But 7 walks have NO full Hamilton sub-cycle yet have non-zero n_q:
  - hyperons: τ⁻, Λ, Σ* (n_q=3)
  - mesons:   ω, ρ, π⁺ (n_q=2)
  - lepton:   e⁻ (n_q=0)

NWT proposes σ-orbit linkage provides the carrier_base for these.

Phase F-3 enumerates structural invariants per compendium walk and
identifies the strongest predictor of n_q across all 25 particles —
especially the non-Hamilton sub-set.

Invariants computed per walk:
  - σ-orbit composition (count of edges in each of 7 σ-orbits)
  - Vertex visit multiplicity (each of 7 vertices)
  - Distinct σ-orbits visited
  - σ-orbit pair coincidence (σ_i + σ_j edges traversed together)
  - Walk decomposition into sub-loops (visits to v_0)
  - σ-orbit COMPLETE coverage (all 3 edges of an orbit used)
  - Bridging-vs-non-bridging ratio

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_f3_sigma_linkage.py
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS, orbit_invariants
from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_f3_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def edge_to_orbit(a: int, b: int) -> int:
    """Find which σ-orbit contains edge (min(a,b), max(a,b)). Returns -1 if not found."""
    e_norm = (min(a, b), max(a, b))
    for oid, orbit in SIGMA_ORBITS.items():
        if e_norm in orbit["edges"]:
            return oid
    return -1


def walk_invariants(walk: list[int]) -> dict:
    """Compute a comprehensive set of K_7 invariants for a closed walk."""
    L = len(walk) - 1
    orbit_counts = Counter()
    vertex_counts = Counter()
    orbit_seq = []
    orbit_pair_coincidence = Counter()
    edges_used = Counter()

    for i in range(L):
        a, b = walk[i], walk[i + 1]
        oid = edge_to_orbit(a, b)
        orbit_counts[oid] += 1
        orbit_seq.append(oid)
        edges_used[(min(a, b), max(a, b))] += 1
        vertex_counts[a] += 1
        # σ-orbit transition coincidence (consecutive pairs)
        if i > 0:
            prev_oid = orbit_seq[i - 1]
            orbit_pair_coincidence[(min(prev_oid, oid), max(prev_oid, oid))] += 1
    vertex_counts[walk[-1]] += 1  # closing visit

    # σ-orbit complete coverage: which orbits had ALL 3 of their edges used?
    sigma_complete = []
    for oid in range(7):
        orbit_edges = set(SIGMA_ORBITS[oid]["edges"])
        used_in_orbit = orbit_edges & set(edges_used.keys())
        if len(used_in_orbit) == 3:
            sigma_complete.append(oid)

    # Distinct σ-orbits visited
    n_distinct_orbits = len([o for o, c in orbit_counts.items() if c > 0])

    # Sub-loop count: number of times we visit v_0 (= return events).
    # Equivalently: walk decomposes into N sub-loops if it visits v_0 N times.
    n_return_zero = sum(1 for v in walk[1:-1] if v == 0) + 1  # +1 for final closure

    # Bridging fraction
    bridging_count = orbit_counts[0] + orbit_counts[1]
    bridging_frac = bridging_count / L if L > 0 else 0.0

    # Equatorial-triangle count (σ_3 + σ_4)
    triangle_count = orbit_counts[3] + orbit_counts[4]
    # Cross count (σ_2 + σ_5 + σ_6)
    cross_count = orbit_counts[2] + orbit_counts[5] + orbit_counts[6]

    return {
        "L": L,
        "orbit_counts": dict(orbit_counts),
        "vertex_counts": dict(vertex_counts),
        "n_distinct_edges": len(edges_used),
        "n_distinct_orbits": n_distinct_orbits,
        "n_sigma_complete": len(sigma_complete),
        "sigma_complete": sigma_complete,
        "n_return_zero": n_return_zero,
        "bridging_count": bridging_count,
        "bridging_frac": bridging_frac,
        "triangle_count": triangle_count,
        "cross_count": cross_count,
        "edges_used": dict(edges_used),
    }


def load_shortest_walks(max_length: int = 25) -> dict:
    from collections import deque
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


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE F-3 — σ-orbit linkage for n_q closure")
    print("=" * 78)
    print()

    print("Step 1 — Load shortest K_7 walks (L ≤ 25)…")
    walks = load_shortest_walks(25)
    print(f"  {len(walks)} (|p|, |q|) classes")
    print()

    # ---- Build per-particle invariants table ---------------------------
    print("Step 2 — Compute invariants per compendium particle.")
    print()
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        inv = walk_invariants(walk)
        rows.append({**entry, **inv, "walk": walk})

    # ---- Print table: σ-orbit composition + key invariants -------------
    print(f"  {'particle':<12} {'n_q':<4} {'L':<3} "
          f"{'σ0':<3} {'σ1':<3} {'σ2':<3} {'σ3':<3} {'σ4':<3} {'σ5':<3} {'σ6':<3} "
          f"{'#σ_complete':<11} {'#dist_σ':<8} {'#ret_0'}")
    print("  " + "-" * 95)
    for r in sorted(rows, key=lambda x: (x['n_q'], x['L'])):
        oc = r['orbit_counts']
        print(f"  {r['name']:<12} {r['n_q']:<4} {r['L']:<3} "
              f"{oc.get(0,0):<3} {oc.get(1,0):<3} {oc.get(2,0):<3} "
              f"{oc.get(3,0):<3} {oc.get(4,0):<3} {oc.get(5,0):<3} {oc.get(6,0):<3} "
              f"{r['n_sigma_complete']:<11} {r['n_distinct_orbits']:<8} "
              f"{r['n_return_zero']}")
    print()

    # ---- Correlation analysis: which invariant best predicts n_q? ------
    print("=" * 78)
    print("CORRELATION ANALYSIS: which K_7 invariant predicts n_q?")
    print("=" * 78)
    print()
    n_q_obs = np.array([r['n_q'] for r in rows], dtype=float)

    candidates = {
        "L (walk length)": np.array([r['L'] for r in rows]),
        "n_distinct_edges": np.array([r['n_distinct_edges'] for r in rows]),
        "n_distinct_orbits": np.array([r['n_distinct_orbits'] for r in rows]),
        "n_sigma_complete": np.array([r['n_sigma_complete'] for r in rows]),
        "n_return_zero": np.array([r['n_return_zero'] for r in rows]),
        "bridging_count (σ0+σ1)": np.array([r['bridging_count'] for r in rows]),
        "bridging_frac": np.array([r['bridging_frac'] for r in rows]),
        "triangle_count (σ3+σ4)": np.array([r['triangle_count'] for r in rows]),
        "cross_count (σ2+σ5+σ6)": np.array([r['cross_count'] for r in rows]),
        "σ_2 count": np.array([r['orbit_counts'].get(2, 0) for r in rows]),
        "σ_3 count": np.array([r['orbit_counts'].get(3, 0) for r in rows]),
        "σ_4 count": np.array([r['orbit_counts'].get(4, 0) for r in rows]),
        "σ_5 count": np.array([r['orbit_counts'].get(5, 0) for r in rows]),
        "σ_3 + σ_4 - σ_2": (np.array([r['orbit_counts'].get(3, 0) for r in rows])
                              + np.array([r['orbit_counts'].get(4, 0) for r in rows])
                              - np.array([r['orbit_counts'].get(2, 0) for r in rows])),
        "(σ_3 × σ_4)/L": np.array([r['orbit_counts'].get(3, 0)
                                    * r['orbit_counts'].get(4, 0) / max(r['L'], 1)
                                    for r in rows]),
        "σ_3 / (σ_2 + 1)": np.array([r['orbit_counts'].get(3, 0)
                                      / (r['orbit_counts'].get(2, 0) + 1)
                                      for r in rows]),
    }
    print(f"  {'invariant':<32} {'Pearson r':>10}")
    print("  " + "-" * 45)
    correlations = {}
    for name, vals in candidates.items():
        if len(set(vals)) > 1:
            r_val = np.corrcoef(n_q_obs, vals)[0, 1]
            correlations[name] = r_val
            print(f"  {name:<32} {r_val:>+10.4f}")
    print()

    # ---- Best candidate detail -----------------------------------------
    best_name = max(correlations.keys(), key=lambda k: abs(correlations[k]))
    print(f"  ★ Best predictor: {best_name}  (|r| = {abs(correlations[best_name]):.4f})")
    print()

    # ---- σ-orbit complete-coverage analysis (NWT linkage hypothesis) ---
    print("=" * 78)
    print("σ-ORBIT COMPLETE-COVERAGE per particle (NWT linkage hypothesis)")
    print("=" * 78)
    print()
    print("  Which σ-orbits have ALL 3 edges traversed at least once?")
    print()
    print(f"  {'particle':<12} {'n_q':<4} {'sigma_complete':<25} "
          f"{'#sigma_complete':<15}")
    print("  " + "-" * 65)
    for r in sorted(rows, key=lambda x: (x['n_q'], x['L'])):
        sc_str = ', '.join(f'σ_{o}' for o in r['sigma_complete'])
        print(f"  {r['name']:<12} {r['n_q']:<4} {sc_str:<25} "
              f"{r['n_sigma_complete']}")
    print()

    # ---- Stratify by sector: lepton, meson, hyperon, nucleon -----------
    print("=" * 78)
    print("STRATIFIED ANALYSIS by n_q sector")
    print("=" * 78)
    print()
    by_nq = defaultdict(list)
    for r in rows:
        by_nq[r['n_q']].append(r)

    for n_q_val in sorted(by_nq.keys()):
        members = by_nq[n_q_val]
        sector = {0: 'leptons', 2: 'mesons', 3: 'hyperons',
                   5: 'nucleons'}.get(n_q_val, '?')
        names = [m['name'] for m in members]
        avg_sigma_complete = np.mean([m['n_sigma_complete'] for m in members])
        avg_distinct = np.mean([m['n_distinct_orbits'] for m in members])
        avg_ret_zero = np.mean([m['n_return_zero'] for m in members])
        avg_bridging = np.mean([m['bridging_count'] for m in members])
        avg_triangle = np.mean([m['triangle_count'] for m in members])
        avg_cross = np.mean([m['cross_count'] for m in members])
        print(f"  {sector} (n_q={n_q_val}, {len(members)} particles): "
              f"{', '.join(names[:6])}{'...' if len(names) > 6 else ''}")
        print(f"    avg #σ_complete={avg_sigma_complete:.2f}  "
              f"#dist_σ={avg_distinct:.2f}  ")
        print(f"    avg #ret_0={avg_ret_zero:.2f}  "
              f"avg #bridging={avg_bridging:.2f}  "
              f"avg #triangle={avg_triangle:.2f}  "
              f"avg #cross={avg_cross:.2f}")
        print()

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) Best predictor scatter vs n_q
    ax = axes[0, 0]
    best_x = candidates[best_name]
    ax.scatter(best_x, n_q_obs, s=120, c='C0', alpha=0.7, edgecolor='k')
    for r, x in zip(rows, best_x):
        ax.annotate(r['name'], (x, r['n_q']), xytext=(4, 4),
                    textcoords='offset points', fontsize=7, alpha=0.7)
    ax.set_xlabel(f"{best_name}")
    ax.set_ylabel("n_q (Paper 11)")
    ax.set_title(f"Best predictor: {best_name}\nPearson r = {correlations[best_name]:+.4f}")
    ax.grid(alpha=0.3)

    # (b) #σ_complete per particle, colored by n_q sector
    ax = axes[0, 1]
    sector_colors = {0: 'C2', 2: 'C0', 3: 'C1', 5: 'C3'}
    sectors_plotted = set()
    for r in rows:
        c = sector_colors.get(r['n_q'], 'gray')
        lbl = f'n_q={r["n_q"]}' if r['n_q'] not in sectors_plotted else None
        sectors_plotted.add(r['n_q'])
        ax.scatter(r['L'], r['n_sigma_complete'], c=c, s=100, alpha=0.7,
                    edgecolor='k', label=lbl)
        ax.annotate(r['name'], (r['L'], r['n_sigma_complete']),
                    xytext=(3, 3), textcoords='offset points', fontsize=7)
    ax.set_xlabel('L (walk length)')
    ax.set_ylabel('# σ-orbits with COMPLETE coverage')
    ax.set_title('σ-orbit complete-coverage by particle + sector')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # (c) σ-orbit count heatmap (rows=particles, cols=σ_0..σ_6)
    ax = axes[1, 0]
    M = np.zeros((len(rows), 7))
    for i, r in enumerate(rows):
        for o in range(7):
            M[i, o] = r['orbit_counts'].get(o, 0)
    im = ax.imshow(M, aspect='auto', cmap='YlOrRd', origin='lower')
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r['name'] for r in rows], fontsize=7)
    ax.set_xticks(range(7))
    ax.set_xticklabels([f'σ_{i}' for i in range(7)])
    for i in range(len(rows)):
        for j in range(7):
            v = int(M[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha='center', va='center',
                        fontsize=6,
                        color='white' if v > M.max() / 2 else 'black')
    ax.set_title('σ-orbit edge counts per particle\n(higher = more uses of that σ-orbit)')
    plt.colorbar(im, ax=ax, label='# edges in σ-orbit')

    # (d) n_q vs σ_3 + σ_4 (triangle count) — possible n_q signature
    ax = axes[1, 1]
    triangle_arr = candidates["triangle_count (σ3+σ4)"]
    sector_color_list = [sector_colors.get(r['n_q'], 'gray') for r in rows]
    ax.scatter(triangle_arr, n_q_obs, s=120, c=sector_color_list, alpha=0.7,
                edgecolor='k')
    for r, x in zip(rows, triangle_arr):
        ax.annotate(r['name'], (x, r['n_q']), xytext=(4, 4),
                    textcoords='offset points', fontsize=7, alpha=0.7)
    ax.set_xlabel('triangle_count = σ_3 + σ_4 (E + F equatorial triangles)')
    ax.set_ylabel('n_q (Paper 11)')
    r_val = correlations.get('triangle_count (σ3+σ4)', 0)
    ax.set_title(f'n_q vs equatorial triangle count\n(Pearson r = {r_val:+.4f})')
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"Phase F-3 — σ-orbit linkage and n_q closure\n"
        f"best predictor: {best_name} (|r| = {abs(correlations[best_name]):.3f})",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_f3_sigma_linkage.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_f3_sigma_linkage.npz",
             names=np.array([r['name'] for r in rows]),
             n_q=n_q_obs,
             L=np.array([r['L'] for r in rows]),
             n_distinct_orbits=np.array([r['n_distinct_orbits'] for r in rows]),
             n_sigma_complete=np.array([r['n_sigma_complete'] for r in rows]),
             n_return_zero=np.array([r['n_return_zero'] for r in rows]),
             bridging_count=np.array([r['bridging_count'] for r in rows]),
             triangle_count=np.array([r['triangle_count'] for r in rows]),
             cross_count=np.array([r['cross_count'] for r in rows]),
             sigma_0=np.array([r['orbit_counts'].get(0, 0) for r in rows]),
             sigma_1=np.array([r['orbit_counts'].get(1, 0) for r in rows]),
             sigma_2=np.array([r['orbit_counts'].get(2, 0) for r in rows]),
             sigma_3=np.array([r['orbit_counts'].get(3, 0) for r in rows]),
             sigma_4=np.array([r['orbit_counts'].get(4, 0) for r in rows]),
             sigma_5=np.array([r['orbit_counts'].get(5, 0) for r in rows]),
             sigma_6=np.array([r['orbit_counts'].get(6, 0) for r in rows]))
    print(f"  data saved {OUT_DIR / 'phase_f3_sigma_linkage.npz'}")


if __name__ == "__main__":
    main()
