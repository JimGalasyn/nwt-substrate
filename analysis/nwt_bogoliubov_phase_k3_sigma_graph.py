"""Bogoliubov Phase K-3 — σ-orbit transition graph topology for n_q.

Test multiple invariants of the σ-orbit transition graph:
- # distinct σ-orbit unordered pairs
- # double-backs (σ_i↔σ_j visited both directions)
- # vertices visited even/odd times (Gauss-code interpretation)
- Cyclomatic complexity (E - V + 1)
- Primary cycle length (= n_q for proton observation)

Triggered by the observation that proton's σ-orbit sequence collapsed
(removing self-loops) is exactly a 5-cycle in the σ-graph:
  σ_0 → σ_3 → σ_5 → σ_4 → σ_1 → σ_0
which matches n_q = 5 exactly.

Plus partial result from earlier: #ω + #ω² = 5 matches ALL 5 nucleons.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_k3_sigma_graph.py
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate import bfs_shortest_walks, walk_sigma_sequence
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_k3_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def sigma_graph_invariants(walk: list[int]) -> dict:
    """Compute multiple σ-orbit transition graph invariants."""
    sig_seq = walk_sigma_sequence(walk)
    n = len(sig_seq)
    # COLLAPSED σ-sequence (remove consecutive duplicates)
    collapsed = []
    for s in sig_seq:
        if not collapsed or collapsed[-1] != s:
            collapsed.append(s)
    # If first and last are same after collapsing, treat as wrap
    if len(collapsed) > 1 and collapsed[0] == collapsed[-1]:
        collapsed = collapsed[:-1]
    # transitions (non-self-loop, with wrap)
    transitions = []
    for i in range(len(collapsed)):
        cur = collapsed[i]
        nxt = collapsed[(i + 1) % len(collapsed)]
        transitions.append((cur, nxt))
    # distinct vertices
    vertices = set(collapsed)
    n_V = len(vertices)
    n_E_total = len(transitions)
    # distinct unordered pairs
    unordered_pairs = set()
    pair_counts = Counter()
    for (a, b) in transitions:
        unordered_pairs.add((min(a, b), max(a, b)))
        pair_counts[(min(a, b), max(a, b))] += 1
    n_unique_pairs = len(unordered_pairs)
    # double-backs: pairs visited 2+ times
    double_backs = sum(1 for p, c in pair_counts.items() if c >= 2)
    # cyclomatic complexity (assume connected — collapsed cycle is)
    # # cycles = E_unique - V + (# connected components)
    # for a single connected component: cycles = E_unique - V + 1
    cyclomatic = n_unique_pairs - n_V + 1
    # vertex visit multiplicity (from collapsed sequence)
    visit_counts = Counter(collapsed)
    n_visited_once = sum(1 for v, c in visit_counts.items() if c == 1)
    n_visited_multiple = sum(1 for v, c in visit_counts.items() if c > 1)
    n_visited_even = sum(1 for v, c in visit_counts.items() if c % 2 == 0)
    n_visited_odd = sum(1 for v, c in visit_counts.items() if c % 2 == 1)
    # primary cycle length: if all vertices visited exactly once,
    # = len(collapsed); else =  some smaller cycle length
    if all(c == 1 for c in visit_counts.values()):
        primary_cycle = len(collapsed)
    else:
        # length of the longest run that visits each vertex once
        # heuristic: # distinct vertices
        primary_cycle = n_V
    return {
        "L": n,
        "collapsed_len": len(collapsed),
        "n_V": n_V,
        "n_unique_pairs": n_unique_pairs,
        "double_backs": double_backs,
        "cyclomatic": cyclomatic,
        "n_visited_once": n_visited_once,
        "n_visited_multiple": n_visited_multiple,
        "n_visited_even": n_visited_even,
        "n_visited_odd": n_visited_odd,
        "primary_cycle": primary_cycle,
        "collapsed": collapsed,
        "transitions": transitions,
        "pair_counts": dict(pair_counts),
    }


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE K-3 — σ-orbit transition graph topology for n_q")
    print("=" * 78)
    print()

    walks = bfs_shortest_walks(max_length=25)
    print(f"Loaded {len(walks)} (|p|, |q|) shortest walks")
    print()

    print(f"{'particle':<12} {'n_q':<4} {'collapsed':<7} {'n_V':<4} "
          f"{'unique':<7} {'dbacks':<7} {'cyclom':<7} {'prim_cyc':<9}")
    print("-" * 80)
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        inv = sigma_graph_invariants(walk)
        rows.append({**entry, "walk": walk, **inv})
        print(f"{entry['name']:<12} {entry['n_q']:<4} {inv['collapsed_len']:<7} "
              f"{inv['n_V']:<4} {inv['n_unique_pairs']:<7} "
              f"{inv['double_backs']:<7} {inv['cyclomatic']:<7} "
              f"{inv['primary_cycle']:<9}")
    print()

    # Test each invariant as n_q predictor
    print("=" * 78)
    print("CORRELATION ANALYSIS")
    print("=" * 78)
    print()
    candidates = {
        "collapsed_len (# σ-orbit visits after collapsing self-loops)": "collapsed_len",
        "n_V (# distinct σ-orbits)": "n_V",
        "n_unique_pairs (# distinct transition pairs)": "n_unique_pairs",
        "double_backs (# pairs visited twice)": "double_backs",
        "cyclomatic complexity (E-V+1)": "cyclomatic",
        "n_visited_once": "n_visited_once",
        "n_visited_multiple": "n_visited_multiple",
        "n_visited_even": "n_visited_even",
        "n_visited_odd": "n_visited_odd",
        "primary_cycle": "primary_cycle",
    }
    n_q_obs = np.array([r["n_q"] for r in rows])
    results = {}
    print(f"{'invariant':<55} {'exact':<8} {'Pearson r'}")
    print("-" * 75)
    for name, attr in candidates.items():
        vals = np.array([r[attr] for r in rows])
        exact = int(sum(1 for r in rows if r[attr] == r["n_q"]))
        r_p = np.corrcoef(n_q_obs, vals)[0, 1] if len(set(vals)) > 1 else float('nan')
        results[name] = (exact, r_p)
        print(f"  {name:<55} {exact}/25    {r_p:+.3f}")
    print()

    # Combined predictors: try simple combinations
    print("=" * 78)
    print("COMBINED PREDICTORS")
    print("=" * 78)
    print()
    combos = {
        "collapsed_len − 2·double_backs": lambda r: r["collapsed_len"] - 2 * r["double_backs"],
        "n_V − 2·double_backs": lambda r: r["n_V"] - 2 * r["double_backs"],
        "primary_cycle − 2·double_backs": lambda r: r["primary_cycle"] - 2 * r["double_backs"],
        "n_unique_pairs − double_backs": lambda r: r["n_unique_pairs"] - r["double_backs"],
        "cyclomatic + 1": lambda r: r["cyclomatic"] + 1,
        "2 × cyclomatic + (1 if cyclomatic > 0 else 0)": lambda r: 2 * r["cyclomatic"] + (1 if r["cyclomatic"] > 0 else 0),
    }
    print(f"{'combination':<55} {'exact':<8} {'Pearson r'}")
    print("-" * 75)
    for name, func in combos.items():
        try:
            vals = np.array([func(r) for r in rows])
            exact = int(sum(1 for r in rows if func(r) == r["n_q"]))
            r_p = np.corrcoef(n_q_obs, vals)[0, 1] if len(set(vals)) > 1 else float('nan')
            print(f"  {name:<55} {exact}/25    {r_p:+.3f}")
        except Exception as e:
            print(f"  {name:<55} ERROR: {e}")
    print()

    # ---- Detailed per-particle σ-graph structure ----------------------
    print("=" * 78)
    print("DETAILED σ-GRAPH STRUCTURE — selected particles")
    print("=" * 78)
    for name in ["e-", "mu-", "p", "K+", "pi+", "Lambda", "Omega-"]:
        r = next((x for x in rows if x["name"] == name), None)
        if r is None: continue
        print()
        print(f"  {name} (n_q={r['n_q']}):")
        print(f"    walk: {'→'.join(str(v) for v in r['walk'])}")
        print(f"    σ-collapsed: {r['collapsed']}")
        print(f"    transitions: {r['transitions']}")
        print(f"    pair counts: {r['pair_counts']}")
        print(f"    n_V={r['n_V']}, unique_pairs={r['n_unique_pairs']}, "
              f"double_backs={r['double_backs']}")
        print(f"    cyclomatic = E-V+1 = {r['n_unique_pairs']}-{r['n_V']}+1 = "
              f"{r['cyclomatic']}")

    # ---- Best predictor summary ---------------------------------------
    best_single = max(results, key=lambda k: results[k][0])
    print()
    print("=" * 78)
    print("HEADLINE — Phase K-3 σ-graph topology")
    print("=" * 78)
    print()
    print(f"  Best single invariant: '{best_single}'")
    print(f"    {results[best_single][0]}/25 exact matches, r = {results[best_single][1]:+.3f}")
    print()
    # Show sector-by-sector breakdown for best
    best_attr = candidates[best_single]
    by_sector = defaultdict(list)
    for r in rows:
        by_sector[r["n_q"]].append(r[best_attr])
    for n_q in sorted(by_sector.keys()):
        vals = sorted(set(by_sector[n_q]))
        sector = {0:'leptons', 2:'mesons', 3:'hyperons', 5:'nucleons'}.get(n_q, '?')
        print(f"    n_q={n_q} ({sector}): values = {vals}")

    # Save
    np.savez(OUT_DIR / "phase_k3_sigma_graph.npz",
             names=np.array([r["name"] for r in rows]),
             n_q=np.array([r["n_q"] for r in rows]),
             collapsed_len=np.array([r["collapsed_len"] for r in rows]),
             n_V=np.array([r["n_V"] for r in rows]),
             n_unique_pairs=np.array([r["n_unique_pairs"] for r in rows]),
             double_backs=np.array([r["double_backs"] for r in rows]),
             cyclomatic=np.array([r["cyclomatic"] for r in rows]),
             primary_cycle=np.array([r["primary_cycle"] for r in rows]))

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, (name, attr) in zip(axes.flatten(), list(candidates.items())[:6]):
        vals = [r[attr] for r in rows]
        sector_colors = {0:'C2', 2:'C0', 3:'C1', 5:'C3'}
        for r in rows:
            c = sector_colors.get(r["n_q"], 'gray')
            ax.scatter(r[attr], r["n_q"], c=c, s=80, alpha=0.6, edgecolor='k')
        lims = [-0.5, max(max(vals), 6) + 1]
        ax.plot(lims, lims, 'k--', alpha=0.4, label='1:1')
        r_p = results[name][1]
        exact = results[name][0]
        ax.set_xlabel(attr)
        ax.set_ylabel('n_q')
        ax.set_title(f'{attr}: {exact}/25, r={r_p:+.2f}', fontsize=10)
        ax.grid(alpha=0.3)
    fig.suptitle("Phase K-3 — σ-orbit transition graph topology vs n_q",
                  fontsize=11, y=1.005)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "phase_k3_sigma_graph.png", dpi=130, bbox_inches='tight')
    plt.close(fig)
    print()
    print(f"  saved figure + npz to {OUT_DIR}")


if __name__ == "__main__":
    main()
