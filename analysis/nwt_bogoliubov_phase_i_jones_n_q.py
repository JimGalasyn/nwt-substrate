"""Bogoliubov Phase I — Jones polynomial closure of n_q via pyknotid.

NWT Q5 hypothesis: n_q is the crossing number of the carrier knot
realized by the walk's σ-orbit traversal — computable via Jones
polynomial / knot identification.

Implementation:
  1. For each K_7 walk, lift the discrete graph walk to a 3D space curve
     on the Heffter torus in R³ (using torus radii R_major, r_minor).
  2. Pass to pyknotid.spacecurves.Knot.identify() to determine the
     carrier knot type.
  3. The identified knot's .min_crossings IS the predicted n_q.
  4. Compare to compendium n_q values for all 25 particles.
  5. Predict n_q for the 4 substrate predictions ((2,2), (2,3), (3,1), (3,3)).

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_i_jones_n_q.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt   # imported BEFORE numpy patches

# numpy compat patches for pyknotid (uses deprecated aliases).
# Only patch the deprecated aliases pyknotid needs; avoid np.bool because
# numpy.ma uses bool internally and patching breaks np.ma.masked.
np.float = float
np.int = int
np.complex = complex
np.object = object
np.str = str
np.long = int

from pyknotid.spacecurves import Knot

from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.condensate.orbit_winding import (
    edge_winding_class, HEFFTER_VERT_UV,
)


OUT_DIR = Path(__file__).parent / "phase_i_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)

# Torus radii for the 3D lift
R_MAJOR = 2.5     # major radius
R_MINOR = 1.0     # minor radius (must be < R_MAJOR for embedded torus)


def walk_to_3d_curve(walk: list[int], n_per_edge: int = 30,
                     R: float = R_MAJOR, r: float = R_MINOR) -> np.ndarray:
    """Lift a K_7 walk (sequence of vertex indices) to a 3D space curve
    on the Heffter torus embedded in R³.

    Each vertex k sits at (u_k, v_k) = (k/7, 3k mod 7 /7) on the unit
    torus.  Each edge is traversed along the SHORTEST geodesic on the
    universal cover; the cumulative (u, v) position is unwrapped to keep
    the curve continuous, then mapped to (x, y, z) via:

      x = (R + r cos(2π v)) cos(2π u)
      y = (R + r cos(2π v)) sin(2π u)
      z = r sin(2π v)

    Returns array of shape (N, 3) where N = (L · n_per_edge).
    """
    pts = []
    # Start at first vertex
    u_curr = HEFFTER_VERT_UV[walk[0]][0]
    v_curr = HEFFTER_VERT_UV[walk[0]][1]
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        # Compute short-path displacement (Δu, Δv) on universal cover
        nu, nv = edge_winding_class(a, b)
        du = nu / 7.0
        dv = nv / 7.0
        # Interpolate from (u_curr, v_curr) to (u_curr + du, v_curr + dv)
        u_targ = u_curr + du
        v_targ = v_curr + dv
        ts = np.linspace(0, 1, n_per_edge, endpoint=False)
        for t in ts:
            u = u_curr + t * du
            v = v_curr + t * dv
            x = (R + r * np.cos(2 * np.pi * v)) * np.cos(2 * np.pi * u)
            y = (R + r * np.cos(2 * np.pi * v)) * np.sin(2 * np.pi * u)
            z = r * np.sin(2 * np.pi * v)
            pts.append((x, y, z))
        u_curr, v_curr = u_targ, v_targ
    # close the curve by appending the start point
    u = HEFFTER_VERT_UV[walk[0]][0]
    v = HEFFTER_VERT_UV[walk[0]][1]
    # Note: don't close exactly — pyknotid handles closed curves
    return np.array(pts)


def identify_walk_knot(walk: list[int], verbose: bool = False) -> dict:
    """Identify the carrier knot of a K_7 walk.

    Returns a dict with:
      - identifier (e.g., '3_1' for trefoil)
      - min_crossings (n_q candidate)
      - alexander, jones, determinant
      - is_unknot (whether the walk is topologically trivial)
    """
    curve = walk_to_3d_curve(walk)
    try:
        k = Knot(curve, verbose=verbose)
        # Identify against catalogue
        candidates = k.identify()
        if not candidates:
            return {"identifier": "?", "min_crossings": None,
                    "is_unknot": False, "n_crossings_diagram": None}
        # Prefer SIMPLEST candidate (smallest crossing number)
        best = min(candidates, key=lambda c: c.min_crossings
                    if hasattr(c, 'min_crossings') else 999)
        result = {
            "identifier": best.identifier,
            "min_crossings": best.min_crossings,
            "is_unknot": best.identifier == "0_1",
            "alexander": str(best.alexander)[:60] if best.alexander else "?",
            "jones": str(best.jones)[:60] if best.jones else "?",
            "determinant": best.determinant if best.determinant else None,
            "n_candidates": len(candidates),
        }
        return result
    except Exception as e:
        return {"identifier": f"ERROR: {e}", "min_crossings": None,
                "is_unknot": False, "error": str(e)}


def load_shortest_walks(max_length: int = 25) -> dict:
    """Phase E-3 BFS reload."""
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
        pp, qq = m_u // 7, m_v // 7
        key = (abs(pp), abs(qq))
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
    print("BOGOLIUBOV PHASE I — Jones polynomial closure of n_q via pyknotid")
    print("=" * 78)
    print()

    print("Step 1 — load shortest K_7 walks for compendium…")
    walks = load_shortest_walks(25)
    print(f"  {len(walks)} (|p|, |q|) classes")
    print()

    # ---- Identify carrier knot for each compendium particle's walk -----
    print("Step 2 — identify carrier knot via pyknotid for each compendium particle.")
    print(f"  This may take a few minutes…")
    print()

    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        knot_info = identify_walk_knot(walk, verbose=False)
        n_q_obs = entry["n_q"]
        n_q_pred = knot_info["min_crossings"]
        # n_q in Paper 11 conventions: 0/1 = unknot (lepton), 2 = Hopf (meson),
        # 3 = trefoil, 5 = cinquefoil (nucleon)
        # The carrier knot's min_crossings should map to n_q (note that n_q=2
        # Paper 11 uses for Hopf link which is a 2-component link not a knot)
        rows.append({
            **entry,
            "walk": walk,
            "knot_id": knot_info["identifier"],
            "knot_min_crossings": n_q_pred,
            "knot_det": knot_info.get("determinant"),
            "n_q_obs": n_q_obs,
            "n_q_pred": n_q_pred,
            "match": n_q_pred == n_q_obs if n_q_pred is not None else False,
        })
        print(f"  {entry['name']:<12} (p,q)=({entry['p']},{entry['q']}) "
              f"n_q_obs={n_q_obs:<2} L={len(walk)-1:<2} "
              f"knot={knot_info['identifier']:<8} "
              f"min_crossings={n_q_pred} "
              f"{'★ MATCH' if n_q_pred == n_q_obs else ''}")
    print()

    # ---- Summary table -------------------------------------------------
    print("=" * 78)
    print("SUMMARY — pyknotid carrier-knot identification vs Paper 11 n_q")
    print("=" * 78)
    print()
    matches = sum(1 for r in rows if r["match"])
    print(f"  Direct n_q matches: {matches}/{len(rows)} "
          f"({matches/len(rows)*100:.0f}%)")
    print()
    # Group by carrier knot
    from collections import defaultdict
    by_knot = defaultdict(list)
    for r in rows:
        by_knot[r["knot_id"]].append(r)
    print(f"  Walks grouped by carrier knot:")
    for kid, ps in sorted(by_knot.items()):
        names = [p["name"] for p in ps]
        n_q_set = sorted(set(p["n_q_obs"] for p in ps))
        cr = ps[0]["knot_min_crossings"]
        print(f"    {kid:<10} min_crossings={cr} → "
              f"compendium n_q values: {n_q_set}  "
              f"({len(ps)} particles): {', '.join(names[:5])}"
              f"{' ...' if len(names) > 5 else ''}")
    print()

    # ---- Predict n_q for the 4 substrate predictions --------------------
    print("=" * 78)
    print("STEP 3 — predict n_q for NWT's 4 substrate predictions")
    print("=" * 78)
    print()
    PREDICTIONS = {
        (2, 2): [0, 1, 2, 5, 1, 4, 0],
        (2, 3): [0, 1, 2, 3, 4, 0, 1, 4, 0],
        (3, 1): [0, 2, 4, 0, 2, 5, 1, 4, 0],
        (3, 3): [0, 1, 2, 3, 6, 2, 5, 1, 4, 0],
    }
    prediction_results = {}
    for pq, walk in PREDICTIONS.items():
        knot_info = identify_walk_knot(walk, verbose=False)
        prediction_results[pq] = knot_info
        print(f"  ({pq[0]},{pq[1]}) walk {'→'.join(str(v) for v in walk)} (L={len(walk)-1})")
        print(f"    carrier knot: {knot_info['identifier']}, "
              f"min_crossings = {knot_info['min_crossings']}")
        print(f"    → predicted n_q = {knot_info['min_crossings']}")
        print()

    # ---- Plot ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) n_q observed vs predicted
    ax = axes[0]
    obs = [r["n_q_obs"] for r in rows]
    pred = [r["n_q_pred"] for r in rows if r["n_q_pred"] is not None]
    obs_clean = [r["n_q_obs"] for r in rows if r["n_q_pred"] is not None]
    ax.scatter(obs_clean, pred, s=150, c=['C2' if o == p else 'C3'
                                            for o, p in zip(obs_clean, pred)],
               alpha=0.7, edgecolor='k')
    for r in rows:
        if r["n_q_pred"] is not None:
            ax.annotate(r["name"], (r["n_q_obs"], r["n_q_pred"]),
                        xytext=(4, 4), textcoords='offset points',
                        fontsize=7, alpha=0.7)
    ax.plot([-0.5, 6], [-0.5, 6], 'k--', alpha=0.4, label='1:1')
    ax.set_xlabel('n_q observed (Paper 11)')
    ax.set_ylabel('n_q predicted (pyknotid carrier crossings)')
    ax.set_title(f'Jones polynomial closure of n_q\n'
                 f'{matches}/{len(rows)} ({matches/len(rows)*100:.0f}%) direct match')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xticks(range(7))
    ax.set_yticks(range(7))

    # (b) Knot distribution
    ax = axes[1]
    knot_counts = {kid: len(ps) for kid, ps in by_knot.items()}
    sorted_knots = sorted(knot_counts.items(),
                          key=lambda x: -x[1])
    labels = [k for k, _ in sorted_knots]
    counts = [c for _, c in sorted_knots]
    ax.barh(range(len(labels)), counts, color='C0')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel('# compendium particles')
    ax.set_title('Carrier-knot distribution\n(pyknotid identification)')
    ax.grid(alpha=0.3, axis='x')
    ax.invert_yaxis()

    fig.suptitle(
        f"Phase I — n_q from pyknotid carrier-knot identification\n"
        f"Test of NWT's hypothesis: n_q = crossing number of carrier",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_i_jones_n_q.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    # Save
    np.savez(OUT_DIR / "phase_i_jones_n_q.npz",
             names=np.array([r["name"] for r in rows]),
             p=np.array([r["p"] for r in rows]),
             q=np.array([r["q"] for r in rows]),
             n_q_obs=np.array([r["n_q_obs"] for r in rows]),
             n_q_pred=np.array([r["n_q_pred"] if r["n_q_pred"] is not None
                                  else -1 for r in rows]),
             knot_id=np.array([r["knot_id"] for r in rows]),
             matches=np.array([r["match"] for r in rows]),
             predictions_pq=np.array(list(PREDICTIONS.keys())),
             predictions_knot_id=np.array(
                [prediction_results[pq]["identifier"]
                 for pq in PREDICTIONS]),
             predictions_n_q=np.array(
                [prediction_results[pq]["min_crossings"]
                 if prediction_results[pq]["min_crossings"] is not None
                 else -1 for pq in PREDICTIONS]),
             match_count=matches, total=len(rows))
    print(f"  data saved {OUT_DIR / 'phase_i_jones_n_q.npz'}")

    print()
    print("=" * 78)
    print("HEADLINE")
    print("=" * 78)
    print()
    print(f"  pyknotid identifies carrier knots from K_7 walks lifted to")
    print(f"  3D space curves on the Heffter torus.")
    print(f"  Direct n_q matches: {matches}/{len(rows)} "
          f"({matches/len(rows)*100:.0f}%)")
    print()
    if matches == len(rows):
        print(f"  ★★★ COMPLETE CLOSURE: n_q = carrier-knot crossing number")
        print(f"  for ALL compendium particles.")
    elif matches > len(rows) * 0.7:
        print(f"  ★ Strong match — NWT's Q5 hypothesis largely confirmed.")
    else:
        print(f"  Partial match — carrier-knot identification needs refinement.")


if __name__ == "__main__":
    main()
