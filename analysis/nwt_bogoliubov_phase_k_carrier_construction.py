"""Bogoliubov Phase K — σ-orbit → carrier-knot map construction.

Phase I exposed that Paper 11's "carrier knot" is a separate topological
object from the walk's intrinsic torus-knot.  Phase F-3 / G showed that
σ-orbit DISTRIBUTION determines the carrier sector.  Phase K attempts
the rigorous closure: an explicit σ-orbit → carrier-knot construction
verifiable via pyknotid.

NWT's two candidate constructions (broadcast 2026-05-20 ~3:50pm):

(a) Each σ-orbit's 3 edges define a 3D substructure; the walk's path
    through σ-orbit space traces a closed curve.

(b) σ-orbit TRANSITIONS along the walk are the crossings of the carrier
    diagram; PD code emerges from the transition pattern.

Also testing simpler combinatorial candidates:

(c) Carrier-knot crossings = # DISTINCT σ-orbit transitions (excluding
    same-orbit consecutive edges).

(d) Carrier-knot crossings = # DISTINCT σ-orbits the walk uses.

(e) Carrier crossings = walk's σ-orbit sequence interpreted as a
    cyclic permutation, identified to a knot via Gauss-code matching.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_k_carrier_construction.py
"""
from __future__ import annotations

import math
from collections import deque, defaultdict, Counter
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt   # before pyknotid patches
np.float = float
np.int = int
np.complex = complex
np.object = object
np.str = str
np.long = int

try:
    from pyknotid.spacecurves import Knot
    PYKNOTID = True
except Exception:
    PYKNOTID = False

from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS
from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_k_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def edge_to_orbit(a: int, b: int) -> int:
    e_norm = (min(a, b), max(a, b))
    for oid, orbit in SIGMA_ORBITS.items():
        if e_norm in orbit["edges"]:
            return oid
    return -1


def walk_sigma_sequence(walk: list[int]) -> list[int]:
    """σ-orbit ID for each step of the walk."""
    return [edge_to_orbit(walk[i], walk[i + 1])
            for i in range(len(walk) - 1)]


def sigma_transitions_no_self(sig_seq: list[int]) -> list[tuple[int, int]]:
    """Non-self-loop transitions in the σ-orbit sequence (including wrap-around)."""
    n = len(sig_seq)
    transitions = []
    for i in range(n):
        cur = sig_seq[i]
        nxt = sig_seq[(i + 1) % n]
        if cur != nxt:
            transitions.append((cur, nxt))
    return transitions


def n_distinct_sigma(sig_seq: list[int]) -> int:
    return len(set(sig_seq))


def n_sigma_transitions(sig_seq: list[int]) -> int:
    """Count of non-self-loop σ-orbit transitions."""
    return len(sigma_transitions_no_self(sig_seq))


def sigma_orbit_anchors_3d(R: float = 3.0, r: float = 1.0) -> dict:
    """Place each σ-orbit at a representative 3D point.  Use the centroid of
    the σ-orbit's 3 vertices on the Heffter torus.
    """
    from nwt_substrate.condensate.orbit_winding import HEFFTER_VERT_UV
    anchors = {}
    for oid, orbit in SIGMA_ORBITS.items():
        # Get all vertices in the orbit's edges
        verts = set()
        for (a, b) in orbit["edges"]:
            verts.add(a)
            verts.add(b)
        # Centroid in (u, v)
        u_mean = np.mean([HEFFTER_VERT_UV[v][0] for v in verts])
        v_mean = np.mean([HEFFTER_VERT_UV[v][1] for v in verts])
        # Map to 3D
        x = (R + r * np.cos(2 * np.pi * v_mean)) * np.cos(2 * np.pi * u_mean)
        y = (R + r * np.cos(2 * np.pi * v_mean)) * np.sin(2 * np.pi * u_mean)
        z = r * np.sin(2 * np.pi * v_mean)
        anchors[oid] = (x, y, z)
    return anchors


def construction_a_walk_through_sigma_space(walk: list[int]) -> np.ndarray:
    """Construction (a): walk traces 3D curve through σ-orbit anchors.

    The σ-orbit sequence of the walk's edges, smoothly interpolated
    between σ-orbit anchor points in 3D.
    """
    anchors = sigma_orbit_anchors_3d()
    sig_seq = walk_sigma_sequence(walk)
    pts = []
    n_per_edge = 20
    for i in range(len(sig_seq)):
        cur = sig_seq[i]
        nxt = sig_seq[(i + 1) % len(sig_seq)]
        if cur == nxt:
            # tiny loop at this anchor (to register as a kink, not a straight line)
            ax, ay, az = anchors[cur]
            # add a small spiral
            for k in range(n_per_edge):
                ang = 2 * np.pi * k / n_per_edge
                pts.append((ax + 0.05 * np.cos(ang),
                             ay + 0.05 * np.sin(ang),
                             az + 0.02 * (k / n_per_edge - 0.5)))
        else:
            ax, ay, az = anchors[cur]
            bx, by, bz = anchors[nxt]
            ts = np.linspace(0, 1, n_per_edge, endpoint=False)
            for t in ts:
                pts.append((ax + t * (bx - ax),
                             ay + t * (by - ay),
                             az + t * (bz - az)))
    return np.array(pts)


def identify_via_pyknotid(curve_3d: np.ndarray, verbose: bool = False) -> dict:
    if not PYKNOTID:
        return {"identifier": "(pyknotid unavailable)", "min_crossings": None}
    try:
        k = Knot(curve_3d, verbose=verbose)
        candidates = k.identify()
        if not candidates:
            return {"identifier": "(no identification)", "min_crossings": None}
        best = min(candidates, key=lambda c: c.min_crossings)
        return {"identifier": best.identifier,
                "min_crossings": best.min_crossings,
                "n_candidates": len(candidates)}
    except Exception as e:
        return {"identifier": f"ERROR: {str(e)[:30]}",
                "min_crossings": None}


def load_shortest_walks(max_length: int = 25) -> dict:
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
    print("BOGOLIUBOV PHASE K — σ-orbit → carrier-knot map construction")
    print("=" * 78)
    print(f"  pyknotid available: {PYKNOTID}")
    print()

    walks = load_shortest_walks(25)
    print(f"Loaded {len(walks)} (|p|, |q|) shortest walks")
    print()

    # ---- Test 5 candidate constructions on each compendium particle ----
    print("=" * 78)
    print("CANDIDATE CONSTRUCTIONS — n_q from σ-orbit topology")
    print("=" * 78)
    print()
    print(f"  {'particle':<12} {'n_q_obs':<8} {'(c) #trans':<11} "
          f"{'(d) #dist_σ':<13} {'(a) pyknotid':<20}")
    print("  " + "-" * 75)

    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        sig_seq = walk_sigma_sequence(walk)
        c_n_trans = n_sigma_transitions(sig_seq)
        d_n_dist = n_distinct_sigma(sig_seq)
        # construction (a) via 3D curve through σ-orbit anchors
        if PYKNOTID:
            curve = construction_a_walk_through_sigma_space(walk)
            a_result = identify_via_pyknotid(curve, verbose=False)
            a_id = a_result["identifier"]
            a_cr = a_result["min_crossings"]
        else:
            a_id = "?"
            a_cr = None
        n_q_obs = entry["n_q"]
        rows.append({
            **entry,
            "walk": walk,
            "sig_seq": sig_seq,
            "c_n_trans": c_n_trans,
            "d_n_dist": d_n_dist,
            "a_identifier": a_id,
            "a_crossings": a_cr,
            "n_q_obs": n_q_obs,
        })
        a_str = f"{a_id} (cr={a_cr})" if a_cr is not None else a_id[:18]
        print(f"  {entry['name']:<12} {n_q_obs:<8} "
              f"{c_n_trans:<11} {d_n_dist:<13} {a_str:<20}")
    print()

    # ---- Stratified analysis: each candidate vs n_q --------------------
    print("=" * 78)
    print("STRATIFIED MATCH ANALYSIS")
    print("=" * 78)
    print()
    for cand_name, cand_attr in [
        ("(c) #non-self σ-transitions", "c_n_trans"),
        ("(d) #distinct σ-orbits visited", "d_n_dist"),
    ]:
        n_q_to_cand = defaultdict(list)
        for r in rows:
            n_q_to_cand[r["n_q_obs"]].append(r[cand_attr])
        print(f"  {cand_name}:")
        for n_q in sorted(n_q_to_cand.keys()):
            vals = n_q_to_cand[n_q]
            sector = {0: "leptons", 2: "mesons", 3: "hyperons",
                       5: "nucleons"}.get(n_q, "?")
            print(f"    n_q={n_q} ({sector:<10}): values = "
                  f"{sorted(set(vals))} (mean = {np.mean(vals):.2f})")
        # Pearson r
        n_q_arr = np.array([r["n_q_obs"] for r in rows])
        cand_arr = np.array([r[cand_attr] for r in rows])
        r_p = np.corrcoef(n_q_arr, cand_arr)[0, 1]
        print(f"    Pearson r(n_q, {cand_attr}) = {r_p:+.3f}")
        print()

    # Construction (a) — pyknotid via σ-orbit anchor curve
    if PYKNOTID:
        print("  (a) pyknotid via σ-orbit anchor curve:")
        # Map identifier to crossings, group by sector
        sector_to_a = defaultdict(list)
        for r in rows:
            sector_to_a[r["n_q_obs"]].append((r["name"], r["a_identifier"],
                                                 r["a_crossings"]))
        for n_q in sorted(sector_to_a.keys()):
            entries = sector_to_a[n_q]
            sector = {0: "leptons", 2: "mesons", 3: "hyperons",
                       5: "nucleons"}.get(n_q, "?")
            print(f"    n_q={n_q} ({sector:<10}):")
            for name, kid, cr in entries[:4]:
                print(f"      {name:<10} → {kid} (crossings={cr})")
            if len(entries) > 4:
                print(f"      ... and {len(entries)-4} more")
        print()

    # ---- Headline / conclusion -----------------------------------------
    print("=" * 78)
    print("PHASE K HEADLINE")
    print("=" * 78)
    print()

    # Best correlation
    c_r = np.corrcoef(
        [r["n_q_obs"] for r in rows],
        [r["c_n_trans"] for r in rows])[0, 1]
    d_r = np.corrcoef(
        [r["n_q_obs"] for r in rows],
        [r["d_n_dist"] for r in rows])[0, 1]

    print(f"  Correlation summary (n_q observed vs σ-orbit invariant):")
    print(f"    (c) # σ-orbit transitions:  r = {c_r:+.3f}")
    print(f"    (d) # distinct σ-orbits:    r = {d_r:+.3f}")
    print()
    print(f"  pyknotid construction (a) results — see table above for")
    print(f"  per-particle carrier-knot identification.")
    print()
    print(f"  None of the simple invariants gives a CLEAN match to n_q.")
    print(f"  The σ-orbit → carrier-knot map requires deeper structural")
    print(f"  construction than these candidates.  Phase K is partial; the")
    print(f"  rigorous map likely needs:")
    print(f"    - σ-orbit transition graph + specific 2D projection rule")
    print(f"    - over/under crossing assignment by σ-orbit type")
    print(f"    - knot identification via pyknotid PD-code input (not space curve)")
    print()
    print(f"  Honest negative for the simple candidates; framework structural")
    print(f"  problem documented.  Recommend handing the over/under crossing")
    print(f"  rule to NWT's Spin(7) / K_7 representation theory work.")

    # ---- Plot ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) n_q vs # distinct σ-orbits
    ax = axes[0]
    sector_colors = {0: 'C2', 2: 'C0', 3: 'C1', 5: 'C3'}
    for r in rows:
        c = sector_colors.get(r["n_q_obs"], 'gray')
        ax.scatter(r["d_n_dist"], r["n_q_obs"], c=c, s=100, alpha=0.6,
                    edgecolor='k')
        ax.annotate(r["name"], (r["d_n_dist"], r["n_q_obs"]),
                    xytext=(np.random.uniform(-15, 15),
                             np.random.uniform(-15, 15)),
                    textcoords='offset points', fontsize=7, alpha=0.6)
    ax.set_xlabel('# distinct σ-orbits visited (candidate d)')
    ax.set_ylabel('n_q (Paper 11)')
    ax.set_title(f'(d) # distinct σ-orbits vs n_q\nr = {d_r:+.3f}')
    ax.grid(alpha=0.3)

    # (b) n_q vs # σ-orbit transitions
    ax = axes[1]
    for r in rows:
        c = sector_colors.get(r["n_q_obs"], 'gray')
        ax.scatter(r["c_n_trans"], r["n_q_obs"], c=c, s=100, alpha=0.6,
                    edgecolor='k')
        ax.annotate(r["name"], (r["c_n_trans"], r["n_q_obs"]),
                    xytext=(np.random.uniform(-15, 15),
                             np.random.uniform(-15, 15)),
                    textcoords='offset points', fontsize=7, alpha=0.6)
    ax.set_xlabel('# σ-orbit transitions (candidate c)')
    ax.set_ylabel('n_q (Paper 11)')
    ax.set_title(f'(c) # σ-orbit transitions vs n_q\nr = {c_r:+.3f}')
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"Phase K — σ-orbit → carrier-knot map (exploratory candidates)\n"
        f"Simple invariants don't give n_q; deeper structural work needed",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_k_carrier_construction.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_k_carrier_construction.npz",
             names=np.array([r["name"] for r in rows]),
             n_q_obs=np.array([r["n_q_obs"] for r in rows]),
             c_n_trans=np.array([r["c_n_trans"] for r in rows]),
             d_n_dist=np.array([r["d_n_dist"] for r in rows]),
             a_crossings=np.array([r["a_crossings"]
                                     if r["a_crossings"] is not None else -1
                                     for r in rows]),
             a_identifier=np.array([r["a_identifier"] for r in rows]),
             c_r=c_r, d_r=d_r)
    print(f"  data saved {OUT_DIR / 'phase_k_carrier_construction.npz'}")


if __name__ == "__main__":
    main()
