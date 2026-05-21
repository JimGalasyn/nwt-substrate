"""Q4 — bridge from σ-orbit signature to mass.

Q1 established that σ-orbit signature is the injective substrate
invariant on (|p|,|q|) classes (16/16 distinct). Paper 6's mass
formula
    m_particle / m_e = (p² + q²)/5 · (β/β_e · ln(8β)/ln(8β_e)) · n_q^q
uses phenomenological (p, q, m, n_q). This script tests:

  (A) Can σ-sig features recover (p, q, m_min, n_q) — the inputs to
      Paper 6?
  (B) Does a clean closed-form mass prediction exist in σ-sig
      coordinates directly, replacing the (p,q,m,n_q) labels?
  (C) Within-multiplet residual — when multiple particles share a
      σ-sig (i.e. share a walk), how does m vary with carrier
      excitation?

Honest framing: σ-sig is a property of the *walk*; particles within
a multiplet (same walk) have different m. So σ-sig → m cannot be a
function on multiplet members; it can at best fix m_min for each
(|p|, |q|) class and leave the ladder Δm = {1, 2, 3} structure to a
separate carrier-excitation invariant.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q4_sigma_sig_to_mass.py
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.walks import edge_to_orbit
from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.particles.mass import paper6_mass_ratio, ME_MEV


def sigma_orbit_signature(walk: list[int]) -> tuple[int, ...]:
    counts = [0] * 7
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        oid = edge_to_orbit(a, b)
        if 0 <= oid < 7:
            counts[oid] += 1
    return tuple(counts)


def derived_features(sig: tuple[int, ...]) -> dict[str, float]:
    n0, n1, n2, n3, n4, n5, n6 = sig
    L = sum(sig)
    return {
        "L":              L,
        "n_polar":        n0 + n1,        # bridging edges (orbits 0, 1)
        "n_polar_avg":    (n0 + n1) / 2,
        "n_polar_asym":   n0 - n1,        # Z_2 polarity asymmetry
        "n_cross":        n2,             # τ-fixed cross-edges
        "n_triangle":     n3 + n4,        # equatorial triangles
        "n_triangle_asym": n3 - n4,       # Z_2 triangle asymmetry
        "n_twist":        n5 + n6,        # cross-twist total
        "n_twist_asym":   n5 - n6,        # Z_3 chirality
        "n_equatorial":   n3 + n4 + n2,   # all-equatorial (no P)
        "n_bridging":     n0 + n1 + n5 + n6,  # all edges through P or cross-twist
    }


def regress(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """Least-squares fit y = X β + intercept. Returns (β, R², residuals)."""
    A = np.hstack([X, np.ones((X.shape[0], 1))])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return beta, r2, y - y_pred


def main():
    print("=" * 78)
    print("Q4 — σ-sig → mass closed-form bridge")
    print("=" * 78)
    print()

    walks = bfs_shortest_walks(max_length=25)

    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        sig = sigma_orbit_signature(walk)
        features = derived_features(sig)
        paper6 = paper6_mass_ratio(entry["p"], entry["q"], entry["m"], entry["n_q"])
        rows.append({
            **entry,
            "key": key,
            "sig": sig,
            "L": sum(sig),
            "features": features,
            "m_obs_MeV": entry["m_obs"],
            "paper6_MeV": paper6 * ME_MEV if paper6 else None,
        })

    # -- Build per-(|p|,|q|)-class summary ------------------------------------
    by_pq = defaultdict(list)
    for r in rows:
        by_pq[r["key"]].append(r)

    print(f"{'class':<8} {'n':<2} {'L':<3} {'σ-sig':<22} "
          f"{'particles':<28} {'m range':<12} {'m_obs range (MeV)'}")
    print("-" * 110)
    classes = []
    for key in sorted(by_pq.keys()):
        entries = by_pq[key]
        ms = [e["m"] for e in entries]
        mobs = [e["m_obs"] for e in entries]
        sig_str = '[' + ' '.join(str(c) for c in entries[0]["sig"]) + ']'
        names = ','.join(e["name"] for e in entries)
        print(f"({key[0]:>2},{key[1]:>2})  {len(entries):<2} {entries[0]['L']:<3} "
              f"{sig_str:<22} {names:<28} "
              f"{min(ms):>2}–{max(ms):<2}      "
              f"{min(mobs):>7.2f}–{max(mobs):<7.2f}")
        classes.append({
            "key": key,
            "sig": entries[0]["sig"],
            "L": entries[0]["L"],
            "features": entries[0]["features"],
            "n_q": entries[0]["n_q"],
            "m_min": min(ms),
            "m_max": max(ms),
            "m_obs_min": min(mobs),
            "m_obs_max": max(mobs),
            "p": key[0], "q": key[1],
        })
    print()

    # -- Hypothesis A: σ-sig predicts m_min (composition-index baseline) ----
    print("=" * 78)
    print("Hypothesis A: σ-sig predicts m_min per (|p|,|q|) class")
    print("=" * 78)
    print()

    feature_names = list(classes[0]["features"].keys())
    X_full = np.array([[c["features"][f] for f in feature_names] for c in classes])
    y_mmin = np.array([c["m_min"] for c in classes])

    # Try (1) walk length L alone, (2) p+q only, (3) (p, q) only, (4) σ-sig only,
    # (5) σ-sig + n_q.
    cases = [
        ("L only",              ["L"]),
        ("polar+triangle+cross+twist+L",
                                ["L", "n_polar", "n_triangle", "n_cross", "n_twist"]),
        ("Z_3-resolved (with chirality)",
                                ["L", "n_polar", "n_triangle", "n_cross",
                                  "n_twist", "n_twist_asym"]),
        ("σ-sig (all 7 components)",
                                ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]),
        ("σ-sig + n_q",
                                ["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n_q"]),
    ]
    for label, feats in cases:
        if feats[0].startswith("n") and feats[0][1:].isdigit():
            # raw σ-sig
            X = np.array([list(c["sig"]) for c in classes])
            if "n_q" in feats:
                X = np.hstack([X, np.array([[c["n_q"]] for c in classes])])
        else:
            X = np.array([[c["features"][f] for f in feats] for c in classes])
        beta, r2, res = regress(X, y_mmin)
        print(f"  {label:<40}  R² = {r2:.4f}   "
              f"|max-res| = {np.max(np.abs(res)):.3f}   "
              f"RMS = {np.sqrt(np.mean(res**2)):.3f}")
    print()

    # -- Hypothesis A.1: closed-form m_min from L + structural features ----
    print("Compact 2-feature attempt:")
    for feats in [["L"], ["L", "n_q"], ["L", "n_polar"], ["L", "n_triangle"],
                   ["L", "n_polar", "n_q"], ["p+q (sum)"], ["L_minus_pq"]]:
        if feats == ["p+q (sum)"]:
            X = np.array([[c["p"] + c["q"]] for c in classes])
        elif feats == ["L_minus_pq"]:
            X = np.array([[c["L"] - c["p"] - c["q"]] for c in classes])
        else:
            X = np.array([[c["features"][f] if f in c["features"] else c["n_q"]
                            for f in feats] for c in classes])
        beta, r2, res = regress(X, y_mmin)
        print(f"  {str(feats):<40}  R² = {r2:.4f}   RMS = {np.sqrt(np.mean(res**2)):.3f}")
    print()

    # -- Hypothesis B: σ-sig predicts log(m_obs_min/m_e) ---------------------
    print("=" * 78)
    print("Hypothesis B: σ-sig predicts log(m_obs/m_e) per class")
    print("=" * 78)
    print()
    log_mobs = np.array([math.log(c["m_obs_min"] / ME_MEV) for c in classes])

    test_feature_sets = [
        ("L",                ["L"]),
        ("L, n_polar",       ["L", "n_polar"]),
        ("L, n_triangle",    ["L", "n_triangle"]),
        ("L, n_polar, n_triangle, n_cross, n_twist",
                              ["L", "n_polar", "n_triangle", "n_cross", "n_twist"]),
        ("Z_3-resolved",     ["L", "n_polar", "n_triangle", "n_cross",
                               "n_twist", "n_twist_asym"]),
        ("full σ-sig",       None),
        ("full σ-sig + n_q", "with_nq"),
    ]
    for label, feats in test_feature_sets:
        if feats is None:
            X = np.array([list(c["sig"]) for c in classes])
        elif feats == "with_nq":
            X = np.array([list(c["sig"]) + [c["n_q"]] for c in classes])
        else:
            X = np.array([[c["features"][f] for f in feats] for c in classes])
        beta, r2, res = regress(X, log_mobs)
        print(f"  {label:<48}  R² = {r2:.4f}   "
              f"RMS(log) = {np.sqrt(np.mean(res**2)):.3f}")
    print()

    # -- Hypothesis C: m_carrier := m - (|p|+|q|) -- a function of σ-sig? ---
    print("=" * 78)
    print("Hypothesis C: m_carrier = m - (|p|+|q|) closed-form from σ-sig?")
    print("=" * 78)
    print()
    print("(Per [[m-derivation-phase-i-exploration]]: leptons have m_carrier=0;")
    print(" m = m_walk + m_carrier where m_walk = |p|+|q|.)")
    print()
    y_carrier_min = np.array([c["m_min"] - c["p"] - c["q"] for c in classes])
    y_carrier_max = np.array([c["m_max"] - c["p"] - c["q"] for c in classes])

    print(f"{'class':<8} {'L':<3} {'n_q':<3} {'m_min':<5} {'m_walk':<6} "
           f"{'m_carrier_min':<14} {'σ-sig'}")
    print("-" * 78)
    for c in classes:
        sig_str = '[' + ' '.join(str(x) for x in c["sig"]) + ']'
        m_walk = c["p"] + c["q"]
        m_carr = c["m_min"] - m_walk
        print(f"({c['p']:>2},{c['q']:>2})  {c['L']:<3} {c['n_q']:<3} "
              f"{c['m_min']:<5} {m_walk:<6} "
              f"{m_carr:<14} {sig_str}")
    print()

    for label, feats in [("L", ["L"]),
                          ("L + n_q", ["L", "n_q"]),
                          ("L, n_polar, n_triangle, n_cross, n_twist",
                           ["L", "n_polar", "n_triangle", "n_cross", "n_twist"]),
                          ("full σ-sig", None),
                          ("full σ-sig + n_q", "with_nq")]:
        if feats is None:
            X = np.array([list(c["sig"]) for c in classes])
        elif feats == "with_nq":
            X = np.array([list(c["sig"]) + [c["n_q"]] for c in classes])
        else:
            X = np.array([[c["features"][f] if f in c["features"] else c["n_q"]
                            for f in feats] for c in classes])
        beta, r2, res = regress(X, y_carrier_min)
        print(f"  {label:<48}  R²(m_carrier_min) = {r2:.4f}   "
              f"RMS = {np.sqrt(np.mean(res**2)):.3f}")
    print()

    # -- Hypothesis D: σ-sig → n_q closed-form ------------------------------
    print("=" * 78)
    print("Hypothesis D: σ-sig → n_q (carrier-knot crossing count)")
    print("=" * 78)
    print()
    y_nq = np.array([c["n_q"] for c in classes])
    for label, feats in [("L", ["L"]),
                          ("n_polar+n_triangle+n_cross+n_twist",
                           ["n_polar", "n_triangle", "n_cross", "n_twist"]),
                          ("full σ-sig", None)]:
        if feats is None:
            X = np.array([list(c["sig"]) for c in classes])
        else:
            X = np.array([[c["features"][f] for f in feats] for c in classes])
        beta, r2, res = regress(X, y_nq)
        print(f"  {label:<48}  R²(n_q) = {r2:.4f}   "
              f"RMS = {np.sqrt(np.mean(res**2)):.3f}")
    print()

    # -- Hypothesis E: pure σ-sig substitution into Paper 6 -----------------
    print("=" * 78)
    print("Hypothesis E: substitute σ-sig features into Paper 6 mass formula")
    print("=" * 78)
    print()
    print("Paper 6: m/m_e = (p²+q²)/5 · (β/β_e · ln(8β)/ln(8β_e)) · n_q^q")
    print("Test if σ-sig-derived p²+q² and m can reproduce Paper 6 numerics.")
    print()
    print(f"{'particle':<10} {'p²+q²':<6} {'paper6_MeV':<12} {'m_obs_MeV':<11}"
          f" {'% residual'}")
    print("-" * 64)
    n_good = 0
    n_total = 0
    for r in rows:
        if r["paper6_MeV"] is None:
            continue
        resid = (r["paper6_MeV"] - r["m_obs_MeV"]) / r["m_obs_MeV"] * 100
        ok = abs(resid) < 10
        n_good += int(ok)
        n_total += 1
        pq2 = r["p"] * r["p"] + r["q"] * r["q"]
        print(f"{r['name']:<10} {pq2:<6} {r['paper6_MeV']:<12.2f} "
              f"{r['m_obs_MeV']:<11.2f} {resid:+.2f}%")
    print()
    print(f"  Paper 6 reference: {n_good}/{n_total} within 10%")
    print()
    print("σ-sig substitution candidates for (p² + q²)/5 prefactor:")
    print(f"{'class':<10} {'p²+q²':<6} {'L':<3} {'n_polar':<8} "
           f"{'n_triangle':<11} {'L·(n_polar+n_triangle)/?'}")
    print("-" * 72)
    for c in classes:
        pq2 = c["p"]**2 + c["q"]**2
        npol = c["features"]["n_polar"]
        ntri = c["features"]["n_triangle"]
        # candidate proxy: L * (npol + ntri)
        proxy = c["L"] * (npol + ntri)
        print(f"({c['p']:>2},{c['q']:>2})    {pq2:<6} {c['L']:<3} {npol:<8} "
              f"{ntri:<11} {proxy}")
    print()


if __name__ == "__main__":
    main()
