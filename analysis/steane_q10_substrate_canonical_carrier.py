"""Q10 — substrate-canonical replacement for Paper 6's n_q^q carrier factor.

Paper 6 mass formula:
    m/m_e = (p² + q²)/5 · (β/β_e · ln(8β)/ln(8β_e)) · n_q^q

where β = sqrt(m²/p² - 1), and n_q ∈ {0, 2, 3, 5} is the carrier-knot
sector label assigned per particle by Paper 11. Reviewer criticism:
the n_q labels are partly empirical and don't have a clean substrate-
canonical derivation (confirmed by Q9 — so(3)_diag rep theory doesn't
recover them).

This script tests substrate-canonical OBSERVABLES that might replace
n_q^q as the carrier-enhancement factor. The carrier enhancement
isolated empirically per walk is:

    R_carrier := (m_obs / m_e) / ((p² + q²)/5 · β-factor)

We compute R_carrier per compendium walk and compare to candidates:

  Paper 6 baseline:   n_q^q   (=1 for n_q∈{0,1}, =2^q, 3^q, 5^q else)
  Candidate 1:         L^q  (walk length to q)
  Candidate 2:         (σ-orbits visited)^q
  Candidate 3:         (n_σ_polar)^q  (sum of σ_0, σ_1 edges)
  Candidate 4:         (n_σ_cross)^q  (sum of cross-block σ_4..σ_6)
  Candidate 5:         (frac_j=1)^(-1) · ... (continuous Casimir-based)
  Candidate 6:         (p + q)^q or (p²+q²)^(q/2) — walk geodesic
  Candidate 7:         (# QR steps)^q
  Candidate 8:         (longest QR-run)^q

For each candidate, compute residuals relative to the empirical
R_carrier and report fit quality vs. Paper 6 baseline.

The hope is to find a SINGLE candidate that is BOTH substrate-canonical
AND fits the compendium as well as n_q^q.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q10_substrate_canonical_carrier.py
"""
from __future__ import annotations

import math
from collections import defaultdict
from itertools import combinations, groupby

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.walks import edge_to_orbit
from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.particles.mass import (
    paper6_mass_ratio, ME_MEV, BETA_E, LN8_BE, beta_phase,
)


N = 7
QR = {1, 2, 4}
NR = {3, 5, 6}


def so7_basis():
    basis = []
    for i, j in combinations(range(N), 2):
        J = np.zeros((N, N))
        J[i, j] = +1.0
        J[j, i] = -1.0
        basis.append(J)
    return basis


SO7_BASIS = so7_basis()
EDGE_INDEX = {}
for k, (i, j) in enumerate(combinations(range(N), 2)):
    EDGE_INDEX[(i, j)] = k
    EDGE_INDEX[(j, i)] = k


def so3_diag_generators():
    def gen(i, j):
        M = np.zeros((N, N))
        M[i, j] = 1.0
        M[j, i] = -1.0
        return M
    L_x = gen(2, 3)
    L_y = gen(1, 3)
    L_z = gen(1, 2)
    L_xp = gen(5, 6)
    L_yp = gen(4, 6)
    L_zp = gen(4, 5)
    return L_x + L_xp, L_y + L_yp, L_z + L_zp


def adjoint_action(J, basis):
    n = len(basis)
    M = np.zeros((n, n))
    for l in range(n):
        comm = J @ basis[l] - basis[l] @ J
        for k in range(n):
            M[k, l] = -0.5 * np.trace(basis[k] @ comm)
    return M


def walk_to_algebra_coords(walk):
    coords = np.zeros(len(SO7_BASIS))
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        k = EDGE_INDEX[(a, b)]
        sign = +1 if a < b else -1
        coords[k] += sign
    return coords


def sigma_signature(walk):
    counts = [0] * 7
    for i in range(len(walk) - 1):
        oid = edge_to_orbit(walk[i], walk[i + 1])
        if 0 <= oid < 7:
            counts[oid] += 1
    return tuple(counts)


def walk_features(walk):
    L = len(walk) - 1
    sig = sigma_signature(walk)
    d_seq = [(walk[i+1] - walk[i]) % 7 for i in range(L)]
    n_qr = sum(1 for d in d_seq if d in QR)
    n_nr = sum(1 for d in d_seq if d in NR)
    qr_runs = [sum(1 for _ in g) for k, g in groupby(d_seq, key=lambda d: d in QR) if k]
    nr_runs = [sum(1 for _ in g) for k, g in groupby(d_seq, key=lambda d: d in NR) if k]
    return {
        "L": L,
        "sig": sig,
        "n_polar": sig[0] + sig[1],
        "n_intra": sig[2] + sig[3],
        "n_cross": sig[4] + sig[5] + sig[6],
        "n_qr": n_qr,
        "n_nr": n_nr,
        "max_qr_run": max(qr_runs) if qr_runs else 0,
        "max_nr_run": max(nr_runs) if nr_runs else 0,
        "n_orbits_visited": sum(1 for c in sig if c > 0),
        "n_vertices_visited": len(set(walk)),
    }


def compute_j_weights():
    """Build so(3)_diag Casimir projectors and return a function
    j_weights(walk) → (frac_j0, frac_j1, frac_j2)."""
    Jx, Jy, Jz = so3_diag_generators()
    adx = adjoint_action(Jx, SO7_BASIS)
    ady = adjoint_action(Jy, SO7_BASIS)
    adz = adjoint_action(Jz, SO7_BASIS)
    C = adx @ adx + ady @ ady + adz @ adz
    C = 0.5 * (C + C.T)
    eigvals, eigvecs = np.linalg.eigh(C)
    groups = defaultdict(list)
    for k, ev in enumerate(eigvals):
        groups[round(ev, 4)].append(k)
    projectors = {}
    for ev, idx in groups.items():
        P = np.zeros_like(C)
        for k in idx:
            v = eigvecs[:, k:k+1]
            P += v @ v.T
        projectors[ev] = P

    def jw(walk):
        A = walk_to_algebra_coords(walk)
        tot = np.linalg.norm(A) ** 2
        weights = {}
        for ev, P in projectors.items():
            weights[ev] = float(np.linalg.norm(P @ A) ** 2)
        if tot == 0:
            return 0.0, 0.0, 0.0
        # Find which key corresponds to j=0, 1, 2
        j0_key = min(projectors.keys(), key=abs)  # nearest 0
        j2_key = min(projectors.keys(), key=lambda k: abs(k - (-6)))
        j1_key = min(projectors.keys(), key=lambda k: abs(k - (-2)))
        return (weights[j0_key] / tot,
                weights[j1_key] / tot,
                weights[j2_key] / tot)

    return jw


def main():
    print("=" * 78)
    print("Q10 — substrate-canonical replacement for n_q^q carrier factor")
    print("=" * 78)
    print()

    walks_dict = bfs_shortest_walks(max_length=25)
    j_weights = compute_j_weights()

    # ---- Compute R_carrier empirically per particle ---------------------
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks_dict:
            continue
        walk = walks_dict[key]
        p, q, m, nq = entry["p"], entry["q"], entry["m"], entry["n_q"]
        b = beta_phase(p, m)
        if b is None or b <= 0:
            continue
        # Paper 6 factors
        torus_gradient = (p * p + q * q) / 5.0
        beta_factor = (b / BETA_E) * (math.log(8 * b) / LN8_BE)
        # R_carrier = (m_obs/m_e) / (torus_gradient · beta_factor)
        m_e_ratio = entry["m_obs"] / ME_MEV
        R_emp = m_e_ratio / (torus_gradient * beta_factor)
        # Paper 6's prediction: n_q^q (with n_q ≤ 1 → 1)
        R_paper6 = 1.0 if nq <= 1 else nq ** q
        # Walk features
        feat = walk_features(walk)
        j0, j1, j2 = j_weights(walk)
        rows.append({
            **entry, "p": p, "q": q,
            "walk": walk,
            "m_e_ratio": m_e_ratio,
            "R_emp": R_emp,
            "R_paper6": R_paper6,
            "torus_grad": torus_gradient,
            "beta_factor": beta_factor,
            **feat,
            "j0": j0, "j1": j1, "j2": j2,
        })

    # ---- Print baseline table -------------------------------------------
    print("Paper 6 baseline check:")
    print(f"{'particle':<10} {'(p,q,m,nq)':<14} {'R_emp':<10} {'R_paper6':<10} "
          f"{'ratio':<8} {'log resid'}")
    print("-" * 78)
    log_resids = []
    for r in rows:
        rp = r['R_emp'] / r['R_paper6'] if r['R_paper6'] > 0 else 0
        lr = math.log(rp) if rp > 0 else float('-inf')
        log_resids.append((r['name'], lr, r))
        print(f"{r['name']:<10} ({r['p']:>2},{r['q']:>2},{r['m']:>2},{r['n_q']:<2})  "
              f"{r['R_emp']:<10.3f} {r['R_paper6']:<10.3f} {rp:<8.4f} {lr:+.4f}")
    rms_p6 = math.sqrt(sum(lr * lr for _, lr, _ in log_resids) / len(log_resids))
    print(f"\n  Paper 6 RMS(log residual) = {rms_p6:.4f} "
          f"(= log of {math.exp(rms_p6):.2f}x error)")
    print()

    # ---- Candidate substrate-canonical replacements ---------------------
    # For each candidate c(walk), test R_predicted = c(walk)^q vs R_emp.
    print("=" * 78)
    print("Candidate substrate-canonical replacements:  R_pred = c(walk)^q")
    print("=" * 78)
    print()

    candidates = [
        ("L (walk length)",           lambda r: r['L']),
        ("L / 7",                     lambda r: r['L'] / 7.0),
        ("n_orbits_visited",          lambda r: r['n_orbits_visited']),
        ("n_polar+1",                 lambda r: r['n_polar'] + 1),
        ("n_intra+1",                 lambda r: r['n_intra'] + 1),
        ("n_cross+1",                 lambda r: r['n_cross'] + 1),
        ("max_qr_run",                lambda r: r['max_qr_run']),
        ("max_nr_run+1",              lambda r: r['max_nr_run'] + 1),
        ("n_qr+1",                    lambda r: r['n_qr'] + 1),
        ("n_vertices_visited",        lambda r: r['n_vertices_visited']),
        ("1+j0_weight·L",             lambda r: 1 + r['j0'] * r['L']),
        ("1+j2_weight·L",             lambda r: 1 + r['j2'] * r['L']),
        ("L/(p+q)",                   lambda r: r['L'] / (r['p'] + r['q'])),
        ("L/q",                       lambda r: r['L'] / r['q']),
        ("(p+q)/q",                   lambda r: (r['p'] + r['q']) / r['q']),
        ("p+q",                       lambda r: r['p'] + r['q']),
        ("max(p,q)",                  lambda r: max(r['p'], r['q'])),
        ("min(p,q)+2",                lambda r: min(r['p'], r['q']) + 2),
        ("(p²+q²)/q²",                lambda r: (r['p']**2 + r['q']**2) / r['q']**2),
    ]

    results = []
    for name, fn in candidates:
        log_resids_c = []
        for r in rows:
            base = fn(r)
            if base <= 0:
                continue
            R_pred = base ** r['q']
            if R_pred <= 0 or r['R_emp'] <= 0:
                continue
            lr = math.log(r['R_emp'] / R_pred)
            log_resids_c.append(lr)
        if log_resids_c:
            rms = math.sqrt(sum(lr * lr for lr in log_resids_c) / len(log_resids_c))
            mean = sum(log_resids_c) / len(log_resids_c)
            results.append((name, rms, mean, len(log_resids_c)))

    # Sort by RMS
    results.sort(key=lambda x: x[1])
    print(f"  {'candidate':<25} {'RMS(log)':<10} {'mean(log)':<10} "
          f"{'n':<3} {'~ x-error'}")
    print("  " + "-" * 70)
    print(f"  {'PAPER 6 BASELINE (n_q^q)':<25} {rms_p6:<10.4f} {'-':<10} {len(rows):<3} "
          f"{math.exp(rms_p6):.2f}x")
    print("  " + "-" * 70)
    for name, rms, mean, n in results:
        print(f"  {name:<25} {rms:<10.4f} {mean:<+10.4f} {n:<3} "
              f"{math.exp(rms):.2f}x")
    print()

    # ---- Affine fit: log R_emp = α + β · log(candidate^q) + γ · ... -----
    print("=" * 78)
    print("Affine fit: log R_emp ≈ α + β·log(candidate)·q  per candidate")
    print("=" * 78)
    print()
    print(f"  {'candidate':<25} {'α':<8} {'β':<8} {'RMS(log)':<10} {'~ x-error'}")
    print("  " + "-" * 70)
    log_R = np.array([math.log(r['R_emp']) for r in rows])
    for name, fn in candidates:
        x_vals = []
        y_vals = []
        for r in rows:
            base = fn(r)
            if base <= 0 or r['R_emp'] <= 0:
                continue
            x_vals.append(math.log(base) * r['q'])
            y_vals.append(math.log(r['R_emp']))
        if len(x_vals) < 2:
            continue
        x = np.array(x_vals)
        y = np.array(y_vals)
        X = np.column_stack([x, np.ones_like(x)])
        beta_coefs, *_ = np.linalg.lstsq(X, y, rcond=None)
        beta, alpha = beta_coefs
        resid = y - (alpha + beta * x)
        rms = float(np.sqrt(np.mean(resid * resid)))
        print(f"  {name:<25} {alpha:<+8.4f} {beta:<+8.4f} {rms:<10.4f} "
              f"{math.exp(rms):.2f}x")
    print()


if __name__ == "__main__":
    main()
