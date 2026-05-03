"""
Phase 3 of the nuclear-binding refinement: two-zone alpha-cluster
model that handles the breakdown of pure LJ at n_alpha >= 9.

Phase 2 finding (memory:connected-sum-composition-law):
  Pure LJ-cluster minima predict empirical N_nn at <0.3% residual
  for n_alpha = 3..8 (light cluster regime), but progressively
  over-predict for n_alpha >= 9: at 56Ni LJ predicts 45 NN where
  data gives 36; at 80Zr it's 72 vs 43.

Two physical effects suppress LJ-cluster coordination at n_alpha >= 9:

  (1) PAULI BLOCKING in medium clusters: when alpha's pack densely,
      their nucleons overlap; Pauli exclusion limits the average
      coordination per alpha.  Empirical N_nn / n_alpha peaks at
      n_alpha = 14 (56Ni) at 2.59, giving z_max = 5.18 effective
      maximum NN bonds per alpha.

  (2) DOUBLY-MAGIC SHELL CLOSURE at Z = N = 28 (n_alpha = 14, 56Ni):
      additional alphas beyond shell closure populate a new outer
      shell with much-reduced coupling.  Empirical Delta(N_nn) per
      added alpha drops from ~3 (n_alpha = 12..14, intra-core) to
      ~1.1 (n_alpha = 15..21, valence shell).

Phase 3 model:

  N_nn_pred(n_alpha) =
      LJ_min(n_alpha)                              n_alpha <= 8
      min(LJ_min(n_alpha), z_max * n_alpha / 2)    9 <= n_alpha <= 14
      N_nn_pred(14) + (n_alpha - 14) * nu_valence  n_alpha > 14

Inputs:
  E_link*       = 2.42 MeV  (per-NN-bond binding, from 12C/16O)
  n_core        = 14        (doubly-magic 56Ni shell closure;
                             standard nuclear shell model -- not fitted)
  z_max         = 5.14      (max average NN bonds per alpha,
                             extracted from 56Ni: 36.3 / 14 * 2)
  nu_valence    = 1.10      (NN bonds per added valence alpha,
                             extracted from 80Zr regime average)

The shell-closure structure (n_core = 14) is taken from the standard
nuclear shell model -- it is the Z = N = 28 doubly-magic point of the
nuclear chart, not a fitted parameter.  z_max and nu_valence are
empirical refinements that quantify Pauli + shell suppression.
"""
from __future__ import annotations
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize


B_E_ALPHA   = 28.296   # MeV   4He binding
E_LINK_STAR = 2.42     # MeV   per-NN-bond, from 12C/16O
N_CORE      = 14       # 56Ni doubly-magic shell closure (Z=N=28)
Z_MAX       = 5.14     # max avg NN per alpha; from 56Ni
NU_VALENCE  = 1.10     # NN per valence alpha beyond shell closure

R0          = 1.0
EPS         = 1.0
NN_CUTOFF   = 1.3 * R0


# ---------------------------------------------------------------------------
# LJ pair potential and cluster minimisation (from Phase 2)
# ---------------------------------------------------------------------------

def lj_potential(r: np.ndarray) -> np.ndarray:
    sr = R0 / r
    return EPS * (sr ** 12 - 2.0 * sr ** 6)


def cluster_energy(coords_flat: np.ndarray, n: int) -> float:
    coords = coords_flat.reshape(n, 3)
    diff = coords[:, None, :] - coords[None, :, :]
    r = np.sqrt(np.sum(diff ** 2, axis=-1))
    iu = np.triu_indices(n, k=1)
    rij = np.maximum(r[iu], 0.4)
    return float(np.sum(lj_potential(rij)))


def cluster_gradient(coords_flat: np.ndarray, n: int) -> np.ndarray:
    coords = coords_flat.reshape(n, 3)
    grad = np.zeros_like(coords)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            d = coords[i] - coords[j]
            r = max(np.linalg.norm(d), 0.4)
            sr = R0 / r
            dVdr = -12.0 * EPS / r * (sr ** 12 - sr ** 6)
            grad[i] += dVdr * d / r
    return grad.flatten()


def lj_min_nn(n: int, n_starts: int = 80, seed: int = 12345) -> int:
    """LJ-cluster minimum nearest-neighbour count."""
    if n < 2: return 0
    if n == 2: return 1
    rng = np.random.default_rng(seed + n)
    best_E = np.inf
    best = None
    for k in range(n_starts):
        radius = max(0.6, 0.55 * n ** (1.0 / 3.0))
        x0 = rng.normal(0.0, radius, size=(n, 3))
        try:
            res = minimize(cluster_energy, x0.flatten(), args=(n,),
                            jac=cluster_gradient, method="L-BFGS-B",
                            options={"gtol": 1e-7, "maxiter": 800})
            if res.fun < best_E - 1e-6:
                best_E = float(res.fun)
                best = res.x.reshape(n, 3)
        except Exception:
            continue
    if best is None: return 0
    diff = best[:, None, :] - best[None, :, :]
    r = np.sqrt(np.sum(diff ** 2, axis=-1))
    iu = np.triu_indices(n, k=1)
    return int(np.sum(r[iu] < NN_CUTOFF))


# ---------------------------------------------------------------------------
# Phase 3 prediction
# ---------------------------------------------------------------------------

# Cache LJ minima (computed once)
_LJ_CACHE: dict[int, int] = {}


def lj_cached(n: int) -> int:
    if n not in _LJ_CACHE:
        _LJ_CACHE[n] = lj_min_nn(n)
    return _LJ_CACHE[n]


def phase3_N_nn(n_alpha: int) -> float:
    """Two-zone model: capped LJ + valence shell after n_core."""
    if n_alpha < 2:
        return 0.0
    if n_alpha <= 8:
        # Light cluster: pure LJ
        return float(lj_cached(n_alpha))
    if n_alpha <= N_CORE:
        # Medium cluster: cap LJ at z_max * n / 2 (Pauli suppression)
        return float(min(lj_cached(n_alpha), Z_MAX * n_alpha / 2))
    # Beyond shell closure: core + valence
    core = phase3_N_nn(N_CORE)
    return core + (n_alpha - N_CORE) * NU_VALENCE


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

ALPHA_CLUSTER_NUCLEI = [
    ("4He",   1,  28.296),
    ("8Be",   2,  56.500),
    ("12C",   3,  92.162),
    ("16O",   4, 127.619),
    ("20Ne",  5, 160.645),
    ("24Mg",  6, 198.257),
    ("28Si",  7, 236.537),
    ("32S",   8, 271.781),
    ("36Ar",  9, 306.717),
    ("40Ca", 10, 342.052),
    ("44Ti", 11, 375.475),
    ("48Cr", 12, 411.469),
    ("52Fe", 13, 447.708),
    ("56Ni", 14, 483.991),
    ("60Zn", 15, 514.999),
    ("64Ge", 16, 545.948),
    ("68Se", 17, 576.398),
    ("72Kr", 18, 607.087),
    ("76Sr", 19, 638.040),
    ("80Zr", 20, 669.797),
    ("84Mo", 21, 697.180),
]


def main() -> None:
    print("=" * 86)
    print(" Phase 3: two-zone alpha-cluster model")
    print(" (LJ light cluster + Pauli-capped medium + valence shell)")
    print("=" * 86)
    print()
    print(f"  E_link*    = {E_LINK_STAR} MeV   (per-bond, from 12C/16O)")
    print(f"  n_core     = {N_CORE}            (Z=N=28 shell closure, "
            "standard shell model)")
    print(f"  z_max      = {Z_MAX}             (Pauli-limited max "
            "coordination)")
    print(f"  nu_valence = {NU_VALENCE}         (valence-alpha bond budget)")
    print()

    print(f"  {'nucleus':<7s} {'n':>3s} "
          f"{'B_E obs':>9s} "
          f"{'LJ N':>5s} {'P3 N':>5s} {'emp':>6s} "
          f"{'B_E pred':>10s} {'res%':>7s} {'zone'}")
    print("  " + "-" * 78)

    rows = []
    for name, n, BE_obs in ALPHA_CLUSTER_NUCLEI:
        BE_extra_obs = BE_obs - n * B_E_ALPHA
        N_emp = BE_extra_obs / E_LINK_STAR
        N_lj = lj_cached(n) if n >= 2 else 0
        N_p3 = phase3_N_nn(n)
        BE_pred = n * B_E_ALPHA + N_p3 * E_LINK_STAR
        res = (BE_pred - BE_obs) / BE_obs * 100
        if n <= 8: zone = "light"
        elif n <= N_CORE: zone = "medium (capped)"
        else: zone = f"valence (+{n-N_CORE})"
        rows.append({"name": name, "n": n, "BE_obs": BE_obs,
                     "N_lj": N_lj, "N_p3": N_p3, "N_emp": N_emp,
                     "BE_pred": BE_pred, "res": res, "zone": zone})
        print(f"  {name:<7s} {n:>3d} "
              f"{BE_obs:>9.2f} "
              f"{N_lj:>5d} {N_p3:>5.1f} {N_emp:>6.2f} "
              f"{BE_pred:>10.2f} {res:>+6.2f}% {zone}")

    res_arr = np.array([r["res"] for r in rows if r["n"] >= 3])
    rms = float(np.sqrt(np.mean(res_arr ** 2)))
    res_p2_arr = np.array([
        ((r["n"] * B_E_ALPHA + r["N_lj"] * E_LINK_STAR) - r["BE_obs"])
        / r["BE_obs"] * 100
        for r in rows if r["n"] >= 3
    ])
    rms_p2 = float(np.sqrt(np.mean(res_p2_arr ** 2)))

    print()
    print(f"  Phase 3 RMS (n_alpha >= 3, {len(res_arr)} nuclei): {rms:.2f}%")
    print(f"  Phase 2 RMS (same set, for comparison):           {rms_p2:.2f}%")
    print(f"  Improvement: {rms_p2/rms:.1f}x")

    # ---------------------------------------------------------------
    # Plot
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    n_arr = np.array([r["n"] for r in rows])
    N_emp = np.array([r["N_emp"] for r in rows])
    N_p2 = np.array([r["N_lj"] for r in rows], dtype=float)
    N_p3 = np.array([r["N_p3"] for r in rows])

    ax = axes[0]
    ax.plot(n_arr, N_emp, "ko", ms=8, zorder=4, label="empirical")
    ax.plot(n_arr, N_p2, "g^", ms=6, alpha=0.6, label="Phase 2 (pure LJ)")
    ax.plot(n_arr, N_p3, "r-", lw=2.0, label="Phase 3 (two-zone)")
    ax.axvline(N_CORE, color="purple", ls=":", lw=1.0,
                label=f"shell closure n={N_CORE}")
    ax.set_xlabel(r"$n_\alpha$")
    ax.set_ylabel(r"$N_{\rm nn}$")
    ax.set_title("Nearest-neighbour count: Phase 2 vs Phase 3")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    res_p2_full = [
        ((r["n"] * B_E_ALPHA + r["N_lj"] * E_LINK_STAR) - r["BE_obs"])
        / r["BE_obs"] * 100 for r in rows
    ]
    res_p3_full = [r["res"] for r in rows]
    ax.plot(n_arr, res_p2_full, "g^-", ms=6, alpha=0.6,
                label=f"Phase 2 (RMS {rms_p2:.2f}%)")
    ax.plot(n_arr, res_p3_full, "ro-", ms=6, lw=1.0,
                label=f"Phase 3 (RMS {rms:.2f}%)")
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(N_CORE, color="purple", ls=":", lw=1.0)
    ax.set_xlabel(r"$n_\alpha$")
    ax.set_ylabel("residual (%)")
    ax.set_title("Binding-energy prediction residual")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out = Path(__file__).resolve().parent / "nwt_alpha_cluster_phase3.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved figure: {out}")


if __name__ == "__main__":
    main()
