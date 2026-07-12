"""Neighbouring-value test (gauntlet mode 5) — the look-elsewhere volume of the
substrate fitting procedure, measured.

The derivation-audit gauntlet's mode 5 asks of every clean numeric coincidence:
*does it hold at neighbouring values, or only at the clean one?*  This module
answers it operationally for the constants-stack fitting pattern — a prefactor
from the K_7/Spin(7) structural-integer menu, times a half-integer power of α,
times up to two rational-coefficient bracket stages:

    T ≈ (p/q) · α^(n/2) · (1 + c₁·α + c₂·α²),   p, q, c-numerators/denominators ∈ menu

by running the SAME greedy procedure against RANDOM targets (log-uniform around
the real one) and reporting the residual distribution.  If the procedure fits
noise as well as it fits the physical constant, the physical hit carries no
information — the agreement is a property of the menu, not of the universe.

Measured result (2026-07-12, the audit that retired the G-relation rescue):

    full structural menu {1,2,3,5,6,7,8,14,15,21}:
        median residual ≈ 2.5 ppm; ~83% of random targets land inside
        CODATA-2018 G's 22 ppm error bar — BETTER than the published
        m_e/M_Pl chain's actual ~5 ppm hit.
    minimal menu {1,2,3,7,8,21}:
        median ≈ 240 ppm; still ~24% inside the G bar.

Under either menu, a ppm-scale landing on m_e/M_Pl is the EXPECTED outcome of
the procedure for an arbitrary target, before multiplying in the documented
form-search freedom (Paper 17's own "identified four independent ways").  The
verdict this feeds: the (8/7)·α^(21/2)·bracket relation is look-elsewhere, not
structure — quantitatively, by the program's own mode-5 test.

The instrument is deliberately generic: pass any target, menu, or stage count.
Every constants-stack claim (present or future) can and should be run through
it BEFORE its agreement is reported anywhere.

    python -m nwt_substrate.benchmarks.neighbouring_value
"""

from __future__ import annotations

import math
import random

ALPHA_MEASURED = 1.0 / 137.035999084     # CODATA-2018 (witness-layer value)

FULL_MENU = (1, 2, 3, 5, 6, 7, 8, 14, 15, 21)
MINIMAL_MENU = (1, 2, 3, 7, 8, 21)

M_E_OVER_M_PL = 4.18546e-23              # the real target (CODATA-2018 primaries)
G_BAR_REL = 22e-6                         # CODATA-2018 G relative 1σ (the fattest bar)


def fit_residual(target: float, menu: tuple[int, ...] = FULL_MENU,
                 alpha: float = ALPHA_MEASURED, stages: int = 2,
                 exp_lo: float = 7.0, exp_hi: float = 15.0) -> float:
    """Greedy best fit of ``target`` by the substrate pattern; returns the final
    relative residual.  Mirrors the published chain's construction order:
    prefactor × exponent first, then one rational bracket coefficient per
    correction stage (c₁·α, c₂·α², ...)."""
    prefactors = sorted({p / q for p in menu for q in menu if 0.05 <= p / q <= 20})
    cmenu = sorted({s * p / q for p in menu for q in menu for s in (1, -1)} | {0.0})
    exps = [n / 2 for n in range(int(2 * exp_lo), int(2 * exp_hi) + 1)]

    _, pref, n = min((abs(math.log(target / (pf * alpha ** e))), pf, e)
                     for pf in prefactors for e in exps)
    base = pref * alpha ** n
    for k in range(1, stages + 1):
        r = target / base - 1.0
        c = min(cmenu, key=lambda cc: abs(r - cc * alpha ** k))
        base *= (1.0 + c * alpha ** k)
    return abs(target / base - 1.0)


def sweep(menu: tuple[int, ...] = FULL_MENU, n_targets: int = 2000,
          decades: float = 3.0, seed: int = 42,
          center: float = M_E_OVER_M_PL) -> dict:
    """Residual distribution of the procedure over random log-uniform targets.
    Returns median residual and the fraction landing inside the G error bar —
    the look-elsewhere volume the real hit must be discounted by."""
    rng = random.Random(seed)
    res = sorted(
        fit_residual(center * 10 ** rng.uniform(-decades / 2, decades / 2), menu)
        for _ in range(n_targets))
    return {
        "menu": menu,
        "median": res[len(res) // 2],
        "frac_within_G_bar": sum(r <= G_BAR_REL for r in res) / len(res),
        "frac_within_100ppm": sum(r <= 100e-6 for r in res) / len(res),
        "n": n_targets,
    }


def main(argv: list[str] | None = None) -> int:
    print("Neighbouring-value test (gauntlet mode 5): fit RANDOM targets with")
    print("(menu prefactor) x alpha^(n/2) x (1 + c1*a + c2*a^2), menu rationals\n")
    for menu, label in ((FULL_MENU, "full structural menu"),
                        (MINIMAL_MENU, "minimal menu")):
        s = sweep(menu)
        print(f"  {label:22s} {menu}")
        print(f"    median residual        : {s['median']*1e6:8.1f} ppm")
        print(f"    within G bar (22 ppm)  : {s['frac_within_G_bar']*100:5.1f}%")
        print(f"    within 100 ppm         : {s['frac_within_100ppm']*100:5.1f}%\n")
    actual = fit_residual(M_E_OVER_M_PL)
    print(f"  the real m_e/M_Pl, same procedure: {actual*1e6:.1f} ppm "
          "(published chain: ~5 ppm)")
    print("  -> a ppm-scale landing is the EXPECTED outcome for an arbitrary")
    print("     target; the physical hit carries no information beyond the menu.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
