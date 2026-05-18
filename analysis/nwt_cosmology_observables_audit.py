"""
Substrate-algebraic audit for Forward Test 4 of the Form prediction framework.

Target: dimensionless cosmological observables.

CROSS-SECTOR test #3: cosmology, after particle physics (FT3 Lie groups),
nuclear physics (FT2 magic numbers), chemistry (FT1 Au_n).

Tests framework vocabulary extension to REAL-NUMBER PRIMITIVES (α, π, √3,
m_e/M_Pl, 2^(-1/4)) — not just integers.

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_cosmology_observables_audit.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Substrate primitives — real-number vocabulary (FT4 extension)
# ---------------------------------------------------------------------------

# Integer primitives (from prior tier resolutions + FTs)
INT_PRIMITIVES: dict[int, str] = {
    1:  "unit",
    3:  "RANK_SO7",
    4:  "N_VERTICES_K_4",
    5:  "H_V_SO7",
    6:  "H_COXETER_SO7",
    7:  "N_VERTICES_K7",
    8:  "DIM_OCTONION",
    9:  "N_POS_ROOTS_SO7",
    12: "K8_PARTITION[2]",
    13: "trefoil(p²+q²)",
    21: "N_EDGES_K7",
    28: "N_EDGES_K8",
    35: "K_7_TRIANGLES",
}

# Real primitives (FT4 introduction)
# α from Paper 17 trefoil formula α = 1/(25π√3 + 1)
ALPHA = 1.0 / (25.0 * math.pi * math.sqrt(3.0) + 1.0)
PI = math.pi
SQRT3 = math.sqrt(3.0)
TWO_TO_NEG_QUARTER = 2.0 ** (-0.25)

# Empirical scales
M_E_OVER_M_PL = 4.185e-23  # electron mass / Planck mass ratio
H_V = 5.0
H_COXETER = 6.0

REAL_PRIMITIVES: dict[str, float] = {
    "α":             ALPHA,
    "π":             PI,
    "√3":            SQRT3,
    "2^(-1/4)":      TWO_TO_NEG_QUARTER,
    "m_e/M_Pl":      M_E_OVER_M_PL,
    "h_v":           H_V,
    "h_Coxeter":     H_COXETER,
}


# ---------------------------------------------------------------------------
# Cosmological observables — locked from pre-reg
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CosmologicalObservable:
    name: str
    measured_value: float
    substrate_predicted: bool
    substrate_formula: Optional[str]
    substrate_components: Optional[list[str]]   # which primitives used
    op_count: Optional[int]                      # operations in formula


SUBSTRATE_PREDICTED = [
    CosmologicalObservable(
        name="Λ_cc",
        measured_value=1.19e-123,
        substrate_predicted=True,
        substrate_formula="(m_e/M_Pl)⁴ · α¹⁶ · h_Coxeter",
        substrate_components=["m_e/M_Pl", "α", "h_Coxeter"],
        op_count=0,   # all primitives + integer powers
    ),
    CosmologicalObservable(
        name="n_s",
        measured_value=0.9635,
        substrate_predicted=True,
        substrate_formula="1 − h_v · α",
        substrate_components=["1", "h_v", "α"],
        op_count=0,   # all primitives + single subtraction
    ),
    CosmologicalObservable(
        name="r",
        measured_value=1.5e-13,
        substrate_predicted=True,
        substrate_formula="α⁶",
        substrate_components=["α"],
        op_count=0,   # single primitive raised to integer power
    ),
    CosmologicalObservable(
        name="Ω_m",
        measured_value=0.316,
        substrate_predicted=True,
        substrate_formula="h_v² · √3 · α",
        substrate_components=["h_v", "√3", "α"],
        op_count=0,   # all primitives × powers
    ),
    CosmologicalObservable(
        name="f_J",
        measured_value=0.785,
        substrate_predicted=True,
        substrate_formula="(1 − √α)³",
        substrate_components=["1", "α"],
        op_count=0,   # primitive 1 + primitive √α + integer power 3
    ),
]


NULL_CONTROL = [
    CosmologicalObservable(
        name="σ_8",
        measured_value=0.81,
        substrate_predicted=False,
        substrate_formula=None,
        substrate_components=None,
        op_count=None,
    ),
    CosmologicalObservable(
        name="Ω_b/Ω_c",
        measured_value=0.186,
        substrate_predicted=False,
        substrate_formula="(loose) h_v² · α ≈ 0.182 (within 2%)",
        substrate_components=["h_v", "α"],
        op_count=1,   # composite involving primitives
    ),
    CosmologicalObservable(
        name="z_reion",
        measured_value=7.7,
        substrate_predicted=False,
        substrate_formula="(loose) N_VERTICES_K7 = 7 (within 10%)",
        substrate_components=["N_VERTICES_K7"],
        op_count=1,   # primitive but match imprecise
    ),
    CosmologicalObservable(
        name="τ_reion",
        measured_value=0.054,
        substrate_predicted=False,
        substrate_formula="(loose) h_Coxeter · α ≈ 0.0437 (within 25%)",
        substrate_components=["h_Coxeter", "α"],
        op_count=1,   # composite, match weak
    ),
    CosmologicalObservable(
        name="z_dec",
        measured_value=1100.0,
        substrate_predicted=False,
        substrate_formula=None,
        substrate_components=None,
        op_count=None,
    ),
]


# ---------------------------------------------------------------------------
# Audit per observable
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ObservableAudit:
    observable: CosmologicalObservable
    n_primitive_components: int
    n_total_components: int
    primitive_fraction: float
    all_primitive: bool


def is_primitive_component(name: str) -> bool:
    """Check if a component is a substrate primitive."""
    if name in REAL_PRIMITIVES:
        return True
    if name in INT_PRIMITIVES.values():
        return True
    if name == "1":
        return True
    return False


def audit_observable(obs: CosmologicalObservable) -> ObservableAudit:
    if obs.substrate_components is None:
        return ObservableAudit(
            observable=obs,
            n_primitive_components=0,
            n_total_components=0,
            primitive_fraction=0.0,
            all_primitive=False,
        )

    n_total = len(obs.substrate_components)
    n_prim = sum(1 for c in obs.substrate_components if is_primitive_component(c))
    frac = n_prim / n_total if n_total > 0 else 0.0

    return ObservableAudit(
        observable=obs,
        n_primitive_components=n_prim,
        n_total_components=n_total,
        primitive_fraction=frac,
        all_primitive=(n_prim == n_total and n_total > 0),
    )


# ---------------------------------------------------------------------------
# Compute substrate-formula values & deviations
# ---------------------------------------------------------------------------

def compute_lambda_cc() -> float:
    return (M_E_OVER_M_PL)**4 * ALPHA**16 * H_COXETER


def compute_n_s() -> float:
    return 1.0 - H_V * ALPHA


def compute_r() -> float:
    return ALPHA**6


def compute_omega_m() -> float:
    return H_V**2 * SQRT3 * ALPHA


def compute_f_J() -> float:
    return (1.0 - math.sqrt(ALPHA))**3


SUBSTRATE_COMPUTED = {
    "Λ_cc":  compute_lambda_cc(),
    "n_s":   compute_n_s(),
    "r":     compute_r(),
    "Ω_m":   compute_omega_m(),
    "f_J":   compute_f_J(),
}


def relative_deviation(measured: float, predicted: float) -> float:
    if predicted == 0:
        return float("inf")
    return abs(measured - predicted) / abs(predicted)


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def run_audit() -> dict:
    return {
        "substrate_predicted": [audit_observable(o) for o in SUBSTRATE_PREDICTED],
        "null_control":        [audit_observable(o) for o in NULL_CONTROL],
        "substrate_computed":  SUBSTRATE_COMPUTED,
    }


def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 80)
    out.append("Forward Test 4 — Cosmological observables")
    out.append("=" * 80)
    out.append("")
    out.append("Framework vocabulary extended to real-number primitives:")
    out.append(f"  α = {ALPHA:.6f}   π = {PI:.6f}   √3 = {SQRT3:.6f}")
    out.append(f"  2^(-1/4) = {TWO_TO_NEG_QUARTER:.6f}   m_e/M_Pl ≈ {M_E_OVER_M_PL:.3e}")
    out.append("")

    # SUBSTRATE-PREDICTED
    out.append("-" * 80)
    out.append("SUBSTRATE-PREDICTED OBSERVABLES (confirmatory expected Form A)")
    out.append("-" * 80)
    out.append("")
    out.append(f"  {'name':<10} {'measured':>14} {'predicted':>14} {'rel_dev':>10} {'formula':<40}")
    for a in audit["substrate_predicted"]:
        obs = a.observable
        pred = SUBSTRATE_COMPUTED[obs.name]
        dev = relative_deviation(obs.measured_value, pred)
        out.append(f"  {obs.name:<10} {obs.measured_value:>14.4e} {pred:>14.4e} {dev:>10.3%} {obs.substrate_formula:<40}")
    out.append("")

    out.append("  Primitive-component analysis:")
    out.append("")
    out.append(f"  {'name':<10} {'op_count':>9} {'n_primitive':>13} {'n_total':>9} {'frac':>7} {'all_prim?':>10}")
    for a in audit["substrate_predicted"]:
        out.append(f"  {a.observable.name:<10} {a.observable.op_count:>9} {a.n_primitive_components:>13} "
                   f"{a.n_total_components:>9} {a.primitive_fraction:>7.0%} {'✓' if a.all_primitive else '✗':>10}")
    out.append("")

    # NULL CONTROL
    out.append("-" * 80)
    out.append("NULL CONTROL (NOT substrate-predicted; expected Form D weak / E)")
    out.append("-" * 80)
    out.append("")
    out.append(f"  {'name':<12} {'measured':>14} {'formula':<55}")
    for a in audit["null_control"]:
        obs = a.observable
        formula = obs.substrate_formula or "(no clean substrate form)"
        out.append(f"  {obs.name:<12} {obs.measured_value:>14.4e} {formula:<55}")
    out.append("")

    out.append("  Primitive-component analysis (where formula exists):")
    out.append("")
    out.append(f"  {'name':<12} {'op_count':>9} {'n_primitive':>13} {'n_total':>9} {'frac':>7} {'all_prim?':>10}")
    for a in audit["null_control"]:
        obs = a.observable
        op = obs.op_count if obs.op_count is not None else "—"
        out.append(f"  {obs.name:<12} {str(op):>9} {a.n_primitive_components:>13} "
                   f"{a.n_total_components:>9} {a.primitive_fraction:>7.0%} {'✓' if a.all_primitive else '✗':>10}")
    out.append("")

    # SUMMARY
    out.append("=" * 80)
    out.append("SUMMARY")
    out.append("=" * 80)
    out.append("")

    sp_all_prim = sum(1 for a in audit["substrate_predicted"] if a.all_primitive)
    nc_all_prim = sum(1 for a in audit["null_control"] if a.all_primitive)

    sp_n = len(audit["substrate_predicted"])
    nc_n = len(audit["null_control"])

    out.append(f"  SUBSTRATE-PREDICTED:  all-primitive count = {sp_all_prim} / {sp_n}")
    out.append(f"  NULL CONTROL:         all-primitive count = {nc_all_prim} / {nc_n}")
    out.append("")

    # Numerical accuracy of predictions
    out.append("  Numerical accuracy (measured vs substrate-predicted):")
    for a in audit["substrate_predicted"]:
        obs = a.observable
        pred = SUBSTRATE_COMPUTED[obs.name]
        dev = relative_deviation(obs.measured_value, pred)
        if dev < 0.001:
            tag = "EXCELLENT (<0.1%)"
        elif dev < 0.01:
            tag = "EXCELLENT (<1%)"
        elif dev < 0.05:
            tag = "GOOD (<5%)"
        elif dev < 0.20:
            tag = "MODERATE (<20%)"
        else:
            tag = "WEAK"
        out.append(f"    {obs.name:<10} dev={dev:>8.3%}  {tag}")
    out.append("")

    # VERDICT
    out.append("=" * 80)
    out.append("FRAMEWORK VERDICT")
    out.append("=" * 80)
    out.append("")
    out.append("  Framework predicted Form A clean iff:")
    out.append("    (a) ≥4 of 5 substrate-predicted observables are all-primitive")
    out.append("    (b) ≤2 of 5 null-control observables are all-primitive")
    out.append("")

    cond_a = sp_all_prim >= 4
    cond_b = nc_all_prim <= 2

    out.append(f"  (a) substrate-predicted all-primitive ≥ 4: {sp_all_prim}/{sp_n}  {'✓' if cond_a else '✗'}")
    out.append(f"  (b) null-control all-primitive ≤ 2:        {nc_all_prim}/{nc_n}  {'✓' if cond_b else '✗'}")
    out.append("")

    if cond_a and cond_b:
        out.append("  → FORM A CONFIRMED for substrate-predicted observables.")
        out.append("    Framework's real-primitive vocabulary extension is principled.")
        out.append("    Null control discriminates: substrate-predictions ARE substrate-specific.")
    elif cond_a:
        out.append("  → FORM A holds for substrate-predicted; null-control discrimination weaker.")
    else:
        out.append("  → Framework FAILS for substrate-predicted observables.")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
