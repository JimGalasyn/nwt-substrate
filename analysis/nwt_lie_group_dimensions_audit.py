"""
Substrate-algebraic audit for Forward Test 3 of the Form prediction framework.

Target: dimensions of fundamental Lie groups in physics, classified into:
  - Empirically REALIZED (SM gauge groups, Lorentz, NWT substrate algebras)
  - Empirically UNREALIZED (GUT candidates ruled out or unobserved)
  - NULL control (mathematical Lie groups not in fundamental physics)

Framework prediction: P-axis (primitive density) discriminates cleanly.
Realized dimensions are substrate primitives (op=0); unrealized and null-
control are composites (op≥1).

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_lie_group_dimensions_audit.py
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from nwt_substrate.isa.constants import (
    DIM_OCTONION,
    H_COXETER_SO7,
    H_V_SO7,
    K8_PARTITION,
    N_EDGES_K7,
    N_EDGES_K8,
    N_POS_ROOTS_SO7,
    N_VERTICES_K7,
    RANK_SO7,
)


# ---------------------------------------------------------------------------
# Substrate-canonical integer set (locked from prior tier resolutions + FTs)
# ---------------------------------------------------------------------------

CANONICAL_PRIMITIVES: dict[int, str] = {
    1:                          "unit",              # 1 (identity)
    RANK_SO7:                   "RANK_SO7",          # 3
    4:                          "N_VERTICES_K_4",    # 4
    H_V_SO7:                    "H_V_SO7",           # 5
    H_COXETER_SO7:              "H_COXETER_SO7",     # 6
    N_VERTICES_K7:              "N_VERTICES_K7",     # 7
    DIM_OCTONION:               "DIM_OCTONION",      # 8
    N_POS_ROOTS_SO7:            "N_POS_ROOTS_SO7",   # 9
    K8_PARTITION[2]:            "K8_PARTITION[2]",   # 12
    13:                         "trefoil(p²+q²)",    # 13
    N_EDGES_K7:                 "N_EDGES_K7",        # 21
    N_EDGES_K8:                 "N_EDGES_K8",        # 28
    35:                         "K_7_TRIANGLES",     # 35
}

# Derived primitives (single-operation forms used in tier resolutions, but
# treated as primitives for cross-sector transfer per [[nuclear-magic-numbers-resolution]])
CANONICAL_DERIVED: dict[int, str] = {
    N_EDGES_K7 - N_VERTICES_K7: "N_EDGES_K7 − N_VERTICES_K7 (=14, G_2 dim)",
    N_EDGES_K7 - H_V_SO7:       "N_EDGES_K7 − H_V_SO7 (=16)",
    N_EDGES_K7 - RANK_SO7:      "N_EDGES_K7 − RANK_SO7 (=18)",
    35 - RANK_SO7:              "K_7_TRIANGLES − RANK_SO7 (=32)",
    H_V_SO7 - RANK_SO7:         "H_V_SO7 − RANK_SO7 (=2)",
}


def is_primitive(n: int) -> bool:
    return n in CANONICAL_PRIMITIVES


def is_derived_primitive(n: int) -> bool:
    """Derived primitives are single-subtraction forms documented in tier
    resolutions. They count as primitives for P-axis scoring."""
    return n in CANONICAL_DERIVED


def primitive_label(n: int) -> Optional[str]:
    if n in CANONICAL_PRIMITIVES:
        return CANONICAL_PRIMITIVES[n]
    if n in CANONICAL_DERIVED:
        return CANONICAL_DERIVED[n]
    return None


# ---------------------------------------------------------------------------
# Operation enumeration (from FT2 audit, slightly extended)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SubstrateForm:
    operation_count: int
    operation_type: str
    formula: str


def find_substrate_forms(n: int) -> list[SubstrateForm]:
    """Enumerate all substrate forms for N, sorted by operation count."""
    forms = []

    # op=0: primitive (direct canonical OR derived)
    if is_primitive(n):
        forms.append(SubstrateForm(0, "primitive", CANONICAL_PRIMITIVES[n]))
    if is_derived_primitive(n):
        forms.append(SubstrateForm(0, "derived-primitive", CANONICAL_DERIVED[n]))

    canonical = sorted(CANONICAL_PRIMITIVES.keys())

    # op=1: subtraction
    seen = set()
    for a in canonical:
        b = a - n
        if b > 0 and b in CANONICAL_PRIMITIVES and (a, b) not in seen:
            seen.add((a, b))
            forms.append(SubstrateForm(
                1, "subtraction",
                f"{CANONICAL_PRIMITIVES[a]} − {CANONICAL_PRIMITIVES[b]}",
            ))

    # op=1: product (a × b)
    seen_prod = set()
    for a in canonical:
        if n % a == 0 and a > 0:
            b = n // a
            if b in CANONICAL_PRIMITIVES:
                pair = tuple(sorted([a, b]))
                if pair not in seen_prod:
                    seen_prod.add(pair)
                    op = "squared" if a == b else "product"
                    forms.append(SubstrateForm(
                        1, op,
                        f"{CANONICAL_PRIMITIVES[a]} × {CANONICAL_PRIMITIVES[b]}",
                    ))

    # op=1: 2-term additive
    seen_add = set()
    for a in canonical:
        b = n - a
        if b > 0 and b in CANONICAL_PRIMITIVES:
            pair = tuple(sorted([a, b]))
            if pair not in seen_add:
                seen_add.add(pair)
                forms.append(SubstrateForm(
                    1, "additive-2",
                    f"{CANONICAL_PRIMITIVES[a]} + {CANONICAL_PRIMITIVES[b]}",
                ))

    # op=2: scaled additive a + k*b
    for k in (2, 3):
        for a in canonical:
            for b in canonical:
                if a + k * b == n:
                    forms.append(SubstrateForm(
                        2, f"additive-{1+k}",
                        f"{CANONICAL_PRIMITIVES[a]} + {k}·{CANONICAL_PRIMITIVES[b]}",
                    ))

    forms.sort(key=lambda f: f.operation_count)
    return forms


def min_op_count(n: int) -> Optional[int]:
    forms = find_substrate_forms(n)
    return forms[0].operation_count if forms else None


# ---------------------------------------------------------------------------
# Target sets — Lie group dimensions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LieGroupTarget:
    name: str
    dimension: int
    realized: bool        # empirically observed in fundamental physics
    role: str             # description


REALIZED = [
    LieGroupTarget("U(1)",          1,  True,  "electromagnetism, hypercharge"),
    LieGroupTarget("SU(2)",         3,  True,  "weak isospin"),
    LieGroupTarget("Lorentz SO(3,1)", 6,True,  "spacetime symmetry"),
    LieGroupTarget("SU(3)",         8,  True,  "QCD color"),
    LieGroupTarget("G_2",           14, True,  "octonion auto / Spin(7) auto"),
    LieGroupTarget("Spin(7)",       21, True,  "NWT substrate algebra"),
    LieGroupTarget("Spin(8)",       28, True,  "NWT K_8 / triality"),
]

UNREALIZED_GUT = [
    LieGroupTarget("SU(4) Pati-Salam", 15, False, "partial GUT"),
    LieGroupTarget("SU(5) Georgi-Glashow", 24, False, "EXCLUDED by proton decay"),
    LieGroupTarget("SO(10)",        45, False, "unrealized GUT"),
    LieGroupTarget("E_6",           78, False, "heterotic-string motivated"),
    LieGroupTarget("E_7",           133,False, "no empirical signal"),
    LieGroupTarget("E_8",           248,False, "M-theory motivated"),
]

NULL_CONTROL = [
    LieGroupTarget("F_4",           52, False, "Jordan algebra auto, BSM only"),
    LieGroupTarget("SU(7)",         48, False, "mathematical only"),
    LieGroupTarget("Spin(11)",      55, False, "mathematical only"),
    LieGroupTarget("Spin(12)",      66, False, "mathematical only"),
    LieGroupTarget("Spin(16)",      120,False, "mathematical only"),
]


# ---------------------------------------------------------------------------
# Audit per target
# ---------------------------------------------------------------------------

def audit_target(target: LieGroupTarget) -> dict:
    forms = find_substrate_forms(target.dimension)
    return {
        "target": target,
        "min_op_count": forms[0].operation_count if forms else None,
        "n_forms": len(forms),
        "is_primitive": is_primitive(target.dimension) or is_derived_primitive(target.dimension),
        "forms": forms,
    }


def run_audit() -> dict:
    return {
        "realized": [audit_target(t) for t in REALIZED],
        "unrealized_gut": [audit_target(t) for t in UNREALIZED_GUT],
        "null_control": [audit_target(t) for t in NULL_CONTROL],
    }


def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 80)
    out.append("Forward Test 3 — Fundamental Lie group dimensions in physics")
    out.append("=" * 80)
    out.append("")
    out.append("Substrate primitives: {3, 4, 5, 6, 7, 8, 9, 12, 13, 21, 28, 35}")
    out.append("Derived primitives:   {2=5-3, 14=21-7, 16=21-5, 18=21-3, 32=35-3}")
    out.append("")

    def render_section(title: str, audits: list[dict]):
        out.append("-" * 80)
        out.append(title)
        out.append("-" * 80)
        out.append("")
        out.append(f"  {'name':<22} {'dim':>5} {'min_op':>8} {'primitive?':>11} {'n_forms':>9}")
        out.append("  " + "-" * 60)
        for a in audits:
            t = a["target"]
            op = str(a["min_op_count"]) if a["min_op_count"] is not None else "—"
            prim = "✓" if a["is_primitive"] else "✗"
            out.append(f"  {t.name:<22} {t.dimension:>5} {op:>8} {prim:>11} {a['n_forms']:>9}")

        out.append("")
        out.append("  Cleanest substrate form per target:")
        for a in audits:
            t = a["target"]
            if a["forms"]:
                f = a["forms"][0]
                out.append(f"    {t.name:<22} dim={t.dimension:>3}  op={f.operation_count}  [{f.operation_type:<18}]  {f.formula}")
            else:
                out.append(f"    {t.name:<22} dim={t.dimension:>3}  NO SUBSTRATE FORM")
        out.append("")

    render_section("REALIZED IN PHYSICS (substrate-primitive prediction)", audit["realized"])
    render_section("UNREALIZED GUT (substrate-composite prediction)", audit["unrealized_gut"])
    render_section("NULL CONTROL (mathematical only)", audit["null_control"])

    # Summary
    out.append("=" * 80)
    out.append("SUMMARY — Primitive count per class")
    out.append("=" * 80)
    out.append("")

    realized_prim = sum(1 for a in audit["realized"] if a["is_primitive"])
    realized_n = len(audit["realized"])
    unrealized_prim = sum(1 for a in audit["unrealized_gut"] if a["is_primitive"])
    unrealized_n = len(audit["unrealized_gut"])
    null_prim = sum(1 for a in audit["null_control"] if a["is_primitive"])
    null_n = len(audit["null_control"])

    out.append(f"  REALIZED (predicted primitive):       {realized_prim} / {realized_n} primitives ({100*realized_prim/realized_n:.0f}%)")
    out.append(f"  UNREALIZED GUT (predicted composite): {unrealized_prim} / {unrealized_n} primitives ({100*unrealized_prim/unrealized_n:.0f}%)")
    out.append(f"  NULL CONTROL (predicted composite):   {null_prim} / {null_n} primitives ({100*null_prim/null_n:.0f}%)")
    out.append("")

    realized_avg = sum(a["min_op_count"] for a in audit["realized"]
                       if a["min_op_count"] is not None) / realized_n
    unrealized_avg = sum(a["min_op_count"] for a in audit["unrealized_gut"]
                         if a["min_op_count"] is not None) / unrealized_n
    null_avg = sum(a["min_op_count"] for a in audit["null_control"]
                   if a["min_op_count"] is not None) / null_n

    out.append(f"  Average min-op count:")
    out.append(f"    Realized:        {realized_avg:.2f}")
    out.append(f"    Unrealized GUT:  {unrealized_avg:.2f}")
    out.append(f"    Null control:    {null_avg:.2f}")
    out.append("")

    # VERDICT
    out.append("=" * 80)
    out.append("FRAMEWORK VERDICT")
    out.append("=" * 80)
    out.append("")
    out.append("  Framework predicted Form A clean iff:")
    out.append("    (a) ≥ 6 of 7 realized dimensions are primitives")
    out.append("    (b) ≤ 1 of 6 unrealized dimensions is a primitive")
    out.append("    (c) 0 of 5 null-control dimensions is a primitive")
    out.append("")

    cond_a = realized_prim >= 6
    cond_b = unrealized_prim <= 1
    cond_c = null_prim == 0

    out.append(f"  (a) realized primitives ≥ 6:    {realized_prim}/{realized_n}  {'✓' if cond_a else '✗'}")
    out.append(f"  (b) unrealized primitives ≤ 1: {unrealized_prim}/{unrealized_n}  {'✓' if cond_b else '✗'}")
    out.append(f"  (c) null primitives = 0:        {null_prim}/{null_n}  {'✓' if cond_c else '✗'}")
    out.append("")

    if cond_a and cond_b and cond_c:
        out.append("  → FORM A CLEAN: framework PASSES Forward Test 3.")
        out.append("    P-axis discriminates realized vs unrealized cleanly.")
        out.append("    Substrate primitives = realized gauge group dimensions in physics.")
    elif cond_a and (cond_b or cond_c):
        out.append("  → FORM D LOAD-BEARING: framework partially PASSES.")
        out.append("    Realized dimensions are primitives; some composite/null leakage.")
    elif cond_a:
        out.append("  → FORM D WEAK: realized prediction holds; composite/null discrimination fails.")
    else:
        out.append("  → FRAMEWORK FAILS: realized dimensions are not all primitives.")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
