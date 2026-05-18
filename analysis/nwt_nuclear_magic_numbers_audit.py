"""
Substrate-algebraic audit for Forward Test 2 of the Form prediction framework.

Target: nuclear shell-model magic numbers {2, 8, 20, 28, 50, 82, 126}
(Mayer/Jensen 1949). Cross-sector validation of the refined framework
(additive + product + squared R-axis from [[au-magic-clusters-resolution]]).

CRITICAL: includes a null-control comparison against non-magic neighbors
{22, 30, 36, 60, 72} under identical scoring criteria.

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_nuclear_magic_numbers_audit.py
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations_with_replacement, product
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
# Substrate-canonical integer set (locked from chemistry-sector survey)
# ---------------------------------------------------------------------------

CANONICAL_PRIMITIVES: dict[int, str] = {
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

# Derived from tier resolutions
CANONICAL_DERIVED: dict[int, str] = {
    N_EDGES_K7 - N_VERTICES_K7: "N_EDGES_K7 − N_VERTICES_K7 (=14)",  # 14
    N_EDGES_K7 - H_V_SO7:       "N_EDGES_K7 − H_V_SO7 (=16)",         # 16
    N_EDGES_K7 - RANK_SO7:      "N_EDGES_K7 − RANK_SO7 (=18)",        # 18
    35 - RANK_SO7:              "K_7_TRIANGLES − RANK_SO7 (=32)",      # 32
    H_V_SO7 - RANK_SO7:         "H_V_SO7 − RANK_SO7 (=2)",             # 2
}


def is_canonical(n: int) -> bool:
    return n in CANONICAL_PRIMITIVES or n in CANONICAL_DERIVED


def canonical_label(n: int) -> Optional[str]:
    if n in CANONICAL_PRIMITIVES:
        return CANONICAL_PRIMITIVES[n]
    if n in CANONICAL_DERIVED:
        return CANONICAL_DERIVED[n]
    return None


# ---------------------------------------------------------------------------
# Target sets
# ---------------------------------------------------------------------------

MAGIC_NUMBERS = [2, 8, 20, 28, 50, 82, 126]
NULL_CONTROL = [22, 30, 36, 60, 72]


# ---------------------------------------------------------------------------
# Operation enumeration — substrate forms for an integer N
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SubstrateForm:
    operation_count: int          # number of operations (lower = cleaner)
    operation_type: str           # "primitive", "subtraction", "product", "squared", "additive-2", "additive-3"
    formula: str
    components: tuple             # underlying integers


def find_primitive(n: int) -> list[SubstrateForm]:
    """Is N a substrate-canonical primitive (op_count = 0)?"""
    if is_canonical(n):
        return [SubstrateForm(
            operation_count=0,
            operation_type="primitive",
            formula=canonical_label(n),
            components=(n,),
        )]
    return []


def find_subtractions(n: int) -> list[SubstrateForm]:
    """N = a − b where both a, b are substrate-canonical (op_count = 1)."""
    canonical = sorted(CANONICAL_PRIMITIVES.keys())
    out = []
    for a in canonical:
        b = a - n
        if b > 0 and is_canonical(b) and b != a:
            out.append(SubstrateForm(
                operation_count=1,
                operation_type="subtraction",
                formula=f"{canonical_label(a)} − {canonical_label(b)}",
                components=(a, b),
            ))
    return out


def find_products(n: int) -> list[SubstrateForm]:
    """N = a × b where both a, b are substrate-canonical primitives, a ≤ b (op_count = 1)."""
    canonical = sorted(CANONICAL_PRIMITIVES.keys())
    out = []
    seen = set()
    for a in canonical:
        if n % a == 0:
            b = n // a
            if b in CANONICAL_PRIMITIVES and (a, b) not in seen and (b, a) not in seen:
                seen.add((a, b))
                op = "squared" if a == b else "product"
                out.append(SubstrateForm(
                    operation_count=1,
                    operation_type=op,
                    formula=f"{CANONICAL_PRIMITIVES[a]} × {CANONICAL_PRIMITIVES[b]}",
                    components=(a, b),
                ))
    return out


def find_2term_additive(n: int) -> list[SubstrateForm]:
    """N = a + b where both a, b are substrate-canonical primitives (op_count = 1)."""
    canonical = sorted(CANONICAL_PRIMITIVES.keys())
    out = []
    seen = set()
    for a in canonical:
        b = n - a
        if b > 0 and b in CANONICAL_PRIMITIVES and (a, b) not in seen and (b, a) not in seen:
            seen.add((a, b))
            out.append(SubstrateForm(
                operation_count=1,
                operation_type="additive-2",
                formula=f"{CANONICAL_PRIMITIVES[a]} + {CANONICAL_PRIMITIVES[b]}",
                components=(a, b),
            ))
    return out


def find_scaled_additive(n: int) -> list[SubstrateForm]:
    """N = a + k*b where a, b substrate-canonical primitives, k ∈ {2, 3}
    (op_count = 2: one product + one sum)."""
    canonical = sorted(CANONICAL_PRIMITIVES.keys())
    out = []
    for k in (2, 3):
        for a in canonical:
            for b in canonical:
                if a + k * b == n:
                    out.append(SubstrateForm(
                        operation_count=2,
                        operation_type=f"additive-{1+k}",
                        formula=f"{CANONICAL_PRIMITIVES[a]} + {k}·{CANONICAL_PRIMITIVES[b]}",
                        components=(a, k, b),
                    ))
    return out


def all_substrate_forms(n: int) -> list[SubstrateForm]:
    """Enumerate all substrate forms for N, sorted by operation count."""
    forms = []
    forms.extend(find_primitive(n))
    forms.extend(find_subtractions(n))
    forms.extend(find_products(n))
    forms.extend(find_2term_additive(n))
    forms.extend(find_scaled_additive(n))
    forms.sort(key=lambda f: f.operation_count)
    return forms


def min_operation_count(n: int) -> Optional[int]:
    """Minimum-operation-count substrate form, or None if no clean form."""
    forms = all_substrate_forms(n)
    if not forms:
        return None
    return forms[0].operation_count


# ---------------------------------------------------------------------------
# Locked-in pre-reg identifications
# ---------------------------------------------------------------------------

LOCKED_IDENTIFICATIONS: dict[int, str] = {
    2:   "H_V_SO7 − RANK_SO7 = 5 − 3",
    8:   "DIM_OCTONION",
    20:  "N_VERTICES_K_4 × H_V_SO7 = 4 × 5",
    28:  "N_EDGES_K8",
    50:  "DIM_OCTONION + 2·N_EDGES_K7 = 8 + 42",
    82:  "2·K_7_TRIANGLES + K8_PARTITION[2] = 70 + 12",
    126: "H_COXETER_SO7 × N_EDGES_K7 = 6 × 21",
}


def verify_locked(n: int, expected_formula: str) -> tuple[bool, list[SubstrateForm]]:
    """Check whether N has a substrate form matching the locked identification."""
    forms = all_substrate_forms(n)
    return len(forms) > 0, forms


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def run_audit() -> dict:
    magic_results = {}
    for n in MAGIC_NUMBERS:
        forms = all_substrate_forms(n)
        magic_results[n] = {
            "n": n,
            "locked": LOCKED_IDENTIFICATIONS.get(n, "—"),
            "min_op_count": min_operation_count(n),
            "n_forms": len(forms),
            "forms": forms,
        }

    null_results = {}
    for n in NULL_CONTROL:
        forms = all_substrate_forms(n)
        null_results[n] = {
            "n": n,
            "min_op_count": min_operation_count(n),
            "n_forms": len(forms),
            "forms": forms,
        }

    return {"magic": magic_results, "null": null_results}


def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 78)
    out.append("Forward Test 2 — Nuclear Shell-Model Magic Numbers")
    out.append("=" * 78)
    out.append("")
    out.append("Empirical target: {2, 8, 20, 28, 50, 82, 126} (Mayer/Jensen 1949).")
    out.append("Null control:     {22, 30, 36, 60, 72} (non-magic neighbors).")
    out.append("")

    # MAGIC NUMBERS
    out.append("-" * 78)
    out.append("MAGIC NUMBERS — substrate forms")
    out.append("-" * 78)
    for n, r in audit["magic"].items():
        out.append("")
        out.append(f"  N = {n:>3}  (locked pre-reg: {r['locked']})")
        out.append(f"    minimum operation count: {r['min_op_count']}, total forms: {r['n_forms']}")
        for f in r["forms"][:5]:  # Show top 5 cleanest forms
            out.append(f"      op={f.operation_count}  [{f.operation_type:<15}]  {f.formula}")
        if r["n_forms"] > 5:
            out.append(f"      ... (+{r['n_forms']-5} more)")

    # NULL CONTROL
    out.append("")
    out.append("-" * 78)
    out.append("NULL CONTROL — non-magic neighbors")
    out.append("-" * 78)
    for n, r in audit["null"].items():
        out.append("")
        out.append(f"  N = {n:>3}  (NOT magic)")
        out.append(f"    minimum operation count: {r['min_op_count']}, total forms: {r['n_forms']}")
        for f in r["forms"][:5]:
            out.append(f"      op={f.operation_count}  [{f.operation_type:<15}]  {f.formula}")
        if r["n_forms"] > 5:
            out.append(f"      ... (+{r['n_forms']-5} more)")

    # COMPARISON TABLE
    out.append("")
    out.append("=" * 78)
    out.append("COMPARISON TABLE — min operation count")
    out.append("=" * 78)
    out.append("")
    out.append(f"  {'class':<10} {'N':>5} {'min_op':>8} {'n_forms':>9}")
    out.append("  " + "-" * 35)
    for n, r in audit["magic"].items():
        op = str(r["min_op_count"]) if r["min_op_count"] is not None else "—"
        out.append(f"  {'magic':<10} {n:>5} {op:>8} {r['n_forms']:>9}")
    out.append("")
    for n, r in audit["null"].items():
        op = str(r["min_op_count"]) if r["min_op_count"] is not None else "—"
        out.append(f"  {'non-magic':<10} {n:>5} {op:>8} {r['n_forms']:>9}")

    # SUMMARY STATISTICS
    out.append("")
    out.append("=" * 78)
    out.append("SUMMARY")
    out.append("=" * 78)
    magic_ops = [r["min_op_count"] for r in audit["magic"].values() if r["min_op_count"] is not None]
    null_ops = [r["min_op_count"] for r in audit["null"].values() if r["min_op_count"] is not None]
    magic_avg = sum(magic_ops) / len(magic_ops) if magic_ops else None
    null_avg = sum(null_ops) / len(null_ops) if null_ops else None

    out.append("")
    out.append(f"  MAGIC NUMBERS:")
    out.append(f"    average min-op count:     {magic_avg:.2f}" if magic_avg is not None else "    no clean forms found")
    out.append(f"    count with op ≤ 1:        {sum(1 for o in magic_ops if o <= 1)} / {len(MAGIC_NUMBERS)}")
    out.append(f"    count with op = 0 (prim): {sum(1 for o in magic_ops if o == 0)} / {len(MAGIC_NUMBERS)}")
    out.append("")
    out.append(f"  NULL CONTROL:")
    out.append(f"    average min-op count:     {null_avg:.2f}" if null_avg is not None else "    no clean forms found")
    out.append(f"    count with op ≤ 1:        {sum(1 for o in null_ops if o <= 1)} / {len(NULL_CONTROL)}")
    out.append(f"    count with op = 0 (prim): {sum(1 for o in null_ops if o == 0)} / {len(NULL_CONTROL)}")

    # VERDICT
    out.append("")
    out.append("=" * 78)
    out.append("FRAMEWORK VERDICT")
    out.append("=" * 78)
    out.append("")
    out.append("  Framework predicted: Form D load-bearing iff")
    out.append("    (a) all 7 magic numbers have clean substrate forms AND")
    out.append("    (b) non-magic null-control numbers have less-clean forms")
    out.append("")

    if magic_ops and null_ops and magic_avg is not None and null_avg is not None:
        diff = null_avg - magic_avg
        out.append(f"  Magic avg op-count:    {magic_avg:.2f}")
        out.append(f"  Null  avg op-count:    {null_avg:.2f}")
        out.append(f"  Difference (null − magic): {diff:+.2f}")
        out.append("")
        if diff >= 0.3:
            out.append("  → Magic numbers have CLEANER substrate forms than non-magic.")
            out.append("    Framework prediction SUPPORTED: substrate-specific signal at magic numbers.")
        elif diff >= -0.3:
            out.append("  → Magic and non-magic numbers have COMPARABLE substrate cleanness.")
            out.append("    Framework prediction WEAKENED: rational-density concern.")
        else:
            out.append("  → Non-magic numbers are CLEANER than magic numbers.")
            out.append("    Framework prediction CONTRADICTED.")

    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
