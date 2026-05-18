"""
Substrate-algebraic audit for Au_n magic-cluster forward test of the
Form prediction framework.

Filed AFTER the pre-registration memo [[au-magic-clusters-prereg]].
Implements the locked-in substrate identifications + rational-density
audit + framework verdict.

This is NOT a 9th chemistry-sector tier item — it is a FRAMEWORK FORWARD
TEST. The chemistry-sector survey is officially done at 8 items.

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_au_magic_cluster_audit.py
"""
from __future__ import annotations

from dataclasses import dataclass
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
# Substrate-canonical integer set (locked at pre-reg time)
# ---------------------------------------------------------------------------

CANONICAL_INTEGERS: dict[int, str] = {
    1:                          "unit (central atom)",
    RANK_SO7:                   "RANK_SO7",
    4:                          "N_VERTICES_K_4",
    H_V_SO7:                    "H_V_SO7",
    H_COXETER_SO7:              "H_COXETER_SO7",
    N_VERTICES_K7:              "N_VERTICES_K7",
    DIM_OCTONION:               "DIM_OCTONION",
    N_POS_ROOTS_SO7:            "N_POS_ROOTS_SO7",
    K8_PARTITION[2]:            "K8_PARTITION[2]",
    13:                         "trefoil(p²+q²)",
    N_EDGES_K7 - N_VERTICES_K7: "N_EDGES_K7 − N_VERTICES_K7 (=14)",
    N_EDGES_K7 - H_V_SO7:       "N_EDGES_K7 − H_V_SO7 (=16)",
    N_EDGES_K7 - RANK_SO7:      "N_EDGES_K7 − RANK_SO7 (=18)",
    N_EDGES_K7:                 "N_EDGES_K7",
    N_EDGES_K8:                 "N_EDGES_K8",
    32:                         "K_7_TRIANGLES − RANK_SO7 (=32)",
    35:                         "K_7_TRIANGLES",
}


# ---------------------------------------------------------------------------
# Au_n magic-cluster reference set
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AuClusterEntry:
    n_au: int
    structure: str
    reference: str
    is_mackay: bool


REFERENCE_CLUSTERS: list[AuClusterEntry] = [
    AuClusterEntry(13,  "Mackay M_1: 1 + 12",        "theoretical + small-cluster MS", True),
    AuClusterEntry(55,  "Mackay M_2: 1 + 12 + 42",   "Schmid 1981 (PPh_3)_12 Cl_6",   True),
    AuClusterEntry(102, "non-Mackay (icos core + ligand shell)", "Kornberg 2007 thiolate", False),
    AuClusterEntry(144, "60-faceted polyhedron + interior",     "Häkkinen 2008 Au_144(SR)_60", False),
    AuClusterEntry(147, "Mackay M_3: 1 + 12 + 42 + 92", "mass spec + theoretical", True),
    AuClusterEntry(309, "Mackay M_4",                  "mass spec",                     True),
    AuClusterEntry(561, "Mackay M_5",                  "larger nanoparticles",         True),
]


def mackay_magic_number(n: int) -> int:
    """M_n = (10n³ + 15n² + 11n + 3) / 3."""
    return (10 * n**3 + 15 * n**2 + 11 * n + 3) // 3


def mackay_shell_size(n: int) -> int:
    """Number of atoms in the n-th Mackay shell (n >= 1)."""
    return 12 + 30 * (n - 1) + 10 * (n - 1) * (n - 2)


# ---------------------------------------------------------------------------
# Locked-in substrate identifications (pre-reg)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LockedIdentification:
    cluster_n: int
    formula: str
    terms: tuple[tuple[int, str], ...]   # ((integer, label), ...)
    distinctiveness: str


LOCKED_IDENTIFICATIONS: list[LockedIdentification] = [
    LockedIdentification(
        cluster_n=13,
        formula="1 + K8_PARTITION[2] = 1 + 12",
        terms=((1, "unit"), (12, "K8_PARTITION[2]")),
        distinctiveness="DOUBLE ROUTE: also 13 = trefoil(p²+q²) = 2² + 3² (Paper 13)",
    ),
    LockedIdentification(
        cluster_n=55,
        formula="1 + K8_PARTITION[2] + 2·N_EDGES_K7 = 1 + 12 + 42",
        terms=((1, "unit"), (12, "K8_PARTITION[2]"), (21, "N_EDGES_K7"), (21, "N_EDGES_K7")),
        distinctiveness="COMPOSITE: shell decomposition matches Mackay (1 + 12 + 42) AND substrate decomposition exactly",
    ),
    LockedIdentification(
        cluster_n=144,
        formula="K8_PARTITION[2]² = 12² = 144",
        terms=((12, "K8_PARTITION[2]"), (12, "K8_PARTITION[2]")),
        distinctiveness="DIRECT SQUARED: first substrate-integer SQUARED producing a cluster magic number. Distinctive operation new to forward test.",
    ),
]


def verify_locked_identification(loc: LockedIdentification) -> bool:
    """Sum of terms must equal cluster_n (for additive identifications).
    For squared identification (Au_144), product of terms must equal cluster_n."""
    if loc.cluster_n == 144:
        # Squared: product
        product_val = 1
        for (v, _) in loc.terms:
            product_val *= v
        return product_val == loc.cluster_n
    return sum(v for (v, _) in loc.terms) == loc.cluster_n


# ---------------------------------------------------------------------------
# Rational-density audit — enumerate substrate-canonical decompositions
# ---------------------------------------------------------------------------

def enumerate_additive_decompositions(target: int, max_terms: int = 4) -> list[tuple[int, ...]]:
    """Enumerate ALL ways to write `target` as a sum of (at most max_terms)
    substrate-canonical integers (with repetition allowed)."""
    canonical = sorted(k for k in CANONICAL_INTEGERS if k > 0)
    out = []
    for k in range(1, max_terms + 1):
        for combo in combinations_with_replacement(canonical, k):
            if sum(combo) == target:
                out.append(combo)
    return out


def enumerate_squared_decompositions(target: int) -> list[tuple[int, int]]:
    """Enumerate (a, b) with a, b substrate-canonical, a*b = target."""
    canonical = sorted(k for k in CANONICAL_INTEGERS if k > 0)
    out = []
    for a in canonical:
        if target % a == 0:
            b = target // a
            if b in CANONICAL_INTEGERS:
                out.append((a, b))
    return out


def mackay_shell_decomposition(cluster_n: int) -> Optional[list[int]]:
    """If cluster_n is a Mackay magic number, return the shell decomposition
    [1, 12, 42, 92, 162, ...] up to that level."""
    n = 1
    cumulative = 1
    shells = [1]
    while cumulative < cluster_n:
        shell = mackay_shell_size(n)
        shells.append(shell)
        cumulative += shell
        if cumulative == cluster_n:
            return shells
        n += 1
        if n > 20:
            break
    return None


# ---------------------------------------------------------------------------
# Audit per cluster
# ---------------------------------------------------------------------------

def audit_cluster(entry: AuClusterEntry) -> dict:
    """Run the audit on a single cluster."""
    n = entry.n_au
    decomps_short = enumerate_additive_decompositions(n, max_terms=3)
    decomps_long = enumerate_additive_decompositions(n, max_terms=4)
    squared = enumerate_squared_decompositions(n)
    mackay_shells = mackay_shell_decomposition(n) if entry.is_mackay else None

    return {
        "cluster": entry,
        "n_short_decomp_3terms": len(decomps_short),
        "n_long_decomp_4terms": len(decomps_long),
        "squared_decompositions": squared,
        "mackay_shells": mackay_shells,
        "short_decomps_sample": decomps_short[:5],
        "long_decomps_sample": decomps_long[:5],
    }


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def run_audit() -> dict:
    locked_results = []
    for loc in LOCKED_IDENTIFICATIONS:
        locked_results.append({
            "locked": loc,
            "verified": verify_locked_identification(loc),
        })

    per_cluster = [audit_cluster(e) for e in REFERENCE_CLUSTERS]

    return {
        "locked_identifications": locked_results,
        "per_cluster": per_cluster,
    }


def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 76)
    out.append("Framework Forward Test 1 — Au_n magic-cluster audit")
    out.append("=" * 76)
    out.append("")

    # Locked identifications
    out.append("LOCKED-IN SUBSTRATE IDENTIFICATIONS (pre-reg)")
    out.append("")
    for r in audit["locked_identifications"]:
        loc = r["locked"]
        status = "VERIFIED ✓" if r["verified"] else "FAILED ✗"
        out.append(f"  Au_{loc.cluster_n}: {loc.formula}  [{status}]")
        out.append(f"    Terms: {loc.terms}")
        out.append(f"    {loc.distinctiveness}")
        out.append("")

    # Per-cluster rational-density audit
    out.append("-" * 76)
    out.append("RATIONAL-DENSITY AUDIT — substrate decompositions per cluster")
    out.append("-" * 76)
    out.append("")
    out.append(f"  {'n_Au':>5} {'Mackay?':<8} {'3-term decomps':>15} {'4-term decomps':>15} {'squared':>10} {'Mackay shells':<25}")
    for c in audit["per_cluster"]:
        entry = c["cluster"]
        mackay_str = "yes" if entry.is_mackay else "no"
        shells_str = str(c["mackay_shells"]) if c["mackay_shells"] else "—"
        out.append(f"  {entry.n_au:>5} {mackay_str:<8} {c['n_short_decomp_3terms']:>15} "
                   f"{c['n_long_decomp_4terms']:>15} {len(c['squared_decompositions']):>10} "
                   f"{shells_str:<25}")
    out.append("")

    # Detail per locked cluster
    out.append("-" * 76)
    out.append("DETAIL — Au_13, Au_55, Au_144 locked-in identification check")
    out.append("-" * 76)
    for c in audit["per_cluster"]:
        if c["cluster"].n_au not in (13, 55, 144, 147, 309):
            continue
        entry = c["cluster"]
        out.append(f"")
        out.append(f"  Au_{entry.n_au} ({entry.structure}):")
        out.append(f"    {c['n_short_decomp_3terms']} substrate decompositions in ≤3 terms")
        out.append(f"    {c['n_long_decomp_4terms']} substrate decompositions in ≤4 terms")
        if c["squared_decompositions"]:
            out.append(f"    Squared decompositions: {c['squared_decompositions']}")
        if c["mackay_shells"]:
            out.append(f"    Mackay shell structure: {c['mackay_shells']}")
            shells_canonical = [(s, CANONICAL_INTEGERS.get(s, "—")) for s in c["mackay_shells"]]
            out.append(f"    Shell-by-shell canonical labels: {shells_canonical}")
    out.append("")

    # Framework verdict
    out.append("=" * 76)
    out.append("FRAMEWORK FORWARD TEST VERDICT")
    out.append("=" * 76)
    out.append("")

    all_verified = all(r["verified"] for r in audit["locked_identifications"])
    out.append(f"  Locked identifications verified: {all_verified}")
    for r in audit["locked_identifications"]:
        loc = r["locked"]
        status = "✓" if r["verified"] else "✗"
        out.append(f"    [{status}] Au_{loc.cluster_n} = {loc.formula}")
    out.append("")

    # Window narrowness check
    # Au_13, Au_55, Au_144 should have specific Mackay-aligned or squared
    # identifications.  Au_147 and Au_309 should weaken.
    au_13_decomps = next(c for c in audit["per_cluster"] if c["cluster"].n_au == 13)
    au_55_decomps = next(c for c in audit["per_cluster"] if c["cluster"].n_au == 55)
    au_144_decomps = next(c for c in audit["per_cluster"] if c["cluster"].n_au == 144)
    au_147_decomps = next(c for c in audit["per_cluster"] if c["cluster"].n_au == 147)
    au_309_decomps = next(c for c in audit["per_cluster"] if c["cluster"].n_au == 309)

    out.append("  Window narrowness (Form-D subset prediction):")
    out.append(f"    Au_13  decomps (≤3 terms): {au_13_decomps['n_short_decomp_3terms']}")
    out.append(f"    Au_55  decomps (≤3 terms): {au_55_decomps['n_short_decomp_3terms']}")
    out.append(f"    Au_144 decomps (≤3 terms): {au_144_decomps['n_short_decomp_3terms']}")
    out.append(f"    Au_147 decomps (≤3 terms): {au_147_decomps['n_short_decomp_3terms']}")
    out.append(f"    Au_309 decomps (≤3 terms): {au_309_decomps['n_short_decomp_3terms']}")
    out.append("")
    out.append("  Framework predicts: load-bearing at small clusters (Au_13-Au_55-Au_144)")
    out.append("                      with degradation at larger Au_147 / Au_309")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
