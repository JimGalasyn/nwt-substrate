"""
Substrate-algebraic audit for Spectra Tier S.1 — polyatomic vibrational
frequency ratios as cross-scale orbital quantization test.

Filed AFTER pre-registration ([[spectra-s1-vibrational-ratios-prereg]]).
Cross-scale framing: substrate dynamics produces orbital quantization at
atomic scale (hydrogen Rydberg) AND planetary scale (Trappist-1
8:5:3:2:3:4:3 resonance chain on 7 planets); does it also produce
molecular vibrational frequency ratios at substrate-canonical values?

Test set: 11 polyatomic molecules in 4 symmetry classes + 3 null-control
halocarbons. For each molecule, compute all unique frequency ratios and
classify each ratio by closest substrate-canonical value (±5% tolerance,
±1% tight).

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_spectra_s1_vibrational_ratios_audit.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations


# ---------------------------------------------------------------------------
# Substrate-canonical value vocabulary (LOCKED at pre-registration)
# ---------------------------------------------------------------------------

ALPHA = 1.0 / (25.0 * math.pi * math.sqrt(3.0) + 1.0)
SQRT_ALPHA = math.sqrt(ALPHA)
SQRT2 = math.sqrt(2.0)
SQRT3 = math.sqrt(3.0)
SQRT5 = math.sqrt(5.0)
SQRT6 = math.sqrt(6.0)
SQRT7 = math.sqrt(7.0)
PI = math.pi


def _label_for_fraction(num: int, den: int) -> str:
    if den == 1:
        return f"{num}"
    return f"{num}/{den}"


# Substrate primitives
PRIMITIVE_INTEGERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 21, 28, 35]

# Build the substrate-canonical ratio set
SUBSTRATE_RATIOS: dict[float, str] = {}

# Small-integer ratios n/m for n, m in primitives, restricted to reasonable range
for n in PRIMITIVE_INTEGERS:
    for m in PRIMITIVE_INTEGERS:
        ratio = n / m
        if 0.5 <= ratio <= 10.0:
            label = _label_for_fraction(n, m)
            # Keep cleanest label if multiple n/m produce same ratio
            if ratio not in SUBSTRATE_RATIOS or len(label) < len(SUBSTRATE_RATIOS[ratio]):
                SUBSTRATE_RATIOS[ratio] = label

# Substrate-real values
SUBSTRATE_RATIOS[SQRT2] = "√2"
SUBSTRATE_RATIOS[SQRT3] = "√3"
SUBSTRATE_RATIOS[SQRT5] = "√5"
SUBSTRATE_RATIOS[SQRT6] = "√6"
SUBSTRATE_RATIOS[SQRT7] = "√7"
SUBSTRATE_RATIOS[PI] = "π"
SUBSTRATE_RATIOS[PI / 2] = "π/2"
SUBSTRATE_RATIOS[PI / 3] = "π/3"
SUBSTRATE_RATIOS[2 * SQRT3] = "2√3"
SUBSTRATE_RATIOS[SQRT3 / 2] = "√3/2"

SUBSTRATE_RATIOS_SORTED = sorted(SUBSTRATE_RATIOS.items())


def closest_substrate_ratio(r: float, tolerance: float = 0.05) -> tuple[str, float] | None:
    """Find closest substrate-canonical ratio to r within tolerance."""
    best_label = None
    best_dev = float("inf")
    for val, label in SUBSTRATE_RATIOS_SORTED:
        if val == 0:
            continue
        dev = abs(r - val) / val
        if dev < best_dev:
            best_dev = dev
            best_label = label
    if best_label is not None and best_dev <= tolerance:
        return (best_label, best_dev)
    return None


# ---------------------------------------------------------------------------
# Molecule reference set (LOCKED, from NIST WebBook + Herzberg)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Molecule:
    formula: str
    symmetry: str
    sym_class: str       # "linear" / "bent" / "pyramidal" / "tetrahedral" / "null"
    frequencies: tuple[float, ...]   # cm^-1, in numbered order


REFERENCE = [
    # Class 1: linear D_∞h or C_∞v
    Molecule("CO_2",   "D_∞h", "linear",      (1388, 667, 2349)),
    Molecule("CS_2",   "D_∞h", "linear",      (658, 397, 1535)),
    Molecule("C_2H_2", "D_∞h", "linear",      (3372, 1974, 3289, 612, 730)),
    Molecule("HCN",    "C_∞v", "linear",      (3311, 712, 2097)),
    Molecule("OCS",    "C_∞v", "linear",      (859, 520, 2062)),

    # Class 2: bent triatomic C_2v
    Molecule("H_2O",   "C_2v", "bent",        (3657, 1595, 3756)),
    Molecule("H_2S",   "C_2v", "bent",        (2615, 1183, 2626)),
    Molecule("SO_2",   "C_2v", "bent",        (1151, 518, 1362)),

    # Class 3: pyramidal C_3v
    Molecule("NH_3",   "C_3v", "pyramidal",   (3337, 950, 3444, 1627)),
    Molecule("PH_3",   "C_3v", "pyramidal",   (2323, 992, 2328, 1118)),

    # Class 4: tetrahedral T_d
    Molecule("CH_4",   "T_d",  "tetrahedral", (2917, 1534, 3019, 1306)),
    Molecule("SiH_4",  "T_d",  "tetrahedral", (2187, 970, 2191, 914)),

    # Null control: lower symmetry / heavy halocarbons
    Molecule("CH_2Cl_2", "C_2v", "null", (2999, 1467, 1268, 717, 282, 757, 898, 3040, 1153)),
    Molecule("CHCl_3",   "C_3v", "null", (3034, 680, 364, 261, 1219, 770)),
    Molecule("CFCl_3",   "C_3v", "null", (1085, 535, 350, 248, 850, 394)),
]


# ---------------------------------------------------------------------------
# Audit per molecule
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RatioHit:
    ratio: float
    indices: tuple[int, int]
    label: str
    deviation: float
    is_tight: bool      # ≤ 1%


def compute_ratios(freqs: tuple[float, ...]) -> list[RatioHit]:
    """All unique pairwise ratios (i > j, ratio > 1)."""
    out = []
    for i, j in combinations(range(len(freqs)), 2):
        f_i, f_j = freqs[i], freqs[j]
        if f_i == f_j:
            continue
        # Always take ratio ≥ 1 so we don't double-count i,j vs j,i
        if f_i > f_j:
            r = f_i / f_j
            indices = (i + 1, j + 1)
        else:
            r = f_j / f_i
            indices = (j + 1, i + 1)
        match = closest_substrate_ratio(r, tolerance=0.05)
        if match:
            label, dev = match
            out.append(RatioHit(
                ratio=r, indices=indices, label=label,
                deviation=dev, is_tight=(dev <= 0.01),
            ))
        else:
            # Still record the ratio for completeness; mark as no-match
            out.append(RatioHit(
                ratio=r, indices=indices, label="—",
                deviation=float("inf"), is_tight=False,
            ))
    return out


def audit_molecule(mol: Molecule) -> dict:
    ratios = compute_ratios(mol.frequencies)
    n_total = len(ratios)
    n_hit = sum(1 for r in ratios if r.label != "—")
    n_tight = sum(1 for r in ratios if r.is_tight)
    return {
        "molecule": mol,
        "ratios": ratios,
        "n_total": n_total,
        "n_hit": n_hit,
        "n_tight": n_tight,
        "hit_rate": n_hit / n_total if n_total > 0 else 0.0,
        "tight_rate": n_tight / n_total if n_total > 0 else 0.0,
    }


def run_audit() -> dict:
    return {m.formula: audit_molecule(m) for m in REFERENCE}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 84)
    out.append("Spectra Tier S.1 — Polyatomic vibrational frequency ratios")
    out.append("=" * 84)
    out.append("")
    out.append(f"Substrate ratio vocabulary: {len(SUBSTRATE_RATIOS)} canonical values in [0.5, 10]")
    out.append(f"Tolerance: ±5% match (hit); ±1% tight hit")
    out.append("")

    # Per-molecule details
    for sym_class in ("linear", "bent", "pyramidal", "tetrahedral", "null"):
        members = [m for m in REFERENCE if m.sym_class == sym_class]
        if not members:
            continue
        out.append("-" * 84)
        out.append(f"CLASS: {sym_class.upper()}")
        out.append("-" * 84)
        for mol in members:
            a = audit[mol.formula]
            out.append("")
            out.append(f"  {mol.formula} ({mol.symmetry}) — frequencies (cm^-1): {list(mol.frequencies)}")
            out.append(f"    {a['n_hit']}/{a['n_total']} ratios match substrate (≤5%), "
                       f"{a['n_tight']} tight (≤1%)")
            for r in a["ratios"]:
                tight = " *" if r.is_tight else ""
                if r.label != "—":
                    out.append(f"      ν_{r.indices[0]}/ν_{r.indices[1]} = {r.ratio:.4f}  ≈  {r.label}  (dev {100*r.deviation:.2f}%){tight}")
                else:
                    out.append(f"      ν_{r.indices[0]}/ν_{r.indices[1]} = {r.ratio:.4f}  no match in tolerance")
        out.append("")

    # Aggregate by symmetry class
    out.append("=" * 84)
    out.append("AGGREGATE BY SYMMETRY CLASS")
    out.append("=" * 84)
    out.append("")
    out.append(f"  {'class':<14} {'#mol':>5} {'avg hit %':>11} {'avg tight %':>13} {'std hit %':>11}")
    out.append("  " + "-" * 60)

    class_stats = {}
    for sym_class in ("linear", "bent", "pyramidal", "tetrahedral", "null"):
        members = [m for m in REFERENCE if m.sym_class == sym_class]
        if not members:
            continue
        hit_rates = [audit[m.formula]["hit_rate"] for m in members]
        tight_rates = [audit[m.formula]["tight_rate"] for m in members]
        avg_hit = sum(hit_rates) / len(hit_rates)
        avg_tight = sum(tight_rates) / len(tight_rates)
        n = len(hit_rates)
        var = sum((h - avg_hit)**2 for h in hit_rates) / max(1, n - 1)
        std_hit = math.sqrt(var)
        class_stats[sym_class] = {
            "n_mol": n,
            "avg_hit": avg_hit,
            "avg_tight": avg_tight,
            "std_hit": std_hit,
        }
        out.append(f"  {sym_class:<14} {n:>5} {100*avg_hit:>10.1f}% {100*avg_tight:>12.1f}% {100*std_hit:>10.1f}%")
    out.append("")

    # Verdict
    out.append("=" * 84)
    out.append("FRAMEWORK VERDICT")
    out.append("=" * 84)
    out.append("")
    out.append("  Pre-reg locked criterion:")
    out.append("    Framework PASSES iff:")
    out.append("      (a) Linear class avg hit rate ≥ 60% AND")
    out.append("      (b) ≥ 2σ above null control")
    out.append("      (c) Null control hit rate ≤ 30%")
    out.append("")

    if "linear" in class_stats and "null" in class_stats:
        lin = class_stats["linear"]
        null = class_stats["null"]
        cond_a = lin["avg_hit"] >= 0.60
        cond_c = null["avg_hit"] <= 0.30
        if null["std_hit"] > 0:
            sigma_distance = (lin["avg_hit"] - null["avg_hit"]) / null["std_hit"]
        else:
            sigma_distance = float("inf") if lin["avg_hit"] > null["avg_hit"] else 0.0
        cond_b = sigma_distance >= 2.0

        out.append(f"  (a) Linear class hit rate ≥ 60%:    {100*lin['avg_hit']:.1f}%  {'✓' if cond_a else '✗'}")
        out.append(f"  (b) Linear − Null ≥ 2σ:             {sigma_distance:.2f}σ  {'✓' if cond_b else '✗'}")
        out.append(f"  (c) Null control hit rate ≤ 30%:    {100*null['avg_hit']:.1f}%  {'✓' if cond_c else '✗'}")
        out.append("")

        if cond_a and cond_b and cond_c:
            out.append("  → FORM D LOAD-BEARING: cross-scale substrate orbital quantization")
            out.append("    extends to molecular vibrational frequencies in high-symmetry molecules.")
        elif cond_a:
            out.append("  → FORM D PARTIAL: high-symmetry pattern present; null discrimination weak.")
        else:
            out.append("  → FORM E: substrate signal absent; frequencies follow condensate dynamics.")

    # Additional class-by-class comparison
    out.append("")
    out.append("  Symmetry-class hit-rate ranking:")
    sorted_classes = sorted(class_stats.items(), key=lambda kv: -kv[1]["avg_hit"])
    for sym_class, stats in sorted_classes:
        out.append(f"    {sym_class:<14}  {100*stats['avg_hit']:.1f}%  ({stats['n_mol']} molecules)")
    out.append("")

    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
