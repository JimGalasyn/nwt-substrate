"""
Substrate-algebraic audit for Spectra Tier S.4 — chromophore absorption
quantization on the αR substrate grid.

Filed AFTER pre-registration ([[spectra-s4-chromophores-prereg]]).

Tests whether visible-band absorption maxima of substrate-protected
chromophores (porphyrins, fullerenes, magic clusters, polyenes) quantize
on substrate-canonical multiples of αR = α × Rydberg ≈ 0.0993 eV.

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_spectra_s4_chromophores_audit.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Substrate primitives
# ---------------------------------------------------------------------------

ALPHA = 1.0 / (25.0 * math.pi * math.sqrt(3.0) + 1.0)
M_E_C2_EV = 510998.95
R_INF = 0.5 * ALPHA**2 * M_E_C2_EV   # Rydberg in eV ≈ 13.6057

# Substrate visible-light quantum: α × R
ALPHA_R = ALPHA * R_INF              # ≈ 0.0993 eV


# Substrate-canonical multiplier set in visible-spectrum range
# Locked at pre-registration: primitives + derived primitives + simple products
SUBSTRATE_MULTIPLIERS: dict[int, str] = {
    # Primitives in [15, 50] range relevant to visible spectrum
    16: "N_EDGES_K7 − H_V_SO7 (= 21−5; also 16e rule)",
    18: "N_EDGES_K7 − RANK_SO7 (= 21−3; periodic-table A.3 shell-18)",
    21: "N_EDGES_K7",
    28: "N_EDGES_K8",
    32: "K_7_TRIANGLES − RANK_SO7 (= 35−3; A.3 lanthanide shell)",
    35: "K_7_TRIANGLES",
    # Composite products of primitives
    30: "h_v · h_Coxeter (= 5·6)",
    42: "2 · N_EDGES_K7 (= 2·21; Mackay shell)",
    36: "h_Coxeter² (= 6²)",
    49: "N_VERTICES_K7² (= 7²)",
    25: "h_v² (= 5²)",
    20: "N_VERTICES_K_4 · H_V_SO7 (= 4·5; nuclear magic shell 20)",
    14: "N_EDGES_K7 − N_VERTICES_K7 (= 21−7; G_2 dim)",
    24: "N_EDGES_K7 + RANK_SO7 (= 21+3)",
    15: "H_COXETER_SO7 + N_POS_ROOTS_SO7 (= 6+9)",
}


# ---------------------------------------------------------------------------
# Chromophore reference set (LOCKED, from spectroscopy literature)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Chromophore:
    name: str
    symmetry: str
    band: str
    lambda_nm: float           # absorption maximum in nm
    energy_eV: float           # = 1239.84 / lambda_nm
    locked_multiplier: int     # predicted substrate integer
    is_substrate_protected: bool


def make_reference():
    return [
        # Substrate-protected chromophores
        Chromophore("Free-base porphine H_2P", "D_2h",  "Q",          612, 2.026, 21, True),
        Chromophore("Mg-porphyrin",            "D_4h",  "Q",          580, 2.139, 21, True),
        Chromophore("Zn-porphyrin",            "D_4h",  "Q",          580, 2.139, 21, True),
        Chromophore("Fe-porphyrin (deoxyhemoglobin)", "D_4h", "Q",   575, 2.157, 21, True),
        Chromophore("Cu-porphyrin",            "D_4h",  "Q",          540, 2.296, 23, True),
        Chromophore("Chlorophyll a",           "C_2v",  "Q_y",        660, 1.879, 18, True),
        Chromophore("Bacteriochlorophyll a",   "C_2v",  "Q_y",        770, 1.610, 16, True),
        Chromophore("β-carotene",              "C_2h",  "π-π*",       450, 2.755, 28, True),
        Chromophore("Indigo",                  "C_2v",  "π-π*",       605, 2.049, 21, True),
        Chromophore("C_60 (h_u→t_1u)",         "I_h",   "π-π*",       410, 3.024, 30, True),
        Chromophore("C_60 strong UV",          "I_h",   "h_u→h_g",    256, 4.843, 49, True),

        # Null control: random non-substrate-protected dyes
        Chromophore("Methylene blue",          "(planar)", "n-π*",    668, 1.856, 0, False),
        Chromophore("Crystal violet",          "(chiral)", "π-π*",    590, 2.102, 0, False),
        Chromophore("Fluorescein",             "(acidic)", "π-π*",    494, 2.510, 0, False),
        Chromophore("Methyl orange",           "(azo)",    "n-π*",    462, 2.683, 0, False),
        Chromophore("Eosin Y",                 "(halog)",  "π-π*",    524, 2.367, 0, False),
        Chromophore("Brilliant green",         "(chiral)", "π-π*",    625, 1.984, 0, False),
    ]


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ChromophoreResult:
    chromophore: Chromophore
    multiplier_observed: float          # E / αR
    closest_multiplier_int: int         # nearest substrate-canonical integer
    closest_label: str
    deviation_pct: float                # % deviation from substrate prediction
    locked_match: bool                  # does observed nearest match the pre-reg locked integer?


def closest_substrate(multiplier: float) -> tuple[int, str, float]:
    """Find closest substrate-canonical integer to observed multiplier."""
    best_int = None
    best_dist = float("inf")
    for k, label in SUBSTRATE_MULTIPLIERS.items():
        d = abs(multiplier - k)
        if d < best_dist:
            best_dist = d
            best_int = k
    return best_int, SUBSTRATE_MULTIPLIERS[best_int], best_dist


def audit_chromophore(c: Chromophore) -> ChromophoreResult:
    mult = c.energy_eV / ALPHA_R
    closest, label, dist = closest_substrate(mult)
    expected = c.locked_multiplier if c.is_substrate_protected else closest
    if expected == 0:
        expected = closest  # for null control, just find closest
    predicted_eV = expected * ALPHA_R
    dev_pct = 100 * abs(c.energy_eV - predicted_eV) / c.energy_eV
    locked_match = (closest == c.locked_multiplier) if c.is_substrate_protected else False
    return ChromophoreResult(
        chromophore=c,
        multiplier_observed=mult,
        closest_multiplier_int=closest,
        closest_label=label,
        deviation_pct=dev_pct,
        locked_match=locked_match,
    )


def run_audit():
    refs = make_reference()
    results = [audit_chromophore(c) for c in refs]
    return {
        "alpha_R_eV": ALPHA_R,
        "results": results,
    }


def count_hits(results, tolerance_pct, substrate_only=True):
    return sum(
        1 for r in results
        if r.chromophore.is_substrate_protected == substrate_only
        and r.deviation_pct <= tolerance_pct
    )


def render_report(audit):
    out = []
    out.append("=" * 96)
    out.append("Spectra Tier S.4 — chromophore absorption quantization on αR substrate grid")
    out.append("=" * 96)
    out.append("")
    out.append(f"Substrate visible-light quantum:  αR = {audit['alpha_R_eV']:.5f} eV")
    out.append(f"Substrate canonical multipliers in visible range [15, 50]:")
    for k in sorted(SUBSTRATE_MULTIPLIERS):
        out.append(f"    {k:>3}: {SUBSTRATE_MULTIPLIERS[k]}")
    out.append("")

    # Substrate-protected chromophores
    out.append("-" * 96)
    out.append("SUBSTRATE-PROTECTED CHROMOPHORES")
    out.append("-" * 96)
    out.append("")
    out.append(f"  {'chromophore':<32} {'E (eV)':>8} {'E/αR':>7} {'closest':>7} {'locked':>7} {'dev %':>7} {'match':>6}")
    out.append("  " + "-" * 90)
    for r in audit["results"]:
        if r.chromophore.is_substrate_protected:
            match = "✓" if r.locked_match else "✗"
            out.append(
                f"  {r.chromophore.name:<32} {r.chromophore.energy_eV:>8.3f} {r.multiplier_observed:>7.2f} "
                f"{r.closest_multiplier_int:>7} {r.chromophore.locked_multiplier:>7} {r.deviation_pct:>6.2f}% {match:>6}"
            )
    out.append("")

    # Null control
    out.append("-" * 96)
    out.append("NULL CONTROL (random organic dyes)")
    out.append("-" * 96)
    out.append("")
    out.append(f"  {'chromophore':<32} {'E (eV)':>8} {'E/αR':>7} {'closest':>7} {'dev %':>7}")
    out.append("  " + "-" * 75)
    for r in audit["results"]:
        if not r.chromophore.is_substrate_protected:
            out.append(
                f"  {r.chromophore.name:<32} {r.chromophore.energy_eV:>8.3f} {r.multiplier_observed:>7.2f} "
                f"{r.closest_multiplier_int:>7} {r.deviation_pct:>6.2f}%"
            )
    out.append("")

    # Statistics
    out.append("=" * 96)
    out.append("SUMMARY STATISTICS")
    out.append("=" * 96)
    out.append("")

    sp = [r for r in audit["results"] if r.chromophore.is_substrate_protected]
    null = [r for r in audit["results"] if not r.chromophore.is_substrate_protected]

    for tol in (1.0, 2.0, 3.0, 5.0):
        sp_hits = sum(1 for r in sp if r.deviation_pct <= tol)
        null_hits = sum(1 for r in null if r.deviation_pct <= tol)
        out.append(
            f"  Tolerance ≤ {tol}%:  "
            f"substrate-protected = {sp_hits}/{len(sp)} ({100*sp_hits/len(sp):.0f}%)  "
            f"null = {null_hits}/{len(null)} ({100*null_hits/len(null):.0f}%)"
        )
    out.append("")

    locked_matches = sum(1 for r in sp if r.locked_match)
    out.append(f"  Locked-prediction match (closest int = pre-reg integer): "
               f"{locked_matches}/{len(sp)}")
    out.append("")

    # Verdict
    out.append("=" * 96)
    out.append("FRAMEWORK VERDICT")
    out.append("=" * 96)
    out.append("")
    out.append("  Pre-reg locked criterion:")
    out.append("    Framework PASSES iff:")
    out.append("      (a) ≥6 of 11 substrate-protected hit ≤ 2% (tight)")
    out.append("      (b) ≤1 of 6 null-control hits ≤ 2%")
    out.append("      (c) β-carotene + C_60 strong UV both hit ≤ 1%")
    out.append("")
    sp_tight = sum(1 for r in sp if r.deviation_pct <= 2.0)
    null_tight = sum(1 for r in null if r.deviation_pct <= 2.0)
    beta_caro = next(r for r in sp if "carotene" in r.chromophore.name)
    c60_strong = next(r for r in sp if "strong" in r.chromophore.band or "h_g" in r.chromophore.band)

    cond_a = sp_tight >= 6
    cond_b = null_tight <= 1
    cond_c = beta_caro.deviation_pct <= 1.0 and c60_strong.deviation_pct <= 1.0

    out.append(f"  (a) substrate ≥6/11 tight: {sp_tight}/11  {'✓' if cond_a else '✗'}")
    out.append(f"  (b) null ≤1/6 tight:        {null_tight}/6  {'✓' if cond_b else '✗'}")
    out.append(f"  (c) β-carotene ≤1%:         {beta_caro.deviation_pct:.2f}%  {'✓' if beta_caro.deviation_pct <= 1.0 else '✗'}")
    out.append(f"      C_60 strong UV ≤1%:      {c60_strong.deviation_pct:.2f}%  {'✓' if c60_strong.deviation_pct <= 1.0 else '✗'}")
    out.append("")

    if cond_a and cond_b and cond_c:
        out.append("  → FORM A/B: substrate orbital-locking extends to multi-electron")
        out.append("    chromophores. Visible-spectrum chromophores quantize on αR substrate grid.")
    elif cond_a and (cond_b or cond_c):
        out.append("  → FORM D LOAD-BEARING: substrate signal present; some criteria miss.")
    elif sp_tight > null_tight:
        out.append("  → FORM D WEAK: substrate hit rate higher than null but discrimination weak.")
    else:
        out.append("  → FORM E: rational density wins; substrate signal indistinguishable from null.")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
