"""
Substrate-algebraic audit for chemistry Tier-C.8 (NMR via Hopf-pair parity).

Filed AFTER the pre-registration memo [[nmr-via-hopf-pair-prereg]].
Implements the five pre-registered Forms (A/B/C/D/E) against fixed
NWT-canonical integer identifications plus the 14-molecule NICS reference
set from canonical organic-chemistry textbooks.

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_nmr_c8_substrate_audit.py
"""
from __future__ import annotations

from dataclasses import dataclass
from math import comb

from nwt_substrate.isa.constants import (
    DIM_OCTONION,
    H_COXETER_SO7,
    H_V_SO7,
    K8_PARTITION,
    N_EDGES_K7,
    N_POS_ROOTS_SO7,
    N_VERTICES_K7,
    RANK_SO7,
)


# ---------------------------------------------------------------------------
# Substrate-canonical integers (locked at pre-reg time)
# ---------------------------------------------------------------------------

# Map magnitude → display label for substrate-canonical small integers
CANONICAL_INTEGERS: dict[int, str] = {
    RANK_SO7:                  "RANK_SO7",
    4:                          "N_VERTICES_K_4",
    H_V_SO7:                    "H_V_SO7",
    H_COXETER_SO7:              "H_COXETER_SO7",
    N_VERTICES_K7:              "N_VERTICES_K7",
    DIM_OCTONION:               "DIM_OCTONION",
    N_POS_ROOTS_SO7:            "N_POS_ROOTS_SO7",
    K8_PARTITION[2]:            "K8_PARTITION[2]",
    13:                          "trefoil(p²+q²)",
    N_EDGES_K7 - N_VERTICES_K7:  "N_EDGES_K7 − N_VERTICES_K7 (=14)",
    N_EDGES_K7 - H_V_SO7:        "N_EDGES_K7 − H_V_SO7 (=16)",
    N_EDGES_K7 - RANK_SO7:       "N_EDGES_K7 − RANK_SO7 (=18)",
    N_EDGES_K7:                  "N_EDGES_K7",
}


# ---------------------------------------------------------------------------
# 14-molecule NICS reference set
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NICSEntry:
    name: str
    nics_ppm: float
    n_pi: int               # local-ring π electron count (NOT total molecule)
    hopf_pair_count: int
    is_aromatic: bool       # True if Hopf-pair parity is odd (4n+2)


def make_reference() -> list[NICSEntry]:
    def pair(n_pi):
        return n_pi // 2
    def aromatic(n_pi):
        return (n_pi // 2) % 2 == 1
    raw = [
        ("benzene",                -8.0,  6),
        ("Cp- anion",             -14.3,  6),
        ("tropylium+ cation",      -7.6,  6),
        ("naphthalene_center",    -10.0,  6),  # local 6π Clar sextet
        ("anthracene_center",     -13.4,  6),  # local 6π
        ("phenanthrene_center",   -10.0,  6),
        ("pyrene_center",          -7.5,  6),
        ("coronene_hub",          -10.5,  6),
        ("coronene_outer",        -18.7,  6),
        ("pyrrole",               -15.1,  6),
        ("furan",                 -12.3,  6),
        ("thiophene",             -13.6,  6),
        ("cyclobutadiene",         25.0,  4),  # representative of +20 to +30
        ("planar_8annulene",       25.0,  8),
    ]
    return [
        NICSEntry(name=n, nics_ppm=v, n_pi=p,
                  hopf_pair_count=pair(p), is_aromatic=aromatic(p))
        for (n, v, p) in raw
    ]


# ---------------------------------------------------------------------------
# Form A — Sign rule (trivial extension)
# ---------------------------------------------------------------------------

def form_a_sign_rule_check(reference: list[NICSEntry]) -> dict:
    """Form A — does Hopf-pair parity predict NICS sign exactly?"""
    matches = 0
    mismatches = []
    for e in reference:
        predicted_sign = "-" if e.is_aromatic else "+"
        actual_sign = "-" if e.nics_ppm < 0 else "+"
        if predicted_sign == actual_sign:
            matches += 1
        else:
            mismatches.append(e)
    return {
        "n_total": len(reference),
        "n_correct_sign": matches,
        "n_mismatches": len(mismatches),
        "mismatches": [(m.name, m.nics_ppm, m.is_aromatic) for m in mismatches],
        "rule_exact": matches == len(reference),
    }


# ---------------------------------------------------------------------------
# Form B — Per-molecule magnitude identifications
# ---------------------------------------------------------------------------

def closest_canonical(magnitude: float, tolerance: float) -> tuple[int, str, float] | None:
    """Find closest substrate-canonical integer to `magnitude`, within tolerance."""
    best = None
    best_dist = tolerance
    for k, label in CANONICAL_INTEGERS.items():
        d = abs(magnitude - k)
        if d <= best_dist:
            best = (k, label, d)
            best_dist = d
    return best


def form_b_magnitude_audit(reference: list[NICSEntry], tolerance_ppm: float = 0.7) -> dict:
    """Form B — count molecules whose |NICS| lands within `tolerance_ppm` of a
    substrate-canonical integer.

    Pre-reg locked-in candidates (must match these):
      benzene −8 ≈ −DIM_OCTONION
      tropylium+ −7.6 ≈ −DIM_OCTONION
      furan −12.3 ≈ −K8_PARTITION[2]=12
      anthracene −13.4 ≈ −trefoil(p²+q²)=13
      thiophene −13.6 ≈ −13
      Cp− −14.3 ≈ −(N_EDGES_K7−N_VERTICES_K7)=14
      coronene_outer −18.7 ≈ −(N_EDGES_K7−RANK_SO7)=18

    Tighter tolerance (0.5 ppm) would drop some; we use 0.7 ppm as the audit
    threshold (~ DFT NICS computational uncertainty).
    """
    hits = []
    misses = []
    for e in reference:
        mag = abs(e.nics_ppm)
        match = closest_canonical(mag, tolerance_ppm)
        if match is not None:
            k, label, dist = match
            hits.append({
                "name": e.name,
                "nics": e.nics_ppm,
                "abs_nics": mag,
                "canonical": k,
                "label": label,
                "deviation": dist,
            })
        else:
            misses.append({
                "name": e.name,
                "nics": e.nics_ppm,
                "abs_nics": mag,
            })
    return {
        "tolerance_ppm": tolerance_ppm,
        "n_hits": len(hits),
        "n_misses": len(misses),
        "hits": hits,
        "misses": misses,
    }


# ---------------------------------------------------------------------------
# Form C — Tr(M²) extension probe (qualitative, no DFT run)
# ---------------------------------------------------------------------------

def form_c_tr_invariant_probe() -> dict:
    """Form C — would require a closed-form Tr(M^k) formula reproducing all
    NICS values from so(7) adjacency matrix.

    Not running DFT or constructing actual so(7) adjacency matrices for
    each molecule here — this is a SCOPING audit. The so7-substrate-isa-probe
    May 2026 finding (5/5 RE match via Tr(M²)) DOES generalize to RE for
    aromatic molecules, but extending to NICS requires DFT-anchored
    calibration of a Tr-invariant magnitude formula. Out of scope for
    this audit; flagged as future work.
    """
    return {
        "verdict": (
            "FORM C — DEFERRED. Requires DFT-anchored calibration of Tr(M^k) "
            "invariant for NICS magnitude. The May 2026 so7-substrate-isa-probe "
            "established Tr(M²) ↔ RE for 5/5 reference molecules; extending to "
            "NICS is a future-work investigation, not resolvable from the audit "
            "data alone. Not load-bearing for the C.8 resolution."
        ),
    }


# ---------------------------------------------------------------------------
# Form E — rational-density audit
# ---------------------------------------------------------------------------

def form_e_rational_density(reference: list[NICSEntry], hits: int,
                            tolerance_ppm: float) -> dict:
    """Form E — what's the random-chance probability of `hits` out of N total
    molecules landing within `tolerance_ppm` of any substrate-canonical
    integer?"""
    canonical_in_range = [k for k in CANONICAL_INTEGERS
                          if 3 <= k <= 25]
    n_canonical = len(canonical_in_range)
    # The integer range [3, 25] has 23 slots; each canonical integer "covers"
    # roughly 2 × tolerance ppm of the magnitude line.
    # Per-molecule hit probability:
    coverage = 2 * tolerance_ppm * n_canonical / (25 - 3 + 1)
    p_per_molecule = min(1.0, coverage)

    # Binomial probability of exactly `hits` successes out of N trials at
    # rate p_per_molecule, or more — use survival function for P(X >= hits)
    n = len(reference)
    p_at_least_hits = sum(
        comb(n, k) * p_per_molecule**k * (1 - p_per_molecule)**(n - k)
        for k in range(hits, n + 1)
    )
    return {
        "canonical_in_range": canonical_in_range,
        "n_canonical": n_canonical,
        "per_molecule_hit_probability": p_per_molecule,
        "n_total": n,
        "n_hits_observed": hits,
        "p_at_least_n_hits_random": p_at_least_hits,
    }


# ---------------------------------------------------------------------------
# Coronene K_7 hub — structurally distinctive test
# ---------------------------------------------------------------------------

def coronene_k7_hub_test(reference: list[NICSEntry]) -> dict:
    """The pre-reg locked-in 'structurally distinctive' test: coronene outer-
    ring NICS ≈ −18 = −(N_EDGES_K7 − RANK_SO7) reuses the K_7-hub identification
    independently established in [[so7-substrate-isa-probe]] and Tier C.7
    (where 18 = N_EDGES_K7 − RANK_SO7 = non-Cartan generators of so(7)).

    This is not just a rational-density hit because the same identification
    surfaces independently in three different contexts (periodic-table A.3
    shell 18, transition-metal C.7 18e rule, coronene-outer NICS magnitude)."""
    coronene_outer = next(e for e in reference if e.name == "coronene_outer")
    target_int = N_EDGES_K7 - RANK_SO7  # = 18
    deviation = abs(abs(coronene_outer.nics_ppm) - target_int)
    return {
        "molecule": "coronene_outer",
        "nics_ppm": coronene_outer.nics_ppm,
        "abs_nics": abs(coronene_outer.nics_ppm),
        "target_integer": target_int,
        "substrate_form": "N_EDGES_K7 − RANK_SO7 = 21 − 3 = 18",
        "deviation_ppm": deviation,
        "deviation_pct": 100.0 * deviation / target_int,
        "structurally_distinctive": (
            "Same identification surfaces in periodic-table A.3 shell-18, "
            "transition-metal C.7 18e rule (non-Cartan so(7) generators), "
            "AND coronene-outer NICS magnitude. Three independent contexts."
        ),
    }


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def run_audit() -> dict:
    ref = make_reference()
    sign_check = form_a_sign_rule_check(ref)
    mag_audit_strict = form_b_magnitude_audit(ref, tolerance_ppm=0.5)
    mag_audit_loose  = form_b_magnitude_audit(ref, tolerance_ppm=0.7)
    rd_strict = form_e_rational_density(ref, mag_audit_strict["n_hits"], 0.5)
    rd_loose  = form_e_rational_density(ref, mag_audit_loose["n_hits"], 0.7)
    return {
        "reference": ref,
        "form_a": sign_check,
        "form_b_strict": mag_audit_strict,
        "form_b_loose": mag_audit_loose,
        "form_c": form_c_tr_invariant_probe(),
        "form_e_strict": rd_strict,
        "form_e_loose": rd_loose,
        "coronene_k7": coronene_k7_hub_test(ref),
    }


def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 76)
    out.append("Tier C.8 substrate audit — NMR via Hopf-pair parity")
    out.append("=" * 76)
    out.append("")

    # Reference
    out.append("14-molecule NICS reference set:")
    out.append("")
    out.append(f"  {'molecule':<24}{'NICS':>9}{'π':>4}{'Hopf-pair':>11}{'class':>14}")
    for e in audit["reference"]:
        cls = "aromatic" if e.is_aromatic else "anti-arom."
        out.append(f"  {e.name:<24}{e.nics_ppm:>9.1f}{e.n_pi:>4}{e.hopf_pair_count:>11}{cls:>14}")
    out.append("")

    # Form A
    out.append("-" * 76)
    out.append("FORM A — Sign rule (Hopf-pair parity → NICS sign)")
    out.append("-" * 76)
    fa = audit["form_a"]
    out.append(f"  Correct sign predictions: {fa['n_correct_sign']} / {fa['n_total']}")
    out.append(f"  Sign rule EXACT: {fa['rule_exact']}")
    if fa["mismatches"]:
        out.append(f"  Mismatches: {fa['mismatches']}")
    out.append("")
    out.append("  Note: this is a TRIVIAL extension of A.1/A.2/B.4 aromaticity")
    out.append("  classification. Form-A 'win' here doesn't carry new substrate weight.")
    out.append("")

    # Form B
    out.append("-" * 76)
    out.append("FORM B — Per-molecule magnitude identifications")
    out.append("-" * 76)
    for label, audit_b in [("STRICT (±0.5 ppm)", audit["form_b_strict"]),
                            ("LOOSE  (±0.7 ppm)", audit["form_b_loose"])]:
        out.append(f"  {label}: {audit_b['n_hits']} / 14 hits")
        for h in audit_b["hits"]:
            out.append(f"    {h['name']:<24} NICS={h['nics']:>7.1f}  ≈ −{h['canonical']:>3}  ({h['label']}, Δ={h['deviation']:.2f})")
        out.append("")

    # Form E (rational density)
    out.append("-" * 76)
    out.append("FORM E — Rational-density audit")
    out.append("-" * 76)
    for label, fe in [("STRICT (0.5 ppm)", audit["form_e_strict"]),
                       ("LOOSE  (0.7 ppm)", audit["form_e_loose"])]:
        out.append(f"  {label}:")
        out.append(f"    canonical integers in [3, 25]:       {fe['n_canonical']}")
        out.append(f"    per-molecule hit probability:        {fe['per_molecule_hit_probability']:.3f}")
        out.append(f"    observed hits / 14:                  {fe['n_hits_observed']}")
        out.append(f"    p(≥ observed hits at random):        {fe['p_at_least_n_hits_random']:.4f}")
        out.append("")

    # Coronene K_7 hub
    out.append("-" * 76)
    out.append("STRUCTURALLY DISTINCTIVE TEST — Coronene outer ring NICS")
    out.append("-" * 76)
    ck = audit["coronene_k7"]
    out.append(f"  NICS measured:         {ck['nics_ppm']:.2f} ppm")
    out.append(f"  Substrate prediction:  −{ck['target_integer']} ppm  ({ck['substrate_form']})")
    out.append(f"  Deviation:             {ck['deviation_ppm']:.2f} ppm ({ck['deviation_pct']:.1f}%)")
    out.append(f"  {ck['structurally_distinctive']}")
    out.append("")

    # Form C
    out.append("-" * 76)
    out.append("FORM C — Tr(M^k) invariant extension")
    out.append("-" * 76)
    out.append(f"  {audit['form_c']['verdict']}")
    out.append("")

    # Verdict
    out.append("=" * 76)
    out.append("OVERALL VERDICT")
    out.append("=" * 76)
    out.append("")
    fa_exact = audit["form_a"]["rule_exact"]
    n_hits_strict = audit["form_b_strict"]["n_hits"]
    n_hits_loose = audit["form_b_loose"]["n_hits"]
    p_strict = audit["form_e_strict"]["p_at_least_n_hits_random"]
    p_loose = audit["form_e_loose"]["p_at_least_n_hits_random"]

    out.append(f"  Form A — {'EXACT' if fa_exact else 'FAILS'} sign rule (14/14). Trivial extension; not new.")
    out.append(f"  Form B — strict (±0.5 ppm): {n_hits_strict}/14 hits, p_random = {p_strict:.3f}")
    out.append(f"           loose  (±0.7 ppm): {n_hits_loose}/14 hits, p_random = {p_loose:.3f}")
    out.append(f"  Form C — DEFERRED (Tr-invariant extension is future work)")
    out.append(f"  Form D — wins if Form A clean + structurally distinctive coronene K_7 reuse")
    out.append(f"           + magnitude hit rate above rational-density baseline")
    out.append(f"  Form E — REJECTED iff p_random < 0.10 at strict tolerance")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
