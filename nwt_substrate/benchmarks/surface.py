"""The pre-registration KILL SURFACE — every dimensionless-constant claim, scored
at the experiment's own precision, with provenance disputes wired to the pinned
external audits.

This module is the canonical CLAIM REGISTRY for the dimensionless-constants
pre-registration (null-worldtube-private
``analysis/DIMENSIONLESS_CONSTANTS_PREREG_draft.md``).  Its contract extends
``predict.py``'s O10 standalone rung in three ways the smoke-level readouts
cannot provide:

1. **S-NOW** (precision confrontation): each row's witness carries a 1σ
   experimental uncertainty, and the verdict is computed at that precision —
   EXACT-COMPATIBLE or DEAD-AS-EXACT — not at a blanket cit tolerance.  A row
   can pass cit at 1% and still be dead by tens of thousands of σ (α itself
   is).  A pass is *compatibility*, never confirmation: every witness value
   below was in hand when the closed forms were found (the contamination
   ledger of the prereg), so S-NOW can only KILL, not confirm.

2. **S-FORWARD** (the only confirmation channel): ``WITNESS_UPDATES`` is an
   APPEND-ONLY register of measured values published after the frozen witness
   set (CODATA-2018 / PDG / Planck-2018).  Each update re-scores the FROZEN
   prediction and reports whether the measurement moved toward or away from
   it.  No formula, order, or witness may change in response to an update —
   a row that needs a new term to survive one is DEAD (logged, no rescue).

3. **Provenance disputes** (the self-tag problem): a provenance lint the
   claim's author satisfies by tagging everything DERIVED is theater.  Rows
   whose self-declared provenance is contested by a PINNED external audit
   carry the dispute on their symbolic node; ``dispute_audit`` suspends their
   cit passes as corroboration until the memory-blind Auditor adjudicates
   (null-worldtube-private ``analysis/PROGRAM_CHARTER.md``; ``~/repos/nwt-audit``).

Scope rule (inherited from ``predict.py``): only DIMENSIONLESS quantities that
are pure functions of the substrate structure (α = 1/(25π√3 + 1) + the
K_7/Spin(7)/walk integers).  Nothing measured enters a prediction; measured
values live only in witnesses.  Quantities excluded from the surface are
listed in ``EXCLUSIONS`` with dated reasons — a silent omission of an
inconvenient row is itself a look-elsewhere sin (prereg §3 rule).

Frozen-order rule: each row pins the ORDER of its closed form (``order``
field).  Adding a correction term after freeze to survive data is forbidden;
the ``test_surface.py`` order-pin test turns any change into a dated,
deliberate amendment rather than a silent retrofit.

Run ``python -m nwt_substrate.benchmarks.surface`` for the full readout
(acceptance invariants, provenance disputes, S-NOW verdict table, S-FORWARD
drift table, open anti-edges).  This output IS the prereg's §3 table — no
hand-typed numbers downstream of this module.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .predict import REFERENCE, predictions
from .o10 import (
    DERIVED, FITTED, MOTIVATED, POST_SELECTED, STATUS_DEFERRED_BRIDGE,
    DerivationDAG, Stage,
    _collapse_lines, _dispute_lines, _provenance_lines, _snow_lines,
)


# ===========================================================================
# Witnesses — QUARANTINED measured values (never read by any prediction).
# Frozen set: CODATA-2018 / PDG / Planck-2018 (the predict.py convention: the
# last fully *measured* set; post-SI2019 defined values excluded).  Sources
# quoted per entry; sigmas are 1σ absolute.
# ===========================================================================

@dataclass(frozen=True)
class Witness:
    value: float
    sigma: float
    source: str


# --- headline constants: values reuse predict.REFERENCE verbatim (single
#     source of truth); this table only ADDS the 1σ error budget. ---
HEADLINE_SIGMA: dict[str, tuple[float, str]] = {
    "inv_alpha":      (2.1e-8,  "CODATA-2018 α⁻¹ = 137.035999084(21)"),
    "sin2_theta_W":   (4.0e-5,  "PDG effective leptonic sin²θ_W = 0.23122(4) "
                                "(scheme mismatch vs the LO on-shell form is the "
                                "marked cit defect — see predict.py)"),
    "cabibbo_lambda": (6.7e-4,  "PDG-2022 CKM global fit λ = 0.22500(67)"),
    "eta_B":          (4.0e-12, "Planck-2018 (Ω_b h² = 0.02237(15) → η_B ≈ 6.1(4e-2)e-10)"),
    # m_e/M_Pl uncertainty is G-dominated: CODATA-2018 G = 6.67430(15)e-11
    # (22 ppm) → M_Pl 11 ppm → σ_rel(m_e/M_Pl) = 1.1e-5.
    "m_e_over_M_Pl":  (4.185e-23 * 1.1e-5, "CODATA-2018 via G = 6.67430(15)e-11 (22 ppm)"),
}

WITNESSES: dict[str, Witness] = {
    key: Witness(REFERENCE[key], sig, src)
    for key, (sig, src) in HEADLINE_SIGMA.items()
}

# m_e/M_Pl: REFERENCE's 6-s.f. value quantizes at ±5e-28 — comparable to its own
# σ, so the S-NOW verdict would hang on a rounding artifact.  Recompute the
# witness from the pinned CODATA-2018 primaries instead (m_e c² = 0.51099895000(15)
# MeV; M_Pl c² = 1.220890(14)e19 GeV — the σ is G-dominated, 11 ppm).
_ME_MEV_PDG = 0.51099895000        # CODATA-2018 m_e, MeV (σ = 1.5e-10 MeV: negligible)
_M_PL_GEV_CODATA18 = 1.220890e19
_M_PL_SIGMA_GEV = 0.000014e19
WITNESSES["m_e_over_M_Pl"] = Witness(
    _ME_MEV_PDG * 1e-3 / _M_PL_GEV_CODATA18,
    (_ME_MEV_PDG * 1e-3 / _M_PL_GEV_CODATA18)
    * (_M_PL_SIGMA_GEV / _M_PL_GEV_CODATA18),
    "CODATA-2018 m_e = 0.51099895000(15) MeV / M_Pl = 1.220890(14)e19 GeV "
    "(σ G-dominated, 11 ppm)")

# --- cosmology beyond predict.py ---
WITNESSES["omega_b_c"] = Witness(
    0.02237 / 0.1200,                       # Planck-2018 Ω_b h² / Ω_c h²
    (0.02237 / 0.1200) * math.hypot(0.00015 / 0.02237, 0.0012 / 0.1200),
    "Planck-2018 TT,TE,EE+lowE+lensing: Ω_b h² = 0.02237(15), Ω_c h² = 0.1200(12)")
WITNESSES["rho_lambda"] = Witness(
    1.20e-123, 0.036e-123,
    "Planck-2018 ρ_Λ/M_Pl⁴ = 1.20(4)e-123 (Ω_Λ = 0.6847(73) + H₀; ~3%)")

# --- electroweak scale ratio (dimensionless: v_EW / m_e) ---
_V_EW_MEV_PDG = 246.21965e3        # PDG via G_F = 1.1663787(6)e-5 GeV⁻² (σ_rel 2.6e-7)
WITNESSES["v_over_m_e"] = Witness(
    _V_EW_MEV_PDG / _ME_MEV_PDG,
    (_V_EW_MEV_PDG / _ME_MEV_PDG) * 2.6e-7,
    "PDG G_F = 1.1663787(6)e-5 GeV⁻² → v = 246.21965(6) GeV; CODATA-2018 m_e")

# --- neutrino mixing (PDG-2022 / NuFIT global fit, normal ordering) ---
WITNESSES["sin2_theta_13"] = Witness(0.0220, 0.0007, "PDG-2022 sin²θ₁₃ = 0.0220(7)")
WITNESSES["sin2_theta_12"] = Witness(0.307, 0.013, "PDG-2022 sin²θ₁₂ = 0.307(13)")
WITNESSES["sin2_theta_23"] = Witness(0.546, 0.021, "PDG-2022 sin²θ₂₃ = 0.546(21) (upper octant)")

# --- particle masses (PDG-2022, MeV) for the mass-ratio block.  The witness
#     for a ratio row m_X/m_e is MASS_MEV[X].value / m_e with σ propagated
#     from the numerator (m_e's σ is 8 orders below every entry here). ---
MASS_MEV: dict[str, Witness] = {
    "e-":      Witness(0.51099895000, 1.5e-10, "CODATA-2018"),
    "mu-":     Witness(105.6583755, 2.3e-6, "PDG-2022"),
    "tau-":    Witness(1776.86, 0.12, "PDG-2022"),
    "pi+":     Witness(139.57039, 0.00018, "PDG-2022"),
    "pi0":     Witness(134.9768, 0.0005, "PDG-2022"),
    "K+":      Witness(493.677, 0.016, "PDG-2022"),
    "K0":      Witness(497.611, 0.013, "PDG-2022"),
    "eta":     Witness(547.862, 0.017, "PDG-2022"),
    "rho":     Witness(775.26, 0.23, "PDG-2022 ρ(770)"),
    "omega":   Witness(782.66, 0.13, "PDG-2022 ω(782)"),
    "p":       Witness(938.27208816, 2.9e-7, "CODATA-2018"),
    "n":       Witness(939.56542052, 5.4e-7, "CODATA-2018"),
    "Sigma+":  Witness(1189.37, 0.07, "PDG-2022"),
    "Sigma0":  Witness(1192.642, 0.024, "PDG-2022"),
    "Sigma-":  Witness(1197.449, 0.030, "PDG-2022"),
    "Lambda":  Witness(1115.683, 0.006, "PDG-2022"),
    "Xi0":     Witness(1314.86, 0.20, "PDG-2022"),
    "Xi-":     Witness(1321.71, 0.07, "PDG-2022"),
    "Delta":   Witness(1232.0, 2.0, "PDG-2022 Δ(1232) Breit-Wigner"),
    "Sigma*":  Witness(1383.7, 1.0, "PDG-2022 Σ(1385)⁰ Breit-Wigner"),
    "Omega-":  Witness(1672.45, 0.29, "PDG-2022"),
    "D+":      Witness(1869.66, 0.05, "PDG-2022"),
    "D0":      Witness(1864.84, 0.05, "PDG-2022"),
    "J/psi":   Witness(3096.900, 0.006, "PDG-2022 J/ψ(1S)"),
    "Upsilon": Witness(9460.30, 0.26, "PDG-2022 Υ(1S)"),
}


# ===========================================================================
# Exclusions — every substrate-claimed quantity NOT on the surface, with a
# dated reason.  Silent omission of an inconvenient row is a look-elsewhere
# sin; this dict is the honesty ledger the Auditor sweep (prereg rider R2)
# checks against the library.
# ===========================================================================

EXCLUSIONS: dict[str, str] = {
    "alpha_s(M_Z)": "2026-07-12: no substrate closed form in the library — "
                    "qcd.alpha_s is the PDG measured input 0.1179.  (A 16α "
                    "memory-claim exists in program notes but has no code "
                    "authority; if it lands as code it must enter as a new row, "
                    "not a retrofit.)",
    "G_Newton":     "2026-07-12: dimensionless content identical to m_e_over_M_Pl "
                    "(G m_e²/ħc = (m_e/M_Pl)²) — a duplicate row would double-count "
                    "one prediction.",
    "G_F":          "2026-07-12: fixed by v_over_m_e (G_F = 1/(√2 v²)) — duplicate.",
    "quark masses": "2026-07-12: scheme-dependent (MS-bar) witnesses; no compendium "
                    "entries; deferred to a v2 row set with a pinned scheme.",
    "CKM A, ρ̄, η̄": "2026-07-12: substrate_ckm exists but scheme/convention pinning "
                    "is unresolved; only the Wolfenstein λ row is frozen in v1.",
    "K_8 DM tower": "2026-07-12: forward-only — no measured witness exists (LZ-G3 / "
                    "XRISM watch).  Kept out of S-NOW by definition; belongs to "
                    "S-FORWARD if a mass is ever measured.",
    "N_eff":        "2026-07-12: no substrate closed form in the library.",
    "Koide relations": "2026-07-12: program-note claims only; no code authority.",
}


# ===========================================================================
# Provenance ADJUDICATIONS — the four disputes deposited to the memory-blind
# Auditor were adjudicated 2026-07-12 (VERDICT, nwt-audit
# `audits/2026-07-12-constants-provenance-disputes/`, commit 4032d04):
# self-tag DERIVED **REFUTED on every row**.  Per the o10 dispute semantics,
# adjudication replaces the tag and clears the `disputed` field, citing the
# verdict — which these records do.  The full contest evidence (the two
# exhibits: the gravity/coupling.py alpha=ALPHA_QED default behind the
# public "G at −11 ppm" headline, substrate-pure ≈7σ outside; and the
# neighbouring_value.py look-elsewhere volume, ~83% of random targets inside
# G's bar) was REPRODUCED by the Auditor (CL-3(a),(b) CONFIRMED, with
# independent seed/menu/stage-count robustness probes) and is preserved in
# the deposit's CLAIM.md and the verdict's work/ directory.
# ===========================================================================

_VERDICT = ("Auditor verdict 2026-07-12-constants-provenance-disputes "
            "(nwt-audit 4032d04)")

# Adjudicated provenance per row, with the verdict's per-row ground.
ADJUDICATIONS: dict[str, tuple[str, str]] = {
    "eta_B": (POST_SELECTED,
              f"{_VERDICT}: DERIVED refuted — kept-closest survivor of a "
              "17-formula sweep (two pinned admissions); docstring's own "
              "multiple-readings for 3 and 14 defeat forcing; mode 5: 54.9% "
              "of random targets fit inside this row's 2σ bar.  Open caveat: "
              "pin the primary sweep artifact by dated CLAIM.md appendix."),
    "omega_b_c": (FITTED,
                  f"{_VERDICT}: DERIVED refuted — the LO term alone was "
                  "already compatible at 1.77σ (nothing forced a correction); "
                  "the added 75α² term lands at z = 0.0056σ (≈0.4% luck), the "
                  "fit-to-central-value signature; mode 5: 81.7%."),
    "m_e_over_M_Pl": (MOTIVATED,
                      f"{_VERDICT}: DERIVED refuted — Paper 17 verified "
                      "verbatim (exponent counted, √α matched via the "
                      "non-rigorous CS integral, 8/7 identified, NLO "
                      "empirical); both exhibits reproduced.  The NNLO α² "
                      "coefficient is FITTED by documentation (computed from "
                      "CODATA, nearest named integer selected) — hence the "
                      "CL-2 NLO order pin."),
    "rho_lambda": (MOTIVATED,
                   f"{_VERDICT}: MOTIVATED inherited via (m_e/M_Pl)⁴; the "
                   "row-specific factors h_Cox = 6 and α¹⁶ are "
                   "UNDERSPECIFIED (no derivation AND no selection history "
                   "in pinned evidence; h = 6 is the unique menu integer "
                   "inside the bar) — treated as contested-upheld until "
                   "discharged by dated amendment; mode 5 at this row's 3% "
                   "bar: 100.0%."),
}


# ===========================================================================
# The surface rows
# ===========================================================================

@dataclass(frozen=True)
class SurfaceRow:
    key: str            # unique row key (DAG node names derive from it)
    section: str        # couplings / cosmology / electroweak / neutrino / mass-ratio / defect
    form: str           # the frozen symbolic form (human-readable)
    order: str          # FROZEN order of the closed form ("exact", "LO", "NLO", "NNLO")
    predicted: float    # computed FROM THE LIBRARY — never hand-typed
    witness: Witness
    provenance: str     # self-declared value provenance (explicit — no defaults here)
    disputed: str = ""  # open external-audit dispute; "" after adjudication
    #                     (the four 2026-07-12 disputes are adjudicated — see
    #                     ADJUDICATIONS, whose grounds ride on `note`)
    note: str = ""
    discharge: str = "" # for deferred-bridge rows: the killable test that closes them


# The frozen order pins.  test_surface.py asserts this dict verbatim: changing
# an order is a dated amendment, never a silent retrofit.
ORDER_PINS: dict[str, str] = {
    "inv_alpha": "exact",
    "sin2_theta_W": "LO",
    "cabibbo_lambda": "exact",
    "eta_B": "exact",
    "m_e_over_M_Pl": "NLO",       # PINNED by the Auditor (CL-2, 2026-07-12):
    #                               the L4(a)-audited form; NNLO retired from
    #                               claim status (coefficient computed from the
    #                               CODATA target — documented target-selection).
    "omega_b_c": "NLO",           # 25α + 75α² — the NLO term adjudicated FITTED.
    "rho_lambda": "NLO",          # follows the CL-2 pin via (m_e/M_Pl)⁴
    #                               (lambda_cc.py call site updated to match).
    "v_over_m_e": "NLO",          # (25/α²)(1 + 25α/(4√3)).
    "sin2_theta_13": "LO",
    "sin2_theta_12": "LO",
    "sin2_theta_23": "LO",
    "paper6_mass_ratio": "exact", # the Paper-6 closed form (per-particle integers
    #                               are the POST_SELECTED part, not the order).
}


def _mass_ratio_rows() -> list[SurfaceRow]:
    """One row per compendium particle: m_X/m_e from the Paper-6 walk formula.

    Provenance is POST_SELECTED — the per-particle walk integers (p, q, m, n_q)
    are assignments, and the library's own compendium records that the nucleon
    family was CORRECTED (1,4,5,3)→(1,3,5,5) on 2026-04-30 to reduce the mass
    residual: the integers were chosen with the observed masses in hand.  The
    discharge (carried by the DAG anti-edge) is a forcing chain that pins the
    integers from charge/isospin/baryon predicates WITHOUT mass input."""
    from ..particles import list_particles, particle
    from ..particles.mass import ME_MEV

    rows = []
    for name in list_particles():
        if name == "e-":
            continue                      # the unit — ratio 1 by construction
        p = particle(name)
        w_mass = MASS_MEV[name]
        m_e = MASS_MEV["e-"].value
        rows.append(SurfaceRow(
            key=f"m_{name}_over_m_e", section="mass-ratio",
            form=f"paper6({p.p},{p.q},{p.m},{p.n_q})",
            order=ORDER_PINS["paper6_mass_ratio"],
            predicted=p.mass_pred / ME_MEV,
            witness=Witness(w_mass.value / m_e, w_mass.sigma / m_e,
                            f"{w_mass.source} m_{name} = {w_mass.value} MeV / CODATA m_e"),
            provenance=POST_SELECTED,
            note="walk integers assigned with masses in hand "
                 "(compendium 2026-04-30 correction note)"))
    return rows


def surface_rows() -> list[SurfaceRow]:
    """The full pre-registration surface, computed from the library."""
    pred = predictions()
    from ..cosmology.lambda_cc import lambda_cc
    from ..cosmology.omega_b_c import omega_b_c
    from ..isa.constants import ALPHA_SUBSTRATE as a
    from ..particles import particle
    from ..particles.mass import ME_MEV

    rows: list[SurfaceRow] = []

    # --- couplings & headline constants (predict.py's five, σ-graded) ---
    forms = {
        "inv_alpha": "25π√3 + 1",
        "sin2_theta_W": "(2 + α)/9",
        "cabibbo_lambda": "√(7α)",
        "eta_B": "3α⁴/14",
        "m_e_over_M_Pl": "(8/7)·(1 + α/7)·α^(21/2)",
    }
    sections = {
        "inv_alpha": "couplings", "sin2_theta_W": "couplings",
        "cabibbo_lambda": "couplings", "eta_B": "cosmology",
        "m_e_over_M_Pl": "couplings",
    }
    for key in forms:
        prov, ground = ADJUDICATIONS.get(key, (DERIVED, ""))
        rows.append(SurfaceRow(
            key=key, section=sections[key], form=forms[key],
            order=ORDER_PINS[key], predicted=pred[key],
            witness=WITNESSES[key], provenance=prov,
            note=ground or "predict.py O10 standalone rung"))

    # --- cosmology beyond predict.py ---
    rows.append(SurfaceRow(
        key="omega_b_c", section="cosmology", form="25α + 75α²",
        order=ORDER_PINS["omega_b_c"], predicted=omega_b_c(),
        witness=WITNESSES["omega_b_c"],
        provenance=ADJUDICATIONS["omega_b_c"][0],
        note=ADJUDICATIONS["omega_b_c"][1]))
    rows.append(SurfaceRow(
        key="rho_lambda", section="cosmology",
        form="(m_e/M_Pl)⁴·α¹⁶·h_Cox",
        order=ORDER_PINS["rho_lambda"], predicted=lambda_cc(),
        witness=WITNESSES["rho_lambda"],
        provenance=ADJUDICATIONS["rho_lambda"][0],
        note=ADJUDICATIONS["rho_lambda"][1],
        discharge="dated amendment pinning either a forcing derivation of "
                  "h_Cox = 6 and the α¹⁶ exponent, or the honest dated "
                  "history of their selection (verdict discharge for the "
                  "UNDERSPECIFIED sub-factors)"))

    # --- electroweak scale ratio ---
    rows.append(SurfaceRow(
        key="v_over_m_e", section="electroweak",
        form="(25/α²)·(1 + 25α/(4√3))",
        order=ORDER_PINS["v_over_m_e"],
        predicted=(25.0 / a**2) * (1.0 + 25.0 * a / (4.0 * math.sqrt(3.0))),
        witness=WITNESSES["v_over_m_e"], provenance=DERIVED,
        note="substrate_gf P7b closed form; 25 = q_cinq² and 4√3 = C_A²·√3 are "
             "identified factors — flagged for the Auditor sweep (R2)"))

    # --- neutrino mixing ---
    rows.append(SurfaceRow(
        key="sin2_theta_13", section="neutrino", form="3α (θ₁₃ = asin√(3α))",
        order=ORDER_PINS["sin2_theta_13"], predicted=3.0 * a,
        witness=WITNESSES["sin2_theta_13"], provenance=DERIVED,
        note="RANK_SO7·α reactor angle (Spin(8) triality)"))
    for key, form, val in (
            ("sin2_theta_12", "1/3 (tri-bimaximal LO)", 1.0 / 3.0),
            ("sin2_theta_23", "1/2 (maximal LO)", 0.5)):
        rows.append(SurfaceRow(
            key=key, section="neutrino", form=form,
            order=ORDER_PINS[key], predicted=val,
            witness=WITNESSES[key], provenance=DERIVED,
            note="LEADING-ORDER value; the NLO charged-lepton rotation is an "
                 "OPEN derivation, deliberately not applied (a NuFIT-reproducing "
                 "magnitude would be a fit)",
            discharge="derive the U_ℓ NLO rotation from the substrate (Paper 20 "
                      "§7.6) BEFORE comparing its magnitude to NuFIT"))

    # --- mass ratios (all compendium particles) ---
    rows.extend(_mass_ratio_rows())

    # --- defect row: the n−p splitting the formula cannot produce ---
    p_n, p_p = particle("n"), particle("p")
    split_pred = (p_n.mass_pred - p_p.mass_pred) / ME_MEV      # identically 0
    m_e = MASS_MEV["e-"].value
    split_wit = Witness(
        (MASS_MEV["n"].value - MASS_MEV["p"].value) / m_e,
        math.hypot(MASS_MEV["n"].sigma, MASS_MEV["p"].sigma) / m_e,
        "CODATA-2018 (m_n − m_p)/m_e = 2.53102")
    rows.append(SurfaceRow(
        key="n_minus_p_over_m_e", section="defect",
        form="paper6(n) − paper6(p) — framing f does not enter the mass formula",
        order=ORDER_PINS["paper6_mass_ratio"], predicted=split_pred,
        witness=split_wit, provenance=DERIVED,
        note="the formula genuinely predicts ZERO splitting (n and p share "
             "(1,3,5,5), differing only in f) — an honest structural gap, "
             "reported rather than patched",
        discharge="derive a framing-dependent mass term (e.g. the EM self-energy "
                  "of the framed carrier) from the substrate, not fitted to "
                  "1.293 MeV"))

    return rows


# ===========================================================================
# S-FORWARD — the append-only post-freeze witness channel.
#
# Rules (enforced by test_surface.py + git history):
#   * entries are APPENDED, never edited or removed;
#   * an entry is (date, key, value, sigma, source) for a measurement PUBLISHED
#     after the frozen witness set;
#   * predictions are never re-derived against an update — the frozen value is
#     re-scored, and the drift direction (toward/away from the prediction) is
#     the readout.  A row that needs a formula change to survive an update is
#     DEAD-AS-EXACT permanently (no rescue orders).
# ===========================================================================

@dataclass(frozen=True)
class WitnessUpdate:
    date: str           # ISO date of the published measurement/adjustment
    key: str            # surface row key
    value: float
    sigma: float
    source: str


WITNESS_UPDATES: tuple[WitnessUpdate, ...] = (
    WitnessUpdate("2024-08-30", "inv_alpha", 137.035999177, 2.1e-8,
                  "CODATA-2022 adjustment α⁻¹ = 137.035999177(21)"),
)


def sforward_readout(rows: list[SurfaceRow] | None = None,
                     n_sigma: float = 2.0) -> list[dict]:
    """Re-score the FROZEN predictions against every witness update: new z,
    verdict at the update's precision, and the drift direction — did the
    measurement move toward the prediction (the only thing that can ever count
    as confirmation on this surface) or away from it?"""
    rows = surface_rows() if rows is None else rows
    by_key = {r.key: r for r in rows}
    out = []
    for u in WITNESS_UPDATES:
        r = by_key.get(u.key)
        if r is None:
            out.append({"date": u.date, "key": u.key, "verdict": "NO-SUCH-ROW",
                        "source": u.source})
            continue
        z_new = abs(r.predicted - u.value) / u.sigma
        d_old = abs(r.predicted - r.witness.value)
        d_new = abs(r.predicted - u.value)
        out.append({
            "date": u.date, "key": u.key, "predicted": r.predicted,
            "old": r.witness.value, "new": u.value, "sigma": u.sigma,
            "z_new": z_new,
            "verdict": "EXACT-COMPATIBLE" if z_new <= n_sigma else "DEAD-AS-EXACT",
            "drift": "toward" if d_new < d_old else
                     ("away" if d_new > d_old else "unchanged"),
            "source": u.source,
        })
    return out


# ===========================================================================
# The surface as an O10 DAG
# ===========================================================================

def build_surface_dag() -> DerivationDAG:
    """The full pre-registration surface as one O10 DerivationDAG: α as the
    root, one sym→output→witness chain per row (witnesses σ-graded so
    ``snow_readout`` scores every row), disputes on the contested symbolic
    nodes, and the anti-edges that keep each headline match from retro-
    justifying its own premise."""
    from ..isa.constants import ALPHA_SUBSTRATE

    g = DerivationDAG()
    g.add("form:25π√3+1", Stage.STRUCTURAL, note="α⁻¹ closed form (no measured input)")
    g.add("α", Stage.STRUCTURAL, ALPHA_SUBSTRATE, "isa.ALPHA_SUBSTRATE")
    g.link("form:25π√3+1", "α")

    for r in surface_rows():
        sym, out, wit = f"sym:{r.key}", r.key, f"wit:{r.key}"
        g.add(sym, Stage.SYMBOLIC, note=f"{r.form}  [order: {r.order}]",
              disputed=r.disputed,
              status=STATUS_DEFERRED_BRIDGE if r.discharge else "",
              discharge=r.discharge)
        g.add(out, Stage.OUTPUT, r.predicted, r.note, provenance=r.provenance)
        g.add(wit, Stage.WITNESS, r.witness.value, r.witness.source,
              sigma=r.witness.sigma)
        if r.section == "mass-ratio":
            walk = f"walk:{r.key}"
            g.add(walk, Stage.STRUCTURAL, note=r.form,
                  provenance=POST_SELECTED)
            g.link(walk, sym)
        g.link("α", sym)
        g.link(sym, out)
        g.link(out, wit)

    # --- anti-edges: no headline match may retro-justify its premise ---
    g.forbid("wit:inv_alpha", "form:25π√3+1",
             reason="the CODATA match ⇏ the closed form that hits it",
             discharge="derive 25, √3 and +1 from the trefoil Aharonov-Bohm walk "
                       "at a NEIGHBOURING winding the form cannot be retuned for, "
                       "or from a second observable sharing the same factors")
    g.forbid("wit:eta_B", "sym:eta_B",
             reason="Planck η_B match ⇏ the 3α⁴/14 form (17-formula-sweep survivor)",
             discharge="a second baryon-sector observable hit by the SAME 3 and 14, "
                       "or the Jones/Murasugi chirality forcing chain")
    g.forbid("wit:omega_b_c", "sym:omega_b_c",
             reason="Planck Ω_b/Ω_c match ⇏ the 75α² NLO term (added inside the "
                    "error bar)",
             discharge="derive 75 = 3·25 from the K₇⊗K₈ partition combinatorics "
                       "with an independent observable it also moves")
    g.forbid("wit:m_p_over_m_e", "walk:m_p_over_m_e",
             reason="the PDG m_p match ⇏ the (1,3,5,5) assignment (corrected "
                    "2026-04-30 with the mass in hand)",
             discharge="a forcing chain that pins (p,q,m,n_q) from charge/isospin/"
                       "baryon predicates WITHOUT mass input — then the mass is a "
                       "genuine prediction")
    # Wiring gap closed per the Auditor verdict (rho_lambda adjudication): the
    # row-specific factors had no anti-edge.  h_Cox = 6 is the unique menu
    # integer landing inside the bar and every integer 1–8 carries an available
    # Lie-theory name, so the Planck match must not retro-justify either factor.
    g.forbid("wit:rho_lambda", "sym:rho_lambda",
             reason="Planck ρ_Λ match ⇏ h_Cox = 6 or the α¹⁶ exponent "
                    "(UNDERSPECIFIED per the verdict; mode 5 hits 100% of "
                    "random targets at this row's 3% bar)",
             discharge="pin a forcing derivation of h_Cox and 16 = 2·dim(𝕆), or "
                       "the honest dated history of their selection (dated "
                       "amendment; verdict discharge)")
    return g


# ===========================================================================
# CLI — this output is the prereg's §3 table generator
# ===========================================================================

def _row_lines(rows: list[SurfaceRow]) -> list[str]:
    lines = ["", f"Surface rows ({len(rows)}; witnesses = CODATA-2018/PDG-2022/"
             "Planck-2018, quarantined):"]
    for r in rows:
        lines.append(f"  {r.key:24s} [{r.section:11s}] {r.form}  "
                     f"(order {r.order}; provenance {r.provenance}"
                     f"{'; DISPUTED' if r.disputed else ''})")
    return lines


def _exclusion_lines() -> list[str]:
    lines = ["", "Exclusions (every claimed quantity NOT on the surface, dated — "
             "silent omission is a look-elsewhere sin):"]
    for k, v in EXCLUSIONS.items():
        lines.append(f"  {k:16s} {v}")
    return lines


def _sforward_lines(n_sigma: float = 2.0) -> list[str]:
    out = sforward_readout(n_sigma=n_sigma)
    if not out:
        return []
    lines = ["", "S-FORWARD readout (append-only witness updates vs FROZEN "
             "predictions — the only confirmation channel):"]
    for r in out:
        if r["verdict"] == "NO-SUCH-ROW":
            lines.append(f"  [{r['date']}]  {r['key']}: no such surface row "
                         f"({r['source']})")
            continue
        z = f"{r['z_new']:.3g}σ" if r["z_new"] < 1e6 else f"{r['z_new']:.2e}σ"
        lines.append(f"  [{r['date']}]  {r['key']:16s} drift {r['drift'].upper():7s} "
                     f"{r['old']:.9g} → {r['new']:.9g}  ({r['verdict']} at {z}; "
                     f"{r['source']})")
    return lines


def main(argv: list[str] | None = None) -> int:
    rows = surface_rows()
    g = build_surface_dag()
    lines = [f"Pre-registration kill surface — {len(rows)} rows, "
             f"{len(g.nodes)} DAG nodes, {len(g.edges)} edges", ""]
    lines.append("Acceptance checklist (O10 structural invariants):")
    for check, ok in g.acceptance_checklist().items():
        lines.append(f"  [{'PASS' if ok else 'FAIL'}]  {check}")
    lines += _row_lines(rows)
    lines += _provenance_lines(g)
    lines += _dispute_lines(g)
    lines += _snow_lines(g)
    lines += _sforward_lines()
    lines += _collapse_lines(g)
    lines += _exclusion_lines()
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
