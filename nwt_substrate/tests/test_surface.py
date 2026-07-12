"""Unit tests for nwt_substrate.benchmarks.surface — the pre-registration
kill surface.

Beyond correctness, several tests here are FREEZE GUARDS: they assert the
frozen order pins, the append-only witness-update discipline, and the
explicit-provenance rule verbatim, so that changing any of them is a
deliberate, dated amendment (the test diff is the receipt) rather than a
silent retrofit.
"""

import dataclasses

import pytest

from nwt_substrate.benchmarks import surface
from nwt_substrate.benchmarks.o10 import (
    DERIVED, POST_SELECTED, Stage, _VALID_PROVENANCE,
)
from nwt_substrate.benchmarks.surface import (
    EXCLUSIONS, ORDER_PINS, WITNESS_UPDATES, WITNESSES,
    build_surface_dag, sforward_readout, surface_rows,
)


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

def test_dag_invariants_pass():
    """The surface DAG satisfies every O10 structural invariant."""
    g = build_surface_dag()
    assert all(g.acceptance_checklist().values())


def test_every_row_scored_or_defect():
    """Every surface row reaches an S-NOW verdict — no silent unscored rows."""
    g = build_surface_dag()
    verdicts = {r["output"]: r["verdict"] for r in g.snow_readout()}
    for row in surface_rows():
        assert row.key in verdicts, f"row {row.key} not S-NOW-scored"
        assert verdicts[row.key] != "UNSCORED", f"row {row.key} has no σ"


def test_all_compendium_mass_ratios_present():
    """Every compendium particle except the e⁻ unit appears as a ratio row."""
    from nwt_substrate.particles import list_particles
    keys = {r.key for r in surface_rows()}
    for name in list_particles():
        if name == "e-":
            continue
        assert f"m_{name}_over_m_e" in keys


def test_explicit_provenance_everywhere():
    """Surface rows never rely on the default-DERIVED inference — the self-tag
    problem this module exists to surface starts with implicit tags."""
    for row in surface_rows():
        assert row.provenance in _VALID_PROVENANCE and row.provenance


# ---------------------------------------------------------------------------
# Verdicts (the honest headline results — these assert the CURRENT truth so a
# library change that flips one is loudly visible)
# ---------------------------------------------------------------------------

def _verdict(key: str) -> str:
    g = build_surface_dag()
    return {r["output"]: r["verdict"] for r in g.snow_readout()}[key]


def test_inv_alpha_dead_as_exact():
    """α itself is the surface's hardest kill: 7.6 ppm off a 1.5e-10 witness."""
    assert _verdict("inv_alpha") == "DEAD-AS-EXACT"


def test_m_p_over_m_e_dead_as_exact():
    """m_p/m_e is measured to 6e-11 relative; the 0.11% walk-formula residual
    is ~3.6e6 σ."""
    assert _verdict("m_p_over_m_e") == "DEAD-AS-EXACT"


def test_n_minus_p_split_predicted_zero():
    """n and p share (1,3,5,5); the formula predicts zero splitting — the
    defect row must report exactly that, not a patched value."""
    row = {r.key: r for r in surface_rows()}["n_minus_p_over_m_e"]
    assert row.predicted == 0.0
    assert row.discharge          # the gap carries its killable close-out test
    assert _verdict("n_minus_p_over_m_e") == "DEAD-AS-EXACT"


def test_undisputed_survivors():
    """The undisputed rows currently compatible at 2σ.  If a library change
    grows or shrinks this set, that is a headline event, not noise."""
    g = build_surface_dag()
    disputed_rows = {r.key for r in surface_rows() if r.disputed}
    survivors = {r["output"] for r in g.snow_readout()
                 if r["verdict"] == "EXACT-COMPATIBLE"} - disputed_rows
    assert survivors == {"cabibbo_lambda", "sin2_theta_13", "m_Sigma*_over_m_e"}


# ---------------------------------------------------------------------------
# Disputes (T1: self-tags vs pinned external audits)
# ---------------------------------------------------------------------------

def test_disputes_present_and_suspending():
    """The three pinned-audit disputes exist and suspend their rows' cit
    passes as corroboration."""
    g = build_surface_dag()
    audit = g.dispute_audit()
    disputed = {r["node"] for r in audit["disputes"]}
    assert {"sym:eta_B", "sym:omega_b_c", "sym:m_e_over_M_Pl",
            "sym:rho_lambda"} <= disputed
    assert "eta_B" in audit["suspended_outputs"]
    assert "omega_b_c" in audit["suspended_outputs"]


def test_alpha_qed_leakage_pinned_in_dispute():
    """The measured-α leakage evidence stays in the m_e_over_M_Pl dispute
    record until the Auditor adjudicates: gravity/coupling.py defaults its
    chain to alpha=ALPHA_QED (the measured CODATA value), which is where the
    'G at 11 ppm' headline came from — substrate-pure the chain misses G by
    ~150 ppm ≈ 7σ.  Dropping this evidence is a dispute-record edit, which is
    forbidden short of adjudication."""
    from nwt_substrate.benchmarks.surface import DISPUTES
    record = DISPUTES["m_e_over_M_Pl"]
    assert "ALPHA_QED" in record
    assert "reverse-smuggle" in record
    # and the leakage itself is real: the coupling chain's default arg is the
    # measured α, not the substrate α.
    import inspect
    from nwt_substrate.gravity import coupling
    from nwt_substrate.gravity.constants import ALPHA_QED
    from nwt_substrate.isa.constants import ALPHA_SUBSTRATE
    sig = inspect.signature(coupling.m_e_over_M_Pl_NNLO)
    default = sig.parameters["alpha"].default
    assert default == ALPHA_QED and default != ALPHA_SUBSTRATE, (
        "coupling.py's default α changed — if it now uses the substrate α, "
        "this is an adjudication-relevant fix: update the dispute record "
        "(dated) rather than deleting this test")


def test_neighbouring_value_lookelsewhere_measured():
    """Gauntlet mode 5, pinned: the substrate fitting pattern fits RANDOM
    targets inside G's error bar most of the time (full menu), so the real
    m_e/M_Pl ppm-hit carries no information beyond the menu.  The thresholds
    here are loose floors — if a code change makes the procedure LESS able to
    fit noise, these bounds still hold; if someone weakens the test to rescue
    the G relation, the diff is the receipt."""
    from nwt_substrate.benchmarks.neighbouring_value import (
        FULL_MENU, MINIMAL_MENU, sweep,
    )
    full = sweep(FULL_MENU, n_targets=500, seed=7)
    assert full["frac_within_G_bar"] > 0.5      # measured ~0.83
    minimal = sweep(MINIMAL_MENU, n_targets=500, seed=7)
    assert minimal["frac_within_G_bar"] > 0.10  # measured ~0.24
    # and the evidence is pinned in the dispute record
    from nwt_substrate.benchmarks.surface import DISPUTES
    assert "neighbouring_value" in DISPUTES["m_e_over_M_Pl"]


def test_mass_rows_post_selected():
    """Mass-ratio rows carry POST_SELECTED (the compendium's own 2026-04-30
    correction note is the in-repo evidence) — they may never quietly become
    DERIVED without a forcing chain discharging the anti-edge."""
    for row in surface_rows():
        if row.section == "mass-ratio":
            assert row.provenance == POST_SELECTED, row.key


def test_anti_edges_carry_discharges():
    """Every forbidden collapse names its killable discharge test — an
    anti-edge with no discharge is an un-cashable IOU."""
    g = build_surface_dag()
    assert g.antiedges
    for fc in g.antiedges:
        assert fc.discharge.strip(), f"{fc.src} ⇏ {fc.dst} has no discharge"
    assert not g.collapse_audit()["violations"]


# ---------------------------------------------------------------------------
# Freeze guards
# ---------------------------------------------------------------------------

def test_order_pins_frozen():
    """The frozen order per closed form, verbatim.  Changing an order (e.g.
    adding an NNLO term to survive a data update) must change this test —
    a dated, deliberate amendment, never a silent retrofit."""
    assert ORDER_PINS == {
        "inv_alpha": "exact",
        "sin2_theta_W": "LO",
        "cabibbo_lambda": "exact",
        "eta_B": "exact",
        "m_e_over_M_Pl": "NNLO",
        "omega_b_c": "NLO",
        "rho_lambda": "NNLO",
        "v_over_m_e": "NLO",
        "sin2_theta_13": "LO",
        "sin2_theta_12": "LO",
        "sin2_theta_23": "LO",
        "paper6_mass_ratio": "exact",
    }


def test_witness_updates_append_only_and_dated():
    """S-FORWARD entries are frozen dataclasses, ISO-dated, and sorted —
    appending is the only legal mutation (git history is the enforcement)."""
    dates = [u.date for u in WITNESS_UPDATES]
    assert dates == sorted(dates)
    for u in WITNESS_UPDATES:
        assert len(u.date) == 10 and u.date[4] == "-" and u.date[7] == "-"
        assert dataclasses.is_dataclass(u) and u.__dataclass_params__.frozen
        assert u.sigma > 0 and u.source


def test_sforward_scores_frozen_prediction():
    """The CODATA-2022 α update re-scores the FROZEN prediction and reports
    the drift; the 2022 adjustment moved AWAY from 25π√3 + 1."""
    out = {r["key"]: r for r in sforward_readout()}
    r = out["inv_alpha"]
    assert r["verdict"] == "DEAD-AS-EXACT"
    assert r["drift"] == "away"


def test_exclusions_are_dated():
    """Every exclusion carries its date — the honesty ledger the Auditor sweep
    checks against the library (silent omission = look-elsewhere sin)."""
    assert EXCLUSIONS
    for reason in EXCLUSIONS.values():
        assert reason[:5] == "2026-", reason


def test_witnesses_never_feed_predictions():
    """Quarantine: witness nodes are sinks in the DAG (predictions cannot read
    a measured value), per the O10 one-way rule."""
    g = build_surface_dag()
    assert g.witnesses_are_sinks()
    for n, nd in g.nodes.items():
        if nd.stage == Stage.WITNESS:
            assert not [d for s, d in g.edges if s == n]
