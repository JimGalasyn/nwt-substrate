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


def test_clean_survivors():
    """Rows compatible at 2σ whose provenance is NOT suspect.  After the
    2026-07-12 adjudication (verdict 4032d04) the fitted/post-selected rows
    are excluded by tag rather than by open dispute; what remains is the
    program's entire surviving compatible-and-clean set.  If a library change
    grows or shrinks this set, that is a headline event, not noise."""
    from nwt_substrate.benchmarks.o10 import SUSPECT_PROVENANCE
    g = build_surface_dag()
    suspect = {r.key for r in surface_rows()
               if r.provenance in SUSPECT_PROVENANCE}
    survivors = {r["output"] for r in g.snow_readout()
                 if r["verdict"] == "EXACT-COMPATIBLE"} - suspect
    assert survivors == {"cabibbo_lambda", "sin2_theta_13"}


# ---------------------------------------------------------------------------
# Adjudications (T1 closed: the Auditor verdict of 2026-07-12, nwt-audit
# 4032d04 — DERIVED refuted on all four rows; write-backs below are the
# pre-committed consequences, executed citing that verdict)
# ---------------------------------------------------------------------------

def test_adjudications_written_back():
    """The four adjudicated tags are in force, every ground cites the
    verdict, no dispute remains open, and the suspended-cit machinery is
    quiescent (nothing left to suspend)."""
    from nwt_substrate.benchmarks.surface import ADJUDICATIONS
    expected = {"eta_B": "post_selected", "omega_b_c": "fitted",
                "m_e_over_M_Pl": "motivated", "rho_lambda": "motivated"}
    rows = {r.key: r for r in surface_rows()}
    for key, tag in expected.items():
        assert rows[key].provenance == tag, key
        assert rows[key].disputed == "", key            # cleared by adjudication
        assert "2026-07-12-constants-provenance-disputes" in rows[key].note, key
        assert ADJUDICATIONS[key][0] == tag
    g = build_surface_dag()
    assert g.dispute_audit()["clean"]


def test_adjudicated_rows_flagged_circular_on_cit():
    """With suspect tags in force, any of the four rows that passes cit must
    surface in provenance_defects as CIRCULAR — agreement is not evidence for
    a fitted/post-selected/motivated value.  (This is the verdict's 'S-NOW
    compatibility may never be cited as support' clause, enforced by lint.)"""
    g = build_surface_dag()
    circular = {r["output"] for r in g.provenance_audit()["circular_passes"]}
    # eta_B, omega_b_c and rho_lambda pass cit at 1% — they must be flagged.
    assert {"eta_B", "omega_b_c", "rho_lambda"} <= circular


def test_alpha_qed_leakage_still_pinned():
    """The measured-α leakage (exhibit (a), CONFIRMED by the Auditor) remains
    documented in the adjudication section and the leakage itself remains in
    coupling.py.  If the default is ever flipped to the substrate α, that is
    a post-verdict correction: record it by dated note, do not delete this
    test."""
    import inspect
    import pathlib

    from nwt_substrate.benchmarks import surface as surface_mod
    src = pathlib.Path(surface_mod.__file__).read_text()
    assert "ALPHA_QED" in src and "neighbouring_value" in src
    from nwt_substrate.gravity import coupling
    from nwt_substrate.gravity.constants import ALPHA_QED
    from nwt_substrate.isa.constants import ALPHA_SUBSTRATE
    default = inspect.signature(coupling.m_e_over_M_Pl_NNLO).parameters["alpha"].default
    assert default == ALPHA_QED and default != ALPHA_SUBSTRATE


def test_neighbouring_value_lookelsewhere_measured():
    """Gauntlet mode 5, pinned (exhibit (b), CONFIRMED by the Auditor with
    independent seeds/menus): the substrate fitting pattern fits RANDOM
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
    """The frozen order per closed form, verbatim.  Changing an order must
    change this test — a dated, deliberate amendment, never a silent retrofit.

    AMENDMENT 2026-07-12 (authorized): m_e_over_M_Pl and rho_lambda
    NNLO → NLO per the Auditor verdict CL-2
    (2026-07-12-constants-provenance-disputes, nwt-audit 4032d04): the NNLO
    α² coefficient is documented target-selection (Paper 17 computed it from
    CODATA before choosing the nearest structural integer); NLO is the
    externally audited L4(a) form.  NNLO is retired from claim status and no
    post-freeze order change can revive it."""
    assert ORDER_PINS == {
        "inv_alpha": "exact",
        "sin2_theta_W": "LO",
        "cabibbo_lambda": "exact",
        "eta_B": "exact",
        "m_e_over_M_Pl": "NLO",
        "omega_b_c": "NLO",
        "rho_lambda": "NLO",
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
