"""Tests for nwt_substrate.sensitivity — the ISA structural-integer sweep."""
import os

import pytest

from nwt_substrate import sensitivity as sens


# ---- pure helpers (no subprocess) ----

def test_leaf_integers_found():
    leaves = dict(sens.leaf_integers())
    # a representative sample of the K_7 / Spin(7) leaf integers
    assert leaves["N_VERTICES_K7"] == 7
    assert leaves["DIM_S_SPIN7"] == 8
    assert leaves["N_EDGES_K7"] == 21
    assert leaves["RANK_SO7"] == 3
    assert len(leaves) >= 15  # there are ~21


def test_patched_source_replaces_one_leaf():
    orig = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    patched = sens._patched_source(orig, "DIM_S_SPIN7", 9)
    assert "DIM_S_SPIN7: int = 9" in patched
    assert "DIM_S_SPIN7: int = 8" not in patched
    # exactly one line changed
    assert sum(a != b for a, b in zip(orig.splitlines(), patched.splitlines())) == 1


def test_patched_source_rejects_unknown_name():
    orig = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    with pytest.raises(ValueError):
        sens._patched_source(orig, "NOT_A_CONSTANT", 1)


def test_patched_scalar_source_single_and_multiline():
    """The scalar-knob patch multiplies a derived constant's RHS by (1+eps) and
    stays valid Python — for a single-line def (α) and a parenthesised def that
    spans several lines (κ), which the balanced-paren scan collapses safely."""
    import ast
    orig = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    for name in ("ALPHA_SUBSTRATE", "KAPPA_MACKEN"):
        patched = sens._patched_scalar_source(orig, name, 1e-3)
        ast.parse(patched)                            # valid Python
        assert f"{name}: float = (" in patched        # RHS wrapped in parens
        assert "(1.0 + (0.001))" in patched           # times (1 + eps)
        assert patched != orig
    with pytest.raises(ValueError):
        sens._patched_scalar_source(orig, "NOT_A_SCALAR", 1e-3)


@pytest.mark.skipif(not os.access(sens.CONSTANTS_PATH, os.W_OK),
                    reason="sensitivity sweep needs a writable (editable) install")
def test_scalar_knob_alpha_couples_qed_observables():
    """The α scalar knob closes the coverage gap for α-anchored benchmarks that
    no integer moves: perturbing α moves pure-QED observables (electron anomaly)
    and the now-isa-anchored α-derivation benchmark, and the source is restored."""
    before = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    try:
        rep = sens.integer_sweep(integers=[], scalars={"ALPHA_SUBSTRATE": 1e-3})
    finally:
        after = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    assert after == before
    movers = rep.movers("ALPHA_SUBSTRATE")
    assert "benchmark_electron_anomaly" in movers
    assert "benchmark_alpha_derivation" in movers     # sourced from isa now, not 25π√3+1 literal
    assert "ALPHA_SUBSTRATE" not in rep.inert_integers


# ---- integration: one-integer sweep, then verify the source is restored ----

@pytest.mark.skipif(not os.access(sens.CONSTANTS_PATH, os.W_OK),
                    reason="sensitivity sweep needs a writable (editable) install")
def test_sweep_couples_and_restores():
    before = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    try:
        report = sens.integer_sweep(integers=["RANK_SO7", "DEGREE_K8"])
    finally:
        after = open(sens.CONSTANTS_PATH, encoding="utf-8").read()
    # source restored byte-identically no matter what
    assert after == before

    # RANK_SO7 (= N_c = N_generations) is load-bearing: it moves predictions.
    assert len(report.movers("RANK_SO7")) > 0
    # DEGREE_K8 is pure bookkeeping in the current code: inert.
    assert report.movers("DEGREE_K8") == []
    assert "DEGREE_K8" in report.inert_integers
    assert len(report.benchmarks) == 38


# ---- structural-criticality layer (Marcel Wende, d12rg) ----

def test_structural_load_and_comovement():
    """Load ranking + correlated-cluster diagnostic on a synthetic report."""
    r = sens.SensitivityReport(
        benchmarks=["a", "b", "c"],
        baseline={"a": "1", "b": "2", "c": "3"},
        per_integer={
            "X": {1: {"status": "ok", "moved": ["a", "b"], "value": 7}},
            "Y": {1: {"status": "ok", "moved": ["a", "b", "c"], "value": 3}},
            "Z": {1: {"status": "ok", "moved": [], "value": 5}},
        },
    )
    load = dict(r.structural_load())
    assert load == {"X": 2, "Y": 3, "Z": 0}
    # ranked by load descending: Y (3) then X (2) then Z (0)
    assert [i for i, _ in r.structural_load()] == ["Y", "X", "Z"]
    # a,b co-move under BOTH X and Y -> 2; a,c and b,c only under Y -> 1
    comov = dict(r.comovement(min_shared=1))
    assert comov[("a", "b")] == 2
    assert comov[("a", "c")] == 1 and comov[("b", "c")] == 1
    # min_shared=2 keeps only the (a,b) correlated pair
    assert r.comovement(min_shared=2) == [(("a", "b"), 2)]
    assert "Structural-load ranking" in r.criticality_summary()


@pytest.mark.skipif(not os.access(sens.CONSTANTS_PATH, os.W_OK),
                    reason="sensitivity sweep needs a writable (editable) install")
def test_gf_refactor_couples_fermi_to_h_v():
    """Regression: after wiring substrate_gf to H_V_SO7, perturbing H_V_SO7
    must move the Fermi-constant and Higgs-VEV predictions (they were inert
    when 25/625 were bare literals)."""
    report = sens.integer_sweep(integers=["H_V_SO7"])
    movers = report.movers("H_V_SO7")
    assert "benchmark_fermi_constant" in movers
    assert "benchmark_higgs_vev" in movers
