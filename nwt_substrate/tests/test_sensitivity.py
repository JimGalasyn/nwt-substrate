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
