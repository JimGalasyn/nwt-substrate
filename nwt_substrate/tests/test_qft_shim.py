"""Tests for the nwt.qft Lagrangian shim."""

import pytest

import nwt_substrate.qft as qft
from nwt_substrate.qft.fields import (
    ScalarField, DiracSpinor, GaugeField, FaddeevPopovGhost,
)
from nwt_substrate.qft.lagrangian import Lagrangian, InteractionTerm


# ---------------------------------------------------------------------------
# Field types
# ---------------------------------------------------------------------------

def test_field_types_carry_substrate_primitive_pointers():
    """Every Field has a non-empty .substrate_primitive description."""
    for f in [qft.electron, qft.muon, qft.up_quark, qft.photon, qft.gluon,
              qft.ghost_qcd]:
        assert f.substrate_primitive
        assert len(str(f)) > 5


def test_dirac_spinor_records_8d_substrate_primitive():
    """Substrate Dirac spinor is 8-dim Cl(1,3) (per memory and dirac.py)."""
    assert "8-dimensional" in qft.electron.substrate_primitive
    assert "Cl(1,3)" in qft.electron.substrate_primitive


def test_quarks_carry_color():
    assert qft.up_quark.n_color == 3
    assert qft.down_quark.n_color == 3
    assert qft.electron.n_color == 1


# ---------------------------------------------------------------------------
# QED Lagrangian
# ---------------------------------------------------------------------------

def test_qed_lagrangian_has_textbook_text():
    assert "L_QED" in qft.qed.text
    assert "γ^μ" in qft.qed.text or "γ" in qft.qed.text


def test_qed_has_9_charged_fermions():
    """3 leptons + 6 quarks = 9 charged Dirac fermions in the loop."""
    fermions = [f for f in qft.qed.fields if isinstance(f, DiracSpinor)]
    assert len(fermions) == 9


def test_qed_has_one_photon():
    photons = [f for f in qft.qed.fields if isinstance(f, GaugeField)]
    assert len(photons) == 1
    assert photons[0].gauge_group == "U(1)_em"


def test_qed_has_9_qed_vertices():
    """One QED vertex per charged fermion species."""
    qed_vertices = [v for v in qft.qed.interactions if "QED vertex" in v.name]
    assert len(qed_vertices) == 9


def test_qed_beta_0_matches_pdg():
    """β_0(L_QED) = (2/3) * b_QED with b_QED = 8 → 16/3 ≈ 5.333."""
    assert qft.qed.beta_0() == pytest.approx(16.0 / 3.0, abs=1e-12)


# ---------------------------------------------------------------------------
# QCD Lagrangian
# ---------------------------------------------------------------------------

def test_qcd_lagrangian_has_textbook_text():
    assert "L_QCD" in qft.qcd.text


def test_qcd_has_6_quarks_one_gluon_one_ghost():
    quarks = [f for f in qft.qcd.fields if isinstance(f, DiracSpinor)]
    gluons = [f for f in qft.qcd.fields if isinstance(f, GaugeField)]
    ghosts = [f for f in qft.qcd.fields if isinstance(f, FaddeevPopovGhost)]
    assert len(quarks) == 6
    assert len(gluons) == 1 and gluons[0].gauge_group == "SU(3)_color"
    assert len(ghosts) == 1


def test_qcd_has_3gluon_4gluon_and_ghost_vertices():
    """L_QCD's interaction list includes the standard non-Abelian gauge vertices."""
    vnames = [v.name for v in qft.qcd.interactions]
    assert any("3-gluon" in n for n in vnames)
    assert any("4-gluon" in n for n in vnames)
    assert any("ghost-gluon" in n for n in vnames)


def test_qcd_beta_0_matches_running_couplings_module():
    """β_0(L_QCD, n_f=5) should equal (11 - 2·5/3) = 23/3, matching the
    structurally-derived value in vacuum_polarization."""
    assert qft.qcd.beta_0(n_f_dirac=5) == pytest.approx(23.0 / 3.0, abs=1e-12)
    assert qft.qcd.beta_0(n_f_dirac=6) == pytest.approx(7.0, abs=1e-12)
    assert qft.qcd.beta_0(n_f_dirac=0) == pytest.approx(11.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Generic Yang-Mills
# ---------------------------------------------------------------------------

def test_yang_mills_su_n_has_n_squared_minus_one_generators():
    for N in [2, 3, 5]:
        L = qft.yang_mills(N)
        gauge = [f for f in L.fields if isinstance(f, GaugeField)]
        assert len(gauge) == 1
        assert gauge[0].n_generators == N * N - 1


def test_yang_mills_su_n_has_3_and_4_gauge_vertices():
    L = qft.yang_mills(3)
    vnames = [v.name for v in L.interactions]
    assert any("3-gauge" in n for n in vnames)
    assert any("4-gauge" in n for n in vnames)


def test_yang_mills_beta_0_recognises_generic_su_n():
    """beta_0() must recognise the bare "SU(N)" tag from yang_mills(N), not
    just the named "SU(3)_color" QCD theory, using β_0 = 11·N/3 − 2·n_f/3.

    Regression: yang_mills(N).beta_0(...) previously fell through to None
    because only the literal "SU(3)_color" tag was matched."""
    from nwt_substrate.amplitudes import vacuum_polarization as vp
    # SU(3) + 6 quarks reproduces the QCD value exactly.
    assert qft.yang_mills(3).beta_0(n_f_dirac=6) == pytest.approx(7.0, abs=1e-12)
    assert (qft.yang_mills(3).beta_0(n_f_dirac=6)
            == qft.qcd.beta_0(n_f_dirac=6))
    # No n_f → pure glue (n_f = 0) for a bare Yang-Mills theory with no matter.
    assert qft.yang_mills(3).beta_0() == pytest.approx(11.0, abs=1e-12)
    assert qft.yang_mills(2).beta_0() == pytest.approx(22.0 / 3.0, abs=1e-12)
    # General N matches the structural vacuum_polarization formula.
    for N, n_f in [(2, 0), (3, 6), (5, 2)]:
        assert qft.yang_mills(N).beta_0(n_f_dirac=n_f) == pytest.approx(
            vp.qcd_beta_0(n_f_dirac=n_f, N_c=N), abs=1e-12)
    # Non-gauge theories still return None.
    assert qft.klein_gordon().beta_0() is None


# ---------------------------------------------------------------------------
# Lagrangian composition
# ---------------------------------------------------------------------------

def test_lagrangian_sum_concatenates_fields_and_vertices():
    L_sum = qft.qed + qft.qcd
    assert len(L_sum.fields) == len(qft.qed.fields) + len(qft.qcd.fields)
    assert (len(L_sum.interactions) ==
            len(qft.qed.interactions) + len(qft.qcd.interactions))
    assert "U(1)_em" in L_sum.gauge_symmetries
    assert "SU(3)_color" in L_sum.gauge_symmetries


def test_standard_model_matter_combines_qed_and_qcd():
    L_sm = qft.standard_model_matter()
    # Both QED and QCD content present
    assert any(isinstance(f, GaugeField) and f.gauge_group == "U(1)_em"
               for f in L_sm.fields)
    assert any(isinstance(f, GaugeField) and f.gauge_group == "SU(3)_color"
               for f in L_sm.fields)


# ---------------------------------------------------------------------------
# Substrate-view introspection
# ---------------------------------------------------------------------------

def test_substrate_view_lists_all_fields():
    view = qft.qed.substrate_view()
    assert "Cl(1,3)" in view
    assert "U(1)_em" in view
    for f in qft.qed.fields:
        assert f.name in view


def test_feynman_rules_point_at_substrate_primitives():
    rules = qft.qed.feynman_rules()
    assert "propagators" in rules
    assert "vertices" in rules
    # QED vertex should reference the substrate qed_vertex helper
    qed_vtx_refs = [v for k, v in rules["vertices"].items() if "QED vertex" in k]
    assert len(qed_vtx_refs) == 9
    assert all("qed_vertex" in str(v) for v in qed_vtx_refs)


def test_feynman_rules_qcd_references_three_and_four_gluon_helpers():
    rules = qft.qcd.feynman_rules()
    vtx_strs = " ".join(str(v) for v in rules["vertices"].values())
    assert "three_gluon_vertex_lorentz" in vtx_strs
    assert "four_gluon_vertex_lorentz" in vtx_strs


# ---------------------------------------------------------------------------
# Klein-Gordon free scalar
# ---------------------------------------------------------------------------

def test_klein_gordon_real_vs_complex():
    L_real = qft.klein_gordon("φ", mass=1.0, charge=0.0)
    L_cmplx = qft.klein_gordon("φ", mass=1.0, charge=1.0)
    s_real = [f for f in L_real.fields if isinstance(f, ScalarField)][0]
    s_cmplx = [f for f in L_cmplx.fields if isinstance(f, ScalarField)][0]
    assert s_real.real is True
    assert s_cmplx.real is False


# ---------------------------------------------------------------------------
# Multi-view introspection
# ---------------------------------------------------------------------------

def test_multiview_electron_shows_all_three_lenses():
    view = qft.multiview("electron")
    assert "QFT view" in view
    assert "Substrate primitive" in view
    # Should mention either 8-dim or Cl(1,3)
    assert "8-dim" in view or "Cl(1,3)" in view


def test_multiview_unknown_particle():
    view = qft.multiview("xyzzy")
    assert "Unknown" in view


# ---------------------------------------------------------------------------
# Cross-shim consistency: QFT β-functions match vacuum_polarization
# ---------------------------------------------------------------------------

def test_qft_beta_0_consistent_with_vacuum_polarization():
    """L_QED.beta_0() and L_QCD.beta_0() are not separate calculations from
    vacuum_polarization.qed_beta_0_total / qcd_beta_0; they're the same
    value labelled in QFT-vocabulary."""
    from nwt_substrate.amplitudes import vacuum_polarization as vp
    species = vp.standard_qed_species()
    assert qft.qed.beta_0() == pytest.approx(vp.qed_beta_0_total(species),
                                              abs=1e-12)
    for n_f in [3, 4, 5, 6]:
        assert qft.qcd.beta_0(n_f_dirac=n_f) == pytest.approx(
            vp.qcd_beta_0(n_f_dirac=n_f, N_c=3), abs=1e-12,
        )
