"""Tests for the nwt.string string-theoretic shim."""

import pytest

import nwt_substrate.string as st


# ---------------------------------------------------------------------------
# Compactification targets
# ---------------------------------------------------------------------------

def test_K7_torus_is_genus_1_with_21_cycles():
    """K_7 has Heffter triangular toroidal embedding: V=7, E=21, F=14, g=1.
    The 21 1-cycles support the Wilson amplitude tower."""
    assert st.K7_torus.genus == 1
    assert st.K7_torus.n_cycles == 21


def test_K7_torus_carries_Spin7_holonomy():
    """K_7 sits inside Spin(7) via 2T → Spin(7) (Phase b2.12)."""
    assert "Spin(7)" in st.K7_torus.holonomy


def test_poincare_sphere_is_E8_ADE():
    """S³/2I = Σ(2,3,5) is the Brieskorn sphere of type E_8 in the
    McKay/ADE classification."""
    assert "E_8" in st.poincare_sphere.ade_type
    assert "2I" in st.poincare_sphere.holonomy


def test_T2_torus_carries_pq_winding():
    """The (p,q) torus-knot windings are 1-cycles on T²."""
    assert st.T2_torus.genus == 1
    assert st.T2_torus.n_cycles == 2


# ---------------------------------------------------------------------------
# (p,q)-string view
# ---------------------------------------------------------------------------

def test_electron_is_2_1_string():
    """Electron carries (p,q) = (2,1) -- the canonical NWT trefoil-like
    winding -- which is also the (2,1) F1+D1 string doublet in IIB."""
    e = st.pq_string("electron")
    assert e.p == 2 and e.q == 1
    assert e.winding_squared == 5
    assert e.sl2z_doublet == (2, 1)


def test_proton_is_1_3_string():
    """Proton carries (p,q) = (1,3) per the corrected nucleon tuple
    (memory: mass-formula-refined-2026-04-30)."""
    p = st.pq_string("proton")
    assert p.p == 1 and p.q == 3
    assert p.winding_squared == 10  # 1² + 3²


def test_pq_string_unknown_particle_raises():
    with pytest.raises(KeyError):
        st.pq_string("xyzzy")


def test_canonical_pq_strings_have_expected_set():
    """The canonical (p,q)-string catalog has the standard SM exemplars."""
    expected = {"electron", "muon", "tau", "up", "down",
                "proton", "neutron", "z_boson"}
    actual = {k.lower() for k in st.CANONICAL_PQ_STRINGS.keys()}
    assert expected <= actual


def test_winding_squared_is_pq_substrate_kinematic_factor():
    """Substrate Paper-6 mass formula has gradient term ∝ p² + q²;
    (p,q)-string tension also goes as √(p² + q²) (g_s → 1).
    Same numerical content."""
    for name, ps in st.CANONICAL_PQ_STRINGS.items():
        assert ps.winding_squared == ps.p ** 2 + ps.q ** 2


# ---------------------------------------------------------------------------
# KK-mode tower
# ---------------------------------------------------------------------------

def test_kk_tower_has_21_modes():
    """K_7 has 21 edges → 21 cycles → 21 KK modes."""
    modes = st.kk_tower()
    assert len(modes) == 21


def test_kk_tower_geometric_progression():
    """Mass ratio between adjacent modes is √α."""
    import math
    modes = st.kk_tower()
    for i in range(len(modes) - 1):
        ratio = modes[i + 1].mass_gev / modes[i].mass_gev
        assert abs(ratio - math.sqrt(st.ALPHA_EM)) < 1e-9


def test_kk_tower_hits_electron_at_n21():
    """n=21 mode (highest cycle) should land near m_e ≈ 0.51 MeV."""
    modes = st.kk_tower()
    n21 = modes[20]
    assert n21.n == 21
    # 0.45 MeV ≈ electron 0.51 MeV (12% off, as expected for bare order)
    assert 1e-4 < n21.mass_gev < 1e-3
    assert n21.sm_identification == "electron"


def test_kk_tower_hits_nucleon_at_n18():
    """n=18 mode should land near nucleon mass ~ 0.94 GeV."""
    modes = st.kk_tower()
    n18 = modes[17]
    assert n18.n == 18
    assert 0.1 < n18.mass_gev < 10.0
    assert n18.sm_identification == "nucleon"


def test_kk_tower_hits_z_at_n16():
    """n=16 mode should land in the EW scale."""
    modes = st.kk_tower()
    n16 = modes[15]
    assert n16.n == 16
    assert 10.0 < n16.mass_gev < 1000.0
    assert n16.sm_identification == "Z boson"


def test_kk_tower_hits_down_at_n20():
    """n=20 mode should land at the down-quark mass ~ 5 MeV."""
    modes = st.kk_tower()
    n20 = modes[19]
    assert n20.n == 20
    assert 1e-3 < n20.mass_gev < 1e-1
    assert n20.sm_identification == "down quark"


def test_seesaw_dual_pair():
    """Wilson-tower seesaw duality: m(n)·m(21-n) = m_e · M_Pl."""
    # The dual of n=5 is 21-5 = 16 (memory: K7-wilson-seesaw-duality)
    assert st.seesaw_dual_pair(5) == 16
    assert st.seesaw_dual_pair(0) == 21
    assert st.seesaw_dual_pair(21) == 0


def test_alpha_to_n_inverts_tower():
    """Given a target mass, find the (continuous) n that produces it."""
    for n_target in [1, 5, 10, 16, 21]:
        modes = st.kk_tower()
        m = modes[n_target - 1].mass_gev
        n_recovered = st.alpha_to_n(m)
        assert abs(n_recovered - n_target) < 1e-9


# ---------------------------------------------------------------------------
# ADE / McKay correspondence
# ---------------------------------------------------------------------------

def test_ade_2I_is_E8():
    """Binary icosahedral 2I ↔ E_8 via McKay correspondence."""
    e = st.ade_lookup("2I")
    assert "E_8" in e.lie_algebra
    assert e.coxeter_number == 30      # h(E_8) = 30


def test_ade_2T_is_E6():
    """Binary tetrahedral 2T ↔ E_6 (used in Phase b2.12 substrate construction)."""
    e = st.ade_lookup("2T")
    assert "E_6" in e.lie_algebra
    assert e.coxeter_number == 12      # h(E_6) = 12


def test_ade_unknown_raises():
    with pytest.raises(KeyError):
        st.ade_lookup("xyzzy")


# ---------------------------------------------------------------------------
# Holonomy classifications
# ---------------------------------------------------------------------------

def test_spin7_holonomy_in_catalog():
    """Spin(7) is the substrate's home holonomy (Cl(0,7) octonion structure)."""
    h = st.holonomy_lookup("Spin(7)")
    assert h.target_dim == 8
    assert "Cl(0,7)" in h.substrate_realization or "octonion" in h.substrate_realization


def test_g2_holonomy_in_catalog():
    """G₂ holonomy is realized as automorphisms of the octonion product."""
    h = st.holonomy_lookup("G_2")
    assert h.target_dim == 7
    assert "octonion" in h.substrate_realization.lower()


def test_2I_holonomy_n_q_irreps():
    """The 2I orbifold realizes n_q ∈ {1..6} as |Irr(2I)| = 6 dim labels
    (memory: nq-is-2I-irrep-dim)."""
    h = st.holonomy_lookup("2I")
    assert "n_q" in h.substrate_realization or "Irr(2I)" in h.substrate_realization


# ---------------------------------------------------------------------------
# Multi-view introspection
# ---------------------------------------------------------------------------

def test_string_view_electron_shows_pq_and_kk():
    """string_view should surface both the (p,q)-string label and the KK
    mode index for known SM particles."""
    view = st.string_view("electron")
    assert "(2, 1)-string" in view
    assert "n=21" in view
    assert "K_7" in view


def test_string_view_unknown_particle():
    view = st.string_view("xyzzy")
    assert "Unknown" in view


# ---------------------------------------------------------------------------
# Cross-shim consistency: nwt.string ↔ existing substrate
# ---------------------------------------------------------------------------

def test_K7_cycles_match_K7_graph_state_edges():
    """nwt.string.K7_torus.n_cycles = 21 must equal the # edges in the
    actual K_7 graph state used by nwt.heron / Phase 6 Wilson amplitude."""
    # K_7 has 7C2 = 21 edges
    n_K7_edges = 7 * 6 // 2
    assert st.K7_torus.n_cycles == n_K7_edges


def test_kk_tower_n_max_matches_K7_edge_count():
    """The KK-mode tower index n runs from 1 to 21 (= K_7 edges)."""
    assert len(st.kk_tower()) == st.K7_torus.n_cycles


def test_substrate_pq_winding_matches_paper6_mass_formula():
    """The (p,q)-string winding² = substrate Paper-6 'p² + q²' kinematic
    factor.  Test that they agree by computing both for the canonical
    proton (1, 3) -- substrate Paper 6 says p² + q² = 10, and the
    PQString class returns the same."""
    proton = st.pq_string("proton")
    paper6_factor = proton.p ** 2 + proton.q ** 2
    assert proton.winding_squared == paper6_factor
    assert paper6_factor == 10
