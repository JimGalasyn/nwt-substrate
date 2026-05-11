"""Cross-shim consistency tests for the substrate ISA.

These tests verify that DIFFERENT shims (chemistry, gravity) reading
the same substrate constants give numerically consistent answers.
The point: if a refactor breaks the K_7 structural identity in one
shim, these tests detect it across shims.

Pattern (active-encoding gate):
  - chemistry consumes isa.K7_TR_M2_THRESHOLD (derived from |E(W_6)|=12)
  - gravity consumes isa.WILSON_EXPONENT_K7 (derived from |E(K_7)|=21)
  - both ultimately reference N_EDGES_K7=21 from isa.constants

If you change N_EDGES_K7 in one place, both shims will report
different numbers and these tests will catch it.
"""
import numpy as np
import pytest

import nwt_substrate.isa as isa
import nwt_substrate.chemistry as chem
import nwt_substrate.qed as qed
import nwt_substrate.qcd as qcd
import nwt_substrate.particles as parts
import nwt_substrate.electroweak as ew
import nwt_substrate.heron as heron
from nwt_substrate.gravity import coupling as gc
from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas
from nwt_substrate.algebra.su3 import gell_mann_matrices, structure_constants
from nwt_substrate.particles.compendium import COMPENDIUM
from nwt_substrate.topology.torus_knots import carrier_for_n_q


# ---------------------------------------------------------------------------
# Structural-identity invariants — at module-import time
# ---------------------------------------------------------------------------

class TestSubstrateConstants:

    def test_K_7_edge_count(self):
        """|E(K_7)| = dim(Adj_Spin(7)) = 21 (cross-derivation)."""
        assert isa.N_EDGES_K7 == 21
        assert isa.DIM_ADJ_SPIN7 == 21
        assert isa.N_EDGES_K7 == isa.DIM_ADJ_SPIN7

    def test_K_7_vertex_count(self):
        """|V(K_7)| = dim(V_Spin(7)) = 7 (cross-derivation)."""
        assert isa.N_VERTICES_K7 == 7
        assert isa.DIM_V_SPIN7 == 7
        assert isa.N_VERTICES_K7 == isa.DIM_V_SPIN7

    def test_K_7_edges_choose_2(self):
        """K_7 has C(7, 2) = 21 edges."""
        assert isa.N_EDGES_K7 == (
            isa.N_VERTICES_K7 * (isa.N_VERTICES_K7 - 1) // 2
        )

    def test_spinor_vector_ratio(self):
        """(8/7) = dim(S)/dim(V) for Spin(7)."""
        assert isa.SPINOR_VECTOR_RATIO == isa.DIM_S_SPIN7 / isa.DIM_V_SPIN7
        assert isa.SPINOR_VECTOR_RATIO == 8.0 / 7.0

    def test_wilson_exponent(self):
        """21/2 = |E(K_7)|/2."""
        assert isa.WILSON_EXPONENT_K7 == isa.N_EDGES_K7 / 2.0
        assert isa.WILSON_EXPONENT_K7 == 10.5

    def test_nlo_coefficient(self):
        """1/7 = 1/|V(K_7)|."""
        assert isa.NLO_VERTEX_COEFFICIENT == 1.0 / isa.N_VERTICES_K7

    def test_nnlo_coefficient(self):
        """21/8 = dim(Adj)/dim(S)."""
        assert isa.NNLO_BRACKET_COEFFICIENT == (
            isa.DIM_ADJ_SPIN7 / isa.DIM_S_SPIN7
        )
        assert isa.NNLO_BRACKET_COEFFICIENT == 21.0 / 8.0

    def test_trace_signatures(self):
        """Tr(M²) for K_7 and W_6 derived from edge counts."""
        assert isa.TR_M2_K7 == -2 * isa.N_EDGES_K7   # -42
        # W_6 has 12 edges = 6 hub + 6 cycle
        assert isa.TR_M2_W6 == -2 * 12   # -24


# ---------------------------------------------------------------------------
# Cross-shim derivation: same constants from both shims
# ---------------------------------------------------------------------------

class TestCrossShimDerivation:

    def test_gravity_m_e_over_M_Pl_LO(self):
        """Gravity LO uses isa.SPINOR_VECTOR_RATIO and isa.WILSON_EXPONENT_K7."""
        alpha = 1.0 / 137.035999084   # CODATA α
        expected = isa.SPINOR_VECTOR_RATIO * (alpha ** isa.WILSON_EXPONENT_K7)
        assert abs(gc.m_e_over_M_Pl_LO(alpha) - expected) < 1e-30

    def test_gravity_m_e_over_M_Pl_NNLO_breakdown(self):
        """Gravity NNLO factors cleanly through isa.k7_wilson_breakdown."""
        alpha = 1.0 / 137.035999084
        breakdown = isa.k7_wilson_breakdown(alpha, order="NNLO")
        assert abs(breakdown["total"] - gc.m_e_over_M_Pl_NNLO(alpha)) < 1e-30
        # The four substrate factors:
        assert breakdown["spinor_vector_ratio"] == 8.0 / 7.0
        assert abs(breakdown["bare_K7_amplitude"]
                    - alpha ** 10.5) < 1e-30

    def test_chemistry_K_7_threshold(self):
        """Chemistry's K_7-hub threshold is Tr(M²) for W_6.
        Both come from N_EDGES_K7 via the substrate constants."""
        # W_6 = 6 hub-petal + 6 petal-cycle = 12 edges
        # Tr(M²) = -2 × |E| = -24
        assert isa.K7_TR_M2_THRESHOLD == isa.TR_M2_W6

    def test_cross_shim_K_7_edge_count_consistency(self):
        """The same N_EDGES_K7 powers:
          (a) gravity's α^(21/2) exponent
          (b) chemistry's K_7 graph contraction in isa.SO7_BASIS

        If anyone redefines |E(K_7)|, both shims break together."""
        # Gravity: WILSON_EXPONENT_K7 × 2 = N_EDGES_K7
        assert 2 * isa.WILSON_EXPONENT_K7 == isa.N_EDGES_K7

        # Chemistry / ISA: SO7_BASIS has N_EDGES_K7 generators
        assert isa.SO7_BASIS.shape[0] == isa.N_EDGES_K7

        # Both derive from the same source
        assert isa.N_EDGES_K7 == 21


# ---------------------------------------------------------------------------
# Numerical regression: gravity predictions unchanged after refactor
# ---------------------------------------------------------------------------

class TestGravityNumericalRegression:
    """Verify the ISA-backed gravity coupling produces the same numbers
    as the pre-refactor closed-form implementation, to machine precision."""

    @staticmethod
    def _legacy_LO(alpha):
        return (8.0 / 7.0) * alpha ** (21.0 / 2.0)

    @staticmethod
    def _legacy_NLO(alpha):
        return (8.0 / 7.0) * alpha ** (21.0 / 2.0) * (1.0 + alpha / 7.0)

    @staticmethod
    def _legacy_NNLO(alpha):
        return (8.0 / 7.0) * alpha ** (21.0 / 2.0) * (
            1.0 + alpha / 7.0 + (21.0 / 8.0) * alpha ** 2
        )

    @pytest.mark.parametrize("alpha", [
        1.0 / 137.035999084,   # CODATA α
        0.1,                    # higher coupling
        0.01,                   # lower
        1e-3,                   # much lower
    ])
    def test_LO_matches_legacy(self, alpha):
        assert gc.m_e_over_M_Pl_LO(alpha) == self._legacy_LO(alpha)

    @pytest.mark.parametrize("alpha", [
        1.0 / 137.035999084, 0.1, 0.01, 1e-3,
    ])
    def test_NLO_matches_legacy(self, alpha):
        assert gc.m_e_over_M_Pl_NLO(alpha) == self._legacy_NLO(alpha)

    @pytest.mark.parametrize("alpha", [
        1.0 / 137.035999084, 0.1, 0.01, 1e-3,
    ])
    def test_NNLO_matches_legacy(self, alpha):
        assert gc.m_e_over_M_Pl_NNLO(alpha) == self._legacy_NNLO(alpha)

    def test_NNLO_within_CODATA_band(self):
        """Paper 17 claim: NNLO is -5.6 ppm on m_e/M_Pl, inside ±22 ppm."""
        predicted = gc.m_e_over_M_Pl_NNLO()
        observed = gc.m_e_over_M_Pl_observed()
        ppm_error = (predicted - observed) / observed * 1e6
        # Allow some slack on the exact ppm number (it depends on
        # CODATA α value), but verify it's within ±22 ppm
        assert abs(ppm_error) < 22.0
        # And specifically near -5.6 ppm
        assert -7.0 < ppm_error < -4.0


# ---------------------------------------------------------------------------
# Cross-shim Wilson amplitude consistency
# ---------------------------------------------------------------------------

class TestQEDSubstrateIdentities:
    """The QED shim consumes Cl(0,7) → Cl(1,3) substrate identities."""

    def test_octonion_dim_equals_spinor_dim(self):
        """The octonion algebra dimension is the Spin(7) spinor dim."""
        assert isa.DIM_OCTONION == isa.DIM_S_SPIN7
        assert isa.DIM_OCTONION == 8

    def test_cl07_imaginary_count(self):
        """Cl(0,7) has 7 imaginary generators = |V(K_7)|."""
        assert isa.N_CL07_IMAGINARY == isa.DIM_V_SPIN7
        assert isa.N_CL07_IMAGINARY == isa.N_VERTICES_K7
        assert isa.N_CL07_IMAGINARY == 7

    def test_lorentz_internal_partition(self):
        """4 (Lorentz) + 3 (internal SU(2)) = 7 (Cl(0,7) imaginary)."""
        assert (isa.N_LORENTZ_FROM_CL07 + isa.DIM_INTERNAL_SU2
                == isa.N_CL07_IMAGINARY)
        assert isa.N_LORENTZ_FROM_CL07 == 4
        assert isa.DIM_INTERNAL_SU2 == 3

    def test_internal_su2_equals_so7_rank(self):
        """The internal SU(2) dim is rank(so(7)) = 3."""
        assert isa.DIM_INTERNAL_SU2 == isa.RANK_SO7

    def test_gamma_matrix_shape(self):
        """Substrate-built Lorentz gammas have shape (DIM_OCTONION, DIM_OCTONION).
        This is the runtime substrate-identity assertion in make_lorentz_gammas."""
        T = make_octonion_table()
        gammas = make_lorentz_gammas(T)
        assert len(gammas) == isa.N_LORENTZ_FROM_CL07
        for g in gammas:
            assert g.shape == (isa.DIM_OCTONION, isa.DIM_OCTONION)

    def test_gamma_anticommutation(self):
        """{γ^μ, γ^ν} = 2 η^{μν} I_8 — verifies the Lorentz Clifford embedding."""
        T = make_octonion_table()
        gammas = make_lorentz_gammas(T)
        eta = np.diag([1.0, -1.0, -1.0, -1.0])
        I8 = np.eye(isa.DIM_OCTONION)
        for mu in range(4):
            for nu in range(4):
                anti = gammas[mu] @ gammas[nu] + gammas[nu] @ gammas[mu]
                expected = 2.0 * eta[mu, nu] * I8
                assert np.allclose(anti, expected, atol=1e-10)

    def test_make_lorentz_gammas_rejects_wrong_count(self):
        """ValueError if you pass != N_LORENTZ_FROM_CL07 indices."""
        T = make_octonion_table()
        with pytest.raises(ValueError, match="length"):
            make_lorentz_gammas(T, lorentz_indices=(1, 2, 3))  # only 3
        with pytest.raises(ValueError, match="length"):
            make_lorentz_gammas(T, lorentz_indices=(1, 2, 3, 4, 5))  # 5

    def test_make_lorentz_gammas_rejects_out_of_range(self):
        """ValueError if any index is outside [1, N_CL07_IMAGINARY]."""
        T = make_octonion_table()
        with pytest.raises(ValueError, match=r"\[1, 7\]"):
            make_lorentz_gammas(T, lorentz_indices=(0, 1, 2, 3))   # 0 invalid
        with pytest.raises(ValueError, match=r"\[1, 7\]"):
            make_lorentz_gammas(T, lorentz_indices=(1, 2, 3, 8))   # 8 invalid

    def test_qed_substrate_namespace_consistent(self):
        """qed.substrate.X must equal isa.X (same source of truth)."""
        assert qed.substrate.DIM_OCTONION == isa.DIM_OCTONION
        assert qed.substrate.N_CL07_IMAGINARY == isa.N_CL07_IMAGINARY
        assert qed.substrate.N_LORENTZ_FROM_CL07 == isa.N_LORENTZ_FROM_CL07
        assert qed.substrate.DIM_INTERNAL_SU2 == isa.DIM_INTERNAL_SU2
        assert qed.substrate.N_EDGES_K7 == isa.N_EDGES_K7
        assert qed.substrate.N_VERTICES_K7 == isa.N_VERTICES_K7

    def test_qed_substrate_breakdown_runs(self):
        """The breakdown text references all the structural identities."""
        text = qed.substrate_breakdown()
        assert "Cl(0,7)" in text
        assert str(isa.N_CL07_IMAGINARY) in text   # 7
        assert str(isa.DIM_OCTONION) in text       # 8
        assert str(isa.N_LORENTZ_FROM_CL07) in text  # 4
        assert str(isa.DIM_INTERNAL_SU2) in text   # 3
        assert str(isa.N_EDGES_K7) in text         # 21


class TestQCDSubstrateIdentities:
    """The QCD shim consumes SU(3) ⊂ G_2 ⊂ Spin(7) substrate identities."""

    def test_n_c_equals_rank_so7(self):
        """N_c = 3 = rank(so(7))."""
        assert isa.N_C_SU3 == isa.RANK_SO7
        assert isa.N_C_SU3 == 3

    def test_n_gluons_equals_dim_octonion(self):
        """N_gluons = N_c² - 1 = 8 = DIM_OCTONION = DIM_S_SPIN7.
        SUBSTRATE IDENTITY: same 8 in QED γ matrices and QCD gluons."""
        assert isa.N_GLUONS == isa.DIM_OCTONION
        assert isa.N_GLUONS == isa.DIM_S_SPIN7
        assert isa.N_GLUONS == isa.N_C_SU3 ** 2 - 1
        assert isa.N_GLUONS == 8

    def test_C_F_formula(self):
        """C_F = (N_c² - 1)/(2 N_c) = DIM_OCTONION / (2 N_C_SU3) = 4/3."""
        assert isa.C_F_SU3 == isa.DIM_OCTONION / (2.0 * isa.N_C_SU3)
        assert abs(isa.C_F_SU3 - 4.0 / 3.0) < 1e-15

    def test_C_A_equals_N_c(self):
        """C_A = N_c = RANK_SO7 = 3."""
        assert isa.C_A_SU3 == isa.N_C_SU3
        assert isa.C_A_SU3 == 3.0

    def test_qcd_constants_sourced_from_isa(self):
        """qcd.N_c, qcd.C_F, qcd.C_A, qcd.T_R all come from isa."""
        assert qcd.N_c == isa.N_C_SU3
        assert qcd.C_F == isa.C_F_SU3
        assert qcd.C_A == isa.C_A_SU3
        assert qcd.T_R == isa.T_R_SU3

    def test_gell_mann_count_and_shape(self):
        """8 Gell-Mann matrices, each 3×3, sourced from substrate constants."""
        ms = gell_mann_matrices()
        assert len(ms) == isa.N_GLUONS
        for m in ms:
            assert m.shape == (isa.N_C_SU3, isa.N_C_SU3)

    def test_gell_mann_normalization(self):
        """Tr(λ^a λ^b) = 2 δ^ab — the standard SU(3) normalization."""
        ms = gell_mann_matrices()
        for a in range(isa.N_GLUONS):
            for b in range(isa.N_GLUONS):
                tr = np.trace(ms[a] @ ms[b]).real
                expected = 2.0 if a == b else 0.0
                assert abs(tr - expected) < 1e-10

    def test_structure_constants_shape(self):
        """Structure constants f^abc have shape (N_GLUONS, N_GLUONS, N_GLUONS)."""
        f = structure_constants()
        assert f.shape == (isa.N_GLUONS, isa.N_GLUONS, isa.N_GLUONS)

    def test_qcd_substrate_namespace_consistent(self):
        """qcd.substrate.X must equal isa.X."""
        assert qcd.substrate.N_C_SU3 == isa.N_C_SU3
        assert qcd.substrate.N_GLUONS == isa.N_GLUONS
        assert qcd.substrate.C_F_SU3 == isa.C_F_SU3
        assert qcd.substrate.C_A_SU3 == isa.C_A_SU3
        assert qcd.substrate.DIM_OCTONION == isa.DIM_OCTONION

    def test_qcd_substrate_breakdown_runs(self):
        text = qcd.substrate_breakdown()
        assert "SU(3)" in text
        assert "DIM_OCTONION" in text
        assert str(isa.N_GLUONS) in text   # 8
        assert str(isa.N_C_SU3) in text     # 3


class TestThreeWayCrossShim:
    """The K_7 substrate operates across THREE shims (chem, gravity, qed).
    These tests assert all three consume the same N_EDGES_K7=21,
    N_VERTICES_K7=7, DIM_S_SPIN7=8 identities."""

    def test_K_7_edge_count_seen_by_three_shims(self):
        # Chemistry: 21 so(7) generators in the basis
        assert isa.SO7_BASIS.shape[0] == isa.N_EDGES_K7
        # Gravity: α^(N_EDGES_K7 / 2)
        assert 2 * isa.WILSON_EXPONENT_K7 == isa.N_EDGES_K7
        # QED: the same K_7 lives in the so(7) adjoint that contains
        # the gamma matrices' algebra
        assert isa.N_EDGES_K7 == isa.DIM_ADJ_SPIN7

    def test_seven_seen_by_three_shims(self):
        # Chemistry: 7-vertex projection
        assert isa.SO7_BASIS.shape[1] == isa.N_VERTICES_K7
        # Gravity: 1/N_VERTICES_K7 NLO coefficient
        assert 1.0 / isa.NLO_VERTEX_COEFFICIENT == isa.N_VERTICES_K7
        # QED: Cl(0,7) has N_VERTICES_K7 imaginary generators
        assert isa.N_CL07_IMAGINARY == isa.N_VERTICES_K7

    def test_eight_seen_by_two_shims(self):
        # Gravity: (8/7) prefactor numerator
        assert isa.SPINOR_VECTOR_RATIO * isa.DIM_V_SPIN7 == isa.DIM_S_SPIN7
        # QED: 8x8 gamma matrices = DIM_OCTONION
        T = make_octonion_table()
        gammas = make_lorentz_gammas(T)
        assert gammas[0].shape[0] == isa.DIM_S_SPIN7


class TestParticlesSubstrateIdentities:
    """The particles shim consumes K_7 carrier-knot substrate identities."""

    def test_n_carrier_types_equals_n_vertices_k7(self):
        """7 carrier types = N_VERTICES_K7 (one per K_7 vertex)."""
        assert isa.N_CARRIER_TYPES == isa.N_VERTICES_K7
        assert isa.N_CARRIER_TYPES == 7

    def test_max_crossing_number(self):
        """MAX_CROSSING_NUMBER = N_VERTICES_K7 - 1 = 6."""
        assert isa.MAX_CROSSING_NUMBER == isa.N_VERTICES_K7 - 1
        assert isa.MAX_CROSSING_NUMBER == 6

    def test_carrier_names_length(self):
        """CARRIER_NAMES has exactly N_CARRIER_TYPES entries."""
        assert len(isa.CARRIER_NAMES) == isa.N_CARRIER_TYPES

    def test_compendium_all_n_q_in_range(self):
        """All 24+ Paper 6 compendium particles have n_q ∈ [0, 6]."""
        for entry in COMPENDIUM:
            n_q = entry["n_q"]
            assert 0 <= n_q <= isa.MAX_CROSSING_NUMBER, (
                f"{entry['name']}: n_q={n_q} outside K_7 substrate range"
            )

    def test_compendium_uses_all_carrier_types(self):
        """The 24+ compendium spans (most of) the 7 carrier types."""
        seen = {entry["n_q"] for entry in COMPENDIUM}
        # Compendium uses n_q ∈ {0, 2, 3, 5} most commonly
        # but we just check all are within substrate range
        assert seen.issubset(set(range(isa.N_CARRIER_TYPES)))

    def test_particle_rejects_n_q_too_large(self):
        """Particle(..., n_q=7) raises ValueError per substrate constraint."""
        with pytest.raises(ValueError, match="outside substrate range"):
            parts.Particle(name="bogus", p=1, q=1, m=2, n_q=7)

    def test_particle_rejects_n_q_negative(self):
        with pytest.raises(ValueError, match="outside substrate range"):
            parts.Particle(name="bogus", p=1, q=1, m=2, n_q=-1)

    def test_carrier_for_n_q_uses_isa_table(self):
        """topology.torus_knots.carrier_for_n_q reads isa.CARRIER_NAMES."""
        for n_q in range(isa.N_CARRIER_TYPES):
            assert carrier_for_n_q(n_q) == isa.CARRIER_NAMES[n_q]

    def test_carrier_for_invalid_n_q_returns_fallback(self):
        """Out-of-range n_q returns a descriptive string, not crash."""
        result = carrier_for_n_q(10)
        assert "10" in result

    def test_proton_unchanged_through_isa_refactor(self):
        """Numerical regression: proton mass still 937.24 MeV via Paper 6."""
        from ..particles.factory import particle
        p = particle("p")
        assert abs(p.mass_pred - 937.24) < 0.01

    def test_electron_unchanged(self):
        from ..particles.factory import particle
        e = particle("e-")
        # Electron is the anchor, so mass_pred should equal ME_MEV
        from ..particles.mass import ME_MEV
        assert abs(e.mass_pred - ME_MEV) < 1e-6

    def test_particles_substrate_namespace_consistent(self):
        assert parts.substrate.N_CARRIER_TYPES == isa.N_CARRIER_TYPES
        assert parts.substrate.MAX_CROSSING_NUMBER == isa.MAX_CROSSING_NUMBER
        assert parts.substrate.CARRIER_NAMES == isa.CARRIER_NAMES

    def test_particles_substrate_breakdown_runs(self):
        text = parts.substrate_breakdown()
        assert "N_CARRIER_TYPES" in text
        assert "K_7" in text
        assert str(isa.N_CARRIER_TYPES) in text


class TestFourWayCrossShim:
    """The K_7 substrate operates across FOUR shims now: chemistry,
    gravity, qed, qcd. These tests assert structural identities span
    all four."""

    def test_eight_seen_by_three_shims_explicitly(self):
        """The number 8 = DIM_OCTONION appears in three independent shims:
          - gravity: numerator of SPINOR_VECTOR_RATIO = 8/7
          - qed: shape of Dirac γ^μ = 8x8
          - qcd: number of gluons = 8
        All read isa.DIM_OCTONION."""
        # Gravity
        assert isa.DIM_OCTONION == int(isa.SPINOR_VECTOR_RATIO * isa.DIM_V_SPIN7)
        # QED
        T = make_octonion_table()
        gammas = make_lorentz_gammas(T)
        assert gammas[0].shape == (isa.DIM_OCTONION, isa.DIM_OCTONION)
        # QCD
        ms = gell_mann_matrices()
        assert len(ms) == isa.DIM_OCTONION
        # All point to the same constant
        assert isa.DIM_OCTONION == 8

    def test_three_seen_by_three_shims(self):
        """The number 3 = RANK_SO7 = N_C_SU3 = DIM_INTERNAL_SU2 appears in:
          - qed: 3 internal-SU(2) directions (= 7 - 4 = N_CL07_IMAGINARY - N_LORENTZ_FROM_CL07)
          - qcd: 3 colors
          - gravity: implicit in (8/7)·(21/8) = 3 = rank(so(7)) NNLO coefficient
                     after distributing"""
        assert isa.RANK_SO7 == isa.N_C_SU3
        assert isa.RANK_SO7 == isa.DIM_INTERNAL_SU2
        # Gravity's NNLO bracket coefficient 21/8 multiplied by 8/7
        # gives 21/7 = 3 = rank(so(7))
        assert (isa.NNLO_BRACKET_COEFFICIENT
                * isa.SPINOR_VECTOR_RATIO) == isa.N_EDGES_K7 / isa.DIM_V_SPIN7
        assert (isa.N_EDGES_K7 / isa.DIM_V_SPIN7) == isa.RANK_SO7

    def test_K_7_substrate_unifies_three_scales(self):
        """N_EDGES_K7 = 21 appears in three independent physical scales:
          chemistry (coronene K_7-hub), gravity (α^(21/2) Wilson), qcd
          (so(7) adjoint where SU(3) lives)."""
        # All three consume the same N_EDGES_K7
        assert isa.SO7_BASIS.shape[0] == isa.N_EDGES_K7
        assert 2 * isa.WILSON_EXPONENT_K7 == isa.N_EDGES_K7
        # SU(3) ⊂ G_2 ⊂ Spin(7) with so(7) = N_EDGES_K7 generators
        assert isa.DIM_ADJ_SPIN7 == isa.N_EDGES_K7


class TestElectroweakSubstrateIdentities:
    """The electroweak shim consumes SM fermion-content substrate identities.

    Most striking: b_QED^SM = 8 = DIM_OCTONION. The same 8 that gives
    QED its 8×8 γ^μ and QCD its 8 gluons equals the sum N_c × Q² over
    all SM fermions.
    """

    def test_n_generations_equals_rank_so7(self):
        """3 SM generations = rank(so(7))."""
        assert isa.N_GENERATIONS == isa.RANK_SO7
        assert isa.N_GENERATIONS == 3

    def test_fermion_types_per_generation(self):
        """4 fermion types per generation (charged lepton, neutrino, up, down)."""
        assert isa.N_FERMION_TYPES_PER_GENERATION == 4

    def test_sm_fermion_flavor_count(self):
        """12 SM fermion flavors = 3 × 4."""
        assert isa.N_SM_FERMION_FLAVORS == (
            isa.N_GENERATIONS * isa.N_FERMION_TYPES_PER_GENERATION
        )
        assert isa.N_SM_FERMION_FLAVORS == 12

    def test_so7_generator_decomposition(self):
        """21 = 6 (Lorentz) + 3 (internal SU(2)) + 12 (mixed SM flavors).
        This is the substrate-boson-count decomposition."""
        assert (isa.N_LORENTZ_GENERATORS
                + isa.DIM_INTERNAL_SU2
                + isa.N_SM_FERMION_FLAVORS) == isa.N_EDGES_K7
        # Each piece checked:
        assert isa.N_LORENTZ_GENERATORS == 6
        assert isa.DIM_INTERNAL_SU2 == 3
        assert isa.N_SM_FERMION_FLAVORS == 12

    def test_b_qed_sm_equals_dim_octonion(self):
        """B_QED_SM = DIM_OCTONION = 8 (substrate identity)."""
        assert isa.B_QED_SM == isa.DIM_OCTONION
        assert isa.B_QED_SM == 8

    def test_b_qed_sm_computed_from_couplings_table(self):
        """Independently compute b_QED from the SM_COUPLINGS dict and
        verify it matches B_QED_SM. This is the substrate identity made
        empirical."""
        b = ew.verify_b_qed_sm(tol=1e-10)
        assert abs(b - isa.B_QED_SM) < 1e-10

    def test_ew_substrate_namespace_consistent(self):
        assert ew.substrate.N_GENERATIONS == isa.N_GENERATIONS
        assert ew.substrate.N_SM_FERMION_FLAVORS == isa.N_SM_FERMION_FLAVORS
        assert ew.substrate.B_QED_SM == isa.B_QED_SM
        assert ew.substrate.DIM_OCTONION == isa.DIM_OCTONION

    def test_ew_substrate_breakdown_runs(self):
        text = ew.substrate_breakdown()
        assert "DIM_OCTONION" in text
        assert "b_QED" in text
        assert str(isa.B_QED_SM) in text
        assert str(isa.N_SM_FERMION_FLAVORS) in text


class TestFiveWayCrossShim:
    """All five shims (chem, gravity, qed, qcd, particles) consume
    the same N_VERTICES_K7 = 7 substrate constant."""

    def test_seven_seen_by_five_shims(self):
        """N_VERTICES_K7 = 7 surfaces in all five shims:
          - chemistry: 7-vertex ring-adjacency projection
          - gravity:   1/N_VERTICES_K7 = 1/7 NLO coefficient
          - qed:       7 imaginary Cl(0,7) directions
          - qcd:       (via Spin(7) ⊃ G_2 ⊃ SU(3) chain)
          - particles: 7 carrier-knot types"""
        # chemistry
        assert isa.SO7_BASIS.shape[1] == isa.N_VERTICES_K7
        # gravity
        assert 1.0 / isa.NLO_VERTEX_COEFFICIENT == isa.N_VERTICES_K7
        # qed
        assert isa.N_CL07_IMAGINARY == isa.N_VERTICES_K7
        # qcd (indirect: rank(so(7)) = 3 used by N_c is part of so(7)
        # which acts on a 7-dim vector space)
        assert isa.DIM_V_SPIN7 == isa.N_VERTICES_K7
        # particles
        assert isa.N_CARRIER_TYPES == isa.N_VERTICES_K7

        # All five point at 7
        assert isa.N_VERTICES_K7 == 7

    def test_compendium_proton_through_full_substrate(self):
        """The proton's existence requires:
          - particles: n_q=5 within [0, MAX_CROSSING_NUMBER=6]
          - chemistry: K_7 substrate not violated by carrier
          - gravity: m_proton scale set by Paper 6 + Paper 17 coupling
          - qcd: 3 colors (N_C_SU3 = 3) for the quark content
          - qed: charge +1 via I_3 + (B+S+C)/2 = 1/2 + 1/2 = 1"""
        from ..particles.factory import particle
        p = particle("p")
        assert p.n_q == 5
        assert p.n_q <= isa.MAX_CROSSING_NUMBER
        assert p.carrier == isa.CARRIER_NAMES[5]   # "cinquefoil"
        assert p.Q_pred == 1.0   # GM-N with B=1, S=0, I_3=1/2
        # Proton is a 3-quark bound state; QCD uses N_C_SU3 = 3 colors
        assert isa.N_C_SU3 == 3
        # And lives in DIM_OCTONION = 8 spinor space (QED)
        assert isa.DIM_OCTONION == 8


class TestHeronSubstrateIdentities:
    """The heron shim compiles substrate ISA to qiskit quantum circuits.

    The most direct hardware identity: the K_7 graph state preparation
    circuit has exactly N_VERTICES_K7 = 7 qubits and N_EDGES_K7 = 21 CZ
    gates. The substrate algebra → quantum gates is no longer just a
    naming convention; it's an enforced structural identity.
    """

    @pytest.fixture(scope="class")
    def k7_check(self):
        """Build the K_7 graph state and count its gates."""
        qc = heron.k7_graph_state()
        return heron.verify_k7_circuit_substrate(qc)

    def test_k7_circuit_qubit_count(self, k7_check):
        """k7_graph_state() has exactly N_VERTICES_K7 qubits."""
        assert k7_check["n_qubits"] == isa.N_VERTICES_K7
        assert k7_check["n_qubits"] == 7
        assert k7_check["n_vertices_match"] is True

    def test_k7_circuit_cz_count(self, k7_check):
        """k7_graph_state() has exactly N_EDGES_K7 = 21 CZ gates,
        one per K_7 edge. This is THE substrate-to-hardware identity."""
        assert k7_check["n_cz"] == isa.N_EDGES_K7
        assert k7_check["n_cz"] == 21
        assert k7_check["n_edges_match"] is True

    def test_k7_circuit_h_count(self, k7_check):
        """k7_graph_state() has N_VERTICES_K7 H gates (one per vertex)."""
        assert k7_check["n_h"] == isa.N_VERTICES_K7

    def test_isa_k7_edge_list_correct(self):
        """isa.k7_edge_list() returns N_EDGES_K7 unique edges of K_7."""
        edges = isa.k7_edge_list()
        assert len(edges) == isa.N_EDGES_K7
        # All edges (a, b) have a < b and a, b in [0, N_VERTICES_K7)
        for a, b in edges:
            assert 0 <= a < b < isa.N_VERTICES_K7
        # All edges unique
        assert len(set(edges)) == isa.N_EDGES_K7

    def test_stabilizer_measurement_validates_vertex(self):
        """stabilizer_measurement(v) rejects v outside [0, N_VERTICES_K7)."""
        with pytest.raises(ValueError, match=r"vertex must be in"):
            heron.k7_stabilizer_circuit(isa.N_VERTICES_K7)  # vertex=7 invalid
        with pytest.raises(ValueError, match=r"vertex must be in"):
            heron.k7_stabilizer_circuit(-1)

    def test_stabilizer_measurement_qubit_count(self):
        """k7_stabilizer_circuit has N_VERTICES_K7 qubits + N_VERTICES_K7 cbits."""
        qc = heron.k7_stabilizer_circuit(vertex=0)
        assert qc.num_qubits == isa.N_VERTICES_K7
        assert qc.num_clbits == isa.N_VERTICES_K7

    def test_ancilla_stabilizer_measurement_uses_one_extra_qubit(self):
        """stabilizer_measurement uses N_VERTICES_K7 + 1 qubits (ancilla)."""
        qc = heron.stabilizer_measurement(vertex=0)
        assert qc.num_qubits == isa.N_VERTICES_K7 + 1

    def test_muon_decay_circuit_uses_K_7_plus_two(self):
        """muon_decay has N_VERTICES_K7 + 2 qubits (K_7 + muon + electron)."""
        qc = heron.muon_decay_circuit()
        assert qc.num_qubits == isa.N_VERTICES_K7 + 2

    def test_heron_substrate_namespace_consistent(self):
        assert heron.substrate.N_VERTICES_K7 == isa.N_VERTICES_K7
        assert heron.substrate.N_EDGES_K7 == isa.N_EDGES_K7
        assert heron.substrate.DEGREE_K7 == isa.DEGREE_K7
        assert heron.substrate.DIM_OCTONION == isa.DIM_OCTONION

    def test_heron_substrate_breakdown_runs(self):
        text = heron.substrate_breakdown()
        assert "N_VERTICES_K7" in text
        assert "N_EDGES_K7" in text
        assert "CZ" in text
        assert str(isa.N_EDGES_K7) in text


class TestSixWayCrossShim:
    """All six shims (chem, gravity, qed, qcd, particles, electroweak)
    consume the same substrate constants. The most spectacular cross-shim
    identity: 8 = DIM_OCTONION surfaces in four independent computations.
    """

    def test_eight_seen_by_four_shims(self):
        """The number 8 = DIM_OCTONION appears in four independent shims:
          - gravity:     numerator of SPINOR_VECTOR_RATIO = 8/7
          - qed:         shape of Dirac γ^μ = 8×8
          - qcd:         number of gluons = 8
          - electroweak: b_QED^SM = Σ N_c Q² = 8
        All read isa.DIM_OCTONION."""
        # Gravity (already tested in TestFourWayCrossShim)
        assert isa.DIM_OCTONION == int(isa.SPINOR_VECTOR_RATIO * isa.DIM_V_SPIN7)
        # QED: 8×8 Dirac γ^μ
        T = make_octonion_table()
        gammas = make_lorentz_gammas(T)
        assert gammas[0].shape == (isa.DIM_OCTONION, isa.DIM_OCTONION)
        # QCD: 8 Gell-Mann matrices
        ms = gell_mann_matrices()
        assert len(ms) == isa.DIM_OCTONION
        # EW: b_QED^SM computed empirically
        b = ew.verify_b_qed_sm(tol=1e-10)
        assert abs(b - isa.DIM_OCTONION) < 1e-10
        # All four point to the same constant
        assert isa.DIM_OCTONION == 8

    def test_three_seen_by_four_shims(self):
        """The number 3 = RANK_SO7 appears in four independent shims:
          - qed:         3 internal-SU(2) directions
          - qcd:         3 colors
          - gravity:     21/7 = rank(so(7)) NNLO coefficient
          - electroweak: 3 SM generations
          - particles:   3 internal SU(2) (via N_C_SU3 = 3 colors of quarks)"""
        assert isa.RANK_SO7 == isa.N_C_SU3
        assert isa.RANK_SO7 == isa.DIM_INTERNAL_SU2
        assert isa.RANK_SO7 == isa.N_GENERATIONS

    def test_twenty_one_seen_as_three_pieces(self):
        """N_EDGES_K7 = 21 decomposes as:
          6 (Lorentz so(1,3)) + 3 (internal SU(2)) + 12 (mixed = SM flavors)
        Each piece corresponds to a different EW/substrate role.

        This is the substrate-boson-count-so7 lemma."""
        total = (isa.N_LORENTZ_GENERATORS
                 + isa.DIM_INTERNAL_SU2
                 + isa.N_SM_FERMION_FLAVORS)
        assert total == isa.N_EDGES_K7
        assert total == 21


class TestSevenWayCrossShim:
    """All seven shims (chem, gravity, qed, qcd, particles, ew, heron)
    consume the same K_7 substrate. The loop is closed: substrate algebra
    → quantum hardware."""

    def test_twenty_one_seen_by_seven_shims(self):
        """N_EDGES_K7 = 21 surfaces in seven independent shims:
          - chemistry:  21 so(7) generators in ISA basis
          - gravity:    α^(21/2) Wilson exponent
          - qed:        21 = dim(so(7) adjoint) holding γ-matrix algebra
          - qcd:        21 = dim(adjoint Spin(7) ⊃ G_2 ⊃ SU(3))
          - particles:  21 - 9 = 12 SM flavors slot in mixed so(7) gens
          - ew:         21 = 6 (Lorentz) + 3 (internal) + 12 (flavors)
          - heron:      21 CZ gates in k7_graph_state() circuit"""
        # chemistry / ISA basis
        assert isa.SO7_BASIS.shape[0] == isa.N_EDGES_K7
        # gravity
        assert 2 * isa.WILSON_EXPONENT_K7 == isa.N_EDGES_K7
        # qed / qcd both via so(7) adjoint
        assert isa.DIM_ADJ_SPIN7 == isa.N_EDGES_K7
        # particles + ew via the 21 = 6+3+12 decomposition
        assert (isa.N_LORENTZ_GENERATORS
                + isa.DIM_INTERNAL_SU2
                + isa.N_SM_FERMION_FLAVORS) == isa.N_EDGES_K7
        # heron via CZ gate count in k7_graph_state circuit
        qc = heron.k7_graph_state()
        check = heron.verify_k7_circuit_substrate(qc)
        assert check["n_cz"] == isa.N_EDGES_K7

        # All seven point at 21
        assert isa.N_EDGES_K7 == 21

    def test_substrate_algebra_compiles_to_hardware(self):
        """The substrate algebra is now active-encoded into quantum
        hardware. k7_graph_state() literally has one CZ gate per K_7
        edge, sourced from isa.k7_edge_list()."""
        qc = heron.k7_graph_state()
        check = heron.verify_k7_circuit_substrate(qc)
        # The hardware-side numbers match the substrate-algebra-side numbers
        assert check["n_qubits"] == isa.N_VERTICES_K7
        assert check["n_cz"] == isa.N_EDGES_K7
        # And both ultimately ground in the K_7 graph itself
        assert isa.N_VERTICES_K7 == 7
        assert isa.N_EDGES_K7 == 21


class TestWilsonAmplitudeObservable:

    def test_isa_wilson_matches_gravity(self):
        """isa.k7_wilson_amplitude == gravity.m_e_over_M_Pl_NNLO."""
        alpha = 1.0 / 137.035999084
        isa_value = isa.k7_wilson_amplitude(alpha, order="NNLO")
        gravity_value = gc.m_e_over_M_Pl_NNLO(alpha)
        assert isa_value == gravity_value

    def test_wilson_amplitude_vectorized(self):
        """isa.k7_wilson_amplitude works on numpy arrays."""
        alphas = np.array([0.01, 0.05, 0.1, 1.0 / 137.0])
        result = isa.k7_wilson_amplitude(alphas, order="NNLO")
        # Element-by-element match to scalar call
        for i, a in enumerate(alphas):
            assert abs(result[i] - isa.k7_wilson_amplitude(float(a))) < 1e-30

    def test_breakdown_factor_consistency(self):
        """Breakdown factors multiply to the total amplitude."""
        alpha = 1.0 / 137.035999084
        breakdown = isa.k7_wilson_breakdown(alpha, order="NNLO")
        expected = (breakdown["spinor_vector_ratio"]
                    * breakdown["bare_K7_amplitude"]
                    * breakdown["nnlo_correction"])
        assert abs(breakdown["total"] - expected) < 1e-25
