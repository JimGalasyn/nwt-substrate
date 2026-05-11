"""Tests for nwt_substrate.isa — the substrate-algebraic ISA layer.

Three-layer verification:
  1. so(7) primitives — basis structure, canonical-graph invariants
  2. Batched kernel — numpy and torch backends agree to machine precision
  3. Observable layer — aromaticity_score, hopf_pair_count, k7_indicator
     match the legacy chemistry-shim implementations on the reference set
"""
import numpy as np
import pytest

import nwt_substrate.isa as isa
import nwt_substrate.chemistry as chem


# ---------------------------------------------------------------------------
# Reference set: 5 PAHs + dicoronylene
# ---------------------------------------------------------------------------

REFERENCE = [
    ("benzene",      "c1ccccc1",                                 36.0),
    ("naphthalene",  "c1ccc2ccccc2c1",                           60.0),
    ("anthracene",   "c1ccc2cc3ccccc3cc2c1",                     84.0),
    ("pyrene",       "c1cc2ccc3cccc4ccc(c1)c2c34",               96.0),
    ("coronene",     "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",  200.0),
]

CORONENE_SMILES = "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"
DICORONYLENE_SMILES = f"{CORONENE_SMILES}-{CORONENE_SMILES}"


# ---------------------------------------------------------------------------
# Layer 1: so(7) primitives
# ---------------------------------------------------------------------------

class TestSO7Primitives:

    def test_basis_shape(self):
        assert isa.SO7_BASIS.shape == (21, 7, 7)

    def test_basis_antisymmetric(self):
        for k in range(21):
            T = isa.SO7_BASIS[k]
            assert np.allclose(T, -T.T)

    def test_basis_unit_norm(self):
        """Each generator has 2 nonzero entries (+1, -1), so Tr(T·Tᵀ) = 2."""
        for k in range(21):
            T = isa.SO7_BASIS[k]
            assert np.trace(T @ T.T) == 2.0

    def test_edge_index_round_trip(self):
        idx = 0
        for a in range(7):
            for b in range(a + 1, 7):
                assert isa.edge_index(a, b) == idx
                assert isa.edge_index(b, a) == idx  # order-independent
                idx += 1

    def test_canonical_K_7(self):
        """K_7 (complete): Tr(M²) = -42, |E| = 21."""
        inv = isa.canonical_invariants()["K_7"]
        assert inv["Tr_M2"] == -42.0
        assert inv["Tr_M4"] == 742.0

    def test_canonical_W_6(self):
        """W_6 (coronene ring-adj): Tr(M²) = -24 = -2 × 12 edges."""
        inv = isa.canonical_invariants()["W_6"]
        assert inv["Tr_M2"] == -24.0
        assert inv["Tr_M4"] == 140.0
        assert inv["Tr_M6"] == -948.0

    def test_canonical_empty(self):
        inv = isa.canonical_invariants()["empty"]
        assert all(v == 0.0 for v in inv.values())

    def test_canonical_single_edge(self):
        """Single edge: Tr(M²) = -2 = -2 × 1 edge."""
        inv = isa.canonical_invariants()["single_edge"]
        assert inv["Tr_M2"] == -2.0


# ---------------------------------------------------------------------------
# Layer 2: Batched kernel — backend equivalence
# ---------------------------------------------------------------------------

BACKENDS = isa.available_backends()


class TestBatchedKernel:

    def _build_test_batch(self):
        """Build a batch with all canonical graphs."""
        from nwt_substrate.isa.so7 import CANONICAL_GRAPHS, edge_list_to_so7
        names = list(CANONICAL_GRAPHS.keys())
        N = len(names)
        adj = np.zeros((N, 7, 7), dtype=np.int64)
        for i, name in enumerate(names):
            for a, b in CANONICAL_GRAPHS[name]:
                adj[i, a, b] = adj[i, b, a] = 1
        return adj, names

    @pytest.mark.parametrize("backend", BACKENDS)
    def test_kernel_matches_canonical(self, backend):
        """All backends reproduce canonical-graph invariants exactly."""
        adj, names = self._build_test_batch()
        inv = isa.graphs_to_invariants(adj, backend=backend)
        canon = isa.canonical_invariants()
        for i, name in enumerate(names):
            assert abs(inv["Tr_M2"][i] - canon[name]["Tr_M2"]) < 1e-9
            assert abs(inv["Tr_M4"][i] - canon[name]["Tr_M4"]) < 1e-9
            assert abs(inv["Tr_M6"][i] - canon[name]["Tr_M6"]) < 1e-9

    def test_backends_agree(self):
        """numpy and any available torch backends give identical results."""
        if len(BACKENDS) < 2:
            pytest.skip("Only one backend available")
        adj, _ = self._build_test_batch()
        ref = isa.graphs_to_invariants(adj, backend="numpy")
        for backend in BACKENDS[1:]:
            other = isa.graphs_to_invariants(adj, backend=backend)
            for key in ref:
                assert np.allclose(ref[key], other[key], atol=1e-9)

    def test_empty_graph_zero_invariants(self):
        adj = np.zeros((3, 7, 7), dtype=np.int64)
        inv = isa.graphs_to_invariants(adj)
        assert np.all(inv["Tr_M2"] == 0)
        assert np.all(inv["Tr_M4"] == 0)
        assert np.all(inv["Tr_M6"] == 0)


# ---------------------------------------------------------------------------
# Layer 3: Observable layer — match legacy chemistry-shim implementation
# ---------------------------------------------------------------------------

class TestAromaticityObservable:

    @pytest.mark.parametrize("name,smiles,expected_re", REFERENCE)
    @pytest.mark.parametrize("backend", BACKENDS)
    def test_re_matches_legacy(self, name, smiles, expected_re, backend):
        """ISA-backed RE matches legacy RE on the 5-molecule reference set."""
        re_legacy = chem.smiles_resonance_energy(smiles, calibration_kcal=12.0)
        re_isa = chem.smiles_resonance_energy_batch([smiles], backend=backend)[0]
        assert abs(re_legacy - re_isa) < 1e-9
        assert abs(re_isa - expected_re) < 1e-9

    @pytest.mark.parametrize("backend", BACKENDS)
    def test_dicoronylene_re(self, backend):
        """Dicoronylene = 2 separate aromatic systems × (1 K_7 each).
        ISA correctly sums RE across systems per molecule."""
        re_legacy = chem.smiles_resonance_energy(DICORONYLENE_SMILES)
        re_isa = chem.smiles_resonance_energy_batch(
            [DICORONYLENE_SMILES], backend=backend
        )[0]
        assert re_legacy == 400.0
        assert re_isa == 400.0

    @pytest.mark.parametrize("backend", BACKENDS)
    def test_batched_reference_set(self, backend):
        """Whole reference set in one batched call gives correct RE per molecule."""
        smiles_list = [s for _, s, _ in REFERENCE]
        re_isa = chem.smiles_resonance_energy_batch(smiles_list, backend=backend)
        for i, (_, _, expected) in enumerate(REFERENCE):
            assert abs(re_isa[i] - expected) < 1e-9


class TestK7Indicator:

    def test_coronene_fires(self):
        """Coronene's W_6 ring-adjacency should fire the K_7 indicator."""
        result = chem.smiles_to_aromaticity(CORONENE_SMILES)
        sys = result.aromatic_ring_systems[0]
        adj, n_used = chem.ring_system_to_so7_adjacency(sys)
        indicator = isa.k7_indicator(adj[None], n_vertices_used=np.array([n_used]))
        assert indicator[0] == 1

    @pytest.mark.parametrize("smiles", [
        "c1ccccc1",                          # benzene
        "c1ccc2ccccc2c1",                    # naphthalene
        "c1ccc2cc3ccccc3cc2c1",              # anthracene
        "c1cc2ccc3cccc4ccc(c1)c2c34",        # pyrene
    ])
    def test_smaller_pahs_do_not_fire(self, smiles):
        result = chem.smiles_to_aromaticity(smiles)
        sys = result.aromatic_ring_systems[0]
        adj, n_used = chem.ring_system_to_so7_adjacency(sys)
        indicator = isa.k7_indicator(adj[None], n_vertices_used=np.array([n_used]))
        assert indicator[0] == 0


class TestHopfPairObservable:

    @pytest.mark.parametrize("name,smiles,expected_pairs", [
        ("benzene",      "c1ccccc1",                                  3),
        ("naphthalene",  "c1ccc2ccccc2c1",                            5),
        ("anthracene",   "c1ccc2cc3ccccc3cc2c1",                      7),
        ("pyrene",       "c1cc2ccc3cccc4ccc(c1)c2c34",                8),
        ("coronene",     "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",   12),
    ])
    def test_pi_pair_count(self, name, smiles, expected_pairs):
        """Hopf-pair count via ISA agrees with chemistry-shim π pairs."""
        result = chem.smiles_to_aromaticity(smiles)
        for sys in result.aromatic_ring_systems:
            adj, n_used = chem.ring_system_to_so7_adjacency(sys)
            hopf = isa.hopf_pair_count(
                adj[None],
                np.array([sys.n_pi_electrons]),
            )
            assert hopf["n_pi_pairs"][0] == sys.n_pi_pairs
        # Verify against expected
        assert result.total_pi_pairs == expected_pairs


# ---------------------------------------------------------------------------
# Layer 4: Signature classification
# ---------------------------------------------------------------------------

class TestSignatureClassification:

    def test_coronene_classifies_as_W_6(self):
        result = chem.smiles_to_aromaticity(CORONENE_SMILES)
        sys = result.aromatic_ring_systems[0]
        adj, _ = chem.ring_system_to_so7_adjacency(sys)
        sigs = isa.classify_signature(adj[None])
        assert sigs[0] == "W_6"

    def test_benzene_classifies_as_empty(self):
        """Benzene has 1 ring → ring-adjacency has 0 edges → empty."""
        result = chem.smiles_to_aromaticity("c1ccccc1")
        sys = result.aromatic_ring_systems[0]
        adj, _ = chem.ring_system_to_so7_adjacency(sys)
        sigs = isa.classify_signature(adj[None])
        assert sigs[0] == "empty"

    def test_naphthalene_classifies_as_single_edge(self):
        result = chem.smiles_to_aromaticity("c1ccc2ccccc2c1")
        sys = result.aromatic_ring_systems[0]
        adj, _ = chem.ring_system_to_so7_adjacency(sys)
        sigs = isa.classify_signature(adj[None])
        assert sigs[0] == "single_edge"
