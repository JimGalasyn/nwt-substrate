"""Tests for substrate Woodward-Hoffmann pericyclic selection rules.

The Hopf-pair parity rule (already used for ground-state aromaticity in
[[bond-orders-substrate-resolution]]) extends to transition-state cyclic
electron counts in pericyclic reactions. Same substrate rule, two domains:
ground-state aromaticity AND transition-state thermal allowance.
"""
from __future__ import annotations

import pytest

from nwt_substrate.chemistry import (
    PERICYCLIC_TS_ELECTRON_COUNT,
    PericyclicSelectionRule,
    electrocyclic_rotation_mode,
    reaction_selection_rule,
    selection_rule,
)


# ---------------------------------------------------------------------------
# Core selection rule (4n+2 / 4n parity)
# ---------------------------------------------------------------------------

class TestSelectionRule:
    def test_4n_plus_2_huckel_thermal(self):
        """6 electrons → 4n+2, Hückel TS, thermal suprafacial allowed."""
        r = selection_rule(6)
        assert r.parity == "odd"
        assert r.hopf_pair_count == 3
        assert r.thermal_allowed_suprafacial is True
        assert r.thermal_required_topology == "Huckel"
        assert r.photochemical_required_topology == "Mobius"

    def test_4n_mobius_thermal(self):
        """4 electrons → 4n, Möbius TS, thermal suprafacial forbidden."""
        r = selection_rule(4)
        assert r.parity == "even"
        assert r.hopf_pair_count == 2
        assert r.thermal_allowed_suprafacial is False
        assert r.thermal_required_topology == "Mobius"
        assert r.photochemical_required_topology == "Huckel"

    def test_2_electrons(self):
        """2 electrons: 1 pair, odd parity (4·0+2)."""
        r = selection_rule(2)
        assert r.parity == "odd"
        assert r.hopf_pair_count == 1
        assert r.thermal_allowed_suprafacial is True

    def test_8_electrons(self):
        """8 electrons: 4 pairs, even parity (4·2)."""
        r = selection_rule(8)
        assert r.parity == "even"
        assert r.hopf_pair_count == 4
        assert r.thermal_allowed_suprafacial is False

    def test_10_electrons(self):
        """10 electrons: 5 pairs, odd parity (4·2+2)."""
        r = selection_rule(10)
        assert r.parity == "odd"
        assert r.hopf_pair_count == 5
        assert r.thermal_allowed_suprafacial is True

    def test_zero_electrons(self):
        """0 electrons is even (degenerate case): treated as Möbius / no thermal allowance."""
        r = selection_rule(0)
        assert r.parity == "even"
        assert r.hopf_pair_count == 0

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            selection_rule(-2)

    def test_odd_raises(self):
        """Pericyclic TS electrons must be even (paired)."""
        with pytest.raises(ValueError):
            selection_rule(5)


# ---------------------------------------------------------------------------
# Canonical pericyclic reactions
# ---------------------------------------------------------------------------

class TestCanonicalReactions:
    def test_diels_alder_thermal_allowed(self):
        """[4+2] Diels-Alder: 6 TS electrons, thermal suprafacial-suprafacial."""
        r = reaction_selection_rule("[4+2]_cycloaddition_Diels_Alder")
        assert r.n_electrons == 6
        assert r.parity == "odd"
        assert r.thermal_allowed_suprafacial is True
        assert r.thermal_required_topology == "Huckel"

    def test_2_plus_2_thermal_forbidden(self):
        """[2+2] cycloaddition: 4 TS electrons, thermal forbidden (suprafacial)."""
        r = reaction_selection_rule("[2+2]_cycloaddition")
        assert r.n_electrons == 4
        assert r.parity == "even"
        assert r.thermal_allowed_suprafacial is False

    def test_4_plus_4_thermal_forbidden(self):
        """[4+4] cycloaddition: 8 TS electrons, thermal forbidden."""
        r = reaction_selection_rule("[4+4]_cycloaddition")
        assert r.n_electrons == 8
        assert r.parity == "even"
        assert r.thermal_allowed_suprafacial is False

    def test_6_plus_4_thermal_allowed(self):
        """[6+4] cycloaddition: 10 TS electrons, thermal allowed."""
        r = reaction_selection_rule("[6+4]_cycloaddition")
        assert r.n_electrons == 10
        assert r.parity == "odd"
        assert r.thermal_allowed_suprafacial is True

    def test_cope_claisen_thermal_allowed(self):
        """[3,3]-sigmatropic Cope/Claisen: 6 TS electrons, thermal allowed."""
        r = reaction_selection_rule("sigmatropic_3_3_Cope_Claisen")
        assert r.n_electrons == 6
        assert r.thermal_allowed_suprafacial is True

    def test_1_3_H_sigmatropic_thermal_forbidden(self):
        """[1,3]-H shift: 4 TS electrons, thermal forbidden suprafacial."""
        r = reaction_selection_rule("sigmatropic_1_3_H")
        assert r.n_electrons == 4
        assert r.thermal_allowed_suprafacial is False

    def test_1_5_H_sigmatropic_thermal_allowed(self):
        """[1,5]-H shift: 6 TS electrons, thermal allowed suprafacial."""
        r = reaction_selection_rule("sigmatropic_1_5_H")
        assert r.n_electrons == 6
        assert r.thermal_allowed_suprafacial is True

    def test_1_7_H_sigmatropic_thermal_forbidden(self):
        """[1,7]-H shift: 8 TS electrons, thermal forbidden."""
        r = reaction_selection_rule("sigmatropic_1_7_H")
        assert r.n_electrons == 8
        assert r.thermal_allowed_suprafacial is False

    def test_unknown_reaction_raises(self):
        with pytest.raises(KeyError):
            reaction_selection_rule("unknown_reaction_name")

    def test_all_canonical_keys_resolve(self):
        """Every entry in PERICYCLIC_TS_ELECTRON_COUNT should resolve to a valid rule."""
        for label, n_e in PERICYCLIC_TS_ELECTRON_COUNT.items():
            r = reaction_selection_rule(label)
            assert r.n_electrons == n_e


# ---------------------------------------------------------------------------
# Electrocyclic rotation modes (con/dis)
# ---------------------------------------------------------------------------

class TestElectrocyclicRotation:
    def test_butadiene_thermal_conrotatory(self):
        """4π electron (butadiene → cyclobutene) thermal: conrotatory."""
        mode = electrocyclic_rotation_mode(4, thermal=True)
        assert mode == "conrotatory"

    def test_hexatriene_thermal_disrotatory(self):
        """6π electron (hexatriene → cyclohexadiene) thermal: disrotatory."""
        mode = electrocyclic_rotation_mode(6, thermal=True)
        assert mode == "disrotatory"

    def test_octatetraene_thermal_conrotatory(self):
        """8π electron thermal: conrotatory (4n parity)."""
        mode = electrocyclic_rotation_mode(8, thermal=True)
        assert mode == "conrotatory"

    def test_butadiene_photochemical_disrotatory(self):
        """4π electron photochemical: disrotatory (complementary to thermal)."""
        mode = electrocyclic_rotation_mode(4, thermal=False)
        assert mode == "disrotatory"

    def test_hexatriene_photochemical_conrotatory(self):
        """6π electron photochemical: conrotatory (complementary to thermal)."""
        mode = electrocyclic_rotation_mode(6, thermal=False)
        assert mode == "conrotatory"

    def test_photochemical_always_complementary(self):
        """Photochemical mode is always opposite to thermal mode."""
        for n_e in [4, 6, 8, 10]:
            t_mode = electrocyclic_rotation_mode(n_e, thermal=True)
            p_mode = electrocyclic_rotation_mode(n_e, thermal=False)
            assert t_mode != p_mode


# ---------------------------------------------------------------------------
# Substrate connection (Hopf-pair parity consistency with ground-state aromaticity)
# ---------------------------------------------------------------------------

class TestSubstrateConsistency:
    def test_pericyclic_parity_matches_aromaticity_parity(self):
        """Same substrate Hopf-pair parity rule for ground-state aromaticity
        and pericyclic transition-state thermal allowance.

        Benzene (6 π electrons, aromatic, Hückel) ↔ Diels-Alder TS (6 e, Hückel)
        Cyclobutadiene (4 π electrons, anti-aromatic, Möbius) ↔ [2+2] TS (4 e, Möbius)
        """
        from nwt_substrate.chemistry import aromaticity_class

        # Ground-state benzene-style 6π: aromatic + odd Hopf-pair
        arom_6 = aromaticity_class("benzene")
        assert arom_6.classification == "aromatic"
        assert arom_6.parity == "odd"

        # Pericyclic TS with 6 electrons: same odd parity
        pc_6 = selection_rule(6)
        assert pc_6.parity == "odd"
        assert pc_6.parity == arom_6.parity

        # Ground-state 4π (anti-aromatic, e.g. cyclobutadiene-style): even parity
        arom_4 = aromaticity_class("cyclobutadiene")
        assert arom_4.parity == "even"

        # Pericyclic TS with 4 electrons: same even parity
        pc_4 = selection_rule(4)
        assert pc_4.parity == "even"
        assert pc_4.parity == arom_4.parity
