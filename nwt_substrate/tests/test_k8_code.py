"""Tests for nwt_substrate.algebra.codes.k8."""

import numpy as np

from nwt_substrate.algebra.codes.k8 import (
    pauli_string,
    steane_stabilizers,
    steane_logicals,
    k8_stabilizers,
    g_higgs,
    h_tensor7,
    commutes,
    verify,
)


def test_six_steane_generators_commute():
    mats = [pauli_string(s) for s in steane_stabilizers().values()]
    assert len(mats) == 6
    for i, a in enumerate(mats):
        for b in mats[i + 1:]:
            assert commutes(a, b)


def test_steane_logicals_commute_with_stabilizers_and_anticommute_pairwise():
    stabs = [pauli_string(s) for s in steane_stabilizers().values()]
    logs = steane_logicals()
    for L in (logs["X_L"], logs["Z_L"], logs["Y_L"]):
        for g in stabs:
            assert commutes(pauli_string(L), g)
    XL, ZL = pauli_string(logs["X_L"]), pauli_string(logs["Z_L"])
    assert np.allclose(XL @ ZL, -(ZL @ XL))  # logical X,Z anticommute


def test_gHiggs_is_Y8_and_a_valid_stabilizer():
    assert g_higgs() == "Y" * 8
    gH = pauli_string(g_higgs())
    for name, s in k8_stabilizers().items():
        if name == "gHiggs":
            continue
        assert commutes(gH, pauli_string(s))


def test_single_qubit_Y0_is_NOT_a_valid_stabilizer():
    # the rigor finding: single-qubit Y_0 anticommutes with a Hamming check
    gH_single = pauli_string("Y" + "I" * 6 + "Y")
    steane8 = [pauli_string(s + "I") for s in steane_stabilizers().values()]
    assert not all(commutes(gH_single, g) for g in steane8)


def test_H7_flips_gHiggs_and_swaps_XZ():
    H7 = h_tensor7()
    gH = pauli_string(g_higgs())
    assert np.allclose(H7 @ gH @ H7, -gH)
    stabs = k8_stabilizers()
    for r in range(3):
        assert np.allclose(
            H7 @ pauli_string(stabs[f"gX{r}"]) @ H7, pauli_string(stabs[f"gZ{r}"])
        )


def test_verify_all_pass():
    v = verify()
    assert all(v.values()), v
