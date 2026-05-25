"""
K₈ [[8,1,1]] substrate stabilizer code and the g^Higgs Higgs stabilizer.

The K₇ substrate is the Steane [[7,1,3]] CSS code (qubits q0..q6 ↔ 7 K₇
vertices). The K₈ extension adds an 8th qubit q7 (the Y₇ octonion direction)
and one extra, non-CSS, Y-type stabilizer

    g^Higgs = Y_0 · Y_7 ,

giving the [[8,1,1]] code (8 physical, 7 stabilizers, 1 logical).
⟨g^Higgs⟩ = ±1 is the Higgs-vacuum direction / matter(+)–antimatter(−)
commitment of the cosmogenic bridge.

★ RIGOR FINDING (the reason this module exists as a verified object):
  Y_0 in g^Higgs must be the Steane *logical* Y_L (the weight-7 transversal
  Y^⊗7), **not** a single-qubit Pauli-Y. A single-qubit Y on a Steane data
  qubit anticommutes with the Hamming X/Z checks supported there, so it is
  NOT a valid [[8,1,1]] stabilizer. With Y_0 = Y_L, g^Higgs = Y^⊗7 ⊗ Y_7 =
  Y^⊗8 commutes with all six Steane generators and is independent of them.

The transversal H⊗7 (Hadamard on q0..q6, identity on q7) implements Steane
Check 1b (swaps X- and Z-stabilizers) and sends g^Higgs → −g^Higgs — the
matter↔antimatter Higgs-vacuum flip.

Two physical channels share the codespace (Stage-7 K₇⊗K₈ partition):
  * K₇-Steane-active  (qubits {1..6} support)  → baryonic matter;
  * K₈-extension-only (qubits {0,7} support)   → cold dark matter
    (sterile-ν-like: gravitates, no SM gauge, Y_0Y_7 W-coupling only).

Provides:
  - pauli_string(s)        : tensor-product Pauli operator from a string
  - steane_stabilizers()   : the 6 Steane [[7,1,3]] generators (7-qubit strings)
  - steane_logicals()      : X_L, Z_L, Y_L (transversal, weight-7)
  - k8_stabilizers()       : the 7 [[8,1,1]] generators (8-qubit strings)
  - g_higgs()              : the Higgs stabilizer string (= "Y"*8)
  - h_tensor7()            : the transversal H⊗7 operator (256×256)
  - commutes(A, B)         : Pauli-operator commutation test
  - verify()               : asserts the code's defining properties (bool dict)

All operations verified by tests in tests/test_k8_code.py.
"""
from __future__ import annotations

import numpy as np

from nwt_substrate.isa.constants import N_VERTICES_K7, N_VERTICES_K8

# single-qubit Paulis
_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_H = (_X + _Z) / np.sqrt(2)
_PAULI = {"I": _I, "X": _X, "Y": _Y, "Z": _Z}

# Steane = Hamming [7,4,3] parity-check supports (columns = binary 1..7).
_STEANE_SUPPORTS = ({3, 4, 5, 6}, {1, 2, 5, 6}, {0, 2, 4, 6})

__all__ = [
    "pauli_string",
    "steane_stabilizers",
    "steane_logicals",
    "k8_stabilizers",
    "g_higgs",
    "h_tensor7",
    "commutes",
    "verify",
]


def pauli_string(s: str) -> np.ndarray:
    """Tensor product of single-qubit Paulis; leftmost char = qubit 0."""
    op = np.array([[1.0 + 0j]])
    for ch in s:
        op = np.kron(op, _PAULI[ch])
    return op


def steane_stabilizers() -> dict[str, str]:
    """6 Steane [[7,1,3]] generators as 7-character Pauli strings (q0..q6)."""
    gens: dict[str, str] = {}
    for r, supp in enumerate(_STEANE_SUPPORTS):
        gens[f"gX{r}"] = "".join("X" if q in supp else "I" for q in range(N_VERTICES_K7))
        gens[f"gZ{r}"] = "".join("Z" if q in supp else "I" for q in range(N_VERTICES_K7))
    return gens


def steane_logicals() -> dict[str, str]:
    """Transversal (weight-7) Steane logical operators."""
    n = N_VERTICES_K7
    return {"X_L": "X" * n, "Z_L": "Z" * n, "Y_L": "Y" * n}


def g_higgs() -> str:
    """g^Higgs = Y_0·Y_7 with Y_0 = logical Y_L  ⟹  Y^⊗8 (8-qubit string)."""
    return "Y" * N_VERTICES_K8


def k8_stabilizers() -> dict[str, str]:
    """The 7 [[8,1,1]] generators (8-qubit strings): 6 Steane (q7=I) + g^Higgs."""
    gens = {name: s + "I" for name, s in steane_stabilizers().items()}
    gens["gHiggs"] = g_higgs()
    return gens


def h_tensor7() -> np.ndarray:
    """Transversal H⊗7 = H on q0..q6, identity on q7 (256×256)."""
    op = np.array([[1.0 + 0j]])
    for _ in range(N_VERTICES_K7):
        op = np.kron(op, _H)
    return np.kron(op, _I)


def commutes(a, b, tol: float = 1e-9) -> bool:
    """True iff operators (matrices, or Pauli strings) commute."""
    A = pauli_string(a) if isinstance(a, str) else a
    B = pauli_string(b) if isinstance(b, str) else b
    return np.allclose(A @ B - B @ A, 0.0, atol=tol)


def verify() -> dict[str, bool]:
    """Assert the defining properties of the [[8,1,1]] code + g^Higgs."""
    stabs = k8_stabilizers()
    mats = {name: pauli_string(s) for name, s in stabs.items()}
    names = list(mats)

    # all stabilizers mutually commute
    all_commute = all(
        commutes(mats[a], mats[b])
        for i, a in enumerate(names) for b in names[i + 1:]
    )

    # the single-qubit reading of Y_0 is NOT a valid stabilizer (rigor finding)
    gH_single = pauli_string("Y" + "I" * (N_VERTICES_K7 - 1) + "Y")
    steane8 = {n: pauli_string(s + "I") for n, s in steane_stabilizers().items()}
    single_qubit_valid = all(commutes(gH_single, m) for m in steane8.values())

    # H⊗7 flips g^Higgs sign
    H7 = h_tensor7()
    gH = pauli_string(g_higgs())
    higgs_sign_flip = np.allclose(H7 @ gH @ H7, -gH)

    # H⊗7 swaps X↔Z Steane stabilizers
    xz_swap = all(
        np.allclose(H7 @ pauli_string(stabs[f"gX{r}"]) @ H7, pauli_string(stabs[f"gZ{r}"]))
        for r in range(3)
    )

    return {
        "stabilizers_commute": all_commute,
        "gHiggs_valid_as_logical_Y_L": all_commute,  # gHiggs ∈ stabs, checked above
        "single_qubit_Y0_invalid": not single_qubit_valid,  # expect True
        "H7_flips_gHiggs": higgs_sign_flip,
        "H7_swaps_XZ": xz_swap,
    }
