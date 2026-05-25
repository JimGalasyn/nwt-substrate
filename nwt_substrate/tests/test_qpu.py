"""M1 tests for nwt_substrate.qpu.

Covers: IR round-trip emission, the canonical-counts contract (the structural
fix for the Braket-vs-Qiskit endianness footgun), and end-to-end decode of the
electron walk on the local simulator.
"""
import importlib

import pytest

from nwt_substrate.qpu import spec, decode
from nwt_substrate.qpu.adapters import canonicalize_braket, canonicalize_ibm

# Electron (2,1) walk Pauli word + its prediction (from the compendium lookup).
ELECTRON_X = [0, 1, 0, 0, 0, 0, 0]
ELECTRON_Z = [0, 0, 1, 1, 1, 0, 1]
ELECTRON_PRED = ("F2", "E2", -1)   # (x_fano, z_fano, logical_z)

_have_braket = importlib.util.find_spec("braket") is not None
_have_qiskit = importlib.util.find_spec("qiskit") is not None


# ----------------------------------------------------------------- IR / spec
def test_destructive_builds_two_specs():
    base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
    z, x = spec.destructive_css(base, "electron")
    assert z.n_qubits == x.n_qubits == 7
    assert z.registers == {"c": 7} and x.registers == {"c": 7}
    assert z.measured_layout()["c"] == tuple(range(7))
    # The basis-change is carried by M.basis (emitted as H at to_qiskit/to_braket
    # time), not as explicit G("h") ops in the spec.
    assert {o.basis for o in z.ops if isinstance(o, spec.M)} == {"z"}
    assert {o.basis for o in x.ops if isinstance(o, spec.M)} == {"x"}


def test_ancilla_builds_13q_with_two_registers():
    base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
    a = spec.ancilla_syndrome(base, "electron")
    assert a.n_qubits == 13
    assert a.registers == {"c_syn": 6, "c_lz": 7}
    assert a.measured_layout()["c_lz"] == tuple(range(7))


@pytest.mark.skipif(not _have_qiskit, reason="qiskit not installed")
def test_to_qiskit_roundtrip():
    base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
    z, x = spec.destructive_css(base, "electron")
    qc = z.to_qiskit()
    assert qc.num_qubits == 7
    assert qc.num_clbits == 7
    # X-readout emits 7 extra H (one per data qubit) vs Z-readout
    n_h = lambda c: sum(1 for inst in c.data if inst.operation.name == "h")
    assert n_h(x.to_qiskit()) == n_h(qc) + 7


@pytest.mark.skipif(not _have_braket, reason="braket not installed")
def test_to_braket_roundtrip():
    base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
    z, _ = spec.destructive_css(base, "electron")
    c = z.to_braket()
    assert c.qubit_count == 7


# ------------------------------------------------- canonical-counts contract
def test_canonical_counts_vendor_equivalence():
    """A Braket-style raw key (leftmost=qubit0) and the matching Qiskit-style
    little-endian per-shot string MUST normalize to the same canonical Counts —
    this is the property that makes decode vendor-agnostic."""
    # Physical Z-readout outcome: qubit q -> bit value below (q=0..6)
    qubit_bits = [1, 0, 1, 1, 0, 0, 1]
    braket_key = "".join(str(b) for b in qubit_bits)          # leftmost = qubit 0
    qiskit_le = "".join(str(b) for b in reversed(qubit_bits))  # rightmost = cbit 0

    cb = canonicalize_braket({braket_key: 5}, {"c": tuple(range(7))})
    ci = canonicalize_ibm({"c": [qiskit_le] * 5}, ("c",))

    assert cb.single() == ci.single() == {braket_key: 5}
    # And both decode identically
    assert decode.fano(decode._bits(cb.single().popitem()[0])) == \
           decode.fano(decode._bits(ci.single().popitem()[0]))


def test_canonical_counts_multi_register_ibm():
    """Ancilla scheme: two registers, each independently reversed."""
    syn_le = "011010"      # cbit5..cbit0
    lz_le = "1000000"      # cbit6..cbit0
    c = canonicalize_ibm({"c_syn": [syn_le], "c_lz": [lz_le]}, ("c_syn", "c_lz"))
    assert c.register("c_syn") == {"010110": 1}   # reversed
    assert c.register("c_lz") == {"0000001": 1}


# ------------------------------------------------------- decode correctness
def test_decode_destructive_recovers_prediction_noiseless():
    """Construct the exact noiseless syndrome bitstrings for the electron walk
    and confirm decode yields (F2, E2, -1)."""
    # A single-shot canonical bitstring whose stabilizer parities give F2/E2 is
    # hard to hand-write; instead verify the fano/logical_z primitives directly
    # against the known electron syndrome from the lookup (x_syndrome F2 etc.).
    # Build the data-bit pattern from a codeword + the electron Pauli word.
    pytest.importorskip("braket")
    from nwt_substrate.qpu.adapters import SimulatorBackend
    backend = SimulatorBackend()
    base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
    z, x = spec.destructive_css(base, "electron")
    zc = backend.submit([backend.compile(z)], shots=400).result()[0]
    xc = backend.submit([backend.compile(x)], shots=400).result()[0]
    x_dist, z_dist = decode.destructive_dists(zc.single(), xc.single())
    v = decode.verdict_destructive("electron", x_dist, z_dist,
                                   ELECTRON_PRED[0], ELECTRON_PRED[1], ELECTRON_PRED[2])
    assert v.passed, v
    assert v.x_half.modal == "F2"
    assert v.z_half.modal == ("E2", -1)


def test_fano_identity_for_codeword():
    """All-zero data bits (a valid |0_L> codeword measurement) -> trivial
    syndrome I and logical_z +1."""
    assert decode.fano([0] * 7) == "I"
    assert decode.logical_z([0] * 7) == +1
