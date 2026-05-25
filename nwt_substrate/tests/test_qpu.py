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


# --------------------------------------------------------------- M2: preflight
# IQM Garnet's real 20-qubit grid (0-indexed).
_GARNET_EDGES = tuple(sorted({
    (min(a - 1, b - 1), max(a - 1, b - 1)) for a, b in [
        (1, 2), (1, 4), (2, 5), (3, 4), (3, 8), (4, 5), (4, 9), (5, 6), (5, 10),
        (6, 7), (6, 11), (7, 12), (8, 9), (8, 13), (9, 10), (9, 14), (10, 11),
        (10, 15), (11, 12), (11, 16), (12, 17), (13, 14), (14, 15), (14, 18),
        (15, 16), (15, 19), (16, 17), (16, 20), (18, 19), (19, 20)]}))


def _garnet_caps():
    from nwt_substrate.qpu.adapters.base import Capabilities
    return Capabilities(name="iqm/Garnet", sdk="braket", n_qubits=20, native_2q="cz",
                        coupling_map=_GARNET_EDGES, per_shot_usd=0.00145, task_fee_usd=0.30,
                        median_2q_fidelity=0.98)  # effective all-in (see capabilities.py)


@pytest.mark.skipif(not _have_qiskit, reason="qiskit not installed")
def test_preflight_flags_ancilla_but_passes_destructive_on_garnet():
    """Reproduces the 2026-05-25 diagnosis as a unit test: the 13-qubit ancilla
    circuit blows Garnet's budget; the destructive variant fits."""
    from nwt_substrate.qpu import capabilities as cap
    base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
    anc = spec.ancilla_syndrome(base, "electron")
    z, _ = spec.destructive_css(base, "electron")

    pf_anc = cap.preflight(anc, _garnet_caps(), shots=100)
    pf_des = cap.preflight(z, _garnet_caps(), shots=100)

    assert pf_anc.n_2q > 40 and not pf_anc.fits_budget       # ~58 CZ -> blows budget
    assert pf_des.n_2q < 25 and pf_des.fits_budget           # ~15 CZ -> fits
    assert pf_des.n_2q < pf_anc.n_2q


def test_window_status_aqt():
    from datetime import datetime, timezone
    from nwt_substrate.qpu.capabilities import window_status
    from nwt_substrate.qpu.adapters.base import Window
    aqt = (Window("Monday", "11:00:00", "15:00:00"),
           Window("Tuesday", "08:00:00", "15:00:00"),
           Window("Wednesday", "08:00:00", "15:00:00"),
           Window("Friday", "08:00:00", "15:00:00"))
    # 2026-05-25 is a Monday; 14:00 UTC is inside the Mon window
    in_win, nxt = window_status(aqt, datetime(2026, 5, 25, 14, 0, tzinfo=timezone.utc))
    assert in_win and nxt is None
    # 2026-05-23 is a Saturday -> closed; next opening is Monday 11:00
    in_win, nxt = window_status(aqt, datetime(2026, 5, 23, 12, 0, tzinfo=timezone.utc))
    assert not in_win and nxt.startswith("2026-05-25T11:00")


def test_window_status_none_is_24_7():
    from nwt_substrate.qpu.capabilities import window_status
    assert window_status(None) == (True, None)


def test_default_scheme_is_destructive_everywhere():
    """M4 policy: destructive CSS readout is the default cross-vendor scheme."""
    from nwt_substrate.qpu.runner import choose_scheme
    from nwt_substrate.qpu.adapters.base import Capabilities
    for sdk in ("ibm", "braket", "simulator"):
        assert choose_scheme(Capabilities("x", sdk, 20)) == "destructive"
    assert choose_scheme(None) == "destructive"


def test_window_status_weekdays_aggregate():
    """IQM Garnet reports executionDay='Weekdays' (Mon-Fri) rather than per-day."""
    from datetime import datetime, timezone
    from nwt_substrate.qpu.capabilities import window_status
    from nwt_substrate.qpu.adapters.base import Window
    garnet = (Window("Weekdays", "03:15:00", "15:30:00"),)
    # Monday 2026-05-25 14:00 UTC -> open
    assert window_status(garnet, datetime(2026, 5, 25, 14, 0, tzinfo=timezone.utc))[0]
    # Saturday 2026-05-23 -> closed, next opening Monday 03:15
    in_win, nxt = window_status(garnet, datetime(2026, 5, 23, 12, 0, tzinfo=timezone.utc))
    assert not in_win and nxt.startswith("2026-05-25T03:15")
