"""
Heron Experiment 11 control channels.

The drift-validation analysis (``analysis/nwt_exp11_drift_validation.py``,
2026-05-06) shows that the bare A/B/C linear correction
``sigma_v = (B - A) - (C - A) / 2`` leaves residual fractional drifts
of ~31% on T1, ~39% on T2, and ~93% on readout error per Heron qubit
at the 95th percentile -- well above the 1% per-vertex sensitivity
the protocol claims.

To attribute any observed sigma_v signal to the substrate axis (and
not to thermal / readout / cosmic-ray drift), we run two control
channels alongside the K_7 measurement at each of the three sidereal
slots A, B, C:

  * **Readout-floor control**: prep |0>^7, measure all 7 qubits in Z.
    This is dominated by readout error (prep0 -> meas0 success rate).
    It must NOT show any sidereal modulation -- if it does, the
    apparent K_7 signal is a readout drift artefact.

  * **T1-floor control**: prep |1>^7, idle for the K_7 circuit's
    duration, measure all 7 in Z.  This is dominated by T1 decay.
    It must NOT show any sidereal modulation -- if it does, the
    apparent K_7 signal is a T1 drift artefact.

The substrate-monism prediction is specifically about the K_7 graph
state's PSL(2,7) symmetry coupling to the substrate axis.  Neither
control state has that symmetry: |0>^7 is a product state with full
permutation symmetry S_7, and idled |1>^7 is identical.  Under the
substrate hypothesis, both controls are predicted to be flat across
A/B/C.

Pass / fail logic (per slot):

  * K_7 signal:   drift-corrected sigma_v
  * RO control:   drift-corrected delta_<Z>_RO (per qubit)
  * T1 control:   drift-corrected delta_<Z>_T1 (per qubit)
  * Verdict:
      - If RO control and T1 control are both consistent with zero
        AND K_7 is non-zero: candidate substrate signal.
      - If any control shows a non-zero drift-corrected delta:
        contaminated; re-run or extend folding.

The controls are dirt cheap (no CZ gates) and run in a few hundred
ms per slot -- well below the K_7 cost.
"""

from __future__ import annotations

from typing import Optional

try:
    from qiskit import QuantumCircuit, ClassicalRegister
    HAS_QISKIT = True
except ImportError:
    HAS_QISKIT = False
    QuantumCircuit = None
    ClassicalRegister = None


def readout_control_circuit() -> "QuantumCircuit":
    """Prep |0>^7, measure all 7 qubits in the Z basis.

    Expected ideal outcome: all 7 measurements return 0.
    Observed deviation: readout error (prob_meas1_prep0).
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed; pip install qiskit")
    qc = QuantumCircuit(7, 7)
    # No prep needed: |0>^7 is the default reset state.
    for q in range(7):
        qc.measure(q, q)
    return qc


def t1_idle_control_circuit(
    idle_duration_ns: int = 8000,
    use_delay_instruction: bool = True,
) -> "QuantumCircuit":
    """Prep |1>^7, idle for ``idle_duration_ns`` ns, measure all in Z.

    Expected ideal outcome: all 7 measurements return 1.
    Observed deviation: T1 decay during the idle, plus readout
    error.

    The default 8000 ns matches the typical full K_7 prep +
    stabilizer measurement circuit time on Heron R2 (depth ~50,
    two-qubit gate ~100 ns; depth ~80 with measurement adds another
    ~3 us of measurement time).  Adjust per-device if the K_7
    circuit time differs.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed; pip install qiskit")
    qc = QuantumCircuit(7, 7)
    for q in range(7):
        qc.x(q)                            # |0> -> |1>
    if use_delay_instruction and idle_duration_ns > 0:
        qc.barrier()
        for q in range(7):
            qc.delay(idle_duration_ns, q, unit="ns")
        qc.barrier()
    for q in range(7):
        qc.measure(q, q)
    return qc


def parse_zbasis_expectations(counts: dict) -> list[float]:
    """From a counts dict on 7 classical bits, compute <Z_q> for each
    qubit.  Convention: bit 0 in the rightmost classical position."""
    expectations = [0.0] * 7
    n_total = 0
    for bitstr, count in counts.items():
        bs = bitstr.replace(" ", "")
        n_total += count
        # qiskit prints bitstrings right-to-left; cbit q at position -(q+1).
        for q in range(7):
            outcome = bs[-(q + 1)]
            # |0> -> +1 in Z; |1> -> -1 in Z
            expectations[q] += count if outcome == "0" else -count
    return [e / n_total for e in expectations] if n_total > 0 else [0.0] * 7


def control_drift_corrected(
    A_expectations: list[float],
    B_expectations: list[float],
    C_expectations: list[float],
) -> list[float]:
    """Apply the same linear A/B/C correction as the K_7 channel:
    delta = (B - A) - (C - A) / 2 per qubit."""
    return [
        (b - a) - (c - a) / 2.0
        for a, b, c in zip(A_expectations, B_expectations, C_expectations)
    ]


__all__ = [
    "readout_control_circuit",
    "t1_idle_control_circuit",
    "parse_zbasis_expectations",
    "control_drift_corrected",
]
