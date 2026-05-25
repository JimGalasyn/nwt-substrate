"""Local simulator backend for free dry-runs.

Uses the Braket LocalSimulator (the canonical-order, no-ancilla destructive path
is Braket-validated) and the same `canonicalize_braket` normalization the real
Braket adapter uses — so dry-run and hardware decode through identical code.
"""
from __future__ import annotations

from .base import Counts, Capabilities, Compiled
from .braket import canonicalize_braket


class _ImmediateJob:
    def __init__(self, results: list[Counts]):
        self._results = results

    @property
    def handles(self) -> list[str]:
        return []

    def result(self) -> list[Counts]:
        return self._results


class SimulatorBackend:
    def __init__(self, n_qubits: int = 32):
        from braket.devices import LocalSimulator
        self._dev = LocalSimulator()
        self.capabilities = Capabilities(
            name="local-simulator", sdk="simulator", n_qubits=n_qubits)

    def compile(self, spec, opt_level: int = 3) -> Compiled:
        return Compiled(spec.to_braket(), spec.measured_layout(), spec.name, spec.count_2q())

    def submit(self, compiled: list[Compiled], shots: int) -> _ImmediateJob:
        results = []
        for c in compiled:
            raw = dict(self._dev.run(c.sdk_circuit, shots=shots).result().measurement_counts)
            results.append(canonicalize_braket(raw, c.register_layout))
        return _ImmediateJob(results)

    def resume(self, handles, compiled):  # pragma: no cover
        raise NotImplementedError("simulator runs are immediate; nothing to resume")
