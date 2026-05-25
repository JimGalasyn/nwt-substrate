"""AWS Braket adapter.

Braket `measurement_counts` keys are already in canonical order (character q ==
qubit q, leftmost = qubit 0), so normalization is a pure slice into named
registers using the circuit's measured layout. `canonicalize_braket` is SDK-free
and unit-tested against a fixture.
"""
from __future__ import annotations

from collections import Counter

from .base import Counts, Capabilities, Compiled

# (region, per_shot_usd, task_fee_usd) — verified 2026-05-23 device inventory
_BRAKET_INFO = {
    "rigetti/Cepheus-1-108Q": ("us-west-1", 0.000425, 0.30),
    "iqm/Garnet": ("eu-north-1", 0.00145, 0.30),
    "iqm/Emerald": ("eu-north-1", 0.0016, 0.30),
    "aqt/Ibex-Q1": ("eu-north-1", 0.0235, 0.30),
}


def canonicalize_braket(raw_counts: dict[str, int], register_layout: dict) -> Counts:
    """raw_counts: {qubit-ordered bitstring (leftmost=qubit0): n}.
    Slice into canonical per-register bitstrings via register_layout."""
    register_order = tuple(register_layout.keys())
    data: Counter = Counter()
    for bs, n in raw_counts.items():
        key = tuple("".join(bs[q] for q in register_layout[r]) for r in register_order)
        data[key] += n
    return Counts(register_order, dict(data))


class _BraketJob:
    def __init__(self, tasks, layouts):
        self._tasks = tasks
        self._layouts = layouts

    @property
    def handles(self) -> list[str]:
        return [t.id for t in self._tasks]

    def result(self) -> list[Counts]:
        out = []
        for t, layout in zip(self._tasks, self._layouts):
            raw = dict(t.result().measurement_counts)
            out.append(canonicalize_braket(raw, layout))
        return out


class BraketBackend:
    def __init__(self, backend: str = "iqm/Garnet", region: str | None = None):
        from braket.aws import AwsDevice
        if backend.startswith("arn:aws:braket:"):
            arn = backend
            region = region or arn.split(":")[3]
            per_shot, fee = None, 0.30
        else:
            exp_region, per_shot, fee = _BRAKET_INFO.get(backend, (region, None, 0.30))
            region = region or exp_region
            arn = f"arn:aws:braket:{region}::device/qpu/{backend}"
        self._device = AwsDevice(arn)
        self.arn = arn
        self.capabilities = Capabilities(
            name=backend, sdk="braket", n_qubits=getattr(self._device, "num_qubits", 0) or 0,
            native_2q="cz", per_shot_usd=per_shot or 0.0, task_fee_usd=fee)

    def compile(self, spec, opt_level: int = 3) -> Compiled:
        return Compiled(spec.to_braket(), spec.measured_layout(), spec.name, spec.count_2q())

    def submit(self, compiled: list[Compiled], shots: int) -> _BraketJob:
        tasks = [self._device.run(c.sdk_circuit, shots=shots) for c in compiled]
        return _BraketJob(tasks, [c.register_layout for c in compiled])

    def resume(self, handles: list[str], compiled: list[Compiled]) -> _BraketJob:
        from braket.aws import AwsQuantumTask
        tasks = [AwsQuantumTask(arn=h) for h in handles]
        return _BraketJob(tasks, [c.register_layout for c in compiled])
