"""IBM Qiskit Runtime adapter.

Qiskit `SamplerV2` returns per-register `BitArray`s whose `get_bitstrings()` are
little-endian (rightmost char = cbit 0). `canonicalize_ibm` reverses each
register to the canonical convention (index 0 leftmost) and zips them per shot.
The reversal lives in exactly this one place — the structural fix for the
endianness footgun. `canonicalize_ibm` is SDK-free and unit-tested.

One asymmetry vs Braket: a SamplerV2 job batches all circuits into a single
job, so an IBM `Job.handles` is the single batch job-id (length 1), whereas
Braket exposes one ARN per circuit.
"""
from __future__ import annotations

from collections import Counter

from .base import Counts, Capabilities, Compiled


def canonicalize_ibm(per_register_bitstrings: dict[str, list[str]], register_order) -> Counts:
    """per_register_bitstrings: {regname: [little-endian bitstring, one per shot]}
    as from BitArray.get_bitstrings(). Returns canonical Counts."""
    register_order = tuple(register_order)
    n_shots = len(next(iter(per_register_bitstrings.values())))
    data: Counter = Counter()
    for s in range(n_shots):
        key = tuple(per_register_bitstrings[r][s][::-1] for r in register_order)
        data[key] += 1
    return Counts(register_order, dict(data))


class _IBMJob:
    def __init__(self, runtime_job, register_orders):
        self._job = runtime_job
        self._register_orders = register_orders  # one per circuit

    @property
    def handles(self) -> list[str]:
        return [self._job.job_id()]

    def result(self) -> list[Counts]:
        res = self._job.result()
        out = []
        for i, reg_order in enumerate(self._register_orders):
            data = res[i].data
            per_reg = {r: getattr(data, r).get_bitstrings() for r in reg_order}
            out.append(canonicalize_ibm(per_reg, reg_order))
        return out


class IBMBackend:
    def __init__(self, backend_name: str, service=None, token_path: str | None = None):
        from qiskit_ibm_runtime import QiskitRuntimeService
        if service is None:
            tok = None
            if token_path:
                from pathlib import Path
                tok = Path(token_path).read_text().strip()
            else:
                from pathlib import Path
                p = Path.home() / ".ibm_quantum_api_key"
                tok = p.read_text().strip() if p.exists() else None
            last = None
            for ch in ("ibm_cloud", "ibm_quantum"):
                try:
                    service = QiskitRuntimeService(channel=ch, token=tok)
                    break
                except Exception as e:  # pragma: no cover - network dependent
                    last = e
            if service is None:
                raise RuntimeError(f"could not connect to IBM: {last}")
        self._service = service
        self._backend = service.backend(backend_name)
        self.capabilities = Capabilities(
            name=backend_name, sdk="ibm", n_qubits=self._backend.num_qubits,
            native_2q="cz", per_shot_usd=0.0, task_fee_usd=0.0)

    def compile(self, spec, opt_level: int = 3) -> Compiled:
        from qiskit import transpile
        qc = transpile(spec.to_qiskit(), backend=self._backend, optimization_level=opt_level)
        return Compiled(qc, spec.measured_layout(), spec.name,
                        qc.num_nonlocal_gates(), qc.depth())

    def submit(self, compiled: list[Compiled], shots: int) -> _IBMJob:
        from qiskit_ibm_runtime import SamplerV2 as Sampler
        sampler = Sampler(mode=self._backend)
        job = sampler.run([c.sdk_circuit for c in compiled], shots=shots)
        return _IBMJob(job, [tuple(c.register_layout.keys()) for c in compiled])

    def resume(self, handles: list[str], compiled: list[Compiled]) -> _IBMJob:
        job = self._service.job(handles[0])  # single batch job-id
        return _IBMJob(job, [tuple(c.register_layout.keys()) for c in compiled])
