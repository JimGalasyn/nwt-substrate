"""Canonical types + Backend/Job protocols for the cross-architecture QPU interface.

See `analysis/qpu_interface_design.md` (null-worldtube-private) for the full
proposal. The single most important type here is `Counts` and its convention,
which is what makes vendor-specific bit-order footguns structurally impossible:
every adapter is responsible for normalizing its SDK's raw output into a
`Counts`, after which all decoding is vendor-agnostic.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Protocol, runtime_checkable


@dataclass(frozen=True)
class Counts:
    """Canonical measurement counts.

    `register_order` lists the classical registers in a fixed order. `data`
    maps a tuple of bitstrings (one per register, aligned to `register_order`)
    to a shot count.

    CONVENTION (the contract every adapter normalizes to): within each
    register's bitstring, character index ``i`` is qubit/cbit ``i`` — index 0
    is the LEFTMOST character. (Braket is already this order; the IBM adapter
    reverses Qiskit's little-endian strings.)
    """

    register_order: tuple[str, ...]
    data: dict[tuple[str, ...], int]

    def total(self) -> int:
        return sum(self.data.values())

    def register(self, name: str) -> dict[str, int]:
        """Marginal counts for one register: {bitstring(index0-left): n}."""
        i = self.register_order.index(name)
        out: Counter[str] = Counter()
        for key, n in self.data.items():
            out[key[i]] += n
        return dict(out)

    def single(self) -> dict[str, int]:
        """Convenience for single-register results."""
        if len(self.register_order) != 1:
            raise ValueError(
                f"single() requires exactly one register, have {self.register_order}")
        return self.register(self.register_order[0])

    @staticmethod
    def from_canonical(register_order, data) -> "Counts":
        return Counts(tuple(register_order), {tuple(k): int(v) for k, v in data.items()})


@dataclass(frozen=True)
class Capabilities:
    """Static-ish description of a backend (subset needed for M1; windows +
    cost guardrails are fleshed out in capabilities.py during M2)."""

    name: str
    sdk: str                       # "ibm" | "braket" | "simulator"
    n_qubits: int
    native_2q: str = "cz"
    coupling_map: tuple[tuple[int, int], ...] = ()
    per_shot_usd: float = 0.0
    task_fee_usd: float = 0.0
    median_2q_fidelity: float = 0.99


@dataclass
class Compiled:
    """An SDK-specific compiled circuit plus the register layout needed to
    normalize its raw counts. `register_layout` maps each register name to the
    qubit indices it measured, in position order."""

    sdk_circuit: Any
    register_layout: dict[str, tuple[int, ...]]
    name: str = "circ"
    n_2q: int = 0
    depth: int = 0


@runtime_checkable
class Job(Protocol):
    def result(self) -> list[Counts]:
        """Block until done; return one canonical Counts per submitted circuit."""
        ...

    @property
    def handles(self) -> list[str]:
        """Per-circuit resume handles (Braket ARNs or IBM job-ids)."""
        ...


@runtime_checkable
class Backend(Protocol):
    capabilities: Capabilities

    def compile(self, spec, opt_level: int = 3) -> Compiled:
        ...

    def submit(self, compiled: list[Compiled], shots: int) -> Job:
        ...

    def resume(self, handles: list[str], compiled: list[Compiled]) -> Job:
        ...
