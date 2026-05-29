"""nwt_substrate.sensitivity — ISA structural-integer sensitivity sweep.

Pasquale Kaboth's question (d12rg, 2026-05-29): *how sensitive are the
predictions to small modifications of the structural identities?* This module
answers it operationally.

For each leaf integer constant ``NAME: int = V`` in :mod:`nwt_substrate.isa.constants`,
it perturbs the value by ±1, re-imports the whole library in a fresh
``python -O`` subprocess (``-O`` strips the import-time structural-identity
``assert``\\s so the perturbed — deliberately inconsistent — ISA still loads),
recomputes every ``benchmark_*`` prediction, and records which predicted values
move. The output is a dependency / robustness map:

  * which predictions are **forced** by a given structural integer (they move
    when it moves), versus
  * which integers are **inert** (perturbing them changes no prediction —
    either pure structural bookkeeping, or the prediction re-encodes the
    integer as a bare literal rather than importing the ISA constant).

That distinction is the look-elsewhere control: an integer that moves many
independent observables at once cannot be a per-observable tuning knob.

Mechanism note
--------------
The sweep works by briefly rewriting the *source* of ``isa/constants.py`` (so
the perturbation propagates through every derived constant and every
``from ..isa import X`` binding in a fresh interpreter), then restoring it. The
original text is held in memory and rewritten in a ``finally`` block, so the
file is always restored even on error/interrupt. Because it writes to the
installed source, this is a **development / editable-install tool** — it
refuses to run if the constants file is not writable.

Usage
-----
::

    python -m nwt_substrate.sensitivity            # full sweep, pretty table

    import nwt_substrate.sensitivity as sens
    report = sens.integer_sweep()
    print(report)
    report.movers("DIM_S_SPIN7")                   # benchmarks moved by 8±1
    report.inert_integers                          # integers that move nothing
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field

from .isa import constants as _constants

CONSTANTS_PATH = _constants.__file__

# leaf integer literal:  NAME: int = 7   (optional trailing comment)
_LEAF_RE = re.compile(r"^([A-Z][A-Z0-9_]*): int = (\d+)\s*(?:#.*)?$", re.M)


# ---------------------------------------------------------------------------
# Child entry point: run every benchmark and emit {func_name: predicted_value}
# ---------------------------------------------------------------------------

def _run_all_predictions() -> dict[str, str]:
    """Run every ``benchmark_*`` callable, returning {func_name: substrate_value}.

    Each benchmark is isolated in its own ``try`` so one failure under a
    perturbed integer does not abort the rest. Invoked in a fresh subprocess
    by the sweep; also callable directly for the un-perturbed baseline.
    """
    from nwt_substrate.benchmarks import compute_speed as cs

    out: dict[str, str] = {}
    for name in sorted(dir(cs)):
        if not name.startswith("benchmark_"):
            continue
        fn = getattr(cs, name)
        if not callable(fn):
            continue
        try:
            result = fn()
            out[name] = str(getattr(result, "substrate_value", result))
        except Exception as exc:  # perturbation may make a benchmark raise
            out[name] = f"__ERROR__: {type(exc).__name__}: {exc}"
    return out


def _child_cmd() -> list[str]:
    return [
        sys.executable, "-O", "-B", "-c",
        "import json, nwt_substrate.sensitivity as s; "
        "print(json.dumps(s._run_all_predictions()))",
    ]


def _run_child(env_extra: dict | None = None, timeout: int = 300) -> dict[str, str]:
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    proc = subprocess.run(
        _child_cmd(), capture_output=True, text=True, timeout=timeout, env=env,
    )
    if proc.returncode != 0:
        tail = (proc.stderr.strip().splitlines() or ["<no stderr>"])[-1]
        return {"__IMPORT_FAILED__": tail}
    try:
        return json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception as exc:
        return {"__PARSE_FAILED__": f"{exc}"}


# ---------------------------------------------------------------------------
# Source patching (always restored)
# ---------------------------------------------------------------------------

def leaf_integers() -> list[tuple[str, int]]:
    """The ``NAME: int = V`` leaf integer constants available to perturb."""
    src = open(CONSTANTS_PATH, encoding="utf-8").read()
    return [(m.group(1), int(m.group(2))) for m in _LEAF_RE.finditer(src)]


def _patched_source(orig: str, name: str, new_value: int) -> str:
    patched, n = re.subn(
        rf"^{re.escape(name)}: int = \d+",
        f"{name}: int = {new_value}", orig, count=1, flags=re.M,
    )
    if n != 1:
        raise ValueError(f"could not uniquely patch {name!r} in {CONSTANTS_PATH}")
    return patched


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

@dataclass
class SensitivityReport:
    """Result of an :func:`integer_sweep`."""
    benchmarks: list[str]                          # all benchmark func names
    baseline: dict[str, str]                       # name -> predicted value
    # integer -> {delta -> {"moved": [...], "status": "ok"|...}}
    per_integer: dict[str, dict[int, dict]] = field(default_factory=dict)

    def movers(self, integer: str) -> list[str]:
        """Benchmarks whose predicted value moves under ``integer`` ±1."""
        moved: set[str] = set()
        for entry in self.per_integer.get(integer, {}).values():
            if entry.get("status") == "ok":
                moved |= set(entry["moved"])
        return sorted(moved)

    @property
    def inert_integers(self) -> list[str]:
        """Integers whose ±1 perturbation moves no prediction at all."""
        return sorted(i for i in self.per_integer if not self.movers(i))

    @property
    def never_moved(self) -> list[str]:
        """Benchmarks not moved by *any* integer (literal-encoded / α-anchored)."""
        moved: set[str] = set()
        for integer in self.per_integer:
            moved |= set(self.movers(integer))
        return sorted(set(self.benchmarks) - moved)

    def coupling(self) -> dict[str, int]:
        """benchmark -> number of distinct integers that move it."""
        counts = {b: 0 for b in self.benchmarks}
        for integer in self.per_integer:
            for b in self.movers(integer):
                counts[b] += 1
        return counts

    # ------------------------------------------------------------------
    # Structural-criticality layer (M. Wende, d12rg v0.3.1 thread):
    # beyond *sensitivity* (how strongly an observable responds), ask which
    # structural integers carry the most *load* (organise many observables)
    # and which observables move together (correlated structural dependence).
    # ------------------------------------------------------------------

    def structural_load(self) -> list[tuple[str, int]]:
        """Integers ranked by structural load = number of benchmarks each moves
        (descending, ties broken by name).  The top entries are the framework's
        load-bearing integers; a load of 0 is inert.  A high load is *critical*
        (a single structural quantity organising many observables) rather than
        merely sensitive."""
        load = [(i, len(self.movers(i))) for i in self.per_integer]
        return sorted(load, key=lambda kv: (-kv[1], kv[0]))

    def comovement(self, min_shared: int = 2) -> list[tuple[tuple[str, str], int]]:
        """Benchmark pairs that move *together*, ranked by how many integers move
        BOTH.  High co-movement reveals a correlated structural dependency — two
        observables governed by the same part of the architecture — even when the
        individual numerical shifts are modest.  Returns ``[((a, b), n), ...]``."""
        pair_counts: dict[tuple[str, str], int] = {}
        for integer in self.per_integer:
            movers = self.movers(integer)
            for i in range(len(movers)):
                for j in range(i + 1, len(movers)):
                    key = tuple(sorted((movers[i], movers[j])))
                    pair_counts[key] = pair_counts.get(key, 0) + 1
        pairs = [(k, n) for k, n in pair_counts.items() if n >= min_shared]
        return sorted(pairs, key=lambda kv: (-kv[1], kv[0]))

    def criticality_summary(self, top: int = 8) -> str:
        """Structural-criticality view: load ranking + correlated benchmark
        clusters, layered on top of the raw sensitivity counts."""
        lines = ["Structural-load ranking (integers by # observables they move):", ""]
        for i, n in self.structural_load():
            if n > 0:
                lines.append(f"  {i:30} moves {n:>2} benchmark(s)")
        lines += ["", "Correlated benchmark pairs (move together under ≥2 integers):", ""]
        comov = self.comovement(min_shared=2)[:top]
        if comov:
            for (a, b), n in comov:
                lines.append(f"  [{n}×]  {a}  <->  {b}")
        else:
            lines.append("  (no benchmark pair shares ≥2 moving integers)")
        return "\n".join(lines)

    def __str__(self) -> str:
        lines = [
            f"ISA sensitivity sweep — {len(self.per_integer)} integers × "
            f"{len(self.benchmarks)} benchmarks",
            "",
            f"  {'integer':30} {'val':>4}  {'# moved':>7}   note",
            "  " + "-" * 64,
        ]
        for name, deltas in self.per_integer.items():
            val = next(iter(deltas.values())).get("value", "?")
            crash = any(e.get("status") != "ok" for e in deltas.values())
            n = len(self.movers(name))
            note = ("load-bearing: breaks executability" if crash
                    else "INERT (no prediction moves)" if n == 0 else "")
            lines.append(f"  {name:30} {str(val):>4}  {n:>7}   {note}")
        nm = self.never_moved
        lines += [
            "",
            f"  {len(self.benchmarks) - len(nm)}/{len(self.benchmarks)} benchmarks "
            f"move under ≥1 integer; {len(nm)} never move:",
            "  " + ", ".join(nm) if nm else "  (all benchmarks couple to some integer)",
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------

def integer_sweep(integers: list[str] | None = None,
                  deltas: tuple[int, ...] = (+1, -1)) -> SensitivityReport:
    """Perturb each leaf ISA integer by ``deltas`` and map which predictions move.

    Parameters
    ----------
    integers
        Names to perturb (default: all leaf ``NAME: int = V`` constants).
    deltas
        Perturbations to apply to each integer (default ±1).

    The ``isa/constants.py`` source is temporarily rewritten per perturbation
    and always restored. Requires a writable (editable / dev) install.
    """
    if not os.access(CONSTANTS_PATH, os.W_OK):
        raise RuntimeError(
            f"sensitivity sweep needs a writable constants source ({CONSTANTS_PATH}); "
            "install nwt-substrate in editable/dev mode (`pip install -e .`)."
        )

    orig = open(CONSTANTS_PATH, encoding="utf-8").read()
    all_leaves = {m.group(1): int(m.group(2)) for m in _LEAF_RE.finditer(orig)}
    targets = integers if integers is not None else list(all_leaves)

    baseline = _run_child()
    if "__IMPORT_FAILED__" in baseline:
        raise RuntimeError(f"baseline import failed: {baseline['__IMPORT_FAILED__']}")

    report = SensitivityReport(benchmarks=sorted(baseline), baseline=baseline)
    try:
        for name in targets:
            if name not in all_leaves:
                raise KeyError(f"{name!r} is not a leaf integer constant")
            val = all_leaves[name]
            report.per_integer[name] = {}
            for delta in deltas:
                try:
                    open(CONSTANTS_PATH, "w", encoding="utf-8").write(
                        _patched_source(orig, name, val + delta))
                    pert = _run_child()
                finally:
                    open(CONSTANTS_PATH, "w", encoding="utf-8").write(orig)
                if "__IMPORT_FAILED__" in pert or "__PARSE_FAILED__" in pert:
                    report.per_integer[name][delta] = {
                        "value": val, "status": "crash",
                        "detail": next(iter(pert.values())), "moved": [],
                    }
                    continue
                moved = sorted(k for k in baseline if pert.get(k) != baseline.get(k))
                report.per_integer[name][delta] = {
                    "value": val, "status": "ok", "moved": moved,
                }
    finally:
        open(CONSTANTS_PATH, "w", encoding="utf-8").write(orig)
    return report


def main(argv: list[str] | None = None) -> int:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    report = integer_sweep()
    print(report)
    if "--criticality" in argv or "--load" in argv:
        print()
        print(report.criticality_summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
