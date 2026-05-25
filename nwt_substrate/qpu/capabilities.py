"""Per-backend capability registry + the `preflight()` guardrail.

`preflight()` estimates routed 2-qubit-gate count, signal fidelity, cost, and
execution-window timing BEFORE anything is submitted — the guardrail that turns
the 2026-05-25 IQM Garnet "$3.56 on a 58-CZ circuit that couldn't fit" lesson
into a refuse-to-submit-blind. Routing is estimated vendor-neutrally by
transpiling onto the backend's real coupling map (the same method that
diagnosed Garnet: 58 CZ ancilla vs 15 CZ destructive).
"""
from __future__ import annotations

import json
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

from .adapters.base import Capabilities, Window

# Static per-backend info not worth a live query. coupling_map / windows are
# fetched live and merged in.
#
# NOTE on median_2q_fidelity: this is an EFFECTIVE all-in per-2q-gate figure
# (folding in routing crosstalk, idle decoherence over depth, and readout), NOT
# the vendor's headline calibration number — headline figures (~0.995 for
# Garnet) badly overpredict deep-circuit success. It is calibrated conservatively
# to observed flat-vs-signal behavior: Garnet went FLAT at ~67 CZ but PASSED 8/8
# at ~15 CZ on 2026-05-25, which brackets the effective f to (0.955, 0.99) at
# FIT_THRESHOLD=0.5; we take 0.98. This is the coarse v1 WARNING heuristic of the
# design doc (§11), to be replaced by a noise-aware estimate later.
REGISTRY: dict[str, dict] = {
    "iqm/Garnet":              dict(sdk="braket", region="eu-north-1", native_2q="cz",  per_shot_usd=0.00145, task_fee_usd=0.30, median_2q_fidelity=0.98),
    "iqm/Emerald":             dict(sdk="braket", region="eu-north-1", native_2q="cz",  per_shot_usd=0.0016,  task_fee_usd=0.30, median_2q_fidelity=0.98),
    "aqt/Ibex-Q1":             dict(sdk="braket", region="eu-north-1", native_2q="xx",  per_shot_usd=0.0235,  task_fee_usd=0.30, median_2q_fidelity=0.99),
    "rigetti/Cepheus-1-108Q":  dict(sdk="braket", region="us-west-1",  native_2q="cz",  per_shot_usd=0.000425, task_fee_usd=0.30, median_2q_fidelity=0.97),
}

CACHE_DIR = Path.home() / ".cache" / "nwt_qpu"
CACHE_TTL_S = 24 * 3600
FIT_THRESHOLD = 0.5
_DAYS = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
# Braket executionDay can be an aggregate (IQM Garnet reports "Weekdays") as well
# as a specific weekday (AQT reports "Monday", "Tuesday", ...).
_DAY_GROUPS = {
    "Everyday": set(range(7)), "Daily": set(range(7)),
    "Weekdays": {0, 1, 2, 3, 4}, "Weekend": {5, 6},
}


def _weekday_indices(day: str) -> set[int]:
    return _DAY_GROUPS[day] if day in _DAY_GROUPS else {_DAYS.index(day)}


@dataclass
class Preflight:
    backend: str
    n_2q: int
    depth: int
    est_signal_fidelity: float
    est_cost_usd: float
    fits_budget: bool
    in_window: bool
    next_window_utc: str | None

    def summary(self) -> str:
        flag = "OK" if self.fits_budget else "!! BLOWS BUDGET"
        win = "open" if self.in_window else f"closed (next {self.next_window_utc})"
        return (f"{self.backend}: 2q={self.n_2q} depth={self.depth} "
                f"fid~{self.est_signal_fidelity:.2f} [{flag}]  "
                f"${self.est_cost_usd:.2f}  window={win}")


# ----------------------------------------------------------------- fetching
def _hms_to_s(hms: str) -> int:
    h, m, s = (int(x) for x in hms.split(":"))
    return h * 3600 + m * 60 + s


def _fetch_braket(backend: str, region: str) -> tuple[tuple, tuple]:
    """Return (coupling_map edges 0-indexed, execution windows) from AWS."""
    arn = f"arn:aws:braket:{region}::device/qpu/{backend}"
    out = subprocess.run(
        ["aws", "braket", "get-device", "--region", region,
         "--device-arn", arn, "--query", "deviceCapabilities", "--output", "text"],
        capture_output=True, text=True, timeout=60).stdout
    d = json.loads(out)
    graph = d.get("paradigm", {}).get("connectivity", {}).get("connectivityGraph", {})
    edges = set()
    for k, vs in graph.items():
        for v in vs:
            a, b = int(k) - 1, int(v) - 1  # IQM labels are 1-indexed
            edges.add((min(a, b), max(a, b)))
    windows = tuple(
        Window(w["executionDay"], w["windowStartHour"], w["windowEndHour"])
        for w in d.get("service", {}).get("executionWindows", []))
    return tuple(sorted(edges)), (windows or None)


def fetch_capabilities(backend: str, region: str | None = None,
                       service=None, use_cache: bool = True) -> Capabilities:
    """Build a fully-populated Capabilities (coupling map + windows) for a backend."""
    if backend.startswith("ibm_"):
        # IBM: free, 24/7, coupling map from the live backend target.
        from .adapters.ibm import IBMBackend
        b = IBMBackend(backend, service=service)
        cm = tuple(sorted({(min(a, b_), max(a, b_)) for a, b_ in b._backend.coupling_map}))
        return Capabilities(name=backend, sdk="ibm", n_qubits=b.capabilities.n_qubits,
                            native_2q="cz", coupling_map=cm, median_2q_fidelity=0.99,
                            execution_windows=None)

    info = REGISTRY.get(backend)
    if info is None:
        raise KeyError(f"no registry entry for {backend!r}")
    region = region or info["region"]

    cache = CACHE_DIR / f"{backend.replace('/', '_')}.json"
    cm = win = None
    if use_cache and cache.exists() and (time.time() - cache.stat().st_mtime) < CACHE_TTL_S:
        c = json.loads(cache.read_text())
        cm = tuple(tuple(e) for e in c["coupling_map"])
        win = tuple(Window(*w) for w in c["windows"]) if c["windows"] else None
    else:
        cm, win = _fetch_braket(backend, region)
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps({
            "coupling_map": [list(e) for e in cm],
            "windows": [[w.day, w.start_hour, w.end_hour] for w in win] if win else None,
            "fetched_at": time.time()}))

    return Capabilities(
        name=backend, sdk=info["sdk"], n_qubits=max(max(e) for e in cm) + 1 if cm else 0,
        native_2q=info["native_2q"], coupling_map=cm,
        per_shot_usd=info["per_shot_usd"], task_fee_usd=info["task_fee_usd"],
        median_2q_fidelity=info["median_2q_fidelity"], execution_windows=win)


# ------------------------------------------------------------- window logic
def window_status(windows, now: datetime | None = None) -> tuple[bool, str | None]:
    """Return (currently_in_window, next_open_iso_utc). None windows = 24/7."""
    if not windows:
        return True, None
    now = now or datetime.now(timezone.utc)
    now_wd, now_s = now.weekday(), now.hour * 3600 + now.minute * 60 + now.second

    spans = [(wd, _hms_to_s(w.start_hour), _hms_to_s(w.end_hour))
             for w in windows for wd in _weekday_indices(w.day)]
    for wd, s0, s1 in spans:
        if wd == now_wd and s0 <= now_s <= s1:
            return True, None

    # soonest future opening within the next 7 days
    best = None
    for wd, s0, _ in spans:
        days_ahead = (wd - now_wd) % 7
        if days_ahead == 0 and s0 <= now_s:
            days_ahead = 7  # today's window already started/passed
        cand = (now + timedelta(days=days_ahead)).replace(
            hour=s0 // 3600, minute=(s0 % 3600) // 60, second=s0 % 60, microsecond=0)
        if best is None or cand < best:
            best = cand
    return False, (best.isoformat() if best else None)


# ----------------------------------------------------------------- preflight
def _basis_for(native_2q: str) -> list[str]:
    two = {"cz": "cz", "ecr": "ecr", "xx": "rxx", "rzz": "rzz"}.get(native_2q, "cz")
    return [two, "rz", "sx", "x"]


def preflight(spec, caps: Capabilities, shots: int, now: datetime | None = None) -> Preflight:
    """Estimate routed 2q-count / fidelity / cost / window for `spec` on `caps`."""
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap

    cm = CouplingMap([[a, b] for a, b in caps.coupling_map]
                     + [[b, a] for a, b in caps.coupling_map])
    t = transpile(spec.to_qiskit(), coupling_map=cm, basis_gates=_basis_for(caps.native_2q),
                  optimization_level=3, seed_transpiler=11)
    n2q = t.num_nonlocal_gates()
    fid = caps.median_2q_fidelity ** n2q
    cost = shots * caps.per_shot_usd + caps.task_fee_usd
    in_window, nxt = window_status(caps.execution_windows, now)
    return Preflight(caps.name, n2q, t.depth(), fid, cost, fid >= FIT_THRESHOLD, in_window, nxt)
