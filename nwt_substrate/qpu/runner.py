"""Experiment lifecycle: build -> (preflight/stage) -> submit -> poll -> decode -> JSON.

A `Backend` adapter does the vendor plumbing; this orchestrates one experiment
across all modes (dry-run / stage / submit / resume) and emits standardized
results. Submit writes resume handles to disk BEFORE polling, so a window-close
or shell timeout mid-run is always recoverable (the orphan-safety that saved the
IQM + AQT weekend runs).
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

from . import spec as _spec
from . import decode as _decode


@dataclass
class Item:
    label: str
    base_ops: list                 # prep + Pauli word (or control) — list[G]
    predicted: tuple               # (x_fano, z_fano, logical_z)


def build_specs(item: Item, scheme: str) -> list[tuple[str, object]]:
    """Return [(tag, CircuitSpec), ...] for an item under a measurement scheme."""
    if scheme == "destructive":
        z, x = _spec.destructive_css(item.base_ops, item.label)
        return [("z", z), ("x", x)]
    if scheme == "ancilla":
        return [("a", _spec.ancilla_syndrome(item.base_ops, item.label))]
    raise ValueError(f"unknown scheme {scheme!r}")


def choose_scheme(caps) -> str:
    """Default measurement-scheme policy (set 2026-05-25): **destructive CSS
    readout everywhere**.

    The ancilla scheme decohered to flat output on IQM Garnet's sparse grid
    (58-67 CZ after routing), while destructive (15 CZ) passed 8/8 there and
    matched it on IBM heavy-hex (8/8, modal 0.76-0.83). Destructive is strictly
    shallower on every topology and trades only an extra circuit per walk, so it
    is the portable default. Pass `scheme='ancilla'` explicitly to override
    (e.g. to reproduce a legacy heavy-hex run)."""
    return "destructive"


def _verdict_for(item: Item, scheme: str, counts_by_tag: dict) -> _decode.Verdict:
    px, pz, plz = item.predicted
    if scheme == "destructive":
        x_dist, z_dist = _decode.destructive_dists(
            counts_by_tag["z"].single(), counts_by_tag["x"].single())
        return _decode.verdict_destructive(item.label, x_dist, z_dist, px, pz, plz)
    return _decode.verdict_ancilla(item.label, _decode.ancilla_dist(counts_by_tag["a"]), px, pz, plz)


def _print_verdict(v):
    if v.scheme == "destructive":
        xh, zh = v.x_half, v.z_half
        print(f"  {v.label:<18} "
              f"X[{xh.modal:>2} p={xh.modal_prob:.2f} pred={xh.predicted:>2} {'ok' if xh.passed else 'XX'}]  "
              f"Z[{str(zh.modal):<9} p={zh.modal_prob:.2f} pred={str(zh.predicted):<9} {'ok' if zh.passed else 'XX'}]  "
              f"-> {'PASS' if v.passed else 'FAIL'}")
    else:
        j = v.joint
        print(f"  {v.label:<18} modal={j.modal} p={j.modal_prob:.2f} pred={j.predicted} "
              f"-> {'PASS' if v.passed else 'FAIL'}")


def _save_handles(out_dir, exp, backend_label, flat_specs, handles):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    p = out_dir / f"{exp}_handles_{backend_label.replace('/', '_')}_{stamp}.json"
    p.write_text(json.dumps(
        [{"order": k, "item_idx": flat_specs[k][0], "tag": flat_specs[k][1],
          "label": flat_specs[k][2].name, "handle": handles[k] if k < len(handles) else None}
         for k in range(len(flat_specs))], indent=2))
    print(f"  handles saved -> {p}")
    return p


def run(items, backend, shots: int = 100, scheme: str = "auto", mode: str = "dry-run",
        caps=None, out_dir=None, exp: str = "exp", resume_handles=None, verbose: bool = True):
    """Run `items` on `backend`. Returns list[Verdict] (None for stage mode)."""
    if scheme == "auto":
        scheme = choose_scheme(caps)
    flat = [(i, tag, sp) for i, it in enumerate(items) for tag, sp in build_specs(it, scheme)]

    if mode == "stage":
        from . import capabilities as _cap
        if caps is None:
            raise ValueError("stage mode requires caps=")
        print(f"  Preflight ({scheme}, {shots} shots) on {caps.name}:")
        any_blown = False
        for i, tag, sp in flat:
            pf = _cap.preflight(sp, caps, shots)
            any_blown |= not pf.fits_budget
            print(f"    {sp.name:<20} {pf.summary()}")
        print(f"  {'!! at least one circuit blows the budget' if any_blown else 'all fit'}")
        return None

    compiled = [backend.compile(sp) for _, _, sp in flat]
    if mode == "resume":
        job = backend.resume(resume_handles, compiled)
    else:  # dry-run / submit
        job = backend.submit(compiled, shots)
        if out_dir and mode == "submit":
            _save_handles(out_dir, exp, backend.capabilities.name, flat, job.handles)

    all_counts = job.result()
    by_item: dict[int, dict] = {}
    for k, (i, tag, _) in enumerate(flat):
        by_item.setdefault(i, {})[tag] = all_counts[k]

    verdicts = []
    if verbose:
        print()
    for i, it in enumerate(items):
        v = _verdict_for(it, scheme, by_item[i])
        verdicts.append(v)
        if verbose:
            _print_verdict(v)

    if out_dir:
        _save_results(out_dir, exp, backend.capabilities.name, mode, shots, scheme, verdicts)
    if verbose:
        print(f"\n  Overall: {sum(v.passed for v in verdicts)}/{len(verdicts)} PASS")
    return verdicts


def _save_results(out_dir, exp, backend_label, mode, shots, scheme, verdicts):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    p = out_dir / f"{exp}_{scheme}_{mode}_{backend_label.replace('/', '_')}_{stamp}.json"
    payload = {
        "experiment": exp, "backend": backend_label, "mode": mode, "scheme": scheme,
        "shots_per_circuit": shots, "timestamp": stamp,
        "n_pass": sum(v.passed for v in verdicts), "n_total": len(verdicts),
        "circuits": [asdict(v) for v in verdicts],
    }
    p.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\n  Saved: {p}")
    return p
