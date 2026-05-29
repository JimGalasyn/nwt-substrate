"""Experiment definitions expressed against the qpu interface.

exp12 Phase 2 is the M3 acceptance consumer: the same 4 walks + 4 controls that
passed 8/8 (ancilla) on IBM Heron and 8/8 (destructive) on ibm_fez + IQM Garnet
on 2026-05-25 — now built once, vendor-neutrally.
"""
from __future__ import annotations

import json
from pathlib import Path

from .runner import Item
from .spec import steane_base_ops, control_base_ops, STEANE_LZ

_LOOKUP = Path(__file__).parents[2] / "analysis" / "steane_exp12_walk_pauli_table.json"

PHASE2_TARGETS = [(1, 3), (2, 1), (3, 5), (1, 8)]   # proton, e-, pi+, mu-
_CONTROL_PRED = {                                    # (x_fano, z_fano, logical_z)
    "identity": ("I", "I", +1),
    "x7": ("I", "I", -1),
    "z7": ("I", "I", +1),
    "h7": ("I", "I", 0),                             # 50/50 logical-Z split
}


def exp12_phase2_items(include_controls: bool = True) -> list[Item]:
    lk = json.loads(_LOOKUP.read_text())
    pq = {(e["p"], e["q"]): e for e in lk["particles"]}
    items: list[Item] = []
    for p, q in PHASE2_TARGETS:
        fp = pq[(p, q)]["forward_pauli"]
        xb, zb = fp["x_part_bits"], fp["z_part_bits"]
        xf, zf = fp["predicted_syndrome_fano_pair"]
        lz = +1 if sum(xb[i] for i in STEANE_LZ) % 2 == 0 else -1
        items.append(Item(f"class_p{p}_q{q}", steane_base_ops(xb, zb), (xf, zf, lz)))
    if include_controls:
        for kind in ("identity", "x7", "z7", "h7"):
            items.append(Item(f"control_{kind}", control_base_ops(kind), _CONTROL_PRED[kind]))
    return items
