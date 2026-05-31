#!/usr/bin/env python3
"""Render the particle-portrait gallery from the library's portraits engine.

Each NWT particle is a BPS vortex defect whose core follows its carrier-knot
curve; the portraits module bends the exact BPS scalar profile along that knot
and renders the abelian-Higgs phase field as a glowing, phase-coloured
filament (see nwt_substrate.portraits).

Usage:
    python3 diagrams/particle_portraits.py [OUTPUT_DIR]

Default output dir: docs/assets/portraits/
Writes: carrier_gallery.png  + a few high-resolution single portraits.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nwt_substrate.portraits as pp


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).resolve().parent.parent / "docs" / "assets" / "portraits"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    # Contact sheet of the whole carrier zoo.
    fig = pp.gallery(N=150, ncols=3,
                     save_to=out_dir / "carrier_gallery.png")
    plt.close(fig)
    print(f"wrote {out_dir/'carrier_gallery.png'}")

    # A few high-resolution single portraits.
    singles = {
        "electron_unknot": dict(p=2, q=1),
        "proton_trefoil": dict(p=2, q=3),
        "nucleon_cinquefoil": dict(p=2, q=5),
        "meson_hopf": dict(hopf=True),
    }
    for name, spec in singles.items():
        fig = pp.portrait(N=220, save_to=out_dir / f"{name}.png", **spec)
        plt.close(fig)
        print(f"wrote {out_dir/(name + '.png')}")


if __name__ == "__main__":
    main()
