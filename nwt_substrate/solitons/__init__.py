"""Topological solitons on a periodic lattice.

Umbrella for NWT's soliton sectors and their shared coordinate grid:

  - `grid.BoxGrid`           periodic cubic grid (shared)
  - `faddeev`                Faddeev-Skyrme (Hopf) solitons -- S^2 field, Hopf
                             charge, knotted/linked hopfions (the L_3 sector)
  - `skyrme`                 (future) Skyrme solitons -- S^3 field, baryon
                             degree, localised lumps

The `faddeev` public API is re-exported here for convenience.  Everything is
pure-numpy and rendering-free; the energy-minimising relaxers live in the
separate JAX/GPU engine.
"""
from .grid import BoxGrid
from .faddeev import (
    n_from_Z,
    solid_angle_triangle,
    plaquette_area_form,
    rational_hopfion,
    faddeev_energy,
    energy_density,
    soliton_size,
    whitehead_hopf_charge,
)

__all__ = [
    "BoxGrid",
    "n_from_Z",
    "solid_angle_triangle",
    "plaquette_area_form",
    "rational_hopfion",
    "faddeev_energy",
    "energy_density",
    "soliton_size",
    "whitehead_hopf_charge",
]
