"""Faddeev-Skyrme (Hopf) soliton primitives -- the construct/measure layer.

Pure-numpy, rendering-free reference implementations of the CP^1 field map, the
faithful (Berg-Luscher) area form, the geometric Whitehead Hopf charge, the
rational-map hopfion initial condition, and the Faddeev energy / soliton-size
diagnostics.  These are the forward, diagnostic half of the soliton work; the
JAX/GPU energy-minimising relaxers live in the separate engine and import these
as their reference (and test oracle).

Quick start
-----------
    from nwt_substrate.solitons import BoxGrid
    from nwt_substrate.solitons.faddeev import rational_hopfion, n_from_Z
    from nwt_substrate.solitons.faddeev import whitehead_hopf_charge

    grid = BoxGrid(N=64, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)        # smooth Q_H = 1 hopfion
    n1, n2, n3 = n_from_Z(Z1, Z2)
    whitehead_hopf_charge(n1, n2, n3, grid)       # ~ +1

References
----------
Faddeev & Niemi, "Stable knot-like structures in classical field theory"
(Nature 1997).  Battye & Sutcliffe (1998); Hietarinta & Salo (2000) for the
rational-map / area-form construction.  NWT Paper 16 (master Lagrangian, L_3
Skyrme-Faddeev sector).
"""
from .cp1 import n_from_Z
from .area_form import solid_angle_triangle, plaquette_area_form
from .seeds import rational_hopfion
from .energy import faddeev_energy, energy_density, soliton_size
from .charge import whitehead_hopf_charge

__all__ = [
    "n_from_Z",
    "solid_angle_triangle",
    "plaquette_area_form",
    "rational_hopfion",
    "faddeev_energy",
    "energy_density",
    "soliton_size",
    "whitehead_hopf_charge",
]
