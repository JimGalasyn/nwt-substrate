"""Classical EM fields from substrate sources (charge density, current density).

FFT–Poisson solvers for the electric and magnetic fields a particle's carrier
radiates, with an optional Euler–Heisenberg nonlinear vacuum correction, plus
multipole moments (the bridge to electromagnetic form factors / scattering) and
field-line tracing.  See :mod:`nwt_substrate.em.fields`.

    >>> import nwt_substrate.em as em
    >>> g, dx, XYZ = em.grid(64, 24.0)
    >>> rho, j = em.deposit_sources([curve], XYZ, dx, Q=-1.0, I=1.0)   # doctest: +SKIP
    >>> E, B = em.maxwell_eh(rho, j, dx, eh_xi=0.0)                    # doctest: +SKIP
"""
from nwt_substrate.em.fields import (
    grid,
    electric_field,
    magnetic_field,
    divergence,
    curl,
    maxwell_eh,
    deposit_sources,
    multipole_moments,
    trace_field_lines,
)

__all__ = [
    "grid", "electric_field", "magnetic_field", "divergence", "curl",
    "maxwell_eh", "deposit_sources", "multipole_moments", "trace_field_lines",
]
