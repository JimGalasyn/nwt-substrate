"""Classical electromagnetic fields from substrate sources.

A particle's carrier carries a charge density ρ and a supercurrent density j;
this module solves Maxwell's static equations for the fields they radiate, on a
periodic grid via FFT–Poisson (natural units ε₀ = μ₀ = 1):

    electric:  ∇²φ = −ρ/ε₀,   E = −∇φ              (Coulomb / near field)
    magnetic:  ∇²A = −μ₀ j,   B = ∇×A              (Biot–Savart)

with an optional **Euler–Heisenberg** nonlinear vacuum correction
(:func:`maxwell_eh`), plus source deposition on a carrier curve, multipole
moments (the bridge to electromagnetic form factors / elastic scattering), and
field-line tracing for visualisation.

The k = 0 (uniform) Fourier mode is dropped, so E is the field of (ρ − ⟨ρ⟩) on
the periodic cell — exact near the source (the region of interest); the far
field on a finite periodic cell carries image effects.  ``j`` is assumed
divergence-free (a closed current), so the vector potential is well defined.

References: Jackson (classical EM); Euler & Heisenberg (1936) / Heisenberg &
Euler nonlinear vacuum; NWT Paper 8a (Helmholtz modes on the vortex-tube
surface), Paper 6 (carrier knots).
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "grid", "electric_field", "magnetic_field", "divergence", "curl",
    "maxwell_eh", "deposit_sources", "multipole_moments", "trace_field_lines",
    "form_factor", "mean_square_radius",
]


def _k_grids(N: int, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX**2 + KY**2 + KZ**2
    K2[0, 0, 0] = 1.0          # avoid /0; the k=0 mode is zeroed in each solve
    return KX, KY, KZ, K2


def grid(N: int, L: float) -> tuple[np.ndarray, float, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Return ``(g, dx, (X, Y, Z))`` for a cubic periodic grid of N points / side L."""
    g = np.linspace(-L / 2, L / 2, N, endpoint=False)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    return g, float(g[1] - g[0]), (X, Y, Z)


def _poisson_periodic(src: np.ndarray, dx: float) -> np.ndarray:
    """Solve ∇²u = −src on a periodic cell (FFT); the k=0 mode is dropped."""
    _, _, _, K2 = _k_grids(src.shape[0], dx)
    sk = np.fft.fftn(src); sk[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(sk / K2))


def _poisson_open(src: np.ndarray, dx: float) -> np.ndarray:
    """Solve ∇²u = −src with open (free-space) BCs, via a zero-padded
    convolution with the Green's function 1/(4π r) (Hockney's method).  No
    periodic images, no neutralizing background — the isolated-source field."""
    N = src.shape[0]; M = 2 * N
    o = np.fft.fftfreq(M) * M                       # offsets [0..N-1, -N..-1]
    OX, OY, OZ = np.meshgrid(o, o, o, indexing="ij")
    r = dx * np.sqrt(OX**2 + OY**2 + OZ**2)
    G = np.zeros((M, M, M)); nz = r > 0
    G[nz] = 1.0 / (4 * np.pi * r[nz])
    G[0, 0, 0] = 1.0 / (4 * np.pi * 0.4 * dx)       # regularized self-cell
    sp = np.zeros((M, M, M)); sp[:N, :N, :N] = src
    u = np.real(np.fft.ifftn(np.fft.fftn(sp) * np.fft.fftn(G))) * dx**3
    return u[:N, :N, :N]


def _grad(u: np.ndarray, dx: float, bc: str) -> list[np.ndarray]:
    if bc == "periodic":
        KX, KY, KZ, _ = _k_grids(u.shape[0], dx)
        uk = np.fft.fftn(u)
        return [np.real(np.fft.ifftn(1j * K * uk)) for K in (KX, KY, KZ)]
    d = np.gradient(u, dx)                           # open: finite-difference
    return [d[0], d[1], d[2]]


def electric_field(rho: np.ndarray, dx: float, eps0: float = 1.0, bc: str = "periodic") -> list[np.ndarray]:
    """E from a charge density via ∇²φ = −ρ/ε₀, E = −∇φ.

    bc : "periodic" (FFT, neutralizing background — fine near the source and
         for Gauss-law checks) or "open" (free-space Green's function, no image
         charges — clean Coulomb field lines at all latitudes).
    Returns ``[Ex, Ey, Ez]``.
    """
    phi = (_poisson_open(rho / eps0, dx) if bc == "open"
           else _poisson_periodic(rho / eps0, dx))
    return [-g for g in _grad(phi, dx, bc)]


def magnetic_field(j: list[np.ndarray], dx: float, mu0: float = 1.0, bc: str = "periodic") -> list[np.ndarray]:
    """B from a current density ``j=[jx,jy,jz]`` via ∇²A = −μ₀ j, B = ∇×A.

    bc : "periodic" or "open" (free-space, no image currents).
    """
    solve = _poisson_open if bc == "open" else _poisson_periodic
    A = [mu0 * solve(jc, dx) for jc in j]
    dAx, dAy, dAz = (_grad(Ac, dx, bc) for Ac in A)
    return [dAz[1] - dAy[2], dAx[2] - dAz[0], dAy[0] - dAx[1]]


def divergence(F: list[np.ndarray], dx: float) -> np.ndarray:
    """Spectral divergence of a vector field ``F=[Fx,Fy,Fz]``."""
    KX, KY, KZ, _ = _k_grids(F[0].shape[0], dx)
    return np.real(np.fft.ifftn(1j * (KX * np.fft.fftn(F[0])
                                      + KY * np.fft.fftn(F[1])
                                      + KZ * np.fft.fftn(F[2]))))


def curl(F: list[np.ndarray], dx: float) -> list[np.ndarray]:
    """Spectral curl of a vector field ``F=[Fx,Fy,Fz]``."""
    KX, KY, KZ, _ = _k_grids(F[0].shape[0], dx)
    Fx, Fy, Fz = (np.fft.fftn(c) for c in F)
    return [np.real(np.fft.ifftn(1j * (KY * Fz - KZ * Fy))),
            np.real(np.fft.ifftn(1j * (KZ * Fx - KX * Fz))),
            np.real(np.fft.ifftn(1j * (KX * Fy - KY * Fx)))]


def maxwell_eh(rho, j, dx, *, eh_xi=0.0, n_iter=6, tol=1e-5, eps0=1.0, mu0=1.0,
               bc="periodic"):
    """Static Maxwell fields with the optional Euler–Heisenberg nonlinear
    vacuum correction.

    The EH effective Lagrangian ``L_EH = ξ[(E²−B²)² + 7(E·B)²]`` makes the
    vacuum a nonlinear medium with polarization ``P = ∂L_EH/∂E`` and
    magnetization ``M = −∂L_EH/∂B``, so Gauss/Ampère read ``∇·(E+P)=ρ`` and
    ``∇×(B−M)=j``.  Solved by fixed-point iteration: linear solve → vacuum
    ``P, M`` from ``(E, B)`` → resolve with ``ρ_eff = ρ−∇·P``, ``j_eff = j+∇×M``.

    Parameters
    ----------
    eh_xi : float
        EH coupling in solver units.  The physical value ``2α²/45m⁴`` is ``~α³``
        at the Compton scale (negligible); boost it to probe the strong-field
        regime.  ``eh_xi = 0`` returns the exact linear Maxwell fields.
    n_iter : int
        Max fixed-point iterations (ignored when ``eh_xi == 0``).
    tol : float
        Convergence tolerance on the relative change in E between iterations;
        the loop stops early once below it.  If ``n_iter`` is exhausted without
        converging (e.g. ``eh_xi`` pushed into the non-convergent strong-field
        regime), a ``RuntimeWarning`` is issued.

    Returns
    -------
    (E, B) : each ``[Fx, Fy, Fz]``.
    """
    E = electric_field(rho, dx, eps0, bc=bc)
    B = magnetic_field(j, dx, mu0, bc=bc)
    if eh_xi == 0.0:
        return E, B
    converged = False
    for _ in range(n_iter):
        inv = sum(c * c for c in E) - sum(c * c for c in B)        # E²−B²
        EdotB = sum(e * b for e, b in zip(E, B))
        # P = ∂L_EH/∂E, M = −∂L_EH/∂B for L_EH = ξ[(E²−B²)² + 7(E·B)²]:
        # the 4 = 2·1 (diagonal term) and 14 = 2·7 (mixed term) are the chain-rule
        # factors off the 1:7 invariant ratio documented in the function header.
        P = [4 * eh_xi * inv * E[i] + 14 * eh_xi * EdotB * B[i] for i in range(3)]
        M = [4 * eh_xi * inv * B[i] - 14 * eh_xi * EdotB * E[i] for i in range(3)]
        rho_eff = rho - divergence(P, dx)
        j_eff = [j[i] + c for i, c in enumerate(curl(M, dx))]
        E_new = electric_field(rho_eff, dx, eps0, bc=bc)
        B = magnetic_field(j_eff, dx, mu0, bc=bc)
        num = sum(np.sum((En - Eo) ** 2) for En, Eo in zip(E_new, E))
        den = sum(np.sum(En ** 2) for En in E_new) + 1e-30
        E = E_new
        if np.sqrt(num / den) < tol:
            converged = True
            break
    if not converged:
        import warnings
        warnings.warn("maxwell_eh did not converge within n_iter; the EH "
                      "fixed point may be in the non-convergent (strong-field) "
                      "regime — increase n_iter or reduce eh_xi.", RuntimeWarning)
    return E, B


def deposit_sources(curves, XYZ, dx, *, Q=0.0, I=1.0, width=0.6):
    """Smear total charge ``Q`` and a tangential loop current ``I`` onto a
    Gaussian shell around the carrier curve(s).

    Parameters
    ----------
    curves : list of (M, 3) arrays
        Carrier curve(s); current flows along each curve's local tangent.
    XYZ : (X, Y, Z) grids (from :func:`grid`).
    Q : total charge (cell integral of ρ); 0 → neutral (no E source).
    I : loop-current scale (the current field is normalised to peak |j| = I).

    Returns ``(rho, j)`` with ``j = [jx, jy, jz]``.
    """
    X, Y, Z = XYZ
    rho = np.zeros_like(X)
    jx = np.zeros_like(X); jy = np.zeros_like(X); jz = np.zeros_like(X)
    for C in curves:
        T = np.gradient(C, axis=0)
        T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
        for (px, py, pz), (tx, ty, tz) in zip(C, T):
            g = np.exp(-((X - px) ** 2 + (Y - py) ** 2 + (Z - pz) ** 2) / (2 * width ** 2))
            rho += g
            jx += g * tx; jy += g * ty; jz += g * tz
    tot = rho.sum() * dx ** 3
    rho = rho * (Q / tot) if (tot > 0 and Q != 0.0) else rho * 0.0
    jmax = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2).max() + 1e-12
    return rho, [I * jx / jmax, I * jy / jmax, I * jz / jmax]


def multipole_moments(rho, j, XYZ, dx):
    """Leading EM multipole moments — the bridge to form factors / scattering.

    Returns a dict with:
      * ``charge``          = ∫ ρ dV
      * ``electric_dipole`` = ∫ r ρ dV
      * ``magnetic_dipole`` = ½ ∫ (r × j) dV
    """
    X, Y, Z = XYZ; dV = dx ** 3
    Q = float(rho.sum() * dV)
    p = np.array([float((X * rho).sum() * dV),
                  float((Y * rho).sum() * dV),
                  float((Z * rho).sum() * dV)])
    jx, jy, jz = j
    m = 0.5 * np.array([float((Y * jz - Z * jy).sum() * dV),
                        float((Z * jx - X * jz).sum() * dV),
                        float((X * jy - Y * jx).sum() * dV)])
    return {"charge": Q, "electric_dipole": p, "magnetic_dipole": m}


def trace_field_lines(field, g, seeds, *, n_steps=400, ds=0.15, both_ways=True):
    """Integrate field lines (RK2) from ``seeds`` along the normalised field.

    Returns a list of ``(k, 3)`` polylines.  Requires SciPy.
    """
    from scipy.interpolate import RegularGridInterpolator
    interp = [RegularGridInterpolator((g, g, g), F, bounds_error=False, fill_value=0.0)
              for F in field]
    lo, hi = g[0], g[-1]
    lines = []
    for seed in seeds:
        for direction in ((1, -1) if both_ways else (1,)):
            p = np.array(seed, float); pts = [p.copy()]
            for _ in range(n_steps):
                v = np.array([f([p])[0] for f in interp]); n = np.linalg.norm(v)
                if n < 1e-9:
                    break
                pm = p + 0.5 * ds * direction * v / n
                vm = np.array([f([pm])[0] for f in interp]); nm = np.linalg.norm(vm)
                if nm < 1e-9:
                    break
                p = p + ds * direction * vm / nm
                if not (lo < p[0] < hi and lo < p[1] < hi and lo < p[2] < hi):
                    break
                pts.append(p.copy())
            if len(pts) > 3:
                lines.append(np.array(pts))
    return lines


def mean_square_radius(rho, XYZ, dx):
    """⟨r²⟩ of a charge density (signed; negative is physical for a net-neutral
    distribution with negative charge at larger r, e.g. the neutron).

    For a **net-charged** body, returns the per-unit-charge ⟨r²⟩ = ∫r²ρ / ∫ρ
    (``sqrt`` it for the rms charge radius).  For a **net-neutral** body
    (|∫ρ| negligible vs ∫|ρ|), the per-charge ratio is ill-defined, so the bare
    second moment ∫r²ρ dV is returned (the neutron-⟨r²⟩ convention).  The
    neutral test is a tolerance relative to the charge scale, not ``!= 0`` —
    a floating-point ``∫ρ ≈ 1e-17`` must not divide.
    """
    X, Y, Z = XYZ
    r2rho = ((X**2 + Y**2 + Z**2) * rho).sum()
    Q, scale = rho.sum(), np.abs(rho).sum()
    if abs(Q) <= 1e-9 * scale:                 # net-neutral → bare 2nd moment
        return float(r2rho * dx**3)
    return float(r2rho / Q)                     # net-charged → per-unit-charge ⟨r²⟩


def form_factor(rho, dx, *, bc="open", pad=2, nbins=60, qmax=None):
    """Spherically-averaged elastic form factor |F(q)| of a charge density,
    normalised to F(0) = 1.

    The form factor is the Fourier transform of the charge distribution; it
    enters elastic scattering as ``dσ/dΩ = (dσ/dΩ)_point · |F(q²)|²``, and its
    low-q slope gives the charge radius (``F ≈ 1 − ⟨r²⟩ q²/6``).

    bc : "open" zero-pads the density (isolated source — clean low-q, no
         periodic-image contamination; recommended) or "periodic" (bare FFT).
    pad : zero-padding factor for ``bc="open"`` (padded box = pad × the grid).
    qmax : restrict the binned q range (in 1/length); ``None`` → full grid.

    Returns ``(q, F)`` arrays (q in inverse length units of ``dx``).
    """
    N = rho.shape[0]
    if bc == "open":
        M = pad * N
        rp = np.zeros((M, M, M)); rp[:N, :N, :N] = rho
    else:
        M, rp = N, rho
    Fk = np.abs(np.fft.fftn(rp))
    F = Fk / Fk.flat[0]                              # |FFT|, normalised to F(0)=1
    k = 2 * np.pi * np.fft.fftfreq(M, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K = np.sqrt(KX**2 + KY**2 + KZ**2)
    hi = qmax if qmax is not None else K.max() * 0.6
    edges = np.linspace(0, hi, nbins)
    q, Fq = [], []
    for i in range(len(edges) - 1):
        m = (K >= edges[i]) & (K < edges[i + 1])
        if m.sum() > 4:
            q.append(K[m].mean()); Fq.append(F[m].mean())
    return np.array(q), np.array(Fq)
