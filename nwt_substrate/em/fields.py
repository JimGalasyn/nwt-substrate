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
]


def _k_grids(N, dx):
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX**2 + KY**2 + KZ**2
    K2[0, 0, 0] = 1.0          # avoid /0; the k=0 mode is zeroed in each solve
    return KX, KY, KZ, K2


def grid(N, L):
    """Return ``(g, dx, (X, Y, Z))`` for a cubic periodic grid of N points / side L."""
    g = np.linspace(-L / 2, L / 2, N, endpoint=False)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    return g, float(g[1] - g[0]), (X, Y, Z)


def electric_field(rho, dx, eps0=1.0):
    """E from a charge density via ∇²φ = −ρ/ε₀, E = −∇φ.  Returns ``[Ex, Ey, Ez]``."""
    KX, KY, KZ, K2 = _k_grids(rho.shape[0], dx)
    rho_k = np.fft.fftn(rho); rho_k[0, 0, 0] = 0.0
    phi_k = rho_k / (eps0 * K2)
    return [np.real(np.fft.ifftn(-1j * K * phi_k)) for K in (KX, KY, KZ)]


def magnetic_field(j, dx, mu0=1.0):
    """B from a current density ``j=[jx,jy,jz]`` via ∇²A = −μ₀ j, B = ∇×A."""
    KX, KY, KZ, K2 = _k_grids(j[0].shape[0], dx)
    Ak = []
    for jc in j:
        jc_k = np.fft.fftn(jc); jc_k[0, 0, 0] = 0.0
        Ak.append(jc_k * mu0 / K2)
    Ax, Ay, Az = Ak
    return [np.real(np.fft.ifftn(1j * (KY * Az - KZ * Ay))),
            np.real(np.fft.ifftn(1j * (KZ * Ax - KX * Az))),
            np.real(np.fft.ifftn(1j * (KX * Ay - KY * Ax)))]


def divergence(F, dx):
    """Spectral divergence of a vector field ``F=[Fx,Fy,Fz]``."""
    KX, KY, KZ, _ = _k_grids(F[0].shape[0], dx)
    return np.real(np.fft.ifftn(1j * (KX * np.fft.fftn(F[0])
                                      + KY * np.fft.fftn(F[1])
                                      + KZ * np.fft.fftn(F[2]))))


def curl(F, dx):
    """Spectral curl of a vector field ``F=[Fx,Fy,Fz]``."""
    KX, KY, KZ, _ = _k_grids(F[0].shape[0], dx)
    Fx, Fy, Fz = (np.fft.fftn(c) for c in F)
    return [np.real(np.fft.ifftn(1j * (KY * Fz - KZ * Fy))),
            np.real(np.fft.ifftn(1j * (KZ * Fx - KX * Fz))),
            np.real(np.fft.ifftn(1j * (KX * Fy - KY * Fx)))]


def maxwell_eh(rho, j, dx, *, eh_xi=0.0, n_iter=6, eps0=1.0, mu0=1.0):
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
        Fixed-point iterations (ignored when ``eh_xi == 0``).

    Returns
    -------
    (E, B) : each ``[Fx, Fy, Fz]``.
    """
    E = electric_field(rho, dx, eps0)
    B = magnetic_field(j, dx, mu0)
    if eh_xi == 0.0:
        return E, B
    for _ in range(n_iter):
        inv = sum(c * c for c in E) - sum(c * c for c in B)        # E²−B²
        EdotB = sum(e * b for e, b in zip(E, B))
        P = [4 * eh_xi * inv * E[i] + 14 * eh_xi * EdotB * B[i] for i in range(3)]
        M = [4 * eh_xi * inv * B[i] - 14 * eh_xi * EdotB * E[i] for i in range(3)]
        rho_eff = rho - divergence(P, dx)
        cM = curl(M, dx)
        j_eff = [j[i] + cM[i] for i in range(3)]
        E = electric_field(rho_eff, dx, eps0)
        B = magnetic_field(j_eff, dx, mu0)
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
