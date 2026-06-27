"""
Dynamical-framing Chern-Simons modification of L_NWT -- the Chiral Magnetic
Effect (CME) that switches on the baryon-asymmetry sign.

The gap (Paper-16 audit)
------------------------
Paper-16's L_NWT carries a Hopf term ``theta * Q_H[n]`` with **theta = 0 or pi,
constant** (the CPT / 0-pi choice).  A constant ``theta * Q_H`` is a topological
total derivative: it contributes no force, no link-protection, and -- crucially
for cosmology -- no dynamics that could *select* a sign for the linking.

The promotion
-------------
Promote the constant Hopf angle to a dynamical framing field ``sigma(x)`` (the
meridional twist / framing phase of the carrier knot, = weak isospin = the
B3-center full twist).  The added Lagrangian is

    dL = 1/2 f^2 (d_mu sigma)^2 + kappa * sigma * eps^{mu nu rho s} F_{mu nu} F_{rho s}
       = 1/2 f^2 (d sigma)^2 + kappa * sigma * (F wedge F),

with ``sigma -> const`` recovering Paper 16 exactly (a strict extension).  Term
(2) is the axion-type Chern-Simons coupling with ``sigma`` in the axion's role.
Integrating it by parts gives a CS 3-form ``-2 kappa (d_mu sigma) eps A d A``;
the gauge equation of motion then gains a current

    j = 2 * mu5 * B ,        mu5 = d_t sigma ,

i.e. the **Chiral Magnetic Effect**.  When the parent black hole imprints a slow
background framing twist rate ``d_t sigma = mu5 ~ sgn(J_parent)``, this is a
chiral chemical potential, and the CME converts that chirality imbalance into a
**net magnetic helicity / net linking** -- the one mechanism that beats the
Rivers-Volovik null (PRL 127, 115702, 2021), which forbids a sign imbalance from
a mere *local* winding-sign bias.  This realizes the postulate
``sgn(eta_B) = sgn(J_parent)`` dynamically (see `cosmology.eta_B_sign`).

Sign discipline (real-time vs imaginary-time)
---------------------------------------------
The CME helicity-amplification sign is a REAL-TIME (induction-time) effect:
``sgn(net helicity) = + sgn(mu5)``.  This is the OPPOSITE of imaginary-time
energy descent, which minimises ``+ mu5 * H`` and so drives the helicity to
``- sgn(mu5)`` (the energetic / link-protection effect).  The two jobs (a:
protection, b: sign selection) of the one ``sigma`` field live in the two time
directions; this module encodes the real-time (baryogenesis) sign.

Validation
----------
Validated against the real-time chiral-MHD induction equation
``d_t A = eta (lap A + 2 mu5 curl A)`` in
``null-worldtube-private/simulations/gauged_relaxer/chiral_mhd_helicity.py``:
the analytic helical-mode growth rate ``eta(-k^2 +/- 2 mu5 k)`` was reproduced
to ~3 decimals, and the achiral-Kibble-Zurek ensemble returned
``sgn<H> = + sgn(mu5)`` (8/8 seeds, ~33 sigma; baseline <H> ~ 0 at mu5 = 0).
The full sign chain was then closed end-to-end in a self-consistent gauged-GPE
quench (`gauged_gpe_cme.py quench`): OUTCOME 3, 7.8 sigma, odd in mu5
(`analysis/cs_term_sign_preregistration.md`).

References
----------
Joyce & Shaposhnikov, PRL 79, 1193 (1997) -- primordial hypermagnetic helicity
from a chiral asymmetry, and the chiral-charge back-reaction.
Boyarsky, Frohlich & Ruchayskiy, PRL 108, 031301 (2012) -- chiral-MHD inverse
cascade / chiral plasma instability.
Gudnason & Nitta, PRD 101, 065011 (2020) -- baryon number = vortex linking.
Rivers & Volovik, PRL 127, 115702 (2021) -- the no-net-linking null this beats.
NWT Paper 16 (L_NWT); `analysis/dynamical_framing_cs_term.md`.
"""
from __future__ import annotations

#: Validated source simulation (null-worldtube-private). Real-time chiral-MHD
#: induction; the ensemble result was sgn<H> = + sgn(mu5) (OUTCOME 3, ~33 sigma).
SOURCE_SIM: str = (
    "null-worldtube-private/simulations/gauged_relaxer/chiral_mhd_helicity.py"
)

#: Joyce-Shaposhnikov back-reaction note. A fixed mu5 pumps the gauge field
#: without bound (the CME instability never saturates with constant chiral
#: charge). Physically, helicity production is sourced by *finite* chiral charge:
#: as helicity H is generated the chiral chemical potential is depleted so that
#:
#:     mu5 + Gamma * H = const
#:
#: (Gamma a positive transfer coefficient; the anomaly trades chiral charge for
#: magnetic helicity). This conservation law saturates the instability at finite
#: amplitude -- it is the correct finite-chiral-charge physics, not a numerical
#: cutoff. Encoded so downstream lattice runs apply it rather than blowing up.
JOYCE_SHAPOSHNIKOV_NOTE: str = (
    "Finite chiral charge: mu5 + Gamma*H = const -> the CME instability "
    "saturates as produced magnetic helicity depletes mu5 "
    "(Joyce-Shaposhnikov 1997)."
)


def cme_growth_rate(
    k: float,
    mu5: float,
    eta: float = 1.0,
    helicity: int = +1,
) -> float:
    """Closed-form growth rate of a helical mode under the chiral-MHD induction.

    For the real-time induction equation ``d_t A = eta (lap A + 2 mu5 curl A)``,
    a helical eigenmode with ``curl A = helicity * k * A`` grows as
    ``A ~ exp(lambda t)`` with

        lambda(k) = eta * (-k^2 + helicity * 2 * mu5 * k).

    The diffusive term ``-eta k^2`` always damps; the CME curl term
    ``helicity * 2 * eta * mu5 * k`` destabilises one helicity sector.  For
    ``mu5 > 0`` the positive-helicity modes with ``k < 2 mu5`` are unstable
    (``lambda > 0``) and grow out of an achiral seed, so the ensemble develops a
    net magnetic helicity of sign ``+ sgn(mu5)`` -- the dynamical realization of
    the eta_B sign.  The negative-helicity sector is the mirror (it decays).

    NOTE on sign: this is the REAL-TIME (induction) growth rate.  It is the
    OPPOSITE sign to imaginary-time energy descent, which minimises ``+ mu5 H``
    and drives the helicity to ``- sgn(mu5)`` (the link-protection effect).

    Parameters
    ----------
    k : float
        Mode wavenumber (>= 0).
    mu5 : float
        Chiral chemical potential, ``mu5 = d_t sigma ~ sgn(J_parent)``.
    eta : float, optional
        Magnetic resistivity / diffusivity (> 0). Default 1.0.
    helicity : int, optional
        Sign of the mode's helicity, +1 (default) or -1.

    Returns
    -------
    float
        The growth rate ``lambda(k)`` of the mode amplitude (the helicity
        ``H ~ |A|^2`` grows at ``2 * lambda``).
    """
    return eta * (-(k ** 2) + helicity * 2.0 * mu5 * k)


def cme_current(B, mu5: float):
    """The Chiral Magnetic Effect current ``j = 2 * mu5 * B``.

    The current that the dynamical-framing Chern-Simons term adds to the gauge
    equation of motion (the lattice CME term).  It is parallel to the magnetic
    field with coefficient ``2 * mu5``; it sources magnetic helicity and is what
    biases vortex linking toward ``sgn(net linking) = + sgn(mu5)``.

    Parameters
    ----------
    B : array_like or float
        Magnetic field ``B = curl A`` (any shape; numpy-broadcastable).
    mu5 : float
        Chiral chemical potential ``mu5 = d_t sigma``.

    Returns
    -------
    Same type/shape as ``B``
        The CME current ``2 * mu5 * B``.
    """
    return 2.0 * mu5 * B
