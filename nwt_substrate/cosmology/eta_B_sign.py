"""
The DYNAMICAL SIGN of the baryon asymmetry -- companion to `eta_B.py`.

`eta_B.py` fixes the MAGNITUDE ``eta_B = 3 alpha^4 / 14`` at zero free
parameters and then merely POSTULATES the sign:

    "The sign sgn(eta_B) = sgn(J_parent) (matter over antimatter) is fixed by
     the parent black hole's spin ... this module gives the magnitude."

This companion module records the MECHANISM behind that postulated sign -- the
chain that makes ``sgn(eta_B) = sgn(J_parent)`` a derived consequence rather
than an assumption:

    J_parent  ->  mu5 = d_t sigma          (parent-BH spin imprints a framing
                                            twist rate on the dynamical Hopf
                                            field sigma; the dynamical-framing
                                            Chern-Simons term)
              ->  net magnetic helicity, sgn = + sgn(mu5)
                                            (Chiral Magnetic Effect; real-time
                                            chiral-MHD; see
                                            condensate.chiral_magnetic)
              ->  net vortex linking, sgn = + sgn(<H>)
                                            (baryon # = vortex linking #,
                                            Gudnason-Nitta PRD 101, 065011, 2020)
              ->  B > 0  =>  eta_B > 0  =  MATTER over antimatter.

It was validated DYNAMICALLY: the real-time chiral-MHD helicity test plus the
self-consistent gauged-GPE quench both returned the pre-registered SUCCESS
(OUTCOME 3) -- ``net linking`` odd in mu5 with the predicted sign, 7.8 sigma
(`analysis/cs_term_sign_preregistration.md`,
`analysis/dynamical_framing_cs_term.md`).  This is the channel that beats the
Rivers-Volovik null (PRL 127, 115702, 2021): the bias is a *bulk* chiral
chemical potential (CME), not a forbidden local winding-sign preference.

The MAGNITUDE is NOT recomputed here -- it stays the substrate-algebra /
Murasugi-Jones result ``3 alpha^4 / 14`` in `eta_B.py`.  This module only
encodes the sign and its mechanism.

References
----------
NWT `cosmology.eta_B` (magnitude); `condensate.chiral_magnetic` (the CME term).
Gudnason & Nitta, PRD 101, 065011 (2020) -- baryon # = vortex linking #.
Joyce & Shaposhnikov, PRL 79, 1193 (1997) -- chiral asymmetry -> net helicity.
Rivers & Volovik, PRL 127, 115702 (2021) -- the no-net-linking null beaten here.
"""
from __future__ import annotations

__all__ = [
    "SIGN_CONVENTION",
    "SIGN_CHAIN",
    "eta_B_sign",
]

#: The fixed sign convention: J_parent > 0 is right-handed parent spin, and
#: maps to matter (eta_B > 0). This is the one external input; everything
#: downstream is mechanism (see SIGN_CHAIN).
SIGN_CONVENTION: str = (
    "J_parent > 0  ==  right-handed parent-BH spin  ==  matter (eta_B > 0)."
)

#: The derived sign chain (each arrow is sign-preserving), validated dynamically.
SIGN_CHAIN: str = (
    "J_parent -> mu5 = d_t sigma -> (CME) net helicity sgn = +sgn(mu5) "
    "-> (Gudnason-Nitta) net vortex linking -> baryon# -> eta_B>0 = matter."
)


def eta_B_sign(j_parent_sign: int = +1) -> int:
    """Sign of the baryon asymmetry eta_B given the sign of the parent-BH spin.

    Each step of the chain ``J_parent -> mu5 -> net helicity -> net linking
    -> B -> eta_B`` is sign-preserving (CME: ``sgn<H> = +sgn(mu5)``;
    baryon # = vortex linking #), so

        sgn(eta_B) = sgn(J_parent).

    The magnitude is given separately by `eta_B.eta_B()` (= ``3 alpha^4 / 14``);
    this function returns ONLY the sign.

    Parameters
    ----------
    j_parent_sign : int, optional
        Sign of the parent black hole's spin J_parent. Positive (default) is the
        right-handed / matter convention (SIGN_CONVENTION); negative flips to
        antimatter. Zero (an unspun parent) gives 0 -- no sign selection,
        a parity-symmetric (achiral) universe with no net asymmetry.

    Returns
    -------
    int
        ``+1`` (matter), ``-1`` (antimatter), or ``0`` (no net asymmetry).
    """
    return (j_parent_sign > 0) - (j_parent_sign < 0)
