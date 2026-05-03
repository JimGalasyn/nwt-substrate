"""
Charge derivation: extended Gell-Mann-Nishijima formula.

  Q = I_3 + (B + S + C + B' + T) / 2

where I_3 = f/2 (framing-derived isospin projection) and
(B, S, C, B', T) are baryon number, strangeness, charm, bottom, top.

For NWT particles:
  - hadrons: Q derived from extended GM-N (38/38 verified, 2026-04-30)
  - leptons: Q from contact-geometry mechanism (Paper 7, separate rule)
"""

from __future__ import annotations


def gell_mann_nishijima(
    I3: float,
    B: int = 0,
    S: int = 0,
    C: int = 0,
    Bottom: int = 0,
    Top: int = 0,
) -> float:
    """
    Extended Gell-Mann-Nishijima charge formula.

    Q = I_3 + (B + S + C + B' + T) / 2

    Parameters
    ----------
    I3 : float
        Isospin projection (I_3 = f/2 in NWT framing).
    B : int
        Baryon number.
    S : int
        Strangeness.
    C : int
        Charm.
    Bottom : int
        Bottomness (B' in standard notation).
    Top : int
        Topness.

    Returns
    -------
    float
        Predicted electric charge in units of e.
    """
    return I3 + (B + S + C + Bottom + Top) / 2.0


def framing_from_I3(I3: float) -> int:
    """Inverse: framing f = 2 * I_3.  Returns rounded integer."""
    return int(round(2 * I3))
