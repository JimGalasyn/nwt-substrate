"""
QPO torus-eigenmode predictions for black-hole X-ray binaries.

Premise: the accretion disk of a Kerr BH is a substrate-condensate vortex
with a Macken-universal aspect ratio κ. Its quasi-periodic oscillations
(QPOs) are (p, q) torus eigenmodes:

    f(p, q) = √((p/κ)² + q²) × c / r_disk     [Hz]

with κ = KAPPA_MACKEN ≈ 9.844 (same aspect ratio as the electron, the
trefoil-carrier particle vortex, and the daughter universe's substrate
torus — universal substrate constant).

The classical "3:2 HFQPO twin peaks" then have the substrate
identification q_high = 3, q_low = 2 (with p = 0 for both), and:

    f(0, q) = q × c / r_disk
    f(0, 3) / f(0, 2) = 3/2 exactly.

Higher modes (p ≥ 1, q ≥ 4) predict a specific hierarchy testable in
RXTE / NICER archival data.

Reference: vortex-vision/analysis/qpo_torus_eigenmode_test.py (2026-05-15).
XTE J1550-564's mysterious 92 Hz HFQPO is the predicted (0, 1) mode at
91.9 Hz (deviation 0.2 %).
"""
from __future__ import annotations

import math

from ..isa import (
    ALPHA_NWT,
    KAPPA_MACKEN,
)


# Speed of light (SI) — only constant we need locally for QPO formulas.
C_LIGHT_M_S = 2.99792458e8


# ---------------------------------------------------------------------------
# Section 1 — Torus eigenmode frequency
# ---------------------------------------------------------------------------

def torus_eigenmode_freq(p: int, q: int,
                          r_disk_m: float,
                          kappa: float = KAPPA_MACKEN) -> float:
    """QPO frequency for (p, q) torus eigenmode at disk radius r_disk_m.

    f(p, q) = √((p/κ)² + q²) × c / r_disk     [Hz]

    Parameters
    ----------
    p, q : int
        Torus mode quantum numbers (poloidal, toroidal).
    r_disk_m : float
        Disk radius in metres.
    kappa : float, default KAPPA_MACKEN
        Disk aspect ratio. Default = κ_Macken ≈ 9.844 (universal).
    """
    f_norm = math.sqrt((p / kappa) ** 2 + q ** 2)
    return f_norm * C_LIGHT_M_S / r_disk_m


def predict_r_disk_from_qpo(f_obs_hz: float,
                            p: int, q: int,
                            kappa: float = KAPPA_MACKEN) -> float:
    """Invert torus eigenmode for disk radius given observed QPO frequency
    and a guessed (p, q) mode."""
    f_norm = math.sqrt((p / kappa) ** 2 + q ** 2)
    return f_norm * C_LIGHT_M_S / f_obs_hz


# ---------------------------------------------------------------------------
# Section 2 — 3:2 twin-peak identification
# ---------------------------------------------------------------------------

def twin_peak_ratio(q_low: int = 2, q_high: int = 3) -> float:
    """Exact 3:2 ratio for the canonical HFQPO twin peaks under p = 0.

    f(0, q_high) / f(0, q_low) = q_high / q_low = 3/2 (exact in the
    substrate-vortex picture; not a tuned resonance condition).
    """
    return q_high / q_low


# ---------------------------------------------------------------------------
# Section 3 — r_disk / r_g and substrate disk-vortex compactness
# ---------------------------------------------------------------------------

def r_disk_over_r_g(f_obs_hz: float, m_bh_kg: float,
                    p: int = 0, q: int = 2,
                    kappa: float = KAPPA_MACKEN,
                    G_si: float = 6.67430e-11) -> float:
    """Disk-to-gravitational-radius ratio implied by an observed QPO.

    r_g = G M / c². The ratio r_disk / r_g is approximately constant
    across Eddington thin-disk sources (~205-308 in the framework's
    QPO survey), supporting the substrate-vortex picture in which the
    disk's vortex aspect ratio is set by κ_Macken rather than by
    accretion-rate parameters.
    """
    r_g = G_si * m_bh_kg / C_LIGHT_M_S ** 2
    r_disk = predict_r_disk_from_qpo(f_obs_hz, p, q, kappa)
    return r_disk / r_g


# ---------------------------------------------------------------------------
# Section 4 — Mode-by-mode predicted frequencies
# ---------------------------------------------------------------------------

def qpo_mode_table(r_disk_m: float,
                    p_max: int = 3, q_max: int = 5,
                    kappa: float = KAPPA_MACKEN) -> dict:
    """Tabulate predicted QPO frequencies for (p, q) ∈ [0, p_max] × [0, q_max].

    Useful for archival searches in RXTE / NICER data: given a source
    with a known mass and a measured 3:2 twin-peak pair (which pins
    r_disk), the table predicts where additional QPO modes should
    appear.
    """
    table = {}
    for p in range(p_max + 1):
        for q in range(q_max + 1):
            if p == 0 and q == 0:
                continue
            table[(p, q)] = torus_eigenmode_freq(p, q, r_disk_m, kappa)
    return table
