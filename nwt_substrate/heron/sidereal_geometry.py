"""
Geometry layer for Heron Experiment 11 -- the sidereal A/B/C
substrate-anisotropy probe.

The bare A/B/C drift correction in
``analysis/nwt_sidereal_AB_C_demo.py`` answers a *weak* question:
"is there sidereal modulation in the K_7 stabilizer expectations?"

For the bold Paper-19 claim ("the substrate axis coincides with the
parent-BH spin axis = the CMB axis-of-evil") we need the *strong*
question:  "does the modulation point at the predicted sky direction?"

This module supplies the geometry needed for the strong question:

  * **Observatory** records the lat/long/height of each Heron host
    site (Yorktown NY for ibm_marrakesh, Ehningen DE for ibm_aachen).
  * **Time -> Earth orientation** converts a UTC unix time to the
    rotation matrix that takes lab-frame Cartesian -> ICRS Cartesian
    at that observatory.  Uses astropy under the hood.
  * **Predicted sigma_v pattern** forward-models what each K_7
    stabilizer's drift-corrected sidereal signal *should* look like
    if the substrate axis points at a given ICRS direction (default:
    Axis of Evil).  Slots A and C share inertial orientation (Earth
    has rotated 360 sidereal degrees between them), so the prediction
    reduces to a 180-deg rotation of qubit positions in the inertial
    frame between slots A and B.
  * **Directional match score** compares observed sigma_v to the
    predicted shape via Pearson correlation and the cosine-similarity
    of the inferred sky axis.  This is the Paper-19 F1 pass/fail.

The chip-layout convention defaults to "chip flat in the lab, chip-x
= lab-east, chip-y = lab-north, chip-z = lab-up": clearly an
approximation of IBM's actual cryostat orientation, but the right
convention for a single-site dry run.  The substrate-axis prediction
amplitude does not depend on this convention, only the predicted
*pattern* across vertices does.  Pre-registration documents (see
``analysis/EXP11_PREREG.md``) MUST cite the chip orientation used.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import math
import numpy as np

from ..cosmology.anisotropy_axes import (
    AXIS_OF_EVIL_CONSENSUS,
    HPA_CONSENSUS,
    SkyAxis,
    galactic_to_cartesian,
)


# =====================================================================
# Observatory locations
# =====================================================================
# Sites are HIGH or LOW confidence per `analysis/heron_database/`
# location inference (commit f9d3338).  The Heron Experiment 11
# protocol commits to running at Yorktown via ibm_marrakesh; the
# Ehningen entry is the EU baseline for the 2-site coincidence run.


@dataclass(frozen=True)
class Observatory:
    """An IBM Quantum host site as a geodetic observatory."""
    name: str
    lat_deg: float        # WGS84 latitude (north positive)
    lon_deg: float        # WGS84 longitude (east positive)
    height_m: float = 0.0
    notes: str = ""


YORKTOWN = Observatory(
    name="IBM Yorktown Heights",
    lat_deg=41.2089, lon_deg=-73.8093, height_m=140.0,
    notes="Hosts ibm_marrakesh (HIGH confidence per Heron Exp 10 "
          "muon-decay provenance).",
)

EHNINGEN = Observatory(
    name="IBM Ehningen",
    lat_deg=48.6629, lon_deg=8.9420, height_m=478.0,
    notes="Hosts ibm_aachen (HIGH confidence per eu-de instance "
          "binding).",
)


# =====================================================================
# Earth orientation: lab Cartesian <-> ICRS Cartesian at time t
# =====================================================================


def lab_to_icrs_matrix(unix_time: float, observatory: Observatory) -> np.ndarray:
    """Return a 3x3 rotation matrix R such that R @ lab_vec = icrs_vec.

    Lab-frame convention: x = local east, y = local north, z = local
    up (zenith).  ICRS-Cartesian convention: x = (RA=0, Dec=0),
    y = (RA=90, Dec=0), z = (Dec=90).

    Implementation uses astropy's AltAz <-> ICRS transform.  Three
    canonical lab unit vectors (east, north, up) are mapped to
    AltAz, transformed to ICRS, and stacked column-wise.
    """
    try:
        from astropy.coordinates import AltAz, EarthLocation, ICRS, SkyCoord
        from astropy.time import Time
        import astropy.units as u
    except ImportError as e:
        raise ImportError(
            "lab_to_icrs_matrix requires astropy; install with "
            "`pip install astropy`."
        ) from e

    loc = EarthLocation(
        lat=observatory.lat_deg * u.deg,
        lon=observatory.lon_deg * u.deg,
        height=observatory.height_m * u.m,
    )
    t = Time(unix_time, format="unix", scale="utc")
    altaz_frame = AltAz(obstime=t, location=loc)

    # Three canonical lab basis vectors as AltAz directions.
    # AltAz az is measured east of north; alt is angle above horizon.
    #   lab x (east) -> az=90, alt=0
    #   lab y (north) -> az=0, alt=0
    #   lab z (up)    -> az=0, alt=90
    basis_altaz = SkyCoord(
        az=[90.0, 0.0, 0.0] * u.deg,
        alt=[0.0, 0.0, 90.0] * u.deg,
        frame=altaz_frame,
    )
    icrs = basis_altaz.icrs
    # Each column of R is a lab basis vector expressed in ICRS Cart.
    R = np.column_stack([
        np.array([
            np.cos(icrs[i].dec.rad) * np.cos(icrs[i].ra.rad),
            np.cos(icrs[i].dec.rad) * np.sin(icrs[i].ra.rad),
            np.sin(icrs[i].dec.rad),
        ])
        for i in range(3)
    ])
    return R


def lab_to_icrs(lab_vec: np.ndarray, unix_time: float,
                observatory: Observatory) -> np.ndarray:
    """Convert a lab-frame Cartesian unit vector to ICRS Cartesian
    at the given time and observatory."""
    R = lab_to_icrs_matrix(unix_time, observatory)
    out = R @ np.asarray(lab_vec, dtype=float)
    n = np.linalg.norm(out)
    return out / n if n > 0 else out


# =====================================================================
# Default K_7 lab-frame embedding
# =====================================================================
# K_7 is the complete graph -- it has no preferred geometric
# embedding.  When transpiled to a Heron, the 7 logical qubits land
# at *some* 7 physical-qubit positions on the heavy-hex chip, with
# fixed coupling-graph distances.  This module's default places the
# 7 qubits at the vertices of a regular hexagon plus a center, in
# the chip plane, with the chip plane assumed flat in the lab.
# Override via the ``qubit_lab_positions`` argument to
# ``predicted_sigma_pattern`` if you have the actual transpiled
# layout.


def default_k7_lab_positions(radius_m: float = 0.005) -> np.ndarray:
    """Return a (7, 3) array of unit lab-frame positions for the
    K_7 vertices, in the regular-hexagon-plus-center layout.

    Vertex 0 sits at the center; vertices 1-6 sit at hexagon corners
    rotated 0, 60, 120, 180, 240, 300 degrees in the chip plane,
    starting from chip-east.

    Returns the 7 *positions* normalised to unit vectors.  The
    radius_m argument is informational only -- the directional
    signature depends on direction, not on chip scale.
    """
    pts = np.zeros((7, 3))
    pts[0] = (1e-9, 0.0, 1.0)  # near-zenith for the center qubit
    for i in range(6):
        ang = i * (math.pi / 3.0)
        pts[1 + i] = (radius_m * math.cos(ang),
                      radius_m * math.sin(ang),
                      1.0)
    # Normalise each row.
    for i in range(7):
        pts[i] /= np.linalg.norm(pts[i])
    return pts


# =====================================================================
# Forward model: predicted sigma_v under the substrate-axis hypothesis
# =====================================================================


def predicted_sigma_pattern(
    t_A_unix: float,
    t_B_unix: float,
    t_C_unix: float,
    observatory: Observatory,
    qubit_lab_positions: Optional[np.ndarray] = None,
    predicted_axis: SkyAxis = AXIS_OF_EVIL_CONSENSUS,
    amplitude: float = 1.0,
) -> dict:
    """Forward-model the predicted sigma_v vector under the
    substrate-axis hypothesis.

    Coupling model (simplest non-trivial):
        At time t the inertial-frame orientation of qubit v is
        R(t) p_v, where R(t) is the lab->ICRS rotation matrix.
        The substrate-coupling at qubit v at time t is
            C_v(t) = ((R(t) p_v) . s_hat)^2
        where s_hat is the predicted axis unit vector in ICRS.

        The drift-corrected sidereal signal is then
            sigma_v_pred = amplitude * [
                (B - A) - (C - A) / 2
            ]_v
        evaluated on the C_v(t) coupling.

    The returned dict carries the predicted sigma_v, the cosines
    that drive it, and the predicted axis in ICRS for downstream
    inspection.
    """
    if qubit_lab_positions is None:
        qubit_lab_positions = default_k7_lab_positions()
    qubit_lab_positions = np.asarray(qubit_lab_positions, dtype=float)
    if qubit_lab_positions.shape != (7, 3):
        raise ValueError(
            f"qubit_lab_positions must be shape (7, 3); got "
            f"{qubit_lab_positions.shape}"
        )

    # Predicted axis in galactic Cartesian, then to ICRS Cartesian
    # via the SkyAxis.to_icrs() helper.
    ra_deg, dec_deg = predicted_axis.to_icrs()
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    s_hat_icrs = np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    ])

    # Compute lab->ICRS at each slot.
    R_A = lab_to_icrs_matrix(t_A_unix, observatory)
    R_B = lab_to_icrs_matrix(t_B_unix, observatory)
    R_C = lab_to_icrs_matrix(t_C_unix, observatory)

    # Couplings C_v(t) = ((R(t) p_v) . s_hat)^2
    qubit_lab_unit = qubit_lab_positions / np.linalg.norm(
        qubit_lab_positions, axis=1, keepdims=True
    )
    cosA = np.array([s_hat_icrs @ (R_A @ p) for p in qubit_lab_unit])
    cosB = np.array([s_hat_icrs @ (R_B @ p) for p in qubit_lab_unit])
    cosC = np.array([s_hat_icrs @ (R_C @ p) for p in qubit_lab_unit])
    C_A = cosA ** 2
    C_B = cosB ** 2
    C_C = cosC ** 2

    sigma_pred = amplitude * ((C_B - C_A) - (C_C - C_A) / 2.0)

    return {
        "sigma_v_pred": sigma_pred,
        "C_A": C_A,
        "C_B": C_B,
        "C_C": C_C,
        "predicted_axis_lb": (predicted_axis.l_deg, predicted_axis.b_deg),
        "predicted_axis_radec": (ra_deg, dec_deg),
        "predicted_axis_name": predicted_axis.name,
    }


# =====================================================================
# Directional match score: observed sigma_v vs predicted pattern
# =====================================================================


def directional_match_score(
    observed_sigma_v: Sequence[float],
    predicted_sigma_v: Sequence[float],
    observed_std: Optional[Sequence[float]] = None,
) -> dict:
    """Compare observed sigma_v to a forward-modelled predicted
    sigma_v under the substrate-axis hypothesis.

    Three statistics:

      * ``pearson_r``: correlation between observed and predicted
        across the 7 vertices.  Range [-1, 1]; +1 means the pattern
        matches exactly up to a positive scale.
      * ``cosine_similarity``: cosine of the angle between the two
        7-vectors.  Sign-sensitive; +1 = pointing at the predicted
        axis, -1 = pointing at the antipode, 0 = uncorrelated.
      * ``best_amplitude_chi2``: the chi^2 of (observed - A * predicted)
        minimised over the amplitude scalar A, weighted by
        ``observed_std`` if provided.  Lower is better; 0 = perfect
        match after rescaling.
    """
    obs = np.asarray(observed_sigma_v, dtype=float)
    pred = np.asarray(predicted_sigma_v, dtype=float)
    if obs.shape != pred.shape:
        raise ValueError("observed and predicted sigma_v must have "
                         "the same shape")

    # Pearson correlation (mean-centred).
    obs_c = obs - obs.mean()
    pred_c = pred - pred.mean()
    denom = np.sqrt((obs_c @ obs_c) * (pred_c @ pred_c))
    pearson_r = float(obs_c @ pred_c / denom) if denom > 0 else 0.0

    # Cosine similarity (uncentred).
    n_obs = np.linalg.norm(obs)
    n_pred = np.linalg.norm(pred)
    cos_sim = float(obs @ pred / (n_obs * n_pred)) if n_obs > 0 and n_pred > 0 else 0.0

    # Best-amplitude chi^2.
    if observed_std is None:
        sigma = np.ones_like(obs)
    else:
        sigma = np.asarray(observed_std, dtype=float)
        sigma = np.where(sigma > 0, sigma, 1.0)

    w = 1.0 / sigma ** 2
    A_hat_num = float(np.sum(w * obs * pred))
    A_hat_den = float(np.sum(w * pred ** 2))
    A_hat = A_hat_num / A_hat_den if A_hat_den > 0 else 0.0
    residuals = obs - A_hat * pred
    chi2 = float(np.sum(w * residuals ** 2))

    return {
        "pearson_r": pearson_r,
        "cosine_similarity": cos_sim,
        "best_amplitude": A_hat,
        "best_amplitude_chi2": chi2,
        "dof": int(obs.size - 1),  # one parameter (amplitude) fitted
    }


__all__ = [
    "Observatory",
    "YORKTOWN",
    "EHNINGEN",
    "lab_to_icrs_matrix",
    "lab_to_icrs",
    "default_k7_lab_positions",
    "predicted_sigma_pattern",
    "directional_match_score",
    "lst_hours",
    "next_lst_match_unix",
    "schedule_triplet_at_lst",
]


# =====================================================================
# Local Sidereal Time anchoring (multi-day folding)
# =====================================================================
# LST at observer = (Greenwich Sidereal Time) + (longitude east).  Two
# unix times t1, t2 at the same observatory have the same LST iff
# (t2 - t1) is an integer number of sidereal days.  The drift-model
# validation showed that single-triplet per-vertex sigma drift is
# ~1.3% on pittsburgh (after backend selection); multi-day folding at
# matched LST averages this down by sqrt(N) days.  These helpers let
# the demo schedule its three slots at a fixed target LST so that
# triplets across many days can be coherently combined.


def lst_hours(unix_time: float, observatory: Observatory) -> float:
    """Local Apparent Sidereal Time in hours [0, 24) at the
    observatory at the given UTC unix time."""
    try:
        from astropy.time import Time
        from astropy.coordinates import EarthLocation
        import astropy.units as u
    except ImportError as e:
        raise ImportError(
            "lst_hours requires astropy; install with `pip install astropy`."
        ) from e

    loc = EarthLocation(
        lat=observatory.lat_deg * u.deg,
        lon=observatory.lon_deg * u.deg,
        height=observatory.height_m * u.m,
    )
    t = Time(unix_time, format="unix", scale="utc")
    return float(t.sidereal_time("apparent", longitude=loc.lon).hour)


def next_lst_match_unix(
    target_lst_hours: float,
    observatory: Observatory,
    after_unix: Optional[float] = None,
    tol_seconds: float = 1.0,
    max_iter: int = 4,
) -> float:
    """Return the next unix time >= ``after_unix`` at which LST
    equals ``target_lst_hours`` at the observatory.

    Uses the exact astropy LST as an oracle and Newton-style
    correction at the sidereal rate (1.00273790935 solar sec per
    sidereal sec).  Converges in ~2 iterations to sub-second
    precision.
    """
    import time as _time
    if after_unix is None:
        after_unix = _time.time()

    SIDEREAL_RATE = 1.00273790935  # solar seconds per sidereal second
    SIDEREAL_DAY_SOLAR_S = 86164.0905  # one sidereal day in solar seconds

    t = float(after_unix)
    for _ in range(max_iter):
        lst_now = lst_hours(t, observatory)
        delta_lst = (target_lst_hours - lst_now) % 24.0
        delta_solar = delta_lst * 3600.0 / SIDEREAL_RATE
        if delta_solar > 0:
            t = t + delta_solar
        # If delta_lst is ~0 we've already arrived; check tolerance.
        residual = (target_lst_hours - lst_hours(t, observatory) + 12) % 24 - 12
        if abs(residual) * 3600.0 < tol_seconds:
            return t
    return t


def schedule_triplet_at_lst(
    target_lst_hours: float,
    observatory: Observatory,
    after_unix: Optional[float] = None,
) -> dict:
    """Plan an A/B/C triplet anchored at a target LST.

    Slot A fires at the next time LST hits ``target_lst_hours``;
    slot B at +12 h sidereal (same lab orientation rotated 180 deg
    in the inertial frame); slot C at +24 h sidereal (same lab and
    inertial orientation as A, used for the linear drift baseline).

    Returns a dict with the three unix timestamps and the actual
    LSTs (which should be target / target+12 / target -- that
    last one because LST is mod 24, but in solar wallclock the C
    slot is 23h 56m later, NOT 24 h, so its LST is ALSO the target).
    """
    SIDEREAL_DAY_S = 86164.0905
    t_A = next_lst_match_unix(target_lst_hours, observatory,
                                after_unix=after_unix)
    t_B = t_A + SIDEREAL_DAY_S / 2
    t_C = t_A + SIDEREAL_DAY_S
    return {
        "target_lst_hours": float(target_lst_hours),
        "observatory": observatory.name,
        "A_unix": t_A,
        "B_unix": t_B,
        "C_unix": t_C,
        "A_lst": lst_hours(t_A, observatory),
        "B_lst": lst_hours(t_B, observatory),
        "C_lst": lst_hours(t_C, observatory),
    }
