"""
Canonical CMB anisotropy axes used by both the cosmogenesis program
(Paper 22 -- parent black hole / accretion-disk fossil reading of the
CMB anomalies) and the Heron Experiment 11 directional layer
(Paper 19 substrate-anisotropy F1 falsifier).

The headline finding (cosmogenesis side, captured in
``analysis/nwt_cmb_coplanarity.py``):

  Hemispherical Power Asymmetry direction (= accretion-disk normal)
  is **89.0 degrees** from the Quadrupole-Octupole axis-of-evil
  (= Kerr spin axis): a 1.0 deg deviation from perfect
  perpendicularity, p = 0.017 (~2.4 sigma).  All 12 HPA x AoE
  measurement pairs span 79.6-96.7 degrees, every one within
  10.4 deg of perpendicular.

The substrate side (Paper 19): the SO(8) -> Spin(7) -> 3+1 holonomy
reduction picks an inertial-frame axis as a *structural*
consequence of the algebra.  The substrate-monism prediction is
that this axis coincides with the parent-BH spin axis, i.e. with
the AoE direction.  Heron Experiment 11 tests this via sidereal
modulation of K_7 stabilizer expectations on real superconducting
qubits.

This module collects:

  * Each measured CMB anomaly direction with its (l, b) in
    Galactic coordinates, an uncertainty estimate, and a citation.
    Multi-measurement collections are exposed as lists; consensus
    values (used as the predicted substrate axis and disk-plane
    normal) are exposed as named ``SkyAxis`` instances.
  * Galactic <-> Cartesian conversions, axis-equivalent angular
    separation, and a 3-vector great-circle fit.
  * A helper that returns all HPA x AoE pair separations -- the
    multi-measurement robustness check from
    ``nwt_cmb_coplanarity_pr4.py``.
  * An optional Galactic -> ICRS conversion (lazy astropy import)
    used by the Heron Experiment 11 directional layer to compare
    the predicted sky axis against a lab-frame qubit direction at
    a given time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import List, Tuple

import numpy as np


# =====================================================================
# Data class
# =====================================================================


@dataclass(frozen=True)
class SkyAxis:
    """A directed sky axis in Galactic coordinates.

    Sky axes are treated as axes (not vectors): the geometric
    pass/fail tests use ``angular_separation_undirected`` which
    folds antiparallel directions onto each other.  The signed
    ``l_deg`` / ``b_deg`` are kept for round-tripping to ICRS and
    for re-running the published analyses verbatim.
    """

    name: str
    l_deg: float
    b_deg: float
    sigma_deg: float = 0.0
    citation: str = ""
    notes: str = ""

    @property
    def cartesian(self) -> np.ndarray:
        """Unit Cartesian vector in Galactic coordinates."""
        return galactic_to_cartesian(self.l_deg, self.b_deg)

    def to_icrs(self) -> Tuple[float, float]:
        """Convert to ICRS (RA, Dec) in degrees via astropy."""
        return galactic_to_icrs(self.l_deg, self.b_deg)


# =====================================================================
# Coordinate conversions
# =====================================================================


def galactic_to_cartesian(l_deg: float, b_deg: float) -> np.ndarray:
    """Galactic (l, b) -> unit Cartesian (x, y, z)."""
    l_rad = np.radians(l_deg)
    b_rad = np.radians(b_deg)
    return np.array([
        np.cos(b_rad) * np.cos(l_rad),
        np.cos(b_rad) * np.sin(l_rad),
        np.sin(b_rad),
    ])


def cartesian_to_galactic(v: np.ndarray) -> Tuple[float, float]:
    """Unit Cartesian -> Galactic (l, b) in degrees, l in [0, 360)."""
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v)
    if v[2] < 0:
        v = -v
    b = float(np.degrees(np.arcsin(np.clip(v[2], -1.0, 1.0))))
    l = float(np.degrees(np.arctan2(v[1], v[0])) % 360.0)
    return l, b


def galactic_to_icrs(l_deg: float, b_deg: float) -> Tuple[float, float]:
    """Galactic (l, b) -> ICRS (RA, Dec) in degrees.

    Uses astropy if available.  Astropy is an optional dependency;
    install with ``pip install astropy`` if you need this helper.
    """
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
    except ImportError as e:
        raise ImportError(
            "galactic_to_icrs requires astropy; install with "
            "`pip install astropy`."
        ) from e
    c = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic").icrs
    return float(c.ra.deg), float(c.dec.deg)


# =====================================================================
# Geometric helpers
# =====================================================================


def angular_separation(a, b) -> float:
    """Signed angular separation (0..180 deg) between two SkyAxis or
    Cartesian vectors."""
    va = a.cartesian if isinstance(a, SkyAxis) else np.asarray(a, dtype=float)
    vb = b.cartesian if isinstance(b, SkyAxis) else np.asarray(b, dtype=float)
    va = va / np.linalg.norm(va)
    vb = vb / np.linalg.norm(vb)
    return float(np.degrees(np.arccos(np.clip(np.dot(va, vb), -1.0, 1.0))))


def angular_separation_undirected(a, b) -> float:
    """Axis-equivalent angular separation (0..90 deg) -- treats
    parallel and antiparallel as the same axis."""
    va = a.cartesian if isinstance(a, SkyAxis) else np.asarray(a, dtype=float)
    vb = b.cartesian if isinstance(b, SkyAxis) else np.asarray(b, dtype=float)
    va = va / np.linalg.norm(va)
    vb = vb / np.linalg.norm(vb)
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(va, vb)), 0.0, 1.0))))


def fit_great_circle(axes_or_vectors) -> dict:
    """Fit a great circle (plane through origin) to a set of unit
    vectors via the smallest-eigenvalue eigenvector of the scatter
    matrix.  Returns ``{normal, normal_lb, rms_deg, mean_deg,
    individual_deg}``.
    """
    vectors = []
    for x in axes_or_vectors:
        v = x.cartesian if isinstance(x, SkyAxis) else np.asarray(x, dtype=float)
        vectors.append(v / np.linalg.norm(v))
    M = np.zeros((3, 3))
    for v in vectors:
        M += np.outer(v, v)
    M /= len(vectors)
    _, eigenvectors = np.linalg.eigh(M)
    normal = eigenvectors[:, 0]
    distances = []
    for v in vectors:
        ang = np.degrees(np.arcsin(np.clip(np.dot(v, normal), -1.0, 1.0)))
        distances.append(abs(float(ang)))
    return {
        "normal": normal,
        "normal_lb": cartesian_to_galactic(normal),
        "rms_deg": float(np.sqrt(np.mean(np.asarray(distances) ** 2))),
        "mean_deg": float(np.mean(distances)),
        "individual_deg": distances,
    }


# =====================================================================
# Tabulated CMB anomaly directions
# =====================================================================
# Provenance:
#   AoE  = quadrupole-octupole alignment, "axis of evil"
#   HPA  = hemispherical power asymmetry direction
#   Cold = CMB Cold Spot center


AXIS_OF_EVIL_MEASUREMENTS: List[SkyAxis] = [
    SkyAxis(
        name="Land & Magueijo 2005",
        l_deg=250.0, b_deg=60.0, sigma_deg=2.0,
        citation="Land & Magueijo, Phys. Rev. Lett. 95, 071301 (2005)",
        notes="Original 'axis of evil'; quadrupole-octupole alignment.",
    ),
    SkyAxis(
        name="de Oliveira-Costa+ 2004",
        l_deg=240.0, b_deg=63.0, sigma_deg=2.0,
        citation="de Oliveira-Costa et al., Phys. Rev. D 69, 063516 (2004)",
        notes="Independent extraction of the preferred axis from WMAP.",
    ),
    SkyAxis(
        name="Planck SMICA Level-2 (Galasyn 2026)",
        l_deg=252.7, b_deg=58.4, sigma_deg=2.5,
        citation="claude-memory/scripts/cmb_level2_planck.py "
                 "(MC over 5000 trials)",
        notes="In-house re-extraction from Planck PR3 SMICA; 2.1 deg "
              "from Land & Magueijo 2005.",
    ),
]


HPA_MEASUREMENTS: List[SkyAxis] = [
    SkyAxis(
        name="Planck 2016",
        l_deg=227.0, b_deg=-27.0, sigma_deg=15.0,
        citation="Planck Collaboration 2016, A&A 594, A16",
        notes="Most widely cited HPA dipole direction.",
    ),
    SkyAxis(
        name="Hoftuft+ 2009",
        l_deg=224.0, b_deg=-22.0, sigma_deg=24.0,
        citation="Hoftuft et al., ApJ 699, 985 (2009)",
        notes="Refined direction using WMAP 5-year.",
    ),
    SkyAxis(
        name="Planck PR4 LVE 14deg",
        l_deg=225.0, b_deg=-33.0, sigma_deg=10.0,
        citation="Planck PR4 (2024 reanalysis)",
        notes="Local-variance estimator at 14 deg patch size.",
    ),
    SkyAxis(
        name="Planck PR4 SEVEM Nside64",
        l_deg=208.0, b_deg=-15.0, sigma_deg=10.0,
        citation="Planck PR4 SEVEM, Nside=64",
    ),
    SkyAxis(
        name="Planck PR4 SEVEM Nside2048",
        l_deg=205.0, b_deg=-20.0, sigma_deg=10.0,
        citation="Planck PR4 SEVEM, Nside=2048",
    ),
    SkyAxis(
        name="Planck PR4 MAP",
        l_deg=247.8, b_deg=-19.6, sigma_deg=10.0,
        citation="Planck PR4 MAP estimator",
    ),
]


# Consensus / preferred values.  These are the values the rest of
# the program uses as the predicted substrate axis (AoE) and disk-
# plane normal (HPA).  Both are within their respective
# multi-measurement spreads; AoE is the average of the two
# canonical published values, HPA is the Planck 2016 value (most
# widely cited).
AXIS_OF_EVIL_CONSENSUS: SkyAxis = SkyAxis(
    name="Axis of Evil (consensus)",
    l_deg=245.0, b_deg=61.5, sigma_deg=5.0,
    citation="Average of Land & Magueijo 2005 and de Oliveira-Costa+ 2004",
    notes="Predicted substrate-anisotropy axis (parent-BH spin axis). "
          "Used as the primary direction for Heron Experiment 11.",
)


HPA_CONSENSUS: SkyAxis = SkyAxis(
    name="Hemispherical Power Asymmetry (consensus)",
    l_deg=227.0, b_deg=-27.0, sigma_deg=15.0,
    citation="Planck Collaboration 2016, A&A 594, A16",
    notes="Disk-plane normal in the parent-BH cosmogenesis reading. "
          "89.0 deg from AXIS_OF_EVIL_CONSENSUS (1.0 deg from "
          "perpendicular, p = 0.017).",
)


# Single measurements (cold spot, dipole, references).
COLD_SPOT: SkyAxis = SkyAxis(
    name="CMB Cold Spot",
    l_deg=207.8, b_deg=-56.3, sigma_deg=2.0,
    citation="Vielva et al., ApJ 609, 22 (2004)",
    notes="59 deg from anti-AoE -- too far for a jet imprint. Lies "
          "0.6-9.0 deg from the great-circle plane defined by HPA + AoE; "
          "interpreted as a gap in the parent BH's accretion flow or "
          "the Eridanus supervoid.",
)


CMB_DIPOLE: SkyAxis = SkyAxis(
    name="CMB kinematic dipole",
    l_deg=264.0, b_deg=48.3, sigma_deg=0.1,
    citation="Planck Collaboration 2018, A&A 641, A1",
    notes="Direction of Local Group's motion w.r.t. the CMB rest frame.",
)


GALACTIC_NORTH: SkyAxis = SkyAxis(
    name="Galactic North Pole",
    l_deg=0.0, b_deg=90.0, sigma_deg=0.0,
    citation="By definition (IAU 1958)",
    notes="Approximately the Sgr A* spin axis (used as a structural "
          "reference in the angular-momentum-chain analysis).",
)


ECLIPTIC_NORTH: SkyAxis = SkyAxis(
    name="Ecliptic North Pole",
    l_deg=96.4, b_deg=29.8, sigma_deg=0.1,
    citation="J2000",
    notes="4.7 deg from perpendicular to AoE in the parent-BH model.",
)


BULK_FLOW: SkyAxis = SkyAxis(
    name="CosmicFlows-4 bulk flow (tangential)",
    l_deg=298.0, b_deg=-8.0, sigma_deg=10.0,
    citation="Whitford et al. 2023; CosmicFlows-4 reanalysis",
    notes="264 km/s tangential component, decomposed from 419 km/s "
          "total flow. Interpreted as fossil accretion-disk rotation "
          "in the parent-BH model. 8.2 deg from perpendicular to AoE.",
)


SHAPLEY: SkyAxis = SkyAxis(
    name="Shapley supercluster",
    l_deg=307.0, b_deg=30.0, sigma_deg=2.0,
    citation="Standard catalogue value",
    notes="Direction of dominant gravitational infall (radial "
          "component of bulk flow, 326 km/s).",
)


# =====================================================================
# Multi-measurement check: HPA x AoE pair separations
# =====================================================================


def hpa_aoe_pair_separations() -> List[Tuple[str, str, float]]:
    """Return all 12 HPA x AoE measurement pair separations in degrees.

    Reproduces the headline result from
    ``analysis/nwt_cmb_coplanarity_pr4.py``: every pair lies within
    10.4 deg of perpendicular (range 79.6 - 96.7 deg, median
    deviation from 90 deg = 5.3 deg).
    """
    out: List[Tuple[str, str, float]] = []
    aoe_pool = AXIS_OF_EVIL_MEASUREMENTS[:2]  # the two canonical
    for hpa_meas, aoe_meas in product(HPA_MEASUREMENTS, aoe_pool):
        sep = angular_separation(hpa_meas, aoe_meas)
        out.append((hpa_meas.name, aoe_meas.name, sep))
    return out
