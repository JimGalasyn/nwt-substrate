"""
nwt.cosmology
=============

The cosmological view of the substrate algebra.

Two threads converge here:

1. **Substrate anisotropy.** The substrate's holonomy reduction
   SO(8) -> Spin(7) -> 3+1 picks a preferred axis as a structural
   consequence of the algebra; the axis is not free.  Paper 19's
   F1 falsifier and the Heron Experiment 11 sidereal A/B/C probe
   compare against this predicted axis.

2. **Cosmogenesis.** Paper 22 (in preparation) reads the CMB
   anomalies as fossils of a parent black hole's accretion-disk
   geometry: the axis-of-evil = parent BH spin axis, the
   hemispherical power asymmetry direction = accretion-disk
   normal, the cold spot = a gap in the disk's accretion flow.

The two threads agree on a *single* sky-fixed axis: the substrate
algebra fixes form/magnitude, cosmogenesis fixes direction.

The :mod:`anisotropy_axes` submodule canonically tabulates the
measured directions of the CMB anomalies, with citations and
uncertainty ranges, plus the rotation matrices and great-circle
fits used by both the cosmogenesis program (Paper 22) and the
Heron Experiment 11 directional layer.

Quick start::

    from nwt_substrate.cosmology import anisotropy_axes as ax

    # Canonical predicted substrate axis (= parent BH spin axis)
    aoe = ax.AXIS_OF_EVIL_CONSENSUS
    print(aoe.l_deg, aoe.b_deg)            # ~245, +60

    # Disk-plane normal (HPA direction)
    hpa = ax.HPA_CONSENSUS
    print(ax.angular_separation(aoe, hpa)) # ~89 degrees

    # All 12 HPA x AoE pair separations (Paper 22 multi-measurement check)
    pairs = ax.hpa_aoe_pair_separations()
"""

from . import anisotropy_axes
from .anisotropy_axes import (
    SkyAxis,
    AXIS_OF_EVIL_CONSENSUS,
    HPA_CONSENSUS,
    COLD_SPOT,
    CMB_DIPOLE,
    GALACTIC_NORTH,
    ECLIPTIC_NORTH,
    BULK_FLOW,
    SHAPLEY,
    AXIS_OF_EVIL_MEASUREMENTS,
    HPA_MEASUREMENTS,
    galactic_to_cartesian,
    cartesian_to_galactic,
    galactic_to_icrs,
    angular_separation,
    angular_separation_undirected,
    fit_great_circle,
    hpa_aoe_pair_separations,
)


__all__ = [
    "anisotropy_axes",
    "SkyAxis",
    "AXIS_OF_EVIL_CONSENSUS",
    "HPA_CONSENSUS",
    "COLD_SPOT",
    "CMB_DIPOLE",
    "GALACTIC_NORTH",
    "ECLIPTIC_NORTH",
    "BULK_FLOW",
    "SHAPLEY",
    "AXIS_OF_EVIL_MEASUREMENTS",
    "HPA_MEASUREMENTS",
    "galactic_to_cartesian",
    "cartesian_to_galactic",
    "galactic_to_icrs",
    "angular_separation",
    "angular_separation_undirected",
    "fit_great_circle",
    "hpa_aoe_pair_separations",
]
