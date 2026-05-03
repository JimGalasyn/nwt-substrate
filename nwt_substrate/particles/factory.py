"""
Particle factory: build a Particle by name from the compendium.

Usage:
    >>> import nwt_substrate as nwt
    >>> p = nwt.particle("p")
    >>> p.mass_pred
    937.24...
    >>> nwt.list_particles()
    ['e-', 'mu-', 'tau-', ...]
"""

from __future__ import annotations

from .particle import Particle
from .compendium import COMPENDIUM, list_names, get_compendium_entry


def particle(name: str) -> Particle:
    """
    Build a Particle from the compendium by name.

    Raises KeyError if name not in compendium.
    """
    entry = get_compendium_entry(name)
    if entry is None:
        raise KeyError(
            f"Unknown particle {name!r}. Available: {list_names()}"
        )
    return Particle(**entry)


def list_particles() -> list[str]:
    """Return all particle names available in the compendium."""
    return list_names()


def all_particles() -> list[Particle]:
    """Build Particle objects for every entry in the compendium."""
    return [Particle(**entry) for entry in COMPENDIUM]
