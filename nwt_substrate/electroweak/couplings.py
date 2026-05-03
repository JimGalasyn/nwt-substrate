"""
Z-boson neutral-current couplings g_V, g_A for SM fermions.

The Z couples to a fermion f as

    g_Z γ^μ (g_V^f - g_A^f γ⁵)

with

    g_V^f = T_3^f - 2 Q_f sin² θ_W
    g_A^f = T_3^f

where T_3 = +1/2 for up-type / neutrinos and -1/2 for down-type / charged
leptons, and Q is electric charge in units of |e|.

The chirality-projecting (γ⁵) structure is the canonical SU(2)_L × U(1)_Y
breaking pattern.  The substrate's internal SU(2) (Cl(0,7) bivector
triplet) is vector under Lorentz; the chiral pattern of the SM EW sector
is a SEPARATE structure that lives at the Lagrangian level (memory:
substrate-boson-count-so7, Walk-phase 4a).

Implementation: a small CouplingTable that returns (T_3, Q, g_V, g_A) for
each substrate fermion, plus catalog of the canonical SM 12 fermions.
"""

from __future__ import annotations

from dataclasses import dataclass

from .constants import SIN2_THETA_W


@dataclass(frozen=True)
class WeakCoupling:
    """Z-coupling parameters for a single fermion species."""
    name: str
    T_3: float           # weak isospin third component, +1/2 or -1/2
    Q: float             # electric charge in units of |e|
    g_V: float
    g_A: float
    n_color: int = 1

    @property
    def gV2_plus_gA2(self) -> float:
        """g_V² + g_A²  -- enters the Z partial width and total |M|²."""
        return self.g_V ** 2 + self.g_A ** 2

    def __str__(self) -> str:
        return (f"{self.name}:  T_3={self.T_3:+.1f}, Q={self.Q:+.3f},  "
                f"g_V={self.g_V:+.4f}, g_A={self.g_A:+.4f},  N_c={self.n_color}")


def _weak_coupling(name: str, T_3: float, Q: float, n_color: int = 1) -> WeakCoupling:
    """Construct a WeakCoupling from (T_3, Q) using the SM relation."""
    g_V = T_3 - 2.0 * Q * SIN2_THETA_W
    g_A = T_3
    return WeakCoupling(name, T_3, Q, g_V, g_A, n_color)


# Canonical SM fermion couplings (charged leptons, neutrinos, quarks)
SM_COUPLINGS = {
    # Charged leptons: T_3 = -1/2, Q = -1
    "e":   _weak_coupling("e",  -0.5, -1.0, 1),
    "mu":  _weak_coupling("mu", -0.5, -1.0, 1),
    "tau": _weak_coupling("tau",-0.5, -1.0, 1),
    # Neutrinos: T_3 = +1/2, Q = 0
    "nu_e":   _weak_coupling("nu_e",   +0.5, 0.0, 1),
    "nu_mu":  _weak_coupling("nu_mu",  +0.5, 0.0, 1),
    "nu_tau": _weak_coupling("nu_tau", +0.5, 0.0, 1),
    # Up-type quarks: T_3 = +1/2, Q = +2/3, N_c = 3
    "u": _weak_coupling("u", +0.5, +2.0 / 3.0, 3),
    "c": _weak_coupling("c", +0.5, +2.0 / 3.0, 3),
    "t": _weak_coupling("t", +0.5, +2.0 / 3.0, 3),
    # Down-type quarks: T_3 = -1/2, Q = -1/3, N_c = 3
    "d": _weak_coupling("d", -0.5, -1.0 / 3.0, 3),
    "s": _weak_coupling("s", -0.5, -1.0 / 3.0, 3),
    "b": _weak_coupling("b", -0.5, -1.0 / 3.0, 3),
}


def coupling(name: str) -> WeakCoupling:
    """Look up the Z neutral-current coupling for a fermion species."""
    n = name.lower()
    aliases = {"electron": "e", "muon": "mu", "up": "u", "down": "d",
               "charm": "c", "strange": "s", "bottom": "b", "top": "t"}
    n = aliases.get(n, n)
    if n not in SM_COUPLINGS:
        raise KeyError(
            f"Unknown fermion '{name}'.  Available: {sorted(SM_COUPLINGS.keys())}"
        )
    return SM_COUPLINGS[n]


def coupling_summary() -> str:
    """Pretty-print all SM Z-couplings."""
    lines = ["SM fermion Z-couplings (g_V = T_3 - 2 Q sin²θ_W, g_A = T_3):"]
    lines.append("")
    for name in ["e", "mu", "tau", "nu_e", "u", "c", "d", "s", "b"]:
        lines.append(f"  {SM_COUPLINGS[name]}")
    return "\n".join(lines)
