"""
Particle class: central data structure for substrate-algebraic particle
representation.

A Particle is specified by 5 substrate quantum numbers + L1 sector:
  (p, q, m, n_q, f) + sector

These determine all physical observables (mass, J, Q, B, S, I) via:
  - mass:   Paper 6 formula
  - J:      L1 sector (scalar=0, spinor=1/2 or 3/2, vector=1)
  - Q:      Gell-Mann-Nishijima Q = f/2 + (B+S+C)/2
  - B:      n_q mod 2 (rough rule, can be overridden)
  - I_3:    f / 2

Physical observables (mass, charge, etc.) can also be supplied directly
when known from PDG, with the framework providing consistency checks.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from ..topology.torus_knots import (
    knot_family,
    seifert_genus,
    crossing_number,
    hopf_charge,
    beta_phase,
    carrier_for_n_q,
)


@dataclass
class Particle:
    """
    Substrate-algebraic representation of a particle.

    Mandatory fields
    ----------------
    name : str
        Display name (e.g., "p", "X(3872)", "Pc(4312)").
    p, q : int
        Toroidal and poloidal winding numbers of the carrier knot.
    m : int
        Phase quantum number (sets beta = sqrt(m^2/p^2 - 1)).
    n_q : int
        Carrier knot crossing number = dim of 2I irrep, in {0, 1, 2, 3, 4, 5, 6}.

    Optional fields
    ---------------
    f : int
        Framing number, gives I_3 = f/2 via Paper 7 mechanism.
    sector : str
        L1 sector: "scalar", "spinor", or "vector".  Determines spin J.
    m_obs : float | None
        Observed mass in MeV, when known (PDG).
    Q_pdg : int | float | None
        Observed electric charge.
    B_pdg : int | None
        Observed baryon number.
    S_pdg : int | None
        Observed strangeness.
    C_pdg : int | None
        Observed charm.
    I_pdg : float | None
        Observed total isospin.
    notes : str
        Free-form notes.
    """

    name: str
    p: int
    q: int
    m: int
    n_q: int
    f: int = 0
    sector: str = "spinor"
    m_obs: Optional[float] = None
    Q_pdg: Optional[float] = None
    B_pdg: Optional[int] = None
    S_pdg: Optional[int] = None
    C_pdg: Optional[int] = None
    I_pdg: Optional[float] = None
    notes: str = ""

    # ---------- Derived properties ----------

    @property
    def J(self) -> float:
        """
        Spin from L1 sector.

        scalar -> J=0, vector -> J=1, spinor -> J=1/2 (or 3/2 for higher modes).
        """
        return {
            "scalar": 0.0,
            "vector": 1.0,
            "spinor": 0.5,
        }.get(self.sector, 0.5)

    @property
    def I3(self) -> float:
        """Isospin projection from framing: I_3 = f/2."""
        return self.f / 2.0

    @property
    def carrier(self) -> str:
        """Named carrier knot family from n_q."""
        return carrier_for_n_q(self.n_q)

    @property
    def torus_knot_family(self) -> str:
        """Named torus-knot family from (p, q)."""
        return knot_family(self.p, self.q)

    @property
    def genus(self) -> int:
        """Seifert genus of the (p, q) torus knot."""
        return seifert_genus(self.p, self.q)

    @property
    def crossings(self) -> int:
        """Crossing number of the (p, q) projection."""
        return crossing_number(self.p, self.q)

    @property
    def Q_H(self) -> int:
        """L3 Hopf charge Q_H = p * m."""
        return hopf_charge(self.p, self.m)

    @property
    def beta(self) -> Optional[float]:
        """Aspect-ratio parameter sqrt(m^2/p^2 - 1)."""
        return beta_phase(self.p, self.m)

    @property
    def B_substrate(self) -> int:
        """Baryon number from substrate (n_q mod 2 rule)."""
        return self.n_q % 2

    @property
    def mass_pred(self) -> Optional[float]:
        """
        Predicted mass in MeV via Paper 6 formula.

        Returns None if (p, q, m, n_q) is invalid (e.g., m <= p).
        """
        from .mass import paper6_mass_mev
        return paper6_mass_mev(self.p, self.q, self.m, self.n_q)

    @property
    def Q_pred(self) -> float:
        """
        Charge from extended Gell-Mann-Nishijima:
          Q = I_3 + (B + S + C + B' + T) / 2

        Uses PDG values when available; falls back to substrate-derived B.
        """
        B = self.B_pdg if self.B_pdg is not None else self.B_substrate
        S = self.S_pdg if self.S_pdg is not None else 0
        C = self.C_pdg if self.C_pdg is not None else 0
        return self.I3 + (B + S + C) / 2.0

    @property
    def mass_residual(self) -> Optional[float]:
        """
        Relative residual (predicted - observed) / observed, in percent.

        Returns None if either prediction or observation is missing.
        """
        if self.mass_pred is None or self.m_obs is None:
            return None
        return (self.mass_pred - self.m_obs) / self.m_obs * 100.0

    # ---------- Display ----------

    def __repr__(self) -> str:
        m_str = f", m_obs={self.m_obs}" if self.m_obs is not None else ""
        return (f"Particle({self.name!r}, p={self.p}, q={self.q}, "
                f"m={self.m}, n_q={self.n_q}, f={self.f}, "
                f"sector={self.sector!r}{m_str})")

    def summary(self) -> str:
        """Multi-line description for human reading."""
        lines = [f"Particle: {self.name}"]
        lines.append(f"  Substrate: (p, q, m, n_q, f) = ({self.p}, {self.q}, "
                     f"{self.m}, {self.n_q}, {self.f})")
        lines.append(f"  L1 sector: {self.sector}")
        lines.append(f"  Carrier knot: {self.carrier}")
        lines.append(f"  Torus-knot family: {self.torus_knot_family}, "
                     f"genus={self.genus}")
        lines.append(f"  Q_H = p*m = {self.Q_H}")
        lines.append(f"  Spin J = {self.J}, I_3 = {self.I3}, "
                     f"Q_predicted = {self.Q_pred:+.2f}")
        if self.mass_pred is not None:
            lines.append(f"  Mass predicted: {self.mass_pred:.2f} MeV")
            if self.m_obs is not None:
                lines.append(f"  Mass observed:  {self.m_obs:.2f} MeV "
                             f"(residual {self.mass_residual:+.2f}%)")
        if self.notes:
            lines.append(f"  Notes: {self.notes}")
        return "\n".join(lines)
