"""
Connected sum (#) as the substrate's molecular composition law.

Schubert (1949) proved every knot has a UNIQUE prime factorisation under
connected sum, which makes (Knots, #) a free commutative monoid generated
by prime knots -- the topological analog of the fundamental theorem of
arithmetic.  In NWT this is the substrate-algebraic operation that
combines two prime-carrier particles into a molecular bound state.

  K_A # K_B   = the carrier knot of the bound state of A and B
  m(K_A # K_B) = m(K_A) + m(K_B) + binding   (LO: binding = 0)

Tested empirically on near-threshold molecular states and the alpha-
cluster ladder (analysis/nwt_connected_sum_test.py and friends).
Headline numbers (Paper-6 end-to-end residuals, no binding correction):

  deuteron     -0.06%    Pc(4312)     +0.013%
  Pc(4457)     +0.37%    X(3872)      -0.53%
  Tcc(3875)    -0.62%    Pc(4440)     +0.75%

For strongly-bound exotic states (d*(2380), heavy nuclei) the # rule
fails by exactly the binding energy, isolating those as Hopf-linked
rather than molecular.  Hopf-link is the complementary substrate
operation, with energy scale Lambda_QCD = 313 MeV per crossing.

This module exposes # as a first-class API:

    >>> import nwt_substrate as nwt
    >>> p = nwt.particle("p")
    >>> n = nwt.particle("n")
    >>> d = nwt.compose(p, n, op="#")
    >>> d.mass_pred                              # ~ 1874.5 MeV
    >>> d.factors                                # [Particle("p"), Particle("n")]

Reference:
  memory:connected-sum-composition-law
  analysis/nwt_connected_sum_test.py
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from ..particles.particle import Particle
from ..particles.mass import paper6_mass_mev


# ---------------------------------------------------------------------------
# Composite particle
# ---------------------------------------------------------------------------


@dataclass
class CompositeParticle:
    """
    A bound state formed by connected sum (#) or Hopf-linking of two or more
    prime-carrier particles.  Carries enough information to recover the
    constituents (Schubert prime factorisation for #) and predict total mass.

    Fields:
      name        Display name (e.g. "X(3872)", "Pc(4312)")
      factors     List of constituent Particle objects (the "primes")
      op          Composition operation: "#" (connected sum, molecular) or
                  "hopf" (Hopf-linked, strong binding)
      binding     Binding energy correction in MeV (negative = bound).
                  Defaults to 0 for LO #.  For Hopf-linked states use
                  - LAMBDA_QCD * crossings.
      m_obs       Observed mass in MeV when known (PDG).
      notes       Free-form notes.
    """

    name: str
    factors: list[Particle]
    op: str = "#"
    binding: float = 0.0
    m_obs: Optional[float] = None
    notes: str = ""

    def __post_init__(self) -> None:
        if self.op not in ("#", "hopf"):
            raise ValueError(f"unknown composition op: {self.op!r}; "
                              "use '#' or 'hopf'")
        if len(self.factors) < 2:
            raise ValueError("composite needs at least two factor particles")

    # ---------- Mass prediction ----------

    @property
    def mass_constituent_sum_pred(self) -> Optional[float]:
        """Sum of paper6_mass over all factor (p,q,m,n_q)."""
        masses: list[float] = []
        for f in self.factors:
            mp = paper6_mass_mev(f.p, f.q, f.m, f.n_q)
            if mp is None:
                return None
            masses.append(mp)
        return sum(masses)

    @property
    def mass_constituent_sum_obs(self) -> Optional[float]:
        """Sum of observed (PDG) masses over factors, when all are known."""
        if any(f.m_obs is None for f in self.factors):
            return None
        return sum(f.m_obs for f in self.factors)  # type: ignore[misc]

    @property
    def mass_pred(self) -> Optional[float]:
        """
        LO + binding mass prediction:  sum_i m_paper6(K_i) + binding.

        For a connected-sum molecular state with binding = 0 this is the
        pure '# is mass-additive' prediction; supply a non-zero `binding`
        (negative, MeV) to capture the explicit binding correction.

        For Hopf-linked states, set binding = - LAMBDA_QCD * N_crossings.
        """
        s = self.mass_constituent_sum_pred
        if s is None:
            return None
        return s + self.binding

    @property
    def mass_residual(self) -> Optional[float]:
        """(predicted - observed) / observed  in percent."""
        mp, mo = self.mass_pred, self.m_obs
        if mp is None or mo is None or mo == 0:
            return None
        return (mp - mo) / mo * 100.0

    # ---------- Substrate algebra ----------

    @property
    def n_factors(self) -> int:
        return len(self.factors)

    @property
    def total_genus(self) -> int:
        """Genus is additive under #: g(K_1 # K_2) = g_1 + g_2 (Schubert)."""
        return sum(f.genus for f in self.factors)

    @property
    def total_crossings(self) -> int:
        """Crossing number is additive for alternating knots (#)."""
        return sum(f.crossings for f in self.factors)

    @property
    def total_baryon_number(self) -> int:
        """Baryon number is additive."""
        return sum(f.B_substrate if f.B_pdg is None else f.B_pdg
                    for f in self.factors)

    # ---------- Display ----------

    def __repr__(self) -> str:
        sym = " # " if self.op == "#" else " <hopf> "
        return f"CompositeParticle({self.name!r} = {sym.join(f.name for f in self.factors)})"

    def summary(self) -> str:
        sym = " # " if self.op == "#" else " <hopf> "
        lines = [f"Composite: {self.name}"]
        lines.append(f"  Operation:    {self.op} ({'molecular' if self.op == '#' else 'Hopf-linked'})")
        lines.append(f"  Factors:      {sym.join(f.name for f in self.factors)}")
        lines.append(f"  Genus total:  {self.total_genus}")
        lines.append(f"  Crossings:    {self.total_crossings}")
        if self.mass_pred is not None:
            lines.append(f"  Mass pred:    {self.mass_pred:.2f} MeV  "
                          f"(constituent sum {self.mass_constituent_sum_pred:.2f} "
                          f"+ binding {self.binding:+.2f})")
        if self.m_obs is not None:
            lines.append(f"  Mass obs:     {self.m_obs:.2f} MeV")
        if self.mass_residual is not None:
            lines.append(f"  Residual:     {self.mass_residual:+.3f}%")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Factory function (exposed at top level as nwt.compose)
# ---------------------------------------------------------------------------


def compose(*factors: Particle,
            op: str = "#",
            name: Optional[str] = None,
            binding: float = 0.0,
            m_obs: Optional[float] = None) -> CompositeParticle:
    """
    Combine two or more particles via the substrate's composition law.

    Parameters
    ----------
    *factors : Particle
        The constituent particles (prime carriers).  Two or more.
    op : {"#", "hopf"}
        Composition operation:
          "#"    -- connected sum, molecular bound state, mass-additive at LO.
          "hopf" -- Hopf-linked, strongly-bound exotic / nuclear binding.
    name : str, optional
        Display name; defaults to "(A # B)" or "(A <hopf> B)".
    binding : float
        Binding-energy correction in MeV (negative = bound).  Default 0
        is the LO # prediction; for Hopf-linked states pass
        binding = -LAMBDA_QCD * n_crossings.
    m_obs : float, optional
        Observed mass in MeV when known.

    Returns
    -------
    CompositeParticle

    Examples
    --------
    Deuteron from p # n:

        >>> p, n = nwt.particle("p"), nwt.particle("n")
        >>> d = nwt.compose(p, n, op="#", name="d", m_obs=1875.61)
        >>> d.mass_pred                  # ~ 1874.48 MeV
        >>> d.mass_residual              # ~ -0.06 %

    A Hopf-linked dibaryon (binding from Lambda_QCD):

        >>> from nwt_substrate.compositions import LAMBDA_QCD_MEV
        >>> dstar = nwt.compose(delta, delta, op="hopf",
        ...                      binding=-LAMBDA_QCD_MEV * 0.27)
    """
    if len(factors) < 2:
        raise ValueError("compose() requires at least two factors")
    if op == "#":
        sym = " # "
    elif op == "hopf":
        sym = " <hopf> "
    else:
        raise ValueError(f"unknown op {op!r}; use '#' or 'hopf'")
    if name is None:
        name = "(" + sym.join(f.name for f in factors) + ")"
    return CompositeParticle(
        name=name,
        factors=list(factors),
        op=op,
        binding=binding,
        m_obs=m_obs,
    )


def decompose(composite: CompositeParticle) -> list[Particle]:
    """
    Schubert prime factorisation: recover the prime constituents of a
    composite particle.  For #, this is exact (Schubert 1949).  For
    Hopf-link, the constituent list is conventional (the input factors).
    """
    return list(composite.factors)
