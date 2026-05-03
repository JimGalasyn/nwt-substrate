"""
Quantum field types for nwt.qft.

A "field" in QFT is a substrate-algebra primitive viewed through QFT
vocabulary.  Each Field instance carries an explicit pointer (via the
`substrate_primitive` attribute) to the underlying substrate object that
implements it.  The thesis: QFT is one view of the substrate algebra.

Field types:
  ScalarField      -- complex scalar (φ).  Substrate: 1-d C-valued field.
  DiracSpinor      -- 4-component Dirac (ψ).  Substrate: Cl(1,3) 8-d spinor
                       built from Cl(0,7) octonion left-multiplication.
  MajoranaSpinor   -- self-conjugate Dirac.  Substrate: real 8-d Cl(1,3).
  GaugeField       -- A^a_μ for a Lie-algebra g.  Substrate: g generators
                       (u(1), su(2), su(3), so(7), ...) and propagator.
  Ghost            -- Faddeev-Popov c, c̄.  Substrate: scalar + Grassmann sign.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Optional


@dataclass(frozen=True)
class ScalarField:
    """
    Complex scalar field φ(x) of mass m.

    QFT view:  L_kin = (∂_μ φ)†(∂^μ φ) − m² φ†φ
    Substrate: 1-dimensional C-valued field; propagator = i/(p² − m² + iε).

    Attributes
    ----------
    name : str
    mass : float | None       # in GeV; None = massless
    charge : float            # U(1)_em charge in units of e
    real : bool               # True for real scalar, False for complex
    spin : int = 0
    """
    name: str
    mass: Optional[float] = None
    charge: float = 0.0
    real: bool = False
    spin: int = 0

    @property
    def substrate_primitive(self) -> str:
        return "1-d complex scalar (real if real=True); propagator i/(p² − m²)"

    def __str__(self) -> str:
        kind = "real" if self.real else "complex"
        m = f", m={self.mass}" if self.mass else ""
        q = f", Q={self.charge}" if self.charge else ""
        return f"{kind} scalar {self.name}{m}{q}"


@dataclass(frozen=True)
class DiracSpinor:
    """
    Dirac fermion field ψ(x) of mass m and charge Q (in units of e).

    QFT view:  L_kin = ψ̄ (i γ^μ ∂_μ − m) ψ
    Substrate: Cl(1,3) 8-dimensional spinor (Cl(0,7) octonion left-mult,
                Wick-rotated to Lorentzian signature).  The 8 components
                = 4 Dirac × 2 SU(2)-doublet (only one carries U(1)_em
                charge per the substrate "internal SU(2) is vector"
                identification).  Propagator: i(p̸ + m)/(p² − m²).
                Vertex (QED): −i e γ^μ.

    Attributes
    ----------
    name : str
    mass : float
    charge : float
    n_color : int = 1         # 1 for leptons, 3 for quarks
    """
    name: str
    mass: float
    charge: float = 0.0
    n_color: int = 1
    spin: float = 0.5

    @property
    def substrate_primitive(self) -> str:
        return (
            "Cl(1,3) 8-dimensional spinor (octonion left-mult, "
            "Wick-rotated); 4 pos-energy + 4 neg-energy spinors per "
            "momentum; propagator (p̸ + m)/(p² − m²) on 8x8 space"
        )

    def __str__(self) -> str:
        c = f", N_c={self.n_color}" if self.n_color != 1 else ""
        return f"Dirac spinor {self.name} (m={self.mass}, Q={self.charge}{c})"


@dataclass(frozen=True)
class GaugeField:
    """
    Yang-Mills gauge field A^a_μ(x) for a compact Lie group G.

    QFT view:  L_kin = −(1/4) F^{aμν} F^a_{μν}, with field strength
                       F^a_{μν} = ∂_μ A^a_ν − ∂_ν A^a_μ + g f^{abc} A^b_μ A^c_ν.
    Substrate: g (Lie algebra) generators T^a as N×N matrices on the
                fundamental representation; structure constants f^{abc}
                from [T^a, T^b] = i f^{abc} T^c.  Propagator:
                  Feynman gauge:  −i η^{μν} δ^{ab} / k²
                  Axial gauge:    (P^{μν}_{axial} δ^{ab}) / k²

    Concrete substrate references (for canonical groups):
      U(1):    trivial; one Abelian generator
      SU(2):   Pauli matrices σ^a / 2; 3 generators
      SU(3):   Gell-Mann λ^a / 2; 8 generators (substrate: nwt.qcd primitives)
      Spin(7): so(7) on 7-dim vector V or 8-dim spinor S
                (substrate: 21 generators = 6 Lorentz + 3 internal SU(2)
                + 12 mixed; see Architectural-Lemmas memory)

    Attributes
    ----------
    name : str
    gauge_group : str         # "U(1)", "SU(2)", "SU(3)", "Spin(7)", ...
    coupling : float          # g; for QED this is e
    n_generators : int        # 1 for U(1), 3 for SU(2), 8 for SU(3)
    massless : bool = True
    spin : int = 1
    """
    name: str
    gauge_group: str
    coupling: float
    n_generators: int
    massless: bool = True
    spin: int = 1

    @property
    def substrate_primitive(self) -> str:
        return f"{self.gauge_group} Lie-algebra generators on fundamental rep"

    def __str__(self) -> str:
        return (
            f"{self.gauge_group} gauge field {self.name} "
            f"(g={self.coupling:.4f}, n_gen={self.n_generators})"
        )


@dataclass(frozen=True)
class FaddeevPopovGhost:
    """
    Faddeev-Popov ghost (c) and anti-ghost (c̄) fields associated with
    a non-Abelian gauge field.

    QFT view:  L_FP = c̄^a (−∂^μ D_μ^{ac}) c^c
                    = c̄^a (−∂² δ^{ac} − g_s f^{abc} ∂^μ A_μ^b) c^c
    Substrate: anti-commuting scalar field; propagator i δ^{ab} / k²;
                ghost-gluon vertex −g_s f^{abc} q^μ (q = outgoing-ghost
                momentum, Peskin convention).  Required for non-Abelian
                gauge invariance (BRST identity); see qcd_gg.py Phase B.8.
    """
    name: str
    gauge_group: str
    coupling: float

    @property
    def substrate_primitive(self) -> str:
        return "anti-commuting scalar; propagator i/k², no spin states"

    def __str__(self) -> str:
        return f"FP ghost {self.name} for {self.gauge_group}"
