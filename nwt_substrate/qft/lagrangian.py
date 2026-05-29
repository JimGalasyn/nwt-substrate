"""
Lagrangian dataclass for nwt.qft.

A Lagrangian is the QFT-vocabulary description of a substrate-algebraic
theory.  It carries:

  - a list of fields (kinetic terms inferred from the Field types)
  - a list of interaction vertices (each pointing to its substrate
    vertex factor)
  - gauge / Lorentz / internal symmetry assertions
  - a `feynman_rules()` accessor that delegates to substrate primitives
  - a `substrate_view()` accessor that surfaces the underlying algebraic
    structure

The thesis: same physics, different vocabulary.  The substrate primitives
are the implementation; the Lagrangian is the QFT-community-standard view.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Callable, Optional

from .fields import ScalarField, DiracSpinor, GaugeField, FaddeevPopovGhost


@dataclass
class InteractionTerm:
    """
    One interaction vertex in a Lagrangian.

    name        -- short label, e.g. "QED vertex", "QCD quark-gluon vertex"
    fields      -- tuple of field names participating in the vertex
    expression  -- LaTeX-like string for textbook display
    coupling    -- coupling constant (g, e, g_s, ...)
    substrate_vertex_fn  -- callable returning the substrate vertex factor,
                             OR a string naming the substrate primitive
                             (when the function lives elsewhere)
    """
    name: str
    fields: tuple
    expression: str
    coupling: float
    substrate_vertex_fn: object = None

    def __str__(self) -> str:
        return f"{self.name}: {self.expression}"


@dataclass
class Lagrangian:
    """
    A quantum-field Lagrangian density L(x) as a substrate-algebra view.

    Attributes
    ----------
    name : str                              -- e.g. "QED", "QCD", "Yang-Mills SU(N)"
    fields : list[Field]                    -- field content (Scalars, Spinors, Gauge, Ghosts)
    interactions : list[InteractionTerm]    -- vertex content
    gauge_symmetries : list[str]            -- e.g. ["U(1)_em"]
    lorentz : str                           -- "SO(1,3)" by default
    text : str                              -- textbook expression for display

    Methods
    -------
    feynman_rules()   -- dict {propagators, vertices} pointing at substrate primitives
    substrate_view()  -- multi-line description of underlying substrate primitives
    fields_summary()  -- pretty-print field list
    + (combine)       -- concatenate two Lagrangians (for SM = L_QED + L_QCD + ...)
    """
    name: str
    fields: list = dc_field(default_factory=list)
    interactions: list = dc_field(default_factory=list)
    gauge_symmetries: list = dc_field(default_factory=list)
    lorentz: str = "SO(1,3)"
    text: str = ""

    def __str__(self) -> str:
        return f"Lagrangian {self.name}: {self.text}"

    def __add__(self, other: "Lagrangian") -> "Lagrangian":
        """Combine two Lagrangians.  Field/interaction lists concatenate."""
        if not isinstance(other, Lagrangian):
            return NotImplemented
        return Lagrangian(
            name=f"{self.name} + {other.name}",
            fields=self.fields + other.fields,
            interactions=self.interactions + other.interactions,
            gauge_symmetries=list(set(self.gauge_symmetries +
                                       other.gauge_symmetries)),
            lorentz=self.lorentz,
            text=f"({self.text}) + ({other.text})",
        )

    def fields_summary(self) -> str:
        if not self.fields:
            return "(no fields)"
        return "\n".join(f"  - {f}" for f in self.fields)

    def feynman_rules(self) -> dict:
        """
        Substrate-implemented Feynman rules for this Lagrangian.

        Returns dict with keys:
          'propagators' -> dict mapping field name -> substrate prop function name
          'vertices'    -> dict mapping interaction name -> substrate vertex factor

        These are the *same* primitives used elsewhere in the library, just
        re-exposed in QFT-style.  The point: building a Lagrangian and
        building substrate amplitudes are two views of the same calculation.
        """
        props = {}
        for f in self.fields:
            if isinstance(f, ScalarField):
                props[f.name] = (
                    "nwt_substrate.amplitudes.propagators.scalar_propagator"
                )
            elif isinstance(f, DiracSpinor):
                props[f.name] = (
                    "nwt_substrate.amplitudes.propagators.fermion_propagator"
                )
            elif isinstance(f, GaugeField):
                if f.massless:
                    props[f.name] = (
                        "nwt_substrate.amplitudes.propagators.photon_propagator"
                        " (Feynman gauge -i η^{μν}/k²; "
                        "axial-gauge alt available in qcd_gg)"
                    )
                else:
                    props[f.name] = (
                        "nwt_substrate.amplitudes.propagators.massive_vector_propagator"
                    )
            elif isinstance(f, FaddeevPopovGhost):
                props[f.name] = "i δ^{ab} / k²  (Faddeev-Popov scalar ghost)"

        vtxs = {}
        for iv in self.interactions:
            vtxs[iv.name] = (
                iv.substrate_vertex_fn if isinstance(iv.substrate_vertex_fn, str)
                else f"{iv.expression}  -> {iv.substrate_vertex_fn}"
            )
        return {"propagators": props, "vertices": vtxs}

    def substrate_view(self) -> str:
        """
        Multi-line description of the substrate primitives this Lagrangian
        is built from.  This is the "look behind the QFT vocabulary" view.
        """
        lines = [f"Lagrangian {self.name} as substrate algebra:"]
        lines.append("")
        lines.append("Fields (→ substrate primitives):")
        for f in self.fields:
            lines.append(f"  {f.name:8s}  →  {f.substrate_primitive}")
        if self.interactions:
            lines.append("")
            lines.append("Interaction vertices (→ substrate factors):")
            for iv in self.interactions:
                ref = (iv.substrate_vertex_fn if isinstance(iv.substrate_vertex_fn, str)
                       else getattr(iv.substrate_vertex_fn, "__name__", "(callable)"))
                lines.append(f"  {iv.name:25s} →  {ref}")
        if self.gauge_symmetries:
            lines.append("")
            lines.append(f"Gauge symmetries: {', '.join(self.gauge_symmetries)}")
        lines.append(f"Lorentz: {self.lorentz}")
        return "\n".join(lines)

    def beta_0(self, n_f_dirac: Optional[int] = None) -> Optional[float]:
        """
        1-loop β-function coefficient if this is a known gauge theory.
        Delegates to nwt_substrate.amplitudes.vacuum_polarization.

        Recognises U(1)_em (QED), SU(3)_color (QCD), and any generic
        non-Abelian ``SU(N)`` tag emitted by :func:`yang_mills`. For a
        non-Abelian SU(N) factor, β_0 = (11/3) N − (2/3) n_f, with n_f
        defaulting to 6 for QCD (its quark content) and to 0 (pure glue)
        for a bare Yang-Mills SU(N) that carries no matter.

        Returns None if not a recognised gauge theory.
        """
        import re
        from ..amplitudes import vacuum_polarization as vp
        if "U(1)_em" in self.gauge_symmetries and n_f_dirac is None:
            species = vp.standard_qed_species()
            return vp.qed_beta_0_total(species)
        if "SU(3)_color" in self.gauge_symmetries:
            n_f = 6 if n_f_dirac is None else n_f_dirac
            return vp.qcd_beta_0(n_f_dirac=n_f, N_c=3)
        # Generic non-Abelian SU(N), e.g. the bare "SU(N)" tag from yang_mills(N).
        for sym in self.gauge_symmetries:
            m = re.fullmatch(r"SU\((\d+)\)", sym)
            if m:
                N_c = int(m.group(1))
                n_f = 0 if n_f_dirac is None else n_f_dirac
                return vp.qcd_beta_0(n_f_dirac=n_f, N_c=N_c)
        return None
