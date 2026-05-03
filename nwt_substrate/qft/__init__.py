"""
nwt.qft
=======

The QFT view of the substrate algebra.

The thesis: textbook quantum field theory and the substrate algebra are
not separate things.  The substrate Cl(1,3) spinors are exactly the
8-dimensional Dirac spinors in the natural Cl(0,7)-octonion basis;
the SU(3) Gell-Mann generators are the substrate's color algebra; the
3-gluon and 4-gluon vertices are encoded in standard non-Abelian Yang-Mills
form; the running couplings come from the same 1-loop vacuum polarisation
diagrams.  Different vocabulary, same primitives.

This shim exposes the QFT-community-standard view -- Lagrangian densities,
field types, Feynman-rule extraction -- as a labelling of substrate-
algebraic constructions that already exist in the library.

Quick start::

    import nwt_substrate.qft as qft

    print(qft.qed)                   # textbook L_QED expression
    print(qft.qed.substrate_view())  # substrate primitives behind it

    # Feynman rules pull from the same primitives used by nwt.qed, nwt.qcd
    rules = qft.qed.feynman_rules()
    print(rules['vertices'])          # qed_vertex, gluon_vertex, ...

    # 1-loop β-function derived from the Lagrangian's particle content,
    # using the substrate-derived QED vacuum-polarisation result
    print(qft.qed.beta_0())           # = 16/3 = 5.333... (PDG match)

    # Compose Lagrangians
    L_sm = qft.qed + qft.qcd
    print(L_sm.fields_summary())

    # Generic Yang-Mills SU(N)
    L_ym3 = qft.yang_mills(N=3)
    print(L_ym3.beta_0(n_f_dirac=6))  # pure-glue + 6 quarks => 7

    # Single-particle multi-view inspection
    qft.multiview('electron')         # show e- in QED / Substrate / Particle views
"""

from .fields import (
    ScalarField,
    DiracSpinor,
    GaugeField,
    FaddeevPopovGhost,
)
from .lagrangian import Lagrangian, InteractionTerm
from .lagrangians import (
    # Particles (re-exported for convenience)
    electron, muon, tau,
    up_quark, down_quark, strange, charm, bottom, top,
    photon, gluon, ghost_qcd,
    # Lagrangians
    qed, qcd, yang_mills, klein_gordon, standard_model_matter,
    # Constants
    ALPHA_QED, E_CHARGE, ALPHA_S, G_S,
)


def multiview(particle_name: str) -> str:
    """
    Single-object multi-view inspector.

    Show how one substrate object appears through different lenses:
      QFT view (Lagrangian role + Field type)
      Substrate-algebra view (Cl(0,7) / Cl(1,3) / 2I / K_7 / SU(N) primitive)
      Particle compendium view (mass formula tuple, J^P, B, S, ...)
      Heron view (qiskit qubit / circuit role, if applicable)

    Parameters
    ----------
    particle_name : str
        e.g. "electron", "muon", "u" (up quark), "photon", "gluon"
    """
    out = [f"=== Multi-view of '{particle_name}' ==="]

    name_map = {
        "e": electron, "electron": electron,
        "mu": muon, "muon": muon, "μ": muon,
        "tau": tau, "τ": tau,
        "u": up_quark, "up": up_quark,
        "d": down_quark, "down": down_quark,
        "s": strange, "c": charm, "b": bottom, "t": top,
        "photon": photon, "γ": photon, "A": photon,
        "gluon": gluon, "g": gluon, "A_g": gluon,
    }
    if particle_name not in name_map:
        return f"Unknown particle '{particle_name}'.  Try: {list(name_map.keys())}"

    f = name_map[particle_name]
    out.append("")
    out.append(f"QFT view:")
    out.append(f"  {f}")
    out.append(f"")
    out.append(f"Substrate primitive:")
    out.append(f"  {f.substrate_primitive}")

    # Particle-compendium view (only fermions are in the compendium)
    try:
        from ..particles import particle as get_particle
        p = get_particle(particle_name)
        out.append(f"")
        out.append(f"Particle compendium:")
        out.append(f"  {p}")
        if hasattr(p, "p") and hasattr(p, "q"):
            out.append(f"  substrate tuple: (p={p.p}, q={p.q}, m={p.m}, n_q={p.n_q})")
        if hasattr(p, "mass_pred"):
            out.append(f"  mass_pred = {p.mass_pred:.4g} MeV")
    except Exception:
        pass

    return "\n".join(out)


__all__ = [
    # Fields
    "ScalarField", "DiracSpinor", "GaugeField", "FaddeevPopovGhost",
    # Lagrangian
    "Lagrangian", "InteractionTerm",
    # Particles
    "electron", "muon", "tau",
    "up_quark", "down_quark", "strange", "charm", "bottom", "top",
    "photon", "gluon", "ghost_qcd",
    # Lagrangians
    "qed", "qcd", "yang_mills", "klein_gordon", "standard_model_matter",
    # Constants
    "ALPHA_QED", "E_CHARGE", "ALPHA_S", "G_S",
    # Helpers
    "multiview",
]
