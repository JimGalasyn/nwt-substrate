"""
Canonical Lagrangians: textbook QFT expressions implemented from substrate
primitives.

Each Lagrangian here:
  - displays in standard textbook form via __str__
  - exposes its substrate primitives via .substrate_view()
  - feeds into amplitude calculations via .feynman_rules() (which delegates
    to nwt_substrate.amplitudes vertex/propagator functions)

The thesis: these textbook Lagrangians are not separately implemented from
the rest of the library -- they are *labels* for substrate-algebraic
constructions that already exist.
"""

from __future__ import annotations

from .fields import ScalarField, DiracSpinor, GaugeField, FaddeevPopovGhost
from .lagrangian import Lagrangian, InteractionTerm
from ..isa.constants import N_C_SU3, N_GLUONS


# ===========================================================================
# Constants
# ===========================================================================
# Structural integers (colour count, gluon count) come from the isa anchor.
# The couplings below are the MEASURED inputs (CODATA/PDG), deliberately
# distinct from the substrate-PREDICTED isa.ALPHA_SUBSTRATE (= 1/(25π√3+1)).

ALPHA_QED = 1.0 / 137.035999084   # CODATA measured α (NOT isa.ALPHA_SUBSTRATE)
E_CHARGE = 0.30282212        # = sqrt(4 pi alpha) at q^2 -> 0
ALPHA_S = 0.1179             # at M_Z
G_S = (4.0 * 3.141592653589793 * ALPHA_S) ** 0.5

# Standard fermion masses (GeV, used throughout the library)
M_E = 0.510998928e-3
M_MU = 0.10566
M_TAU = 1.77686
M_U = 2.16e-3
M_D = 4.67e-3
M_S = 0.0934
M_C = 1.27
M_B = 4.18
M_T = 172.7


# ===========================================================================
# Standard fermions (re-used across L_QED, L_QCD, L_EW)
# ===========================================================================

# Charged leptons (Q in units of e; lepton charge is -1 conventionally)
electron = DiracSpinor("e",   M_E,  charge=-1.0, n_color=1)
muon     = DiracSpinor("μ",   M_MU, charge=-1.0, n_color=1)
tau      = DiracSpinor("τ",   M_TAU, charge=-1.0, n_color=1)

# Quarks
up_quark    = DiracSpinor("u", M_U, charge=+2.0/3.0, n_color=N_C_SU3)
down_quark  = DiracSpinor("d", M_D, charge=-1.0/3.0, n_color=N_C_SU3)
strange     = DiracSpinor("s", M_S, charge=-1.0/3.0, n_color=N_C_SU3)
charm       = DiracSpinor("c", M_C, charge=+2.0/3.0, n_color=N_C_SU3)
bottom      = DiracSpinor("b", M_B, charge=-1.0/3.0, n_color=N_C_SU3)
top         = DiracSpinor("t", M_T, charge=+2.0/3.0, n_color=N_C_SU3)


# ===========================================================================
# QED Lagrangian
# ===========================================================================

photon = GaugeField(
    name="A",
    gauge_group="U(1)_em",
    coupling=E_CHARGE,
    n_generators=1,
)


def _build_qed() -> Lagrangian:
    """L_QED = ψ̄ (i γ^μ D_μ − m) ψ − (1/4) F^{μν} F_{μν}, summed over
    charged fermions, with D_μ = ∂_μ + i e Q_f A_μ."""

    fermions = [electron, muon, tau,
                up_quark, down_quark, strange, charm, bottom, top]

    interactions = [
        InteractionTerm(
            name=f"QED vertex ({f.name})",
            fields=(f.name, "A", f.name),
            expression=fr"-i e Q_{{{f.name}}} \bar{{ψ}}_{{{f.name}}} γ^μ ψ_{{{f.name}}} A_μ",
            coupling=E_CHARGE * f.charge,
            substrate_vertex_fn="nwt_substrate.amplitudes.vertices.qed_vertex",
        )
        for f in fermions
        if f.charge != 0.0
    ]

    text = (
        r"L_QED = − (1/4) F^{μν} F_{μν} + Σ_f ψ̄_f (i γ^μ D_μ − m_f) ψ_f"
        ", with D_μ = ∂_μ + i e Q_f A_μ"
    )

    return Lagrangian(
        name="QED",
        fields=fermions + [photon],
        interactions=interactions,
        gauge_symmetries=["U(1)_em"],
        lorentz="SO(1,3)",
        text=text,
    )


qed = _build_qed()


# ===========================================================================
# QCD Lagrangian
# ===========================================================================

gluon = GaugeField(
    name="A_g",
    gauge_group="SU(3)_color",
    coupling=G_S,
    n_generators=N_GLUONS,
)

ghost_qcd = FaddeevPopovGhost(
    name="c",
    gauge_group="SU(3)_color",
    coupling=G_S,
)


def _build_qcd() -> Lagrangian:
    """L_QCD = Σ_q ψ̄_q (i γ^μ D_μ − m_q) ψ_q
              − (1/4) G^{aμν} G^a_{μν}
              + L_FP (Faddeev-Popov ghosts)
       with D_μ = ∂_μ − i g_s T^a A^a_μ on quark fundamentals."""

    quarks = [up_quark, down_quark, strange, charm, bottom, top]

    interactions = [
        InteractionTerm(
            name=f"qqg ({q.name})",
            fields=(q.name, "A_g", q.name),
            expression=fr"-i g_s \bar{{ψ}}_{{{q.name}}} γ^μ T^a ψ_{{{q.name}}} A^a_μ",
            coupling=G_S,
            substrate_vertex_fn="nwt_substrate.amplitudes.vertices.gluon_vertex",
        )
        for q in quarks
    ]

    interactions += [
        InteractionTerm(
            name="ggg (3-gluon)",
            fields=("A_g", "A_g", "A_g"),
            expression=r"-i g_s f^{abc} V^{μνρ}(k_1, k_2, k_3)",
            coupling=G_S,
            substrate_vertex_fn=(
                "nwt_substrate.amplitudes.vertices.three_gluon_vertex_lorentz"
            ),
        ),
        InteractionTerm(
            name="gggg (4-gluon)",
            fields=("A_g", "A_g", "A_g", "A_g"),
            expression=(r"-i g_s² × Σ_χ f-tensor × Lorentz_χ "
                        r"(3 colour structures, see qcd_gg)"),
            coupling=G_S * G_S,
            substrate_vertex_fn=(
                "nwt_substrate.amplitudes.vertices.four_gluon_vertex_lorentz"
            ),
        ),
        InteractionTerm(
            name="cgc (ghost-gluon)",
            fields=("c", "A_g", "c"),
            expression=r"-g_s f^{abc} q^μ  (q = outgoing-ghost momentum)",
            coupling=G_S,
            substrate_vertex_fn=(
                "(no top-level helper yet; see qcd_gg._M_ghost_amplitude_squared_avg "
                "for the BRST-closure use)"
            ),
        ),
    ]

    text = (
        r"L_QCD = Σ_q ψ̄_q (i γ^μ D_μ − m_q) ψ_q "
        r"− (1/4) G^{aμν} G^a_{μν} + L_FP"
        ", with D_μ = ∂_μ − i g_s T^a A^a_μ"
    )

    return Lagrangian(
        name="QCD",
        fields=quarks + [gluon, ghost_qcd],
        interactions=interactions,
        gauge_symmetries=["SU(3)_color"],
        lorentz="SO(1,3)",
        text=text,
    )


qcd = _build_qcd()


# ===========================================================================
# Generic Yang-Mills SU(N)
# ===========================================================================

def yang_mills(N: int, coupling: float = 1.0) -> Lagrangian:
    """
    Pure Yang-Mills L_YM(SU(N)) = -(1/4) F^{aμν} F^a_{μν}.

    Substrate: SU(N) Lie-algebra generators (N²-1 of them) on the
    fundamental representation, structure constants f^{abc} from
    [T^a, T^b] = i f^{abc} T^c.  For N=3, this is the gluon part of
    QCD; for N=2, the W bosons (after EW symmetry breaking, they get mass).
    """
    n_gen = N * N - 1
    A = GaugeField(
        name=f"A_{{SU({N})}}",
        gauge_group=f"SU({N})",
        coupling=coupling,
        n_generators=n_gen,
    )
    interactions = [
        InteractionTerm(
            name=f"AAA (3-gauge, SU({N}))",
            fields=("A", "A", "A"),
            expression=fr"-i g f^{{abc}} V^{{μνρ}}(k_1, k_2, k_3)",
            coupling=coupling,
            substrate_vertex_fn=(
                "nwt_substrate.amplitudes.vertices.three_gluon_vertex_lorentz"
            ),
        ),
        InteractionTerm(
            name=f"AAAA (4-gauge, SU({N}))",
            fields=("A", "A", "A", "A"),
            expression=(r"-i g² × (3 colour-Lorentz structures)"),
            coupling=coupling * coupling,
            substrate_vertex_fn=(
                "nwt_substrate.amplitudes.vertices.four_gluon_vertex_lorentz"
            ),
        ),
    ]
    text = (
        f"L_YM(SU({N})) = − (1/4) F^{{aμν}} F^a_{{μν}}, "
        f"a = 1..{n_gen}, [T^a, T^b] = i f^{{abc}} T^c"
    )
    return Lagrangian(
        name=f"Yang-Mills SU({N})",
        fields=[A],
        interactions=interactions,
        gauge_symmetries=[f"SU({N})"],
        lorentz="SO(1,3)",
        text=text,
    )


# ===========================================================================
# Free scalar / Klein-Gordon
# ===========================================================================

def klein_gordon(name: str = "φ", mass: float = 0.0,
                  charge: float = 0.0) -> Lagrangian:
    """
    Free Klein-Gordon Lagrangian:
       L_KG = (1/2) (∂_μ φ)(∂^μ φ) − (1/2) m² φ²    (real scalar)
    or  L_KG = (∂_μ φ)†(∂^μ φ) − m² φ†φ              (complex scalar).
    """
    real = (charge == 0.0)
    phi = ScalarField(name=name, mass=mass, charge=charge, real=real)
    text = (
        fr"L_KG = (∂_μ {name})† (∂^μ {name}) − m² {name}† {name}"
        if not real else
        fr"L_KG = (1/2) (∂_μ {name})(∂^μ {name}) − (1/2) m² {name}²"
    )
    return Lagrangian(
        name=f"Klein-Gordon ({name}, m={mass})",
        fields=[phi],
        interactions=[],
        gauge_symmetries=[],
        lorentz="SO(1,3)",
        text=text,
    )


# ===========================================================================
# Standard Model assemblage (sum of QED + QCD + ... )
# ===========================================================================

def standard_model_matter() -> Lagrangian:
    """
    Compose the matter+gauge sector L_QED + L_QCD.  Higgs/EW pieces are
    placeholders pending a dedicated nwt.electroweak shim.

    Substrate identification: this Lagrangian sums all standard-model
    matter with their U(1)_em + SU(3)_color gauge content.  The substrate
    primitives (Cl(1,3) spinors, U(1)+SU(3) generators, Feynman propagators
    and vertices) are exactly those used in the library's amplitude module.
    """
    return qed + qcd
