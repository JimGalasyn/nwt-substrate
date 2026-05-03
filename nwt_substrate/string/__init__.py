"""
nwt.string
==========

The string-theoretic view of the substrate algebra.

The thesis: NWT proposes no new string theory.  The substrate's familiar
objects -- (p,q) torus-knot windings, K_7 toroidal embedding, α^(n/2)
Wilson amplitude tower, 2I icosahedral structure, Spin(7) on Cl(0,7),
Σ(2,3,5) Brieskorn sphere -- carry textbook string-theoretic *labels*:

  Substrate                    String-theoretic label
  ---------                    ----------------------
  (p,q) torus-knot winding     (p,q)-string SL(2,Z) doublet
  m_geom² ∝ p² + q²            (p,q)-string tension (weak-coupling g_s → 1)
  K_7 toroidal embedding       T² compactification with 21 1-cycles
  α^(n/2)·M_Pl tower (n=1..21) Kaluza-Klein mode tower on K_7
  2I (binary icosahedral)      E_8 ADE via McKay correspondence
  Σ(2,3,5) = Poincaré sphere   E_8 ADE singularity (heterotic / F-theory)
  Spin(7) on Cl(0,7)           Spin(7) holonomy on 8-manifolds (M-theory)
  2T → Spin(7) (Phase b2.12)   finite Spin(7) orbifold compactification
  Trefoil = T(2,3)             (2,3) torus knot ↔ E_8 fibre

Different vocabulary, same primitives.

Quick start::

    import nwt_substrate.string as st

    # Compactification targets
    print(st.K7_torus)             # K_7 with toroidal embedding
    print(st.poincare_sphere)      # Σ(2,3,5) = E_8 ADE

    # (p,q)-string view
    print(st.pq_string("electron")) # (2, 1)-string on K_7 toroidal
    print(st.pq_string("proton"))   # (1, 3)-string on K_7 toroidal

    # KK-mode tower
    print(st.kk_tower_summary())    # n=1..21 with 4 SM-scale clean hits

    # Holonomy classifications
    print(st.holonomy_summary())    # Spin(7), G₂, 2I, ...

    # ADE / McKay correspondence
    print(st.ade_lookup("2I"))      # 2I ↔ E_8

    # Multi-view from a single particle
    print(st.string_view("electron"))
"""

from .compactifications import (
    CompactificationTarget,
    K7_torus,
    poincare_sphere,
    T2_torus,
    ADEEntry,
    ADE_CATALOG,
    ade_lookup,
)

from .strings import (
    PQString,
    CANONICAL_PQ_STRINGS,
    pq_string,
)

from .kk_tower import (
    KKMode,
    M_PLANCK_GEV,
    ALPHA_EM,
    kk_tower,
    kk_tower_summary,
    seesaw_dual_pair,
    alpha_to_n,
)

from .holonomy import (
    HolonomyEntry,
    HOLONOMY_CATALOG,
    holonomy_lookup,
    holonomy_summary,
)


def string_view(particle_name: str) -> str:
    """
    Substrate-particle multi-view in string-theoretic vocabulary.

    Shows the (p,q)-string identification, KK-mode index in the K_7
    Wilson tower (when known), and which compactification target hosts
    the particle's substrate construction.
    """
    out = [f"=== string_view('{particle_name}') ==="]
    try:
        ps = pq_string(particle_name)
        out.append("")
        out.append(f"(p,q)-string view:  {ps}")
        out.append(f"  winding² (p² + q²)         = {ps.winding_squared}")
        out.append(f"  SL(2,Z) doublet            = {ps.sl2z_doublet}")
        out.append(f"  Compactification target    = {ps.compactification}")
        if ps.notes:
            out.append(f"  Notes: {ps.notes}")
    except KeyError as e:
        return f"Unknown particle '{particle_name}'.  Try one of: " \
               + ", ".join(CANONICAL_PQ_STRINGS.keys())

    # KK mode identification if applicable
    out.append("")
    tower = kk_tower()
    matches = [m for m in tower
               if m.sm_identification and particle_name.lower() in m.sm_identification.lower()]
    if matches:
        out.append("KK-mode (K_7 Wilson tower) identification:")
        for m in matches:
            out.append(f"  {m}")
    return "\n".join(out)


__all__ = [
    # compactifications
    "CompactificationTarget", "K7_torus", "poincare_sphere", "T2_torus",
    "ADEEntry", "ADE_CATALOG", "ade_lookup",
    # (p,q)-strings
    "PQString", "CANONICAL_PQ_STRINGS", "pq_string",
    # KK tower
    "KKMode", "M_PLANCK_GEV", "ALPHA_EM",
    "kk_tower", "kk_tower_summary", "seesaw_dual_pair", "alpha_to_n",
    # holonomy
    "HolonomyEntry", "HOLONOMY_CATALOG", "holonomy_lookup", "holonomy_summary",
    # multi-view helper
    "string_view",
]
