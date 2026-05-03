"""
24-particle compendium (Paper 6 benchmark, with 2026-04-30 nucleon correction).

Each particle specifies the substrate quantum numbers (p, q, m, n_q, f) and
PDG observables (mass, charge, baryon number, strangeness, charm, isospin).

Nucleon family corrected from (1,4,*,3) to (1,3,*,5) preserving Q_H = p*m;
this drops the median mass residual from 1.06% to 0.76% across all 24 particles.

Used by the factory to build particles by name.
"""

from __future__ import annotations

# Compendium entries: name, p, q, m, n_q, f, sector, m_obs, Q_pdg, B, S, C, I, notes
# (f = 2 * I_3; isospin projection encoded as framing)
COMPENDIUM = [
    # ---- Leptons (n_q = 0, framing handles charge separately per Paper 7) ----
    dict(name="e-",     p=2, q=1, m=3,  n_q=0, f=0,  sector="spinor",
         m_obs=0.511,    Q_pdg=-1, B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0,
         notes="anchor: m_e and v sets all scales; charge via Paper 7 contact"),
    dict(name="mu-",    p=1, q=8, m=9,  n_q=0, f=0,  sector="spinor",
         m_obs=105.66,   Q_pdg=-1, B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0),
    dict(name="tau-",   p=3, q=4, m=17, n_q=3, f=0,  sector="spinor",
         m_obs=1776.86,  Q_pdg=-1, B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0,
         notes="stealth baryon: substrate-topologically baryon, gauge lepton"),

    # ---- Light mesons ----
    dict(name="pi+",    p=3, q=5, m=5,  n_q=2, f=2,  sector="scalar",
         m_obs=139.57,   Q_pdg=+1, B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=1.0),
    dict(name="pi0",    p=7, q=3, m=18, n_q=2, f=0,  sector="scalar",
         m_obs=135.00,   Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=1.0),
    dict(name="K+",     p=2, q=5, m=8,  n_q=2, f=1,  sector="scalar",
         m_obs=493.68,   Q_pdg=+1, B_pdg=0, S_pdg=+1, C_pdg=0, I_pdg=0.5),
    dict(name="K0",     p=7, q=5, m=15, n_q=2, f=-1, sector="scalar",
         m_obs=497.61,   Q_pdg=0,  B_pdg=0, S_pdg=+1, C_pdg=0, I_pdg=0.5),
    dict(name="eta",    p=6, q=5, m=15, n_q=2, f=0,  sector="scalar",
         m_obs=547.86,   Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0),

    # ---- Vector mesons ----
    dict(name="rho",    p=5, q=7, m=7,  n_q=2, f=0,  sector="vector",
         m_obs=775.26,   Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=1.0),
    dict(name="omega",  p=4, q=5, m=17, n_q=2, f=0,  sector="vector",
         m_obs=782.66,   Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0),

    # ---- Nucleons (CORRECTED 2026-04-30 to (1,3,5,5)) ----
    dict(name="p",      p=1, q=3, m=5,  n_q=5, f=1,  sector="spinor",
         m_obs=938.27,   Q_pdg=+1, B_pdg=1, S_pdg=0, C_pdg=0, I_pdg=0.5,
         notes="corrected from Paper 6 (1,4,5,3) to (1,3,5,5); Q_H=5 preserved"),
    dict(name="n",      p=1, q=3, m=5,  n_q=5, f=-1, sector="spinor",
         m_obs=939.57,   Q_pdg=0,  B_pdg=1, S_pdg=0, C_pdg=0, I_pdg=0.5),

    # ---- Sigma family (CORRECTED 2026-04-30 to (1,3,6,5)) ----
    dict(name="Sigma+", p=1, q=3, m=6,  n_q=5, f=2,  sector="spinor",
         m_obs=1189.4,   Q_pdg=+1, B_pdg=1, S_pdg=-1, C_pdg=0, I_pdg=1.0),
    dict(name="Sigma0", p=1, q=3, m=6,  n_q=5, f=0,  sector="spinor",
         m_obs=1192.6,   Q_pdg=0,  B_pdg=1, S_pdg=-1, C_pdg=0, I_pdg=1.0),
    dict(name="Sigma-", p=1, q=3, m=6,  n_q=5, f=-2, sector="spinor",
         m_obs=1197.4,   Q_pdg=-1, B_pdg=1, S_pdg=-1, C_pdg=0, I_pdg=1.0),

    # ---- Other J=1/2 baryons ----
    dict(name="Lambda", p=3, q=4, m=12, n_q=3, f=0,  sector="spinor",
         m_obs=1115.7,   Q_pdg=0,  B_pdg=1, S_pdg=-1, C_pdg=0, I_pdg=0.0),
    dict(name="Xi0",    p=5, q=4, m=16, n_q=3, f=1,  sector="spinor",
         m_obs=1314.86,  Q_pdg=0,  B_pdg=1, S_pdg=-2, C_pdg=0, I_pdg=0.5),
    dict(name="Xi-",    p=5, q=4, m=16, n_q=3, f=-1, sector="spinor",
         m_obs=1321.71,  Q_pdg=-1, B_pdg=1, S_pdg=-2, C_pdg=0, I_pdg=0.5),

    # ---- J=3/2 baryons ----
    dict(name="Delta",  p=5, q=4, m=15, n_q=3, f=0,  sector="spinor",
         m_obs=1232.0,   Q_pdg=0,  B_pdg=1, S_pdg=0, C_pdg=0, I_pdg=1.5,
         notes="J=3/2 isoquadruplet"),
    dict(name="Sigma*", p=3, q=4, m=14, n_q=3, f=0,  sector="spinor",
         m_obs=1385.0,   Q_pdg=0,  B_pdg=1, S_pdg=-1, C_pdg=0, I_pdg=1.0,
         notes="J=3/2"),
    dict(name="Omega-", p=7, q=4, m=19, n_q=3, f=0,  sector="spinor",
         m_obs=1672.5,   Q_pdg=-1, B_pdg=1, S_pdg=-3, C_pdg=0, I_pdg=0.0,
         notes="J=3/2"),

    # ---- Charm mesons ----
    dict(name="D+",     p=2, q=7, m=5,  n_q=2, f=1,  sector="scalar",
         m_obs=1869.7,   Q_pdg=+1, B_pdg=0, S_pdg=0, C_pdg=+1, I_pdg=0.5),
    dict(name="D0",     p=3, q=7, m=7,  n_q=2, f=-1, sector="scalar",
         m_obs=1864.84,  Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=+1, I_pdg=0.5),

    # ---- Heavy quarkonium ----
    dict(name="J/psi",  p=2, q=7, m=7,  n_q=2, f=0,  sector="vector",
         m_obs=3096.9,   Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0),
    dict(name="Upsilon",p=4, q=9, m=8,  n_q=2, f=0,  sector="vector",
         m_obs=9460.3,   Q_pdg=0,  B_pdg=0, S_pdg=0, C_pdg=0, I_pdg=0.0),
]


# Index by name
_BY_NAME = {entry["name"]: entry for entry in COMPENDIUM}


def get_compendium_entry(name: str) -> dict | None:
    """Return the compendium dict for a particle, or None if not found."""
    return _BY_NAME.get(name)


def list_names() -> list[str]:
    """Return all compendium particle names."""
    return [entry["name"] for entry in COMPENDIUM]
