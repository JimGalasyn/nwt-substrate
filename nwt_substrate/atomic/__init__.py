"""
nwt.atomic
==========

The atomic-physics view of the substrate algebra: substrate Coulomb
potential V(r) = -α ℏc/r from P8 closure (Galasyn 2026-05-23) plus
Bohr-Rydberg spectrum and precision-EM predictions for hydrogen.

The substrate Coulomb potential is derived from the P1 polar-photon
vertex sustained at frequency ω, with substrate distance r = c/ω
(D2 identification: each polar-photon vertex application is a
"light-cone tick").

All hydrogen-spectroscopy predictions follow from this Coulomb potential
+ substrate α (Paper 17 K_7 Wilson) at zero free EW parameters; m_e is
the dimensional anchor.

Quick start::

    import nwt_substrate.atomic as at

    print(at.rydberg())              # 13.6059 eV   (+15 ppm vs PDG)
    print(at.bohr_radius())          # 52.917 pm    (-8 ppm vs PDG)
    print(at.hyperfine_21cm())       # 5.88 μeV     (+0.1% vs PDG)
    print(at.electron_a_e_one_loop()) # 0.00116142  (+8 ppm vs PDG)

    # Full precision chain
    print(at.precision_chain_summary())

    # Verification
    result = at.verify_substrate_hydrogen()
    assert result['pass']
"""

from .hydrogen import (
    # Physical constants
    M_E_EV, M_P_EV, HBAR_C_EV_PM, G_P,
    # Substrate hydrogen formulas
    rydberg,
    bohr_radius,
    hydrogen_R_H,
    lyman_alpha,
    hyperfine_21cm,
    fine_structure_2P,
    electron_a_e_one_loop,
    lamb_shift_scale,
    # Precision chain + verification
    precision_chain,
    verify_substrate_hydrogen,
    precision_chain_summary,
)


__all__ = [
    "M_E_EV", "M_P_EV", "HBAR_C_EV_PM", "G_P",
    "rydberg", "bohr_radius", "hydrogen_R_H", "lyman_alpha",
    "hyperfine_21cm", "fine_structure_2P",
    "electron_a_e_one_loop", "lamb_shift_scale",
    "precision_chain", "verify_substrate_hydrogen", "precision_chain_summary",
]
