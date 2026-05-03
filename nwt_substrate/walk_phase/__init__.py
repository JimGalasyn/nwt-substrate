"""
Walk-phase: substrate-algebraic scattering and decay computations.

Verified processes (from yesterday's walk-phase work, 2026-04-29):
  - Compton scattering (Klein-Nishina to 4.8e-16)
  - Muon decay (V-A vertex, ratio 4.000000 to 8.9e-15)
  - Neutron beta-decay (full g_A=1.27 vertex, ratio 4.000000)

TODO: refactor existing octonion_dirac_step{3,3b,4a,4b,4c}.py into modules:
  - spinors.py        (spinor solutions to Dirac equation)
  - propagators.py    (S_F = (slash{P}+m)/(P^2-m^2))
  - currents.py       (V-A, EM, strong)
  - compton.py        (gamma + e -> gamma + e)
  - decay.py          (n-body decay rates)
  - phase_space.py    (Monte Carlo n-body phase space)
"""
