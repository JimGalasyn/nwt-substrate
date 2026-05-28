"""
nwt_substrate.benchmarks — substrate algebra vs traditional methods.

Times the substrate algebra's closed-form predictions and compares against
typical computational cost of traditional methods (lattice QCD, CKMfitter,
explicit Lanczos diagonalization, etc.).

Key finding: substrate algebra is 10⁶-10¹⁸× faster for FORWARD prediction
of SM parameters, topological invariants, and DM phenomenology.  The
"speedup" is for prediction at substrate-accuracy (typically 0.1-7%);
traditional methods retain the precision lead for sub-ppm validation.

Quick start:

    >>> from nwt_substrate.benchmarks import run_all
    >>> run_all()        # runs all benchmarks, prints comparison table

For individual benchmarks see compute_speed module.
"""

from .compute_speed import (
    BenchmarkResult,
    benchmark_alpha_derivation,
    benchmark_mass_spectrum,
    benchmark_modular_data,
    benchmark_ckm_cabibbo,
    benchmark_full_ckm,
    benchmark_higgs_vev,
    benchmark_higgs_mass_vs_98gev,
    benchmark_fermi_constant,
    benchmark_z_boson_width,
    benchmark_muon_lifetime,
    benchmark_neutrino_sector,
    benchmark_pmns_angles,
    benchmark_decay_constants,
    benchmark_vector_meson_decay,
    benchmark_atomic_hydrogen,
    benchmark_electron_anomaly,
    benchmark_qcd_constants,
    benchmark_sin2_theta_W,
    benchmark_black_hole_thermodynamics,
    benchmark_chemistry,
    benchmark_k7_face_structure,
    benchmark_wimp_tower,
    benchmark_gravitational_constant,
    benchmark_lambda_cc,
    benchmark_omega_b_c,
    benchmark_eta_B,
    run_all,
)

__all__ = [
    "BenchmarkResult",
    "benchmark_alpha_derivation",
    "benchmark_mass_spectrum",
    "benchmark_modular_data",
    "benchmark_ckm_cabibbo",
    "benchmark_full_ckm",
    "benchmark_higgs_vev",
    "benchmark_higgs_mass_vs_98gev",
    "benchmark_fermi_constant",
    "benchmark_z_boson_width",
    "benchmark_muon_lifetime",
    "benchmark_neutrino_sector",
    "benchmark_pmns_angles",
    "benchmark_decay_constants",
    "benchmark_vector_meson_decay",
    "benchmark_atomic_hydrogen",
    "benchmark_electron_anomaly",
    "benchmark_qcd_constants",
    "benchmark_sin2_theta_W",
    "benchmark_black_hole_thermodynamics",
    "benchmark_chemistry",
    "benchmark_k7_face_structure",
    "benchmark_wimp_tower",
    "benchmark_gravitational_constant",
    "benchmark_lambda_cc",
    "benchmark_omega_b_c",
    "benchmark_eta_B",
    "run_all",
]
