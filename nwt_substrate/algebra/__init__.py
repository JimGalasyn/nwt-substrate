"""Substrate algebra primitives: octonions, Clifford, Dirac, 2I irreps."""

from .octonions import (
    make_octonion_table,
    octo_mul,
    octo_assoc,
    octo_norm_sq,
    basis_vector,
    FANO_LINES,
)

from .clifford import (
    left_mult_matrix,
    right_mult_matrix,
    bivector_matrix,
)

from .dirac import (
    make_lorentz_gammas,
    slash_p,
    make_gamma5,
)

from .icosahedral import (
    sl2_f5_elements,
    conjugacy_classes,
    irrep_dimensions_2I,
)

from .su3 import (
    gell_mann_matrices,
    su3_generators,
    structure_constants,
    d_constants,
    fundamental_casimir,
    adjoint_casimir,
    C_F,
    C_A,
    N_C,
)

from .paley import (
    is_prime,
    quadratic_residues,
    non_residues,
    is_paley_admissible,
    PaleyDesign,
    paley_design,
    AGLBreakingAnalysis,
    agl_breaking_analysis,
    k7_paley_bifold,
)

from .g2_bridge import (
    PALEY_TO_BAEZ_LABELING,
    AGL_Z3_SIGN_LIFT,
    baez_permutation_of_agl_z3,
    signed_permutation_matrix,
    agl_z3_g2_matrix,
    complex_structure_J,
    holomorphic_basis_C3,
    BridgeVerification,
    verify_bridge,
)

__all__ = [
    "make_octonion_table",
    "octo_mul",
    "octo_assoc",
    "octo_norm_sq",
    "basis_vector",
    "FANO_LINES",
    "left_mult_matrix",
    "right_mult_matrix",
    "bivector_matrix",
    "make_lorentz_gammas",
    "slash_p",
    "make_gamma5",
    "sl2_f5_elements",
    "conjugacy_classes",
    "irrep_dimensions_2I",
    "gell_mann_matrices",
    "su3_generators",
    "structure_constants",
    "d_constants",
    "fundamental_casimir",
    "adjoint_casimir",
    "C_F",
    "C_A",
    "N_C",
    "is_prime",
    "quadratic_residues",
    "non_residues",
    "is_paley_admissible",
    "PaleyDesign",
    "paley_design",
    "AGLBreakingAnalysis",
    "agl_breaking_analysis",
    "k7_paley_bifold",
    "PALEY_TO_BAEZ_LABELING",
    "AGL_Z3_SIGN_LIFT",
    "baez_permutation_of_agl_z3",
    "signed_permutation_matrix",
    "agl_z3_g2_matrix",
    "complex_structure_J",
    "holomorphic_basis_C3",
    "BridgeVerification",
    "verify_bridge",
]
