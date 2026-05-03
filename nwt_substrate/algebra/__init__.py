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
]
