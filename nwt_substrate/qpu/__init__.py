"""nwt_substrate.qpu — cross-architecture QPU interface.

Vendor-neutral substrate physics (spec, decode) separated from vendor-specific
plumbing (adapters). See analysis/qpu_interface_design.md in null-worldtube-private.

M1 (this milestone): spec + decode + adapters with the canonical counts contract.
"""
from . import spec, decode
from .adapters import (
    Counts, Capabilities, Compiled, Backend, Job,
    BraketBackend, IBMBackend, SimulatorBackend,
    canonicalize_braket, canonicalize_ibm,
)

__all__ = [
    "spec", "decode",
    "Counts", "Capabilities", "Compiled", "Backend", "Job",
    "BraketBackend", "IBMBackend", "SimulatorBackend",
    "canonicalize_braket", "canonicalize_ibm",
]
