"""QPU backend adapters: vendor-specific plumbing behind a common protocol."""
from .base import Counts, Capabilities, Compiled, Backend, Job
from .braket import BraketBackend, canonicalize_braket
from .ibm import IBMBackend, canonicalize_ibm
from .simulator import SimulatorBackend

__all__ = [
    "Counts", "Capabilities", "Compiled", "Backend", "Job",
    "BraketBackend", "canonicalize_braket",
    "IBMBackend", "canonicalize_ibm",
    "SimulatorBackend",
]
