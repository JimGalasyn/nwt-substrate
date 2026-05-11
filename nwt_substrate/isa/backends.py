"""Backend dispatch for substrate ISA — numpy + torch (CPU/CUDA).

The ISA primitives (so(7) basis, batched trace contractions) are
written in terms of an abstract `array` interface that's satisfied
by both numpy and torch. This module selects the backend and
provides factory functions that return the right module/device.

Default backend is "numpy". Torch backends are opt-in via
`select_backend("torch_cuda")` or by passing `backend=` to ISA calls.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np


BackendName = Literal["numpy", "torch_cpu", "torch_cuda"]


@dataclass(frozen=True)
class Backend:
    """Backend descriptor.

    Attributes
    ----------
    name : "numpy" | "torch_cpu" | "torch_cuda"
    xp : array module (numpy or torch)
    device : str or None (None for numpy, "cpu"/"cuda" for torch)
    """
    name: BackendName
    xp: Any
    device: str | None

    def asarray(self, x, dtype=None):
        """Convert input to backend array on the right device."""
        if self.name == "numpy":
            return np.asarray(x, dtype=dtype)
        else:
            import torch
            t = torch.as_tensor(x, dtype=_torch_dtype(dtype) if dtype else None)
            return t.to(self.device)

    def to_numpy(self, x) -> np.ndarray:
        """Bring result back to host as numpy array."""
        if self.name == "numpy":
            return np.asarray(x)
        else:
            return x.detach().cpu().numpy()


def _torch_dtype(np_dtype):
    """Map numpy dtype to torch dtype."""
    import torch
    mapping = {
        np.float32: torch.float32,
        np.float64: torch.float64,
        np.int32: torch.int32,
        np.int64: torch.int64,
        np.bool_: torch.bool,
    }
    if np_dtype is None:
        return None
    if hasattr(np_dtype, "type"):
        np_dtype = np_dtype.type
    return mapping.get(np_dtype, torch.float64)


_NUMPY_BACKEND = Backend(name="numpy", xp=np, device=None)


def get_backend(name: BackendName = "numpy") -> Backend:
    """Return Backend descriptor for the named backend.

    Raises ImportError if torch is requested but not installed.
    Raises RuntimeError if torch_cuda is requested but CUDA unavailable.
    """
    if name == "numpy":
        return _NUMPY_BACKEND
    if name in ("torch_cpu", "torch_cuda"):
        try:
            import torch
        except ImportError as e:
            raise ImportError(
                f"backend={name!r} requires torch; pip install torch"
            ) from e
        if name == "torch_cuda":
            if not torch.cuda.is_available():
                raise RuntimeError(
                    "backend='torch_cuda' requested but CUDA is unavailable"
                )
            return Backend(name="torch_cuda", xp=torch, device="cuda")
        return Backend(name="torch_cpu", xp=torch, device="cpu")
    raise ValueError(f"Unknown backend: {name!r}")


def available_backends() -> list[BackendName]:
    """List backends currently usable on this machine."""
    avail: list[BackendName] = ["numpy"]
    try:
        import torch
        avail.append("torch_cpu")
        if torch.cuda.is_available():
            avail.append("torch_cuda")
    except ImportError:
        pass
    return avail
