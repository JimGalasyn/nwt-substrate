"""The single, vendor-agnostic Steane decode + verdict.

Consumes only canonical `Counts` (see adapters/base.py) — bitstrings with
character index i == qubit/cbit i (index 0 leftmost). Because bit order is
interpreted exactly once (in each adapter's normalization), this module cannot
be vendor-sensitive: the same `Counts` decodes identically regardless of which
SDK produced it.

Replaces the ~26 hand-rolled decoders scattered across analysis/.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from nwt_substrate.isa.constants import RANK_SO7, N_VERTICES_K7
from .spec import STEANE_STAB, STEANE_LZ

# Hamming(7,4) syndrome position -> Fano point (mirror of the legacy pipeline)
QUBIT_TO_FANO = {0: "E3", 1: "E2", 2: "F1", 3: "E1", 4: "F2", 5: "F3", 6: "P"}


def _bits(bitstring: str) -> list[int]:
    """Canonical bitstring (index 0 leftmost) -> list of ints indexed by qubit."""
    return [int(ch) for ch in bitstring]


def fano(data_bits: list[int]) -> str:
    """Map 7 data-qubit parities over the 3 stabilizer supports to a Fano point.
    
    Parameters
    ----------
    data_bits : list[int]
        List of 7 measurement outcomes (0 or 1).
        
    Returns
    -------
    str
        Fano point label ('I', 'E1', 'E2', 'E3', 'F1', 'F2', 'F3') or 'I' if no error.
        
    Raises
    ------
    ValueError
        If data_bits does not have exactly 7 elements.
    """
    if len(data_bits) != N_VERTICES_K7:
        raise ValueError(f"Expected 7 data bits, got {len(data_bits)}")
    s = [sum(data_bits[i] for i in STEANE_STAB[k]) % 2 for k in range(RANK_SO7)]
    pos = sum(s[k] * _HAMMING_BIT_WEIGHTS[k] for k in range(RANK_SO7))
    return "I" if pos == 0 else QUBIT_TO_FANO[pos - 1]


def logical_z(data_bits: list[int]) -> int:
    return +1 if sum(data_bits[i] for i in STEANE_LZ) % 2 == 0 else -1


# --------------------------------------------------------------- distributions
def destructive_dists(z_counts: dict[str, int], x_counts: dict[str, int]) -> tuple[object, object]:
    """From the two destructive-readout registers (each {bitstring: n}) return
    (x_dist over x_fano, z_dist over (z_fano, logical_z)).
    
    Parameters
    ----------
    z_counts : dict[str, int]
        Dictionary mapping Z-basis readout bitstrings to outcome counts.
    x_counts : dict[str, int]
        Dictionary mapping X-basis readout bitstrings to outcome counts.
        
    Returns
    -------
    tuple
        (x_dist, z_dist) where x_dist is a Counter over x_fano points and
        z_dist is a Counter over (z_fano, logical_z) tuples.
    """
    x_dist: Counter = Counter()
    z_dist: Counter = Counter()
    for bs, n in z_counts.items():
        b = _bits(bs)
        z_dist[(fano(b), logical_z(b))] += n
    for bs, n in x_counts.items():
        x_dist[fano(_bits(bs))] += n
    return x_dist, z_dist


def ancilla_dist(counts: object) -> "Counter":
    """From an ancilla-scheme Counts (registers c_syn[6], c_lz[7]) return a joint
    distribution over ((x_fano, z_fano), logical_z).
    
    Parameters
    ----------
    counts : Counts
        Canonical Counts object with c_syn and c_lz registers.
        
    Returns
    -------
    Counter
        Joint distribution over ((x_fano, z_fano), logical_z) outcomes.
    """
    i_syn = counts.register_order.index("c_syn")
    i_lz = counts.register_order.index("c_lz")
    dist: Counter = Counter()
    for key, n in counts.data.items():
        syn, lz = key[i_syn], _bits(key[i_lz])
        xpos = sum(int(syn[k]) * _HAMMING_BIT_WEIGHTS[k] for k in range(RANK_SO7))
        zpos = sum(int(syn[RANK_SO7 + k]) * _HAMMING_BIT_WEIGHTS[k] for k in range(RANK_SO7))
        xf = "I" if xpos == 0 else QUBIT_TO_FANO[xpos - 1]
        zf = "I" if zpos == 0 else QUBIT_TO_FANO[zpos - 1]
        dist[((xf, zf), logical_z(lz))] += n
    return dist


# ----------------------------------------------------------------- verdicts
@dataclass
class HalfResult:
    modal: object
    modal_prob: float
    predicted: object
    pred_bin_prob: float
    passed: bool


@dataclass
class Verdict:
    label: str
    scheme: str
    passed: bool
    x_half: HalfResult | None = None
    z_half: HalfResult | None = None
    joint: HalfResult | None = None
    detail: dict | None = None


def verdict_destructive(label, x_dist, z_dist, pred_x_fano, pred_z_fano, pred_lz) -> Verdict:
    n_x = sum(x_dist.values()) or 1
    n_z = sum(z_dist.values()) or 1
    modal_x, mxn = x_dist.most_common(1)[0]
    (modal_zf, modal_zlz), mzn = z_dist.most_common(1)[0]
    x_pass = modal_x == pred_x_fano

    detail = None
    if pred_lz == 0:  # h7-style control: require modal syndrome + ~50/50 LZ split
        n_plus = sum(v for (zf, lz), v in z_dist.items() if zf == pred_z_fano and lz == +1)
        n_minus = sum(v for (zf, lz), v in z_dist.items() if zf == pred_z_fano and lz == -1)
        tot = n_plus + n_minus
        pct_plus = (n_plus / tot) if tot else None
        z_pass = (modal_zf == pred_z_fano) and (pct_plus is not None and 0.35 <= pct_plus <= 0.65)
        detail = {"n_plus": n_plus, "n_minus": n_minus, "pct_plus": pct_plus}
    else:
        z_pass = (modal_zf == pred_z_fano) and ((pred_lz not in (+1, -1)) or modal_zlz == pred_lz)

    return Verdict(
        label=label, scheme="destructive", passed=bool(x_pass and z_pass), detail=detail,
        x_half=HalfResult(modal_x, mxn / n_x, pred_x_fano, x_dist.get(pred_x_fano, 0) / n_x, bool(x_pass)),
        z_half=HalfResult((modal_zf, modal_zlz), mzn / n_z, (pred_z_fano, pred_lz),
                          z_dist.get((pred_z_fano, pred_lz), 0) / n_z if pred_lz in (+1, -1) else float("nan"),
                          bool(z_pass)),
    )


def verdict_ancilla(label, dist, pred_x_fano, pred_z_fano, pred_lz) -> Verdict:
    n = sum(dist.values()) or 1
    modal_key, mn = dist.most_common(1)[0]
    pred_key = ((pred_x_fano, pred_z_fano), pred_lz)
    if pred_lz == 0:
        syn = (pred_x_fano, pred_z_fano)
        n_plus = dist.get((syn, +1), 0)
        n_minus = dist.get((syn, -1), 0)
        tot = n_plus + n_minus
        pct_plus = (n_plus / tot) if tot else None
        passed = (modal_key[0] == syn) and (pct_plus is not None and 0.35 <= pct_plus <= 0.65)
        detail = {"n_plus": n_plus, "n_minus": n_minus, "pct_plus": pct_plus}
    else:
        passed = modal_key == pred_key
        detail = None
    return Verdict(
        label=label, scheme="ancilla", passed=bool(passed), detail=detail,
        joint=HalfResult(modal_key, mn / n, pred_key, dist.get(pred_key, 0) / n, bool(passed)),
    )
