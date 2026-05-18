"""
Paley (q, k, λ) bifold construction of K_q for q ≡ 3 (mod 4).

The Paley design is the Steiner system constructed from the quadratic
residues of F_q. For q = 7, this is the (7, 3, 1) Paley design, equal
to the Fano plane PG(2, 2) — and identical (up to vertex relabeling)
to the substrate K_7 structure used throughout NWT.

Substrate framing
-----------------
The Paley (7, 3, 1) design is the **AGL(1, 7)-equivariant unique
bifold organization of K_7's 21 edges** into Fano triangulations. Its
automorphism group AGL(1, 7) = Z_7 ⋊ F_7* has order 42 = 6 × 7.

NWT's Z_3 × Z_2 partition of K_7 vertices (Paper 19 §V1, Closure 2)
is conjugate to the Paley QR / NR partition under S_7; the two
framings describe the **same** discrete substrate-algebraic structure.

The cosmogenic parent BH spin axis spontaneously breaks AGL(1, 7) → Z_3:

  AGL(1, 7) order 42
        ↓ vertex 0 = spin-axis polar direction (Z_7 translation choice)
        ↓ Z_2 polarity = parent BH spin sign (QR vs NR / Type-A vs Type-B)
  Residual symmetry order 42 / 14 = **3 = Z_3** = three-fermion-generation cycle

This is the discrete-substrate-vertex picture of Paper 20's continuous
G_2 → SU(3) breaking. The two layers describe the same cosmogenic
mechanism; both leave the same Z_3 residual.

Coding-theory connection
------------------------
The Paley (7, 3, 1) design is the foundation of the Hamming (7, 4, 3)
code, the Steane [[7, 1, 3]] quantum error-correcting code (used in
NWT's Heron Exp 11 work), and the Golay (24, 12, 8) code family.

NWT's existing K_7 stabilizer-code work IS Paley-structured. This
module makes the equivalence explicit.

The Paley-Hadamard family at q ≡ 3 (mod 4) — q = 7, 11, 19, 23, 31,
... — gives the natural cross-substrate extension. Each q is a
candidate Cl(0, k) substrate with its own AGL(1, q) symmetry and
residual subgroup after spin-axis breaking. See
``paley_design`` for the general construction.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable


def is_prime(n: int) -> bool:
    """True if n is a prime."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    for k in range(3, int(math.isqrt(n)) + 1, 2):
        if n % k == 0:
            return False
    return True


def quadratic_residues(p: int) -> set[int]:
    """Quadratic residues of F_p^*.

    For p = 7: {1, 2, 4}. These are the nonzero squares mod p.
    """
    if not is_prime(p):
        raise ValueError(f"p must be a prime, got {p}")
    return {(k * k) % p for k in range(1, p)}


def non_residues(p: int) -> set[int]:
    """Non-residues of F_p^*.

    For p = 7: {3, 5, 6}. The complement of QR in F_p^*.
    """
    return set(range(1, p)) - quadratic_residues(p)


def is_paley_admissible(q: int) -> bool:
    """True if q is admissible for the Paley (q, (q-1)/2, (q-3)/4) construction.

    Requires q ≡ 3 (mod 4) and q prime. Includes q = 7, 11, 19, 23, 31, ...
    """
    return is_prime(q) and q % 4 == 3


@dataclass(frozen=True)
class PaleyDesign:
    """The Paley (q, k, λ) design constructed from F_q quadratic residues."""

    q: int
    """The field order. Must be prime with q ≡ 3 (mod 4)."""

    k: int
    """Block size = (q - 1) / 2."""

    lambda_: int
    """Every 2-element subset of F_q lies in exactly λ blocks. = (q - 3) / 4."""

    points: tuple[int, ...]
    """F_q = {0, 1, ..., q - 1}."""

    qr: frozenset[int]
    """Quadratic residues of F_q^*."""

    nr: frozenset[int]
    """Non-residues of F_q^*."""

    qr_blocks: tuple[frozenset[int], ...]
    """QR-Paley blocks: translates of QR around F_q. q blocks total."""

    nr_blocks: tuple[frozenset[int], ...]
    """NR-Paley blocks: translates of NR around F_q. q blocks total."""

    def num_blocks(self) -> int:
        """Total number of blocks in QR-Paley design = q."""
        return len(self.qr_blocks)


def paley_design(q: int) -> PaleyDesign:
    """Construct the Paley (q, k, λ) design from F_q quadratic residues.

    Parameters
    ----------
    q : int
        Prime with q ≡ 3 (mod 4). Common substrate-relevant values:
        7, 11, 19, 23, 31.

    Returns
    -------
    PaleyDesign

    Examples
    --------
    For q = 7: the (7, 3, 1) Paley design = PG(2, 2) Fano plane =
    Steiner triple system on 7 points. AGL(1, 7) is its full
    automorphism group, order 42.

        >>> design = paley_design(7)
        >>> design.k
        3
        >>> design.lambda_
        1
        >>> sorted(design.qr)
        [1, 2, 4]
        >>> design.num_blocks()
        7
    """
    if not is_paley_admissible(q):
        raise ValueError(
            f"q = {q} not Paley-admissible; need prime with q ≡ 3 (mod 4)"
        )

    k = (q - 1) // 2
    lambda_ = (q - 3) // 4

    qr = quadratic_residues(q)
    nr = non_residues(q)

    # QR-Paley blocks: translates of QR around F_q. Each block is
    # QR + i = {(r + i) mod q : r in QR} for i in F_q.
    qr_blocks = tuple(
        frozenset((r + i) % q for r in qr) for i in range(q)
    )
    nr_blocks = tuple(
        frozenset((s + i) % q for s in nr) for i in range(q)
    )

    return PaleyDesign(
        q=q,
        k=k,
        lambda_=lambda_,
        points=tuple(range(q)),
        qr=frozenset(qr),
        nr=frozenset(nr),
        qr_blocks=qr_blocks,
        nr_blocks=nr_blocks,
    )


# ---------------------------------------------------------------------------
# AGL(1, q) automorphism group and symmetry-breaking analysis
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AGLBreakingAnalysis:
    """Symmetry-breaking decomposition of AGL(1, q) → residual under
    parent-BH-spin-axis-style breaking.
    """

    q: int
    """Field order."""

    agl_order: int
    """|AGL(1, q)| = q · (q - 1)."""

    vertex_0_choices: int
    """Z_q translation factor — choice of "polar vertex" / spin-axis direction."""

    polarity_choices: int
    """Z_2 polarity factor — choice of QR vs NR (parent BH spin sign)."""

    broken_subgroup_order: int
    """Product of all broken-direction factors."""

    residual_order: int
    """|residual unbroken subgroup| = agl_order / broken_subgroup_order."""

    residual_label: str
    """Name of the residual subgroup, e.g. 'Z_3' for q = 7."""


def agl_breaking_analysis(q: int) -> AGLBreakingAnalysis:
    """Compute the AGL(1, q) → residual breaking under the
    "spin-axis polar + Z_2 polarity" pattern (parent BH cosmogenic
    breaking).

    For q = 7: AGL(1, 7) order 42 → Z_3 residual via 7 × 2 = 14 broken.
    This Z_3 is the three-fermion-generation cycle in NWT's Paper 20.

    For q = 11: AGL(1, 11) order 110 → Z_5 residual via 11 × 2 = 22 broken.
    Suggests a 5-generation-class substrate for Cl(0, k) with q = 11.

    For q = 19: AGL(1, 19) order 342 → Z_9 residual via 19 × 2 = 38 broken.

    Each q gives a specific residual subgroup whose order is the
    "generation count" prediction for that substrate.
    """
    if not is_paley_admissible(q):
        raise ValueError(
            f"q = {q} not Paley-admissible; need prime with q ≡ 3 (mod 4)"
        )

    agl_order = q * (q - 1)
    broken = q * 2
    residual = agl_order // broken

    # Label small residuals
    labels = {3: "Z_3", 5: "Z_5", 9: "Z_9", 11: "Z_11", 15: "Z_15", 21: "Z_21"}
    label = labels.get(residual, f"Z_{residual}")

    return AGLBreakingAnalysis(
        q=q,
        agl_order=agl_order,
        vertex_0_choices=q,
        polarity_choices=2,
        broken_subgroup_order=broken,
        residual_order=residual,
        residual_label=label,
    )


# ---------------------------------------------------------------------------
# Public convenience: K_7 Paley bifold (NWT's substrate case)
# ---------------------------------------------------------------------------

def k7_paley_bifold() -> PaleyDesign:
    """The Paley (7, 3, 1) bifold on K_7 — NWT's substrate case.

    Returns the Paley (7, 3, 1) design. QR = {1, 2, 4}, NR = {3, 5, 6};
    each QR translate is a 3-element block of the Fano plane PG(2, 2).
    AGL(1, 7) is its automorphism group, order 42.

    Maps onto NWT's Z_3 × Z_2 ⊂ Aut(K_7):
    - Z_3 fixes vertex 0 (polar vertex P) and cycles QR / NR orbits
    - Z_2 = additive negation x ↦ -x swaps QR and NR (parent BH spin sign)
    """
    return paley_design(7)


__all__ = [
    "is_prime",
    "quadratic_residues",
    "non_residues",
    "is_paley_admissible",
    "PaleyDesign",
    "paley_design",
    "AGLBreakingAnalysis",
    "agl_breaking_analysis",
    "k7_paley_bifold",
]
