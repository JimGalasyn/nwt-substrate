"""Standalone substrate PREDICTIONS — O10 DAG-cit-readout for the constants stack.

Implements the standalone-output rung of L. Leighton's O10 "Ladder Derivation
Protocol" (the DAG cit-readout specialization), whose proof-order edges are:

    analytic proof -> exact symbolic parent -> exact witness evaluator
                   -> standalone Python output -> CODATA-2018 compatibility witness

The edges are ONE-WAY: a downstream witness may TEST a prediction but may never
become upstream source authority for it (no decimal, fit, or endpoint convention
moves backward across an edge).  Concretely, this module emits ONLY values
*derived* from the substrate structure — α is the substrate closed form
``1/(25π√3 + 1)`` (``isa.ALPHA_SUBSTRATE``), and every quantity below is a pure
function of the K_7 / Spin(7) integers.  Only *dimensionless* quantities appear,
so there is no scale input either (no m_e, no M_Pl): nothing measured can enter.

The measured values are QUARANTINED in ``REFERENCE`` (witness-only) and are never
read by ``predictions()``.  Per the O10 constants rule the witness layer is
CODATA-2018 / PDG — the last fully *measured* set; post-SI2019 hard-defined
constant values are excluded even from witness use (comparing a prediction to a
*defined* value is circular).  Emit the two streams separately and let the diff
be the cit readout (cit = transversal witness-invariance across the path):

    python -m nwt_substrate.benchmarks.predict             > predicted.txt
    python -m nwt_substrate.benchmarks.predict --reference  > measured.txt
    diff -u predicted.txt measured.txt

Each differing line is a place the substrate departs from measurement — with
nowhere for a measured input to have crept into the derivation.  (Conflating a
measured *effective* angle with the predicted *on-shell* angle caused the v0.3.1
sin²θ_W bug; this separation prevents recurrence — and the diff makes plain that
sin²θ_W is a ~3.5% leading-order angle, not a sub-1% result.)
"""

from __future__ import annotations

import math

# Import ONLY the substrate structure.  Nothing here is a measured value:
# ALPHA_SUBSTRATE = 1/(25π√3 + 1) is the substrate-derived fine-structure constant.
from ..isa.constants import (
    ALPHA_SUBSTRATE,
    N_EDGES_K7,
    N_VERTICES_K7,
    NLO_VERTEX_COEFFICIENT,
    SPINOR_VECTOR_RATIO,
)


def predictions() -> dict[str, float]:
    """Dimensionless substrate predictions — pure functions of the structure.

    No argument, no measured input, no scale: a standalone derivation.

    ORDER PIN (Auditor verdict `2026-07-12-constants-provenance-disputes`,
    CL-2): the frozen m_e/M_Pl claim is the NLO form (8/7)(1+α/7)·α^(21/2) —
    the externally audited L4(a) form.  The NNLO α² term is RETIRED from
    claim status: Paper 17 documents its coefficient as computed from the
    CODATA target before the "structural" integer 21/8 was selected
    (target-selection), and substrate-pure it buys nothing (−65 ppm NLO vs
    +75 ppm NNLO on the ratio, both ~6σ dead against the 11.5 ppm witness
    bar).  ``isa.NNLO_BRACKET_COEFFICIENT`` remains as code; it may never be
    cited as the claim, and no post-freeze order change can revive it.
    """
    a = ALPHA_SUBSTRATE                                  # 1/(25π√3 + 1)
    nlo_bracket = 1.0 + NLO_VERTEX_COEFFICIENT * a
    return {
        "inv_alpha":      1.0 / a,                       # 25π√3 + 1
        "sin2_theta_W":   (2.0 + a) / 9.0,               # (2 + α)/9
        "cabibbo_lambda": math.sqrt(N_VERTICES_K7 * a),  # √(7α)
        "eta_B":          3.0 * a ** 4 / 14.0,           # baryon asymmetry 3α⁴/14
        "m_e_over_M_Pl":  SPINOR_VECTOR_RATIO * a ** (N_EDGES_K7 / 2.0) * nlo_bracket,
    }


# Measured / experimental reference values — QUARANTINED.  These are NOT imported
# by predictions(); they exist only so the two streams can be diffed externally.
# Per O10 the witness layer is CODATA-2018 / PDG / Planck (the last fully *measured*
# set); post-SI2019 hard-defined constant values are excluded even as witnesses.
# Sources: CODATA-2018 (α⁻¹ = 137.035999084), PDG (sin²θ_W effective, λ),
# Planck-2018 (η_B), CODATA-2018 (m_e/M_Pl).
REFERENCE: dict[str, float] = {
    "inv_alpha":      137.035999084,
    # sin²θ_W: PDG *effective* leptonic angle. The substrate (2+α)/9 is a
    # LEADING-ORDER angle ~3.5% below it (the on-shell 1−M_W²/M_Z² scheme is
    # much closer) — kept here precisely so the diff shows a genuine deviation,
    # not a curated all-green table. (Conflating these two angles caused the
    # v0.3.1 sin²θ_W benchmark bug; the separation here prevents recurrence.)
    "sin2_theta_W":   0.23122,
    "cabibbo_lambda": 0.22500,
    "eta_B":          6.10e-10,
    "m_e_over_M_Pl":  4.18540e-23,
}


def _emit(values: dict[str, float]) -> None:
    for key in sorted(values):
        print(f"{key:18s} = {values[key]:.9g}")


def main(argv: list[str] | None = None) -> int:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    _emit(REFERENCE if "--reference" in argv else predictions())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
