"""
Multi-particle composite states: molecular and Hopf-linked.

Two distinct binding mechanisms in NWT, formalised here as the
substrate's composition operations:

  MOLECULAR via CONNECTED SUM (#)  -- Yukawa-style, few-MeV binding:
    Schubert prime factorisation makes (Knots, #) a free commutative
    monoid; LO mass rule m(K_A # K_B) = m(K_A) + m(K_B).  Empirical
    test (memory:connected-sum-composition-law) gives sub-percent
    Paper-6 residuals for deuteron, X(3872), Tcc, Pc(4312/4440/4457),
    and the alpha-cluster ladder up to 0.61% RMS via Phase 2/3
    cluster-geometry refinement.  Use compose(..., op="#").

  HOPF-LINKED -- topological, ~Lambda_QCD = 313 MeV per crossing:
    Heavy nuclei (linked trefoils, Paper 15) and predicted exotic
    states (~395 Hopf-linked pentaquarks across 12 flavour sectors).
    Strong binding regime where # fails (e.g. d*(2380) breaks the #
    rule by ~84 MeV, isolating it as Hopf-linked).
    Use compose(..., op="hopf", binding=-LAMBDA_QCD_MEV * n_crossings).

The two operations are LAYERED, not competing: # supplies the carrier-
knot topology (which prime factors are sewn together), Hopf-link
supplies the binding-energy magnitude.  Most nuclei have both; light
near-threshold molecules have only #.
"""

from .connected_sum import (
    CompositeParticle,
    compose,
    decompose,
)

LAMBDA_QCD_MEV = 313.0  # NWT linking-energy scale (Paper 15 nuclear binding)

__all__ = [
    "CompositeParticle",
    "compose",
    "decompose",
    "LAMBDA_QCD_MEV",
]
