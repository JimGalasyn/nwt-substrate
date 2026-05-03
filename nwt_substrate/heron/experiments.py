"""
NWT experiments on IBM Heron — registry of run + proposed experiments.

The substrate's central object is the K_7 graph state.  Each experiment
listed here is a substrate-prediction test that has either been run on
real Heron hardware, or is queued for a future run.

The framework: any Heron experiment is a substrate-algebra computation
realized in qubits.  Substrate predictions about K_7 stabilizer
expectations, entanglement entropy bipartitions, syndrome distributions,
or vortex-defect lifetimes (Experiment 9) become specific quantum-circuit
runs whose outcomes are *direct* observations of substrate physics —
not analog simulations.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class HeronExperiment:
    number: int
    name: str
    status: str               # "run", "proposed", "in-flight"
    description: str
    backend: Optional[str] = None
    paper_section: Optional[str] = None
    notes: str = ""

    def __repr__(self):
        return f"HeronExperiment({self.number}, {self.name!r}, status={self.status!r})"


EXPERIMENTS = [
    HeronExperiment(
        number=1,
        name="K_7 graph-state preparation + stabilizer check",
        status="run",
        description=(
            "Prepare |K_7> via H + CZ over all 21 K_7 edges, then measure "
            "each of 7 stabilizers S_v.  All seven should yield +1 within "
            "device noise.  Substrate's most basic existence test."
        ),
        backend="ibm_marrakesh",
        paper_section="Paper 17 §13 / Paper 19 W1+W2",
    ),
    HeronExperiment(
        number=2,
        name="Z-basis joint distribution + entanglement entropy",
        status="run",
        description=(
            "Prep |K_7>, measure all 7 qubits in Z basis, sample to "
            "reconstruct the bipartite entanglement entropy across "
            "every (k, 7-k) cut for k in {1,2,3}.  Substrate predicts "
            "log(2) entropy per cut for the maximally entangled K_7."
        ),
        backend="ibm_marrakesh",
        paper_section="Paper 19 W1+W2",
    ),
    HeronExperiment(
        number=3,
        name="X-basis tomography",
        status="run",
        description=(
            "Same as Experiment 2 but in X basis (rotated Hadamard).  "
            "Combined with Z and Y runs, full single-qubit tomography."
        ),
        backend="ibm_marrakesh",
        paper_section="Paper 19 W1+W2",
    ),
    HeronExperiment(
        number=4,
        name="3-body <Y_u Y_v Y_w> = 0 null test",
        status="run",
        description=(
            "Substrate predicts <Y_u Y_v Y_w> = 0 for any three "
            "non-Fano-line vertices on |K_7>.  Direct test of the "
            "PSL(2,7) Fano-line structure underpinning Paper 17."
        ),
        backend="ibm_marrakesh",
        paper_section="Paper 17 §13.3",
    ),
    HeronExperiment(
        number=5,
        name="Syndrome distribution test",
        status="run",
        description=(
            "Apply random local errors, measure syndrome, check that "
            "the distribution matches the K_7 stabilizer code prediction.  "
            "Probes the QEC interpretation of the substrate."
        ),
        backend="ibm_marrakesh",
        paper_section="Paper 17 §13.4 / Paper 19 W3",
    ),
    HeronExperiment(
        number=6,
        name="K_7 entanglement tomography (proposed)",
        status="proposed",
        description=(
            "Full 7-qubit tomography across all bipartitions to reconstruct "
            "the entanglement spectrum.  Compare to substrate prediction "
            "of ln(2) per cut in the trivial case, with corrections from "
            "Spin(7) -> SU(3) x SU(2)_int decomposition."
        ),
        paper_section="Paper 19 W6",
        notes="Higher shot count needed (~10^6 shots).",
    ),
    HeronExperiment(
        number=7,
        name="Substrate metrology (proposed)",
        status="proposed",
        description=(
            "Use |K_7> as a metrology resource — measure α via the K_7 "
            "Wilson amplitude on hardware, compare to CODATA.  Validates "
            "the n=21 -> m_e prediction (Paper 17) on a quantum device."
        ),
        paper_section="Paper 18 §X (TBD)",
        notes="Requires gate-level α-extraction protocol.",
    ),
    HeronExperiment(
        number=8,
        name="Two-universe slicing of |G_aug> (proposed)",
        status="proposed",
        description=(
            "Augment the K_7 graph state with extra qubits to represent "
            "two universes; slice the joint state and measure cross-"
            "universe correlations.  Substrate prediction: zero for "
            "decoupled, nonzero for entangled cosmologies."
        ),
        paper_section="Paper 19 W7 / Paper 22",
    ),
    HeronExperiment(
        number=9,
        name="Watch a neutron decay (proposed)",
        status="proposed",
        description=(
            "Encode a neutron's substrate state as a topology of stabilizers, "
            "evolve under a substrate-algebra Hamiltonian for one substrate "
            "lifetime, observe the decay pattern.  At Heron's energy scale "
            "this is a literal neutron decay — not a simulation.  Walk-phase "
            "4c established the substrate-algebraic |M|^2 needed."
        ),
        paper_section="Paper 20 (proposed)",
        notes=(
            "Requires ~16-20 qubits; encoding + evolution + measurement "
            "designed in checkpoint_2026-04-29.md."
        ),
    ),
    HeronExperiment(
        number=10,
        name="Muon decay on K_7 vacuum (run 2026-05-01)",
        status="run",
        backend="ibm_marrakesh",
        description=(
            "Prepare |K_7> on 7 qubits (substrate vacuum) plus 2 ancillae "
            "(muon at qubit 7, electron at qubit 8).  Apply a parameterised "
            "V-A XY-mixer between the ancillae at angle θ, sweep θ ∈ [0, π] "
            "in one Sampler.run call (parameterised batch).  Expected "
            "muon survival: P(muon at θ) = cos²(θ).  Substrate-physics "
            "narrative: muon decays via V-A (substrate's Cl(0,7) chirality "
            "structure); the K_7 vacuum is unchanged."
        ),
        paper_section="Paper 20 / electroweak shim demo",
        notes=(
            "Job d7qgmo4f3ras73b5orqg on ibm_marrakesh, 41.5s wall clock. "
            "9 qubits, 40 gates, depth 13, 4096 shots × 16 angles. "
            "Noise envelope fit: P(θ) = 0.963 × cos²(θ) + 0.020 "
            "(contrast 96%, floor 2%, RMS deviation 1.6% absolute, max 3.3%). "
            "Plateau values P(0)=0.9946, P(π)=0.9924 -- substrate K_7 vacuum "
            "preserved through 40-gate circuit.  Trough P(π/2)=0.027 vs "
            "theory 0 -- decoherence drift toward equilibration.  "
            "Data saved to analysis/output/heron_runs/exp10_muon_decay_2026-05-01.json."
        ),
    ),
    HeronExperiment(
        number=11,
        name="Sidereal A/B/C substrate-anisotropy probe (proposed)",
        status="proposed",
        description=(
            "Three identical K_7 stabilizer-measurement runs at sidereal "
            "times t_s, t_s + 12h sidereal, t_s + 24h sidereal (i.e. solar "
            "wallclock 0 / +11h58m / +23h56m).  A and C share the same lab-"
            "to-inertial-frame orientation (Earth has rotated 360 sidereal "
            "degrees); B is rotated ~180 degrees in the inertial frame from "
            "both.  Drift-corrected sidereal signal is "
            "sigma_v = (B - A) - (C - A)/2 per stabilizer S_v.  Substrate-"
            "monism prediction: if SO(8) -> Spin(7) -> 3+1 reduction picks "
            "an inertial-frame axis (parent-BH-set per Paper 22 / two-layer "
            "unification), B should differ from A,C by a coherent per-"
            "vertex pattern.  Detection threshold ~1% per vertex; expected "
            "magnitude unknown but bounded above by FLRW null-test signal "
            "magnitude (a few percent).  Null result tightens a substrate-"
            "isotropy upper bound; positive result identifies which graph "
            "vertex aligns with the substrate axis and motivates a full "
            "sidereal-modulation scan."
        ),
        paper_section="Paper 19 (substrate monism, anisotropy probe) / Paper 22 (BH cosmogenesis joint test)",
        notes=(
            "Total Heron compute time ~2 minutes (3 runs x 4096 shots).  "
            "Wallclock ~24 hours dominated by sidereal scheduling.  Fits "
            "10 min/month open-instance budget.  Use Experiment 1's "
            "circuit unchanged; only the timing changes.  Multi-site "
            "extension: simultaneous runs at geographically separated "
            "processors (e.g. ibm_marrakesh + ibm_kingston) act as a "
            "VLBI-style array probing the inertial-frame direction of "
            "the substrate anisotropy axis."
        ),
    ),
]


def list_experiments(status: Optional[str] = None) -> list[HeronExperiment]:
    """Return all experiments, or filter by status."""
    if status is None:
        return EXPERIMENTS
    return [e for e in EXPERIMENTS if e.status == status]


def experiment(n: int) -> HeronExperiment:
    """Look up a single experiment by number."""
    for e in EXPERIMENTS:
        if e.number == n:
            return e
    raise KeyError(f"No experiment numbered {n}")


def summary() -> dict:
    """Quick summary of the experiment registry."""
    by_status = {}
    for e in EXPERIMENTS:
        by_status.setdefault(e.status, []).append(e.number)
    return {
        "total": len(EXPERIMENTS),
        "by_status": {k: len(v) for k, v in by_status.items()},
        "run_numbers": by_status.get("run", []),
        "proposed_numbers": by_status.get("proposed", []),
    }
