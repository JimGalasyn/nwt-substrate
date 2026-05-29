"""O10 DAG cit-readout — the substrate constants derivation as a validated graph.

Operationalizes the DAG specialization of L. Leighton's O10 "Ladder Derivation
Protocol".  Where ``predict.py`` emits the standalone-output *rung*, this module
builds the whole directed acyclic graph it belongs to and checks the O10
invariants on it.

The proof-order stages (edges run strictly upward through them) are:

    STRUCTURAL  ->  SYMBOLIC  ->  EVALUATOR  ->  OUTPUT  ->  WITNESS
    (K_7/Spin(7)   (closed       (numeric       (standalone  (CODATA-2018 /
     integers)      form)         evaluator)     output)      PDG, sink-only)

O10 rules enforced here:
  * edges are ONE-WAY — no edge may point to an earlier stage (a witness or
    output sourcing upstream authority); STRUCTURAL->STRUCTURAL identities aside;
  * the graph is ACYCLIC at proof-authority level;
  * WITNESSES are sinks — CODATA-2018 / PDG, never a source (post-SI2019
    *defined* constants are excluded even as witnesses);
  * ``cit`` is the transversal witness-invariance readout: each OUTPUT must agree
    with its WITNESS within tolerance — and a failing edge is *marked as a
    defect*, never silently repaired (O10: "mark defective edges");
  * multi-parent convergence is the Horn-clause readout ``(s_1 ∧ … ∧ s_n) → W``;
  * a node reached by two directed paths is a *commuting diagram* — these are
    exactly isa's import-time structural identities (e.g. 21 = C(7,2) = 3·7).

``python -m nwt_substrate.benchmarks.o10`` prints the acceptance checklist, the
cit readout (with any defect edges marked), the structural-load ranking, and the
commutative-diagram checks.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum

from .predict import REFERENCE, predictions


class Stage(IntEnum):
    """Proof-order stages.  O10 edges run strictly upward (lower -> higher)."""
    STRUCTURAL = 0      # K_7/Spin(7) integers — source authority
    SYMBOLIC = 1        # exact closed-form parent
    EVALUATOR = 2       # numeric evaluator (isa / predict)
    OUTPUT = 3          # standalone Python output
    WITNESS = 4         # CODATA-2018 / PDG measured value — sink only


@dataclass(frozen=True)
class Node:
    name: str
    stage: Stage
    value: float | None = None
    note: str = ""
    commutative: bool = False     # True iff this node is a "two-paths-to-one-value"
    #                               identity (parent values must agree); multi-INPUT
    #                               formula nodes are Horn-clause conjunctions, not this.


@dataclass
class DerivationDAG:
    """A directed acyclic graph of one constants-derivation stack."""
    nodes: dict[str, Node] = field(default_factory=dict)
    edges: list[tuple[str, str]] = field(default_factory=list)

    def add(self, name: str, stage: Stage, value: float | None = None,
            note: str = "", commutative: bool = False) -> str:
        if name not in self.nodes:
            self.nodes[name] = Node(name, stage, value, note, commutative)
        return name

    def link(self, src: str, dst: str) -> None:
        self.edges.append((src, dst))

    def _children(self, name: str) -> list[str]:
        return [d for s, d in self.edges if s == name]

    def _parents(self, name: str) -> list[str]:
        return [s for s, d in self.edges if d == name]

    # ---- O10 structural invariants ----

    def backward_edges(self) -> list[tuple[str, str]]:
        """Edges violating the one-way rule: dst at a strictly earlier stage than
        src (an output/witness sourcing upstream).  Same-stage edges are allowed
        (STRUCTURAL->STRUCTURAL integer identities)."""
        return [(s, d) for s, d in self.edges
                if self.nodes[d].stage < self.nodes[s].stage]

    def is_acyclic(self) -> bool:
        indeg = {n: 0 for n in self.nodes}
        for _, d in self.edges:
            indeg[d] += 1
        queue = [n for n, k in indeg.items() if k == 0]
        seen = 0
        while queue:
            n = queue.pop()
            seen += 1
            for c in self._children(n):
                indeg[c] -= 1
                if indeg[c] == 0:
                    queue.append(c)
        return seen == len(self.nodes)

    def witnesses_are_sinks(self) -> bool:
        """No edge leaves a WITNESS node (witnesses test, never source)."""
        return not any(self.nodes[s].stage == Stage.WITNESS for s, _ in self.edges)

    # ---- readouts ----

    def reachable_outputs(self, name: str) -> set[str]:
        seen: set[str] = set()
        stack, outs = [name], set()
        while stack:
            for c in self._children(stack.pop()):
                if c in seen:
                    continue
                seen.add(c)
                if self.nodes[c].stage == Stage.OUTPUT:
                    outs.add(c)
                stack.append(c)
        return outs

    def load_ranking(self) -> list[tuple[str, int]]:
        """Upstream nodes ranked by # downstream OUTPUT nodes they reach — the
        O10/Horn-clause load-bearing measure (cf. M. Wende's structural
        criticality).  The most load-bearing node organises the most observables."""
        internal = [n for n, nd in self.nodes.items()
                    if nd.stage in (Stage.STRUCTURAL, Stage.SYMBOLIC, Stage.EVALUATOR)]
        load = [(n, len(self.reachable_outputs(n))) for n in internal]
        return sorted((kv for kv in load if kv[1] > 0), key=lambda kv: (-kv[1], kv[0]))

    def horn_frontier(self, target: str) -> list[str]:
        """STRUCTURAL premises {s_i} with a directed path into ``target`` —
        the Horn-clause readout (s_1 ∧ … ∧ s_n) → target."""
        seen: set[str] = set()
        stack, frontier = [target], set()
        while stack:
            for p in self._parents(stack.pop()):
                if p in seen:
                    continue
                seen.add(p)
                if self.nodes[p].stage == Stage.STRUCTURAL:
                    frontier.add(p)
                stack.append(p)
        return sorted(frontier)

    def cit_readout(self, tol: float = 0.01) -> list[dict]:
        """Per OUTPUT->WITNESS edge: derived vs witness + admissibility.
        cit holds on the path iff the output agrees with its witness within tol."""
        rows = []
        for s, d in self.edges:
            if self.nodes[s].stage == Stage.OUTPUT and self.nodes[d].stage == Stage.WITNESS:
                pv, wv = self.nodes[s].value, self.nodes[d].value
                dev = abs(pv - wv) / abs(wv)
                rows.append({"output": s, "witness": d, "predicted": pv,
                             "measured": wv, "rel_dev": dev, "admissible": dev <= tol})
        return rows

    def cit_defects(self, tol: float = 0.01) -> list[str]:
        """Output names whose cit edge fails at ``tol`` (marked, not repaired)."""
        return [r["output"] for r in self.cit_readout(tol) if not r["admissible"]]

    def commutative_checks(self) -> list[dict]:
        """Identity nodes (``commutative=True``) reached by two directed paths:
        the parent values must agree.  These mirror isa's import-time structural
        identities (e.g. 21 = C(7,2) = 3·7).  Multi-*input* formula nodes are NOT
        commuting diagrams — they are Horn-clause conjunctions of distinct
        premises — so they are excluded here."""
        rows = []
        for n, nd in self.nodes.items():
            if not nd.commutative:
                continue
            ps = self._parents(n)
            vals = [self.nodes[p].value for p in ps if self.nodes[p].value is not None]
            rows.append({"node": n, "parents": ps, "values": vals,
                         "commutes": (len({round(v, 9) for v in vals}) == 1) if vals else True})
        return rows

    def acceptance_checklist(self) -> dict[str, bool]:
        """The O10 structural invariants (all must hold for an admissible DAG).
        cit *content* is reported separately by cit_readout/cit_defects, because
        a defect edge is marked rather than failing the structure."""
        return {
            "all_edges_directed_one_way": not self.backward_edges(),
            "acyclic_at_proof_authority": self.is_acyclic(),
            "witnesses_are_sinks": self.witnesses_are_sinks(),
            "commutative_diagrams_agree": all(r["commutes"] for r in self.commutative_checks()),
        }


def build_constants_dag() -> DerivationDAG:
    """The constants stack as an O10 DAG: α as the load-bearing root feeding the
    five dimensionless predictions, plus two commutative-diagram identities."""
    from ..isa.constants import (
        ALPHA_SUBSTRATE, DIM_OCTONION, DIM_S_SPIN7, DIM_V_SPIN7,
        N_EDGES_K7, N_VERTICES_K7, RANK_SO7,
    )
    pred, ref = predictions(), REFERENCE
    g = DerivationDAG()

    # α — derived from the closed form 1/(25π√3+1), and once evaluated it serves
    # as the upstream authority for the other four predictions (O10: the readout
    # of one ladder is the structural input of the next).  Both sit at the
    # STRUCTURAL stage; the form->α edge is a STRUCTURAL->STRUCTURAL derivation.
    g.add("form:25π√3+1", Stage.STRUCTURAL, note="α⁻¹ closed form (no measured input)")
    g.add("α", Stage.STRUCTURAL, ALPHA_SUBSTRATE, "isa.ALPHA_SUBSTRATE")
    g.link("form:25π√3+1", "α")

    # structural integers used by the prediction symbolic forms
    structurals = {
        "N_VERTICES_K7": (float(N_VERTICES_K7), "|V(K_7)| = 7"),
        "N_EDGES_K7":    (float(N_EDGES_K7), "|E(K_7)| = 21 (Wilson exponent ×2)"),
        "8/7":           (DIM_S_SPIN7 / DIM_V_SPIN7, "spinor/vector ratio"),
        "int:2": (2.0, ""), "int:3": (3.0, ""), "int:9": (9.0, ""), "int:14": (14.0, ""),
    }
    for nm, (val, note) in structurals.items():
        g.add(nm, Stage.STRUCTURAL, val, note)

    # the five dimensionless predictions: form -> output -> witness
    specs = {
        "inv_alpha":      ("25π√3+1",            ["α"]),
        "sin2_theta_W":   ("(2+α)/9",            ["α", "int:2", "int:9"]),
        "cabibbo_lambda": ("√(7α)",              ["α", "N_VERTICES_K7"]),
        "eta_B":          ("3α⁴/14",             ["α", "int:3", "int:14"]),
        "m_e_over_M_Pl":  ("(8/7)·α^(21/2)·NNLO", ["α", "8/7", "N_EDGES_K7"]),
    }
    for key, (form, parents) in specs.items():
        sym, out, wit = f"sym:{key}", key, f"wit:{key}"
        g.add(sym, Stage.SYMBOLIC, note=form)
        g.add(out, Stage.OUTPUT, pred[key], "standalone derived")
        g.add(wit, Stage.WITNESS, ref[key], "CODATA-2018 / PDG")
        for p in parents:
            g.link(p, sym)
        g.link(sym, out)
        g.link(out, wit)

    # commutative-diagram identity: 21 reached two ways (isa import-time assert).
    g.add("C(7,2)", Stage.STRUCTURAL, N_VERTICES_K7 * (N_VERTICES_K7 - 1) / 2, "7·6/2")
    g.add("RANK·DIM_V", Stage.STRUCTURAL, RANK_SO7 * DIM_V_SPIN7, "3·7")
    g.add("id:21", Stage.SYMBOLIC, float(N_EDGES_K7), "N_EDGES_K7 — two paths must agree",
          commutative=True)
    g.link("C(7,2)", "id:21")
    g.link("RANK·DIM_V", "id:21")
    # and 8 = DIM_OCTONION = DIM_S_SPIN7
    g.add("DIM_OCTONION", Stage.STRUCTURAL, float(DIM_OCTONION), "octonion dim")
    g.add("DIM_S_SPIN7", Stage.STRUCTURAL, float(DIM_S_SPIN7), "Spin(7) spinor dim")
    g.add("id:8", Stage.SYMBOLIC, float(DIM_OCTONION), "octonion = spinor dim",
          commutative=True)
    g.link("DIM_OCTONION", "id:8")
    g.link("DIM_S_SPIN7", "id:8")
    return g


def main(argv: list[str] | None = None) -> int:
    g = build_constants_dag()
    lines = [f"O10 DAG cit-readout — {len(g.nodes)} nodes, {len(g.edges)} edges", ""]

    lines.append("Acceptance checklist (O10 structural invariants):")
    for check, ok in g.acceptance_checklist().items():
        lines.append(f"  [{'PASS' if ok else 'FAIL'}]  {check}")

    lines += ["", "cit readout (output vs CODATA-2018/PDG witness; tol 1%):"]
    for r in g.cit_readout(tol=0.01):
        mark = "ok " if r["admissible"] else "DEFECT"
        lines.append(f"  [{mark}]  {r['output']:16s} {r['predicted']:.6g}  "
                     f"vs {r['measured']:.6g}  ({r['rel_dev']*100:.2f}%)")
    defects = g.cit_defects(tol=0.01)
    if defects:
        lines.append(f"  marked defect edge(s): {', '.join(defects)} "
                     f"(open, not repaired — e.g. sin²θ_W is a leading-order angle)")

    lines += ["", "Structural-load ranking (nodes by # observables reached):"]
    for n, load in g.load_ranking():
        lines.append(f"  {n:16s} reaches {load} output(s)")

    lines += ["", "Commutative diagrams (two directed paths must agree):"]
    for r in g.commutative_checks():
        lines.append(f"  [{'commutes' if r['commutes'] else 'BREAKS'}]  {r['node']}"
                     f"  via {r['parents']}  = {r['values']}")

    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
