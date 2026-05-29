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
commutative-diagram checks.  Add ``--redundancy`` for M. Wende's
derivation-route-redundancy readout: how many *independent* routes converge on
each answer, the single points of failure that survive removal, and node
*criticality* (outputs ungrounded if a node is removed) beside *load* (outputs
reached) — the leave-one-route-out dual of the load ranking.  ``--suite
--sensitivity --redundancy`` runs it over the 38-benchmark coupling graph.
"""

from __future__ import annotations

import re
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
    dev: float | None = None      # OUTPUT relative deviation from its witness (fraction),
    #                               set directly for suite benchmarks; None -> computed
    #                               from values, or qualitative (no metric).
    kind: str = ""                # provenance of dev: ppm / pct / exact / score / qualitative


@dataclass
class DerivationDAG:
    """A directed acyclic graph of one constants-derivation stack."""
    nodes: dict[str, Node] = field(default_factory=dict)
    edges: list[tuple[str, str]] = field(default_factory=list)

    def add(self, name: str, stage: Stage, value: float | None = None,
            note: str = "", commutative: bool = False,
            dev: float | None = None, kind: str = "") -> str:
        if name not in self.nodes:
            self.nodes[name] = Node(name, stage, value, note, commutative, dev, kind)
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
        """Per OUTPUT->WITNESS edge: relative deviation + admissibility.  cit holds
        iff the output agrees with its witness within ``tol``.  Qualitative outputs
        (no numeric metric) are skipped here and listed by qualitative_outputs()."""
        rows = []
        for s, d in self.edges:
            so, do = self.nodes[s], self.nodes[d]
            if so.stage != Stage.OUTPUT or do.stage != Stage.WITNESS:
                continue
            if so.dev is not None:                       # suite: deviation set directly
                dev, pv, wv = so.dev, so.value, do.value
            elif so.value is not None and do.value is not None:
                pv, wv = so.value, do.value              # headline: compute from values
                dev = abs(pv - wv) / abs(wv)
            else:
                continue                                  # qualitative — not cit-checked
            rows.append({"output": s, "witness": d, "predicted": pv,
                         "measured": wv, "rel_dev": dev, "admissible": dev <= tol})
        return rows

    def qualitative_outputs(self) -> list[str]:
        """OUTPUT nodes with no numeric deviation — O10 marks these as
        qualitative/uncertain rather than fabricating a number."""
        return sorted(n for n, nd in self.nodes.items()
                      if nd.stage == Stage.OUTPUT and nd.dev is None and nd.value is None)

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

    # ---- derivation-route redundancy (M. Wende's "how many independent routes
    #      converge, and what survives if one is removed?") ----

    def _or_node(self, name: str, parents: list[str]) -> bool:
        """Whether ``name`` grounds as an OR of its parents (alternative routes)
        rather than an AND (a conjunction of premises).  Two cases are OR:
          * a ``commutative`` identity — its parents are independent routes to
            one value that must merely agree (21 = C(7,2) = 3·7);
          * an OUTPUT fed *only* by structural integers — the sensitivity DAG's
            influence fan-in, where each integer independently moves the
            benchmark (removing one does not unground it).
        Everything else is an AND: a closed form needs every premise."""
        nd = self.nodes[name]
        if nd.commutative:
            return True
        return (nd.stage == Stage.OUTPUT and bool(parents)
                and all(self.nodes[p].stage == Stage.STRUCTURAL for p in parents))

    def _grounded(self, exclude: frozenset[str] = frozenset()) -> set[str]:
        """Least-fixpoint set of nodes still *derivable from the structural
        axioms* once ``exclude`` is removed.  A STRUCTURAL leaf is an axiom; an
        OR node (see ``_or_node``) needs any parent; every other node needs all
        parents.  Removing one premise of a conjunction ungrounds it; removing
        one route of an OR node does not — which is what makes a "single point of
        failure" well-defined."""
        grounded: set[str] = set()
        changed = True
        while changed:
            changed = False
            for n, nd in self.nodes.items():
                if n in grounded or n in exclude:
                    continue
                ps = [p for p in self._parents(n)]
                if not ps:
                    ok = nd.stage == Stage.STRUCTURAL                  # axiom
                elif self._or_node(n, ps):
                    ok = any(p in grounded for p in ps)                # OR — alt routes
                else:
                    ok = all(p in grounded for p in ps)                # AND — conjunction
                if ok:
                    grounded.add(n)
                    changed = True
        return grounded

    def is_grounded(self, node: str, exclude: frozenset[str] = frozenset()) -> bool:
        return node in self._grounded(exclude)

    def cut_nodes(self, target: str) -> list[str]:
        """Single points of failure for ``target``: internal/structural nodes
        whose removal leaves ``target`` underivable from the axioms.  A node on
        only some of several independent routes is *not* a cut node — the others
        still ground the target."""
        if not self.is_grounded(target):
            return []
        return sorted(v for v, nd in self.nodes.items()
                      if v != target
                      and nd.stage not in (Stage.OUTPUT, Stage.WITNESS)
                      and not self.is_grounded(target, frozenset({v})))

    @staticmethod
    def _max_flow(cap: dict, adj: dict, s, t) -> int:
        """Edmonds-Karp max flow (unit-ish; graphs here are tiny)."""
        from collections import deque
        flow = 0
        while True:
            parent = {s: s}
            q = deque([s])
            while q:
                u = q.popleft()
                if u == t:
                    break
                for w in adj.get(u, ()):
                    if w not in parent and cap.get((u, w), 0) > 0:
                        parent[w] = u
                        q.append(w)
            if t not in parent:
                return flow
            b, v = 1 << 30, t
            while v != s:
                b, v = min(b, cap[(parent[v], v)]), parent[v]
            v = t
            while v != s:
                u = parent[v]
                cap[(u, v)] -= b
                cap[(v, u)] = cap.get((v, u), 0) + b
                v = u
            flow += b

    def independent_routes(self, target: str) -> int:
        """Internally node-disjoint derivation routes reaching ``target`` from
        the structural ground (Menger: = the min node-cut, via unit-node-capacity
        max flow).  Conjunctive premises that funnel through a single closed form
        count as ONE route; genuinely alternative routes — a value reached two
        disjoint ways (a commuting identity) or a benchmark moved by several
        independent integers — count as many.  Routes ≥ 2 is the redundancy."""
        INF = 1 << 30
        SRC = ("<src>", "o")
        cap: dict[tuple, int] = {}
        adj: dict[tuple, set] = {}

        def edge(u, v, c):
            cap[(u, v)] = cap.get((u, v), 0) + c
            adj.setdefault(u, set()).add(v)
            adj.setdefault(v, set()).add(u)
            cap.setdefault((v, u), 0)

        for n, nd in self.nodes.items():
            edge((n, "i"), (n, "o"), INF if n == target else 1)   # unit node capacity
            if nd.stage == Stage.STRUCTURAL:
                edge(SRC, (n, "i"), INF)                          # ground feeds each axiom
        for s, d in self.edges:
            edge((s, "o"), (d, "i"), INF)
        return self._max_flow(cap, adj, SRC, (target, "i"))

    def route_redundancy(self, targets: list[str] | None = None) -> list[dict]:
        """Per target (OUTPUT nodes + commuting identities by default):
        independent route count, single points of failure, and resilience
        (routes ≥ 2)."""
        if targets is None:
            targets = [n for n, nd in self.nodes.items()
                       if nd.stage == Stage.OUTPUT or nd.commutative]
        return [{"target": t, "routes": self.independent_routes(t),
                 "spof": self.cut_nodes(t), "resilient": self.independent_routes(t) >= 2}
                for t in sorted(targets)]

    def criticality_ranking(self) -> list[dict]:
        """Internal/structural nodes ranked by *criticality* — the number of
        OUTPUTs that become ungrounded if the node is removed — beside their
        *load* (# OUTPUTs reached, as in ``load_ranking``).  ``load − critical``
        is the redundancy: a high-load node that is rarely a sole route is well
        backed up; a node whose criticality equals its load is irreplaceable.
        This is the leave-one-route-out dual of ``load_ranking``."""
        outs = [n for n, nd in self.nodes.items()
                if nd.stage == Stage.OUTPUT and self.is_grounded(n)]
        crit: dict[str, int] = {}
        for o in outs:
            for v in self.cut_nodes(o):
                crit[v] = crit.get(v, 0) + 1
        rows = []
        for v, nd in self.nodes.items():
            if nd.stage not in (Stage.STRUCTURAL, Stage.SYMBOLIC, Stage.EVALUATOR):
                continue
            load, c = len(self.reachable_outputs(v)), crit.get(v, 0)
            if load or c:
                rows.append({"node": v, "load": load, "critical": c,
                             "redundancy": load - c})
        return sorted(rows, key=lambda r: (-r["critical"], -r["load"], r["node"]))

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


# ---------------------------------------------------------------------------
# Full benchmark-suite DAG
# ---------------------------------------------------------------------------

def _parse_deviation(accuracy: str) -> tuple[float | None, str]:
    """Parse a BenchmarkResult.substrate_accuracy into (relative deviation
    fraction, provenance kind).  'X%' -> X/100; 'X ppm' -> X/1e6;
    'exact'/'machine'/'100%'/'N/N' -> 0.0 (perfect); otherwise None — O10 marks
    a missing metric as *qualitative* rather than fabricating a number.  A
    leading headline metric wins; ranges and embedded metrics take the upper
    bound (conservative for defect-marking)."""
    s = accuracy.strip()
    low = s.lower()
    if low.startswith("exact") or "machine" in low:
        return 0.0, "exact"
    m = re.match(r"[~<≈]?\s*([\d.]+)\s*(ppm|%)", s)
    if m:
        v = float(m.group(1))
        if m.group(2) == "ppm":
            return v / 1e6, "ppm"
        return (0.0, "score") if v >= 100 else (v / 100.0, "pct")
    m = re.search(r"\b(\d+)\s*/\s*(\d+)\b", s)            # N/M score
    if (m and m.group(1) == m.group(2)) or re.search(r"\b100\s*%", s):
        return 0.0, "score"
    pcts = [float(x) for x in re.findall(r"([\d.]+)\s*%", s) if float(x) < 100]
    ppms = [float(x) for x in re.findall(r"([\d.]+)\s*ppm", s)]
    if pcts:
        return max(pcts) / 100.0, "pct"
    if ppms:
        return max(ppms) / 1e6, "ppm"
    return None, "qualitative"


def benchmark_functions() -> dict:
    """The 38 suite benchmark functions keyed by name (the sensitivity sweep's keys)."""
    from . import compute_speed as cs
    return {n: f for n, f in vars(cs).items()
            if n.startswith("benchmark_") and callable(f)}


def build_suite_dag(report=None) -> DerivationDAG:
    """The whole 38-benchmark suite as one O10 DAG.

    Each benchmark is an OUTPUT node with a WITNESS edge carrying its deviation
    (parsed from ``substrate_accuracy``; qualitative benchmarks left unscored and
    surfaced by ``qualitative_outputs()``).  The STRUCTURAL->OUTPUT edges are the
    *computed* coupling: pass a ``SensitivityReport`` and an edge
    ``integer -> benchmark`` is drawn for every integer that moves that benchmark,
    so the DAG's ``load_ranking()`` equals the sweep's structural load.  Without a
    report, a single 'ISA' structural hub feeds every benchmark."""
    g = DerivationDAG()
    g.add("ISA", Stage.STRUCTURAL, note="K_7/Spin(7) structural integers (isa)")
    for name, fn in sorted(benchmark_functions().items()):
        res = fn()
        dev, kind = _parse_deviation(res.substrate_accuracy)
        g.add(name, Stage.OUTPUT, note=res.name, dev=dev, kind=kind)
        g.add(f"wit:{name}", Stage.WITNESS, note="measured (CODATA-2018 / PDG / Planck)")
        g.link(name, f"wit:{name}")
        if report is None:
            g.link("ISA", name)
    if report is not None:
        for integer in report.per_integer:
            g.add(integer, Stage.STRUCTURAL, note="isa structural knob (integer or derived scalar)")
            for bench in report.movers(integer):
                if bench in g.nodes:
                    g.link(integer, bench)
    return g


def _redundancy_lines(g: DerivationDAG, list_all: bool = True,
                      label: str = "targets") -> list[str]:
    """The M. Wende derivation-route-redundancy readout: per-target independent
    route count + single points of failure, then node criticality vs load.  For
    the suite (``list_all=False``) only the fragile single-route targets are
    enumerated; the constants DAG lists every target."""
    rr = g.route_redundancy()
    resilient = [r for r in rr if r["routes"] >= 2]
    fragile = [r for r in rr if r["routes"] == 1]      # exactly one structural route
    unswept = [r for r in rr if r["routes"] == 0]      # no structural route in this DAG
    lines = ["", "Derivation-route redundancy (M. Wende — independent routes + "
             "single points of failure):",
             f"  {len(rr)} {label}: {len(resilient)} resilient (≥2 independent routes), "
             f"{len(fragile)} single-route (one SPOF), "
             f"{len(unswept)} with no detected structural route (0 routes — exact-combinatorial "
             f"output, measured-input-anchored, or a bare-literal isa-anchoring item; not fragility)"]
    for r in sorted(resilient if list_all else [], key=lambda r: (-r["routes"], r["target"])):
        spof = ", ".join(r["spof"]) if r["spof"] else "—"
        lines.append(f"  [resilient]  {r['target']:26s} {r['routes']} route(s)   SPOF: {spof}")
    for r in sorted(fragile, key=lambda r: r["target"]):
        lines.append(f"  [SPOF     ]  {r['target']:26s} 1 route    SPOF(s): "
                     f"{', '.join(r['spof']) or '—'}")
    if unswept and not list_all:
        names = ", ".join(r["target"].replace("benchmark_", "") for r in unswept)
        lines.append(f"  [unswept  ]  {len(unswept)}: {names}")
    lines += ["", "Node criticality vs load (ungrounded-if-removed | reached | redundancy "
              "= load−critical):"]
    for r in g.criticality_ranking()[:12]:
        lines.append(f"  {r['node']:26s} critical {r['critical']:>2}   load {r['load']:>2}"
                     f"   redundancy {r['redundancy']:>2}")
    return lines


def _main_suite(use_sensitivity: bool = False, show_redundancy: bool = False) -> int:
    report = None
    if use_sensitivity:
        from ..sensitivity import integer_sweep
        report = integer_sweep()
    g = build_suite_dag(report=report)
    cit, defects, qual = g.cit_readout(0.01), g.cit_defects(0.01), g.qualitative_outputs()
    n_bench = len(benchmark_functions())
    lines = [f"O10 suite DAG — {len(g.nodes)} nodes, {len(g.edges)} edges, {n_bench} benchmarks", ""]
    lines.append("Acceptance checklist (O10 structural invariants):")
    for k, v in g.acceptance_checklist().items():
        lines.append(f"  [{'PASS' if v else 'FAIL'}]  {k}")
    lines += ["", f"cit readout: {len(cit)} scored, {len(cit) - len(defects)} admissible, "
              f"{len(defects)} defect edge(s) marked, {len(qual)} qualitative:"]
    for r in sorted(cit, key=lambda r: -r["rel_dev"]):
        if not r["admissible"]:
            lines.append(f"  [DEFECT]  {r['output']:36s} {r['rel_dev'] * 100:6.2f}%")
    if qual:
        lines.append(f"  qualitative (no vs-measurement metric): {', '.join(qual)}")
    if report is not None:
        lines += ["", "Structural-load ranking (isa integers by # benchmarks moved):"]
        for n, load in g.load_ranking()[:10]:
            lines.append(f"  {n:30s} reaches {load} benchmark(s)")
    else:
        lines += ["", "(add --sensitivity for the real structural->benchmark coupling "
                  "edges + load ranking; slow — the sweep patches isa per integer.)"]
    if show_redundancy:
        if report is None:
            lines += ["", "(add --sensitivity for route redundancy — with the single ISA "
                      "hub every benchmark is trivially one route.)"]
        else:
            lines += _redundancy_lines(g, list_all=False, label="benchmarks")
    print("\n".join(lines))
    return 0


def main(argv: list[str] | None = None) -> int:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if "--suite" in argv:
        return _main_suite(use_sensitivity="--sensitivity" in argv,
                           show_redundancy="--redundancy" in argv)
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

    if "--redundancy" in (argv or []):
        lines += _redundancy_lines(g)

    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
