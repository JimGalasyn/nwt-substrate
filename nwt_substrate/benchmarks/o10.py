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
--sensitivity --redundancy`` runs it over the 38-benchmark coupling graph, and
appends M. Wende's derivation-*diversity* layer: route count is multiplicity,
not independence, so it clusters the integers into structural sectors by shared
perturbation response and re-counts *effective* routes per sector — surfacing
benchmarks whose graph-independent routes collapse into one (hidden SPOFs).
Add ``--priority`` for M. Wende's *prioritization* layer: the directional dual of
the diagnostics, ranking the open items by closure gained per unit effort (and
triaging the sweep-unreachable benchmarks whose effort is unknown until read).

A **value-provenance** layer runs in both readouts (it is cheap, so always on).
The structural invariants and ``cit`` check *where* a value sits and *whether* it
matches; they cannot see *how* the value was obtained — and a value that was
fitted, normalised, post-selected, or convention-pinned matches its witness by
construction, so a clean cit edge corroborates the choice, not the theory.  Each
``Node`` carries a ``provenance`` (DERIVED / MEASURED / DEFINITION / NORMALIZED /
FITTED / POST_SELECTED / CONVENTION / ASSERTED; default inferred), and three lints
mark — never repair — what the rest of the DAG misses: ``provenance_defects`` (a
non-derived OUTPUT passing cit = circular), ``tautology_nodes`` (a value defined
to equal its own premise), and ``asserted_operators`` (an operator depended on but
never constructed).  These encode the failure modes a Maxwell-from-D12RG review
and our occasion-inflation κ=6 / H0 audit both turned on.
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


# Value provenance — HOW a node's value was obtained, orthogonal to its proof
# Stage.  Stage says *where in the ladder* a node sits; provenance says *how the
# number got there*, which is what the structural invariants and cit cannot see:
# a value that was fitted / normalised / post-selected / convention-pinned AGREES
# with its witness by construction, so a clean cit edge corroborates nothing.
# This is the failure mode an external Maxwell-from-D12RG derivation review and
# our own occasion-inflation κ=6 / H0 audit both turned on: cit asks "does it
# match?"; the provenance lints ask "would it match no matter what?".  Defects
# here are MARKED, never repaired (the same O10 discipline as cit defects).
DERIVED = "derived"            # forced from upstream structure (the good case)
MEASURED = "measured"          # empirical witness (CODATA / PDG)
DEFINITION = "definition"      # assigned by fiat; value IS a parent -> tautology risk
NORMALIZED = "normalized"      # set to a fixed value (=1, ...) by a normalisation choice
FITTED = "fitted"             # tuned to match the target
POST_SELECTED = "post_selected"  # the closest of a search (look-elsewhere)
CONVENTION = "convention"      # pinned by an arbitrary convention (e.g. radius λ̄ vs λ̄/2)
MOTIVATED = "motivated"        # structural rationale, not forced (the L4(a) audit
#                                vocabulary: "motivated rather than derived" — a counted
#                                exponent, a matched amplitude, an identified prefactor)
ASSERTED = "asserted"          # depended on as an operator but never constructed

# Values whose agreement with a witness is circular, not evidence.  MOTIVATED is
# suspect: a motivated-not-forced value could have been motivated differently had
# the target differed, so its match corroborates the identification, not the theory.
SUSPECT_PROVENANCE = frozenset(
    {DEFINITION, NORMALIZED, FITTED, POST_SELECTED, CONVENTION, MOTIVATED})
_VALID_PROVENANCE = SUSPECT_PROVENANCE | {DERIVED, MEASURED, ASSERTED}


# Result-level CLAIM status — what the HEADLINE asserts, orthogonal to a node's
# value-provenance (provenance says HOW the number got there; status says what
# epistemic standing the result claims).  Adapted from the per-section audit
# blocks of an external D12RG "trace-determinant carrier" paper, whose discipline
# (carrier first, downstream readouts second, NO REVERSE SMUGGLING) is exactly
# what catches a downstream coincidence retro-justifying an upstream premise — a
# smuggle that is usually NOT a cycle (it is a prose inference never added as an
# edge, so is_acyclic() stays green).  Only STATUS_DEFERRED_BRIDGE carries an
# obligation; the rest are descriptive.
STATUS_DEFINITION = "status_definition"        # the headline is a definition
STATUS_THEOREM = "status_theorem"              # forced from upstream (claims proof)
STATUS_MEASURED_MATCH = "status_measured_match"  # agrees with an empirical witness
STATUS_DEFERRED_BRIDGE = "status_deferred_bridge"  # an IOU to future work — MUST
#                              name the killable test that would discharge it.
_VALID_STATUS = frozenset(
    {STATUS_DEFINITION, STATUS_THEOREM, STATUS_MEASURED_MATCH, STATUS_DEFERRED_BRIDGE})


@dataclass(frozen=True)
class ForbiddenCollapse:
    """An *anti-edge*: a declaration that ``dst`` must NOT come to depend on
    ``src`` — the structural form of "X readout ⇏ the theorem that produced X".
    The DAG audits only edges that exist in it, so a reverse-smuggle (using a
    downstream match to prop up an upstream premise) is invisible to is_acyclic()
    until someone wires it up; this records the forbidden wiring so the lint can
    fire the moment it appears.

    ``discharge`` is the killable test that *would* license the collapse.  An
    empty discharge is itself a defect: forbidding a collapse you cannot say how
    to license is an un-cashable IOU (the unfalsifiability the source paper fell
    into — deferral as a permanent shield rather than a tracked debt)."""
    src: str
    dst: str
    reason: str = ""
    discharge: str = ""


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
    provenance: str = ""          # how the VALUE was obtained (module constants above);
    #                               "" -> inferred (MEASURED for WITNESS, else DERIVED).
    status: str = ""              # result-level CLAIM standing (STATUS_* above); "" -> none.
    discharge: str = ""           # for a STATUS_DEFERRED_BRIDGE: the killable test that
    #                               would discharge the IOU; "" -> an un-cashable obligation.
    sigma: float | None = None    # WITNESS 1σ experimental uncertainty (absolute), enabling
    #                               the S-NOW readout; None -> witness has no error budget
    #                               and cannot be S-NOW-scored (cit-only).
    disputed: str = ""            # provenance DISPUTE record: a pinned EXTERNAL audit
    #                               contests this node's self-declared provenance
    #                               ("contested=<tag> per <repo/path@sha>").  A dispute is
    #                               marked, never repaired: cit agreement on any OUTPUT
    #                               downstream of a disputed node is suspended as
    #                               corroboration until the dispute is adjudicated
    #                               (the Auditor's verdict replaces the tag AND clears
    #                               this field, citing the verdict).  Self-tags and
    #                               external audits disagreeing silently is the failure
    #                               mode this field exists to surface.


@dataclass
class DerivationDAG:
    """A directed acyclic graph of one constants-derivation stack."""
    nodes: dict[str, Node] = field(default_factory=dict)
    edges: list[tuple[str, str]] = field(default_factory=list)
    antiedges: list[ForbiddenCollapse] = field(default_factory=list)

    def add(self, name: str, stage: Stage, value: float | None = None,
            note: str = "", commutative: bool = False,
            dev: float | None = None, kind: str = "", provenance: str = "",
            status: str = "", discharge: str = "",
            sigma: float | None = None, disputed: str = "") -> str:
        if provenance and provenance not in _VALID_PROVENANCE:
            raise ValueError(f"unknown provenance {provenance!r} for node {name!r}")
        if status and status not in _VALID_STATUS:
            raise ValueError(f"unknown status {status!r} for node {name!r}")
        if name not in self.nodes:
            self.nodes[name] = Node(name, stage, value, note, commutative,
                                    dev, kind, provenance, status, discharge,
                                    sigma, disputed)
        return name

    def link(self, src: str, dst: str) -> None:
        self.edges.append((src, dst))

    def forbid(self, src: str, dst: str, reason: str = "", discharge: str = "") -> None:
        """Declare an anti-edge: ``dst`` must not come to depend on ``src``.  Like
        ``link``, it does not require the endpoints to exist yet (they may be added
        later).  ``discharge`` names the killable test that would license the
        collapse; leaving it empty is itself flagged by ``collapse_defects``."""
        self.antiedges.append(ForbiddenCollapse(src, dst, reason, discharge))

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

    # ---- value-provenance lints (the "would it match no matter what?" axis cit
    #      cannot see — fitted/normalised/post-selected values, tautologies, and
    #      asserted-but-unconstructed operators) ----

    def provenance(self, name: str) -> str:
        """A node's value provenance, inferring the default from its stage when
        unset: a WITNESS is MEASURED, everything else is DERIVED until declared
        otherwise.  Backwards-compatible: existing DAGs read as fully derived."""
        nd = self.nodes[name]
        if nd.provenance:
            return nd.provenance
        return MEASURED if nd.stage == Stage.WITNESS else DERIVED

    def provenance_defects(self, tol: float = 0.01) -> list[dict]:
        """OUTPUTs whose value was NOT genuinely derived (fitted / normalised /
        post-selected / convention-pinned / defined) yet pass cit anyway.  Such a
        value agrees with its witness *by construction*, so a clean cit edge is
        circular — corroboration of the choice, not the theory.  This is the
        anti-numerology counterpart to ``cit_readout``: cit asks "does it match?",
        this asks "would it match regardless?".  Marked, never repaired."""
        cit = {r["output"]: r for r in self.cit_readout(tol)}
        rows = []
        for n, nd in self.nodes.items():
            if nd.stage != Stage.OUTPUT:
                continue
            prov = self.provenance(n)
            if prov not in SUSPECT_PROVENANCE:
                continue
            r = cit.get(n)
            passes = bool(r and r["admissible"])
            rows.append({"output": n, "provenance": prov, "passes_cit": passes,
                         "rel_dev": (r["rel_dev"] if r else None),
                         # a SUSPECT value that *passes* cit is the circular case;
                         # one that *fails* is at least honestly marked by cit too.
                         "circular": passes, "note": nd.note})
        return sorted(rows, key=lambda r: (not r["circular"], r["output"]))

    def asserted_operators(self) -> list[str]:
        """Nodes that some derivation *depends on* but that are never constructed:
        explicit ASSERTED provenance, or a non-axiom, non-witness node that is used
        as a parent yet has no value and no incoming construction edge.  This is
        the undefined-Hodge-star failure mode — a step that leans on an operator
        doing all the work without ever building it."""
        used = {s for s, _ in self.edges}
        out = set()
        for n, nd in self.nodes.items():
            if n not in used:                 # nothing depends on it -> not a smuggle
                continue
            if self.provenance(n) == ASSERTED:
                out.add(n)
            elif (nd.value is None and not self._parents(n)
                  and nd.stage not in (Stage.STRUCTURAL, Stage.WITNESS)):
                out.add(n)
        return sorted(out)

    def tautology_nodes(self, rel: float = 1e-9) -> list[dict]:
        """SYMBOLIC/OUTPUT nodes whose value equals a parent's to machine precision
        *and* whose provenance is by-fiat (DEFINITION / NORMALIZED / CONVENTION) —
        a "theorem" that is its own premise relabelled (the J:=δF → "derive" d*F=J
        pattern; ChargeCell² defined-then-"found"=1).  Value identity alone is fine
        (a pass-through evaluator is honest); identity *plus* a by-fiat value is the
        tautology."""
        suspect = {DEFINITION, NORMALIZED, CONVENTION}
        out = []
        for n, nd in self.nodes.items():
            if nd.stage not in (Stage.SYMBOLIC, Stage.OUTPUT) or nd.value is None:
                continue
            if self.provenance(n) not in suspect:
                continue
            for p in self._parents(n):
                pv = self.nodes[p].value
                if pv is None:
                    continue
                # absolute floor of 1 so a zero-valued definition/convention is
                # still caught (the lint is value identity, not relative error).
                if abs(nd.value - pv) <= rel * max(abs(pv), 1.0):
                    out.append({"node": n, "equals_parent": p, "value": nd.value,
                                "provenance": self.provenance(n)})
                    break
        return out

    def provenance_audit(self, tol: float = 0.01) -> dict:
        """One-call summary of the three value-provenance lints — the gauntlet axis
        the structural invariants and cit do not cover."""
        suspect = self.provenance_defects(tol)
        circular = [r for r in suspect if r["circular"]]
        asserted = self.asserted_operators()
        tautologies = self.tautology_nodes()
        return {
            "circular_passes": circular,
            "suspect_outputs": suspect,
            "asserted_operators": asserted,
            "tautologies": tautologies,
            "clean": not asserted and not tautologies and not circular,
        }

    # ---- provenance-dispute lint (the T1 axis: the DAG's SELF-declared tags vs
    #      the program's own PINNED external audits.  A provenance lint that the
    #      claim's author satisfies by tagging everything DERIVED is theater; this
    #      lint surfaces every node whose tag a pinned audit contests, and
    #      suspends cit corroboration downstream until adjudication.  Marked,
    #      never repaired — adjudication (the memory-blind Auditor's verdict)
    #      replaces the tag and clears the dispute, citing the verdict.) ----

    def disputed_nodes(self) -> list[str]:
        """Nodes whose self-declared provenance is contested by a pinned external
        audit (``disputed`` non-empty)."""
        return sorted(n for n, nd in self.nodes.items() if nd.disputed)

    def dispute_audit(self, tol: float = 0.01) -> dict:
        """Per disputed node: the dispute record and every OUTPUT it reaches —
        those outputs' cit passes are SUSPENDED as corroboration pending
        adjudication (agreement with a witness cannot corroborate a chain whose
        provenance is itself in dispute)."""
        cit = {r["output"]: r for r in self.cit_readout(tol)}
        rows = []
        suspended: set[str] = set()
        for n in self.disputed_nodes():
            outs = sorted(self.reachable_outputs(n) | (
                {n} if self.nodes[n].stage == Stage.OUTPUT else set()))
            suspended.update(o for o in outs if cit.get(o, {}).get("admissible"))
            rows.append({"node": n, "dispute": self.nodes[n].disputed,
                         "self_tag": self.provenance(n), "reaches": outs})
        return {"disputes": rows,
                "suspended_outputs": sorted(suspended),
                "clean": not rows}

    # ---- S-NOW readout (precision confrontation): cit with a fixed tolerance is
    #      a smoke test — at tol=1% every sub-percent postdiction "passes".  The
    #      kill surface asks the sharper question: is the frozen form compatible
    #      with the measured value AT THE EXPERIMENT'S OWN PRECISION?  A row can
    #      pass cit and still be DEAD-AS-EXACT by thousands of σ (α itself is).
    #      Requires witness nodes to carry ``sigma``; rows without one are
    #      reported UNSCORED rather than silently skipped. ----

    def snow_readout(self, n_sigma: float = 2.0) -> list[dict]:
        """Per OUTPUT->WITNESS edge: z = |predicted − measured| / σ_measured and
        the verdict — EXACT-COMPATIBLE (z ≤ n_sigma), DEAD-AS-EXACT (z > n_sigma:
        the frozen form is excluded at current experimental precision and
        survives only as the leading order of an undeclared series), or UNSCORED
        (witness has no σ).  NOTE: a pass is *compatibility*, not confirmation —
        postdictive values are compatible with their own targets by construction;
        only the S-FORWARD channel can upgrade a row to evidence."""
        rows = []
        for s, d in self.edges:
            so, do = self.nodes[s], self.nodes[d]
            if so.stage != Stage.OUTPUT or do.stage != Stage.WITNESS:
                continue
            if so.value is None or do.value is None:
                continue                                   # qualitative — cit's domain
            if do.sigma is None or do.sigma <= 0:
                rows.append({"output": s, "witness": d, "predicted": so.value,
                             "measured": do.value, "sigma": None, "z": None,
                             "verdict": "UNSCORED"})
                continue
            z = abs(so.value - do.value) / do.sigma
            rows.append({"output": s, "witness": d, "predicted": so.value,
                         "measured": do.value, "sigma": do.sigma, "z": z,
                         "verdict": "EXACT-COMPATIBLE" if z <= n_sigma
                                    else "DEAD-AS-EXACT"})
        return rows

    # ---- forbidden-collapse lints (the *negative-edge* axis: a downstream readout
    #      retro-justifying an upstream premise — "X readout ⇏ the theorem that
    #      produced X".  The provenance lints ask "would this value match no matter
    #      what?"; these ask "is a conclusion being smuggled back into its own
    #      premise?", the reverse-implication the directed DAG cannot see until the
    #      smuggle is wired up).  Marked, never repaired. ----

    def _reaches(self, src: str, dst: str) -> bool:
        """Whether a directed path src -> ... -> dst exists (src itself excluded)."""
        if src not in self.nodes or dst not in self.nodes:
            return False
        stack, seen = list(self._children(src)), set()
        while stack:
            n = stack.pop()
            if n == dst:
                return True
            if n in seen:
                continue
            seen.add(n)
            stack.extend(self._children(n))
        return False

    def deferred_bridges(self) -> list[str]:
        """Nodes whose headline status is a DEFERRED_BRIDGE (an IOU to future work)."""
        return sorted(n for n, nd in self.nodes.items()
                      if nd.status == STATUS_DEFERRED_BRIDGE)

    def collapse_defects(self, tol: float = 0.01) -> list[dict]:
        """Audit the declared anti-edges and deferred bridges.  Four defect kinds:
          * ``violated``     — a directed path src->dst now EXISTS: the forbidden
                               wiring was added (a conclusion feeds its own premise).
          * ``coincidence``  — a MEASURED node whose value equals src's (within
                               ``tol``) reaches dst: the premise is propped up by a
                               numeric match laundered through a different node.
          * ``undischarged`` — a ``forbid`` with no ``discharge``: an un-cashable
                               IOU (forbidding a collapse you can't say how to
                               license — the source paper's unfalsifiability trap).
          * ``bridge_iou``   — a STATUS_DEFERRED_BRIDGE node with no ``discharge``.
        ``violated``/``coincidence`` are hard (a real smuggle); the two IOU kinds
        are marked obligations, surfaced but not failing the structural checklist."""
        rows: list[dict] = []
        for fc in self.antiedges:
            if self._reaches(fc.src, fc.dst):
                rows.append({"kind": "violated", "src": fc.src, "dst": fc.dst,
                             "reason": fc.reason, "discharge": fc.discharge})
            else:
                sv = self.nodes[fc.src].value if fc.src in self.nodes else None
                if sv is not None:
                    for m, nd in self.nodes.items():
                        if (m != fc.src and self.provenance(m) == MEASURED
                                and nd.value is not None
                                and abs(nd.value - sv) <= tol * max(abs(sv), 1.0)
                                and (m == fc.dst or self._reaches(m, fc.dst))):
                            rows.append({"kind": "coincidence", "src": fc.src,
                                         "dst": fc.dst, "via": m, "reason": fc.reason,
                                         "discharge": fc.discharge})
                            break
            if not fc.discharge.strip():
                rows.append({"kind": "undischarged", "src": fc.src, "dst": fc.dst,
                             "reason": fc.reason})
        for n in self.deferred_bridges():
            if not self.nodes[n].discharge.strip():
                rows.append({"kind": "bridge_iou", "node": n,
                             "note": self.nodes[n].note})
        return rows

    def collapse_audit(self, tol: float = 0.01) -> dict:
        """One-call summary of the forbidden-collapse layer.  ``violations`` are the
        hard reverse-smuggles (a forbidden path or a laundered numeric match);
        ``open_obligations`` are the un-cashable IOUs (anti-edges or bridges with no
        discharge).  ``clean`` requires both empty."""
        defects = self.collapse_defects(tol)
        violations = [d for d in defects if d["kind"] in ("violated", "coincidence")]
        obligations = [d for d in defects if d["kind"] in ("undischarged", "bridge_iou")]
        return {
            "violations": violations,
            "open_obligations": obligations,
            "deferred_bridges": self.deferred_bridges(),
            "clean": not violations and not obligations,
        }

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
            # a declared anti-edge that is now actually wired up (or propped by a
            # laundered numeric match) is a structural defect, not just a warning:
            # the DAG derives something it forbade.  Un-cashable IOUs are reported
            # by collapse_audit but kept out of the hard checklist (marked, not failed).
            "no_forbidden_collapse_violated": not self.collapse_audit()["violations"],
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

    # Forbidden collapse (dogfood of the anti-edge layer): η_B = 3α⁴/14 lands on the
    # Planck baryon asymmetry, but int:3 and int:14 are bare structural literals here.
    # The Planck match must NOT be allowed to retro-justify the choice of 3 and 14 —
    # that would be the conclusion (the witness) feeding its own premise.  The path
    # wit:eta_B -> int:14 does not exist (witnesses are sinks), so this is not yet
    # violated; the discharge names the killable test that *would* license it, so the
    # obligation is honestly OPEN rather than an un-cashable IOU.
    g.forbid("wit:eta_B", "int:14",
             reason="Planck η_B match ⇏ the integers 3,14 that hit it",
             discharge="derive 3 & 14 from Jones/Murasugi knot-chirality independently "
                       "of the Planck value — e.g. a second baryon-sector observable, or "
                       "a neighbouring-α test the integers cannot be retuned for")
    g.forbid("wit:eta_B", "int:3",
             reason="Planck η_B match ⇏ the integers 3,14 that hit it",
             discharge="as int:14 — forced from chirality combinatorics, not back-fit")
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


def _diversity_lines(report, threshold: float = 0.5) -> list[str]:
    """M. Wende's derivation-*diversity* readout (d12rg v0.4.1 follow-up):
    route count is multiplicity; this is independence.  Cluster the load-bearing
    integers into structural sectors by shared perturbation response, then flag
    benchmarks whose graph-independent routes collapse into one sector — routes
    that look redundant but fail to the same perturbation.  Needs the sweep."""
    div = report.route_diversity(threshold)
    sectors = report.integer_sectors(threshold)
    n_int = sum(len(s) for s in sectors)
    collapsed = [r for r in div if r["collapsed"] > 0]
    hidden_spof = [r for r in collapsed if r["effective"] == 1]
    lines = ["", "Derivation diversity (M. Wende — route independence beyond route "
             f"count; sectors = integers sharing a perturbation response, Jaccard ≥ {threshold}):",
             f"  {len(sectors)} structural sector(s) among {n_int} load-bearing integers; "
             f"{len(collapsed)} benchmark(s) lose routes to shared response "
             f"({len(hidden_spof)} collapse to a hidden single sector)"]
    multi = [s for s in sectors if len(s) > 1]
    if multi:
        lines.append("  sectors that merge multiple integers (one effective route each):")
        for s in multi:
            lines.append(f"    {{{', '.join(s)}}}")
    for r in sorted(collapsed, key=lambda r: (r["effective"], -r["raw"], r["benchmark"])):
        tag = "hidden SPOF" if r["effective"] == 1 else "diversity<routes"
        names = " | ".join("+".join(g) for g in r["sectors"])
        lines.append(f"  [{tag:16s}]  {r['benchmark'].replace('benchmark_', ''):24s} "
                     f"raw {r['raw']} -> {r['effective']} effective   sectors: {names}")
    if not collapsed:
        lines.append("  (every multi-route benchmark draws on genuinely distinct sectors — "
                     "route count == route independence)")
    return lines


def _provenance_lines(g: DerivationDAG, tol: float = 0.01) -> list[str]:
    """The value-provenance readout: the "would it match no matter what?" axis.
    Marks OUTPUTs that pass cit on a non-derived value (circular), operators that
    are asserted but never built, and definitional tautologies — none of which the
    structural invariants or cit can see.  Clean DAGs print one reassuring line."""
    audit = g.provenance_audit(tol)
    lines = ["", "Value-provenance lints (HOW each value was obtained — the axis "
             "cit cannot see; marked, not repaired):"]
    if audit["clean"]:
        lines.append("  [clean]  every output is DERIVED/MEASURED; no asserted operators, "
                     "no definitional tautologies, no fitted value passing cit.")
        return lines
    for r in audit["circular_passes"]:
        lines.append(f"  [CIRCULAR ]  {r['output']:26s} provenance={r['provenance']} "
                     f"passes cit — agreement is not evidence (value not derived)")
    for r in audit["tautologies"]:
        lines.append(f"  [TAUTOLOGY]  {r['node']:26s} = parent {r['equals_parent']} "
                     f"by {r['provenance']} — a theorem that is its own premise")
    for n in audit["asserted_operators"]:
        lines.append(f"  [ASSERTED ]  {n:26s} depended on but never constructed "
                     "(undefined-operator smuggle)")
    return lines


def _dispute_lines(g: DerivationDAG, tol: float = 0.01) -> list[str]:
    """The provenance-dispute readout: self-declared tags contested by pinned
    external audits.  Every affected output's cit pass is listed as SUSPENDED —
    agreement cannot corroborate a chain whose provenance is in dispute.  A DAG
    with no disputes prints nothing (the section only appears when it bites)."""
    audit = g.dispute_audit(tol)
    if audit["clean"]:
        return []
    lines = ["", "Provenance disputes (self-declared tag vs pinned external audit; "
             "marked, not repaired — adjudication replaces the tag):"]
    for r in audit["disputes"]:
        lines.append(f"  [DISPUTED ]  {r['node']:26s} self={r['self_tag']}   {r['dispute']}")
    if audit["suspended_outputs"]:
        lines.append("  cit passes SUSPENDED as corroboration pending adjudication: "
                     + ", ".join(audit["suspended_outputs"]))
    return lines


def _snow_lines(g: DerivationDAG, n_sigma: float = 2.0) -> list[str]:
    """The S-NOW readout: each scored output confronted at the experiment's own
    precision.  Pass = compatibility (postdictions are compatible with their own
    targets by construction), never confirmation."""
    rows = g.snow_readout(n_sigma)
    if not rows:
        return []
    lines = ["", f"S-NOW readout (|predicted − measured| / σ_experiment; "
             f"verdict at {n_sigma:g}σ — compatibility, NOT confirmation):"]
    for r in sorted(rows, key=lambda r: -(r["z"] if r["z"] is not None else -1)):
        if r["z"] is None:
            lines.append(f"  [unscored      ]  {r['output']:24s} witness has no σ")
            continue
        tag = "ok  " if r["verdict"] == "EXACT-COMPATIBLE" else "DEAD"
        z = f"{r['z']:.3g}σ" if r["z"] < 1e6 else f"{r['z']:.2e}σ"
        lines.append(f"  [{tag} {z:>10s}]  {r['output']:24s} "
                     f"{r['predicted']:.9g}  vs  {r['measured']:.9g} ± {r['sigma']:.2g}")
    return lines


def _collapse_lines(g: DerivationDAG, tol: float = 0.01) -> list[str]:
    """The forbidden-collapse readout: declared anti-edges (no-reverse-smuggling
    rules) and deferred-bridge IOUs.  A VIOLATED line is a conclusion now feeding
    its own premise; an OPEN line is a tracked debt that names its killable test;
    an IOU line is an obligation with no discharge — the unfalsifiability trap."""
    audit = g.collapse_audit(tol)
    if not g.antiedges and not audit["deferred_bridges"]:
        return []
    lines = ["", "Forbidden-collapse lints (no-reverse-smuggling — a downstream "
             "readout propping up an upstream premise; marked, not repaired):"]
    if audit["clean"] and not g.antiedges:
        lines.append("  [clean]  no anti-edges declared.")
        return lines
    for d in audit["violations"]:
        if d["kind"] == "violated":
            lines.append(f"  [VIOLATED ]  {d['src']} ⇒ {d['dst']} now has a directed path "
                         f"— {d['reason'] or 'forbidden collapse wired up'}")
        else:
            lines.append(f"  [LAUNDERED]  {d['src']} ⇒ {d['dst']} via {d['via']} (matching "
                         f"value) — {d['reason'] or 'premise propped by a numeric match'}")
    for d in audit["open_obligations"]:
        if d["kind"] == "undischarged":
            lines.append(f"  [IOU      ]  {d['src']} ⇏ {d['dst']} has NO discharge test "
                         "— un-cashable (name the killable test or drop the anti-edge)")
        else:
            lines.append(f"  [IOU      ]  deferred bridge {d['node']!r} has no discharge "
                         "— name the test that would close it")
    # honestly-open anti-edges (named discharge, not yet met) — the tracked-debt case
    viol = {(d["src"], d["dst"]) for d in audit["violations"]}
    obl = {(d.get("src"), d.get("dst")) for d in audit["open_obligations"]}
    for fc in g.antiedges:
        if (fc.src, fc.dst) in viol or (fc.src, fc.dst) in obl:
            continue
        lines.append(f"  [open     ]  {fc.src} ⇏ {fc.dst}   discharge: {fc.discharge}")
    return lines


def _redundancy_lines(g: DerivationDAG, list_all: bool = True,
                      label: str = "targets", report=None) -> list[str]:
    """The M. Wende derivation-route-redundancy readout: per-target independent
    route count + single points of failure, then node criticality vs load.  For
    the suite (``list_all=False``) only the fragile single-route targets are
    enumerated; the constants DAG lists every target.  When a sweep ``report`` is
    supplied, the derivation-*diversity* layer (route independence) is appended."""
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
    if report is not None:
        lines += _diversity_lines(report)
    return lines


# Benchmarks confirmed by code inspection to derive a quantity from a *measured*
# literal (a real tier-1 fix, not a guess) — the muon-decay rate still reaches for
# a local CODATA-alpha constant, the residue of the v0.4.2 QED alpha-fix.
_VERIFIED_LEAKS = frozenset({"benchmark_muon_decay_rate"})


def closure_priority(g: DerivationDAG, report) -> dict:
    """M. Wende's prioritization layer (d12rg) — rank the open items by closure
    gained per unit effort, the directional dual of the diagnostic layers.

    Each benchmark scores a closure weight: 0.0 ungrounded (no structural route),
    0.5 single effective sector (a hidden SPOF), 1.0 with >= 2 independent sectors
    — minus 0.5 if it is a marked defect edge (derived, but wrong).  Open items
    split into two queues:

      * ``items``  — effort is *known* (a hidden SPOF, a defect edge, a verified
        measured-input leak), ranked by ``roi = gain / effort``;
      * ``triage`` — benchmarks the sweep's leaf-integer perturbations cannot
        reach at all, so the effort is *unknown* until a code read tells a cheap
        measured-input leak (the alpha-fix pattern) from a real derivation gap.
        Ranked by *potential* closure; resolving each is itself the leverage.

    Effort tiers are a declared heuristic keyed to category (1 mechanical /
    2 known-layer / 4 research), the one place judgement enters.  Needs a sweep
    ``report`` (the diversity signatures come from it)."""
    defects = set(g.cit_defects(0.01))
    qual = set(g.qualitative_outputs())
    div = {r["benchmark"]: r for r in report.route_diversity()}
    items, triage, total = [], [], 0.0
    outs = [n for n, nd in g.nodes.items() if nd.stage == Stage.OUTPUT]
    for b in outs:
        raw = sum(1 for i in report.per_integer if b in report.movers(i))
        eff = div[b]["effective"] if b in div else 0
        is_defect = b in defects
        base = 0.0 if raw == 0 else (0.5 if eff <= 1 else 1.0)
        total += max(0.0, base - (0.5 if is_defect else 0.0))
        name = b.replace("benchmark_", "")
        if raw == 0:                                   # sweep cannot reach it
            triage.append({"benchmark": name, "defect": is_defect,
                           "potential": round(0.5 + (0.5 if is_defect else 0.0), 2)})
            continue
        if eff <= 1:                                   # single-sector hidden SPOF
            if b in _VERIFIED_LEAKS:
                gain, effort, cat = 0.5, 1, "hidden SPOF — verified measured-input leak (route via isa)"
            elif is_defect:
                gain, effort, cat = 1.0, 2, "hidden SPOF + defect (route + correction layer)"
            else:
                gain, effort, cat = 0.5, 2, "hidden SPOF (wire an independent route)"
        elif is_defect:                                # resilient but wrong
            gain, effort, cat = 0.5, 2, "defect edge (add a known correction layer)"
        else:
            continue                                   # resilient + admissible: nothing to do
        items.append({"benchmark": name, "category": cat, "gain": round(gain, 2),
                      "effort": effort, "roi": round(gain / effort, 3)})
    items.sort(key=lambda r: (-r["roi"], -r["gain"], r["benchmark"]))
    triage.sort(key=lambda r: (-r["potential"], r["benchmark"]))
    return {"items": items, "triage": triage,
            "closure": round(total, 1), "ceiling": len(outs)}


def _priority_lines(g: DerivationDAG, report) -> list[str]:
    """The closure-priority readout: current closure + the effort-known ROI
    ranking + the triage queue + the directional "do first" line."""
    res = closure_priority(g, report)
    pct = res["closure"] / res["ceiling"] * 100 if res["ceiling"] else 0.0
    lines = ["", "Closure priority (M. Wende — closure gained per unit effort; "
             "the directional dual of the diagnostic layers):",
             f"  structural closure now: {res['closure']} / {res['ceiling']} "
             f"({pct:.0f}% derived-and-resilient)",
             f"  [A] effort-known ({len(res['items'])}) — ROI = gain / effort "
             f"(tiers: 1 mechanical, 2 known-layer, 4 research):"]
    for rank, it in enumerate(res["items"], 1):
        lines.append(f"    {rank:>2}  ROI {it['roi']:.2f}  +{it['gain']:.1f}/eff{it['effort']}  "
                     f"{it['benchmark']:24s} {it['category']}")
    lines.append(f"  [B] triage ({len(res['triage'])}) — high *potential* closure, effort unknown: "
                 f"sweep can't reach these, so each is a leak (cheap) or a gap (dear) until read:")
    for t in res["triage"]:
        lines.append(f"      +{t['potential']:.1f}  {t['benchmark']}"
                     f"{'  [+defect]' if t['defect'] else ''}")
    t1 = [i for i in res["items"] if i["effort"] == 1]
    if t1:
        lines.append(f"  → do first (verified tier-1): {', '.join(i['benchmark'] for i in t1)} "
                     f"(+{sum(i['gain'] for i in t1):.1f} closure, mechanical cost)")
    return lines


def _main_suite(use_sensitivity: bool = False, show_redundancy: bool = False,
                show_priority: bool = False) -> int:
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
    lines += _provenance_lines(g)
    lines += _collapse_lines(g)
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
            lines += _redundancy_lines(g, list_all=False, label="benchmarks", report=report)
    if show_priority:
        if report is None:
            lines += ["", "(add --sensitivity for closure-priority — it ranks the open "
                      "items by closure-per-effort, which needs the sweep's coupling.)"]
        else:
            lines += _priority_lines(g, report)
    print("\n".join(lines))
    return 0


def main(argv: list[str] | None = None) -> int:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if "--suite" in argv:
        return _main_suite(use_sensitivity="--sensitivity" in argv,
                           show_redundancy="--redundancy" in argv,
                           show_priority="--priority" in argv)
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

    lines += _provenance_lines(g)
    lines += _dispute_lines(g)
    lines += _snow_lines(g)
    lines += _collapse_lines(g)

    if "--redundancy" in (argv or []):
        lines += _redundancy_lines(g)

    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
