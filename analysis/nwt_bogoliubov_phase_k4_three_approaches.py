"""Bogoliubov Phase K-4 — three new approaches to σ-orbit → carrier-knot.

Phase K-3 left LEPTONS + HEAVY HADRONS uncovered by the sector-conditional
σ-graph cyclomatic rule. K-4 tests three new constructions to close
them:

1. Gauss linking with reference Hamilton cycles (d=1 QR, d=2 QR, d=3 NR)
2. Spectral invariants of the 7x7 σ-transition matrix
3. Heffter face traversal count (the 14 triangular faces of K_7 embedding)

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_k4_three_approaches.py
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# numpy compat for pyknotid
np.float = float
np.int = int
np.complex = complex
np.object = object
np.str = str
np.long = int

try:
    from pyknotid.spacecurves import Link
    PYKNOTID = True
except Exception as e:
    PYKNOTID = False
    print(f"pyknotid unavailable: {e}")

from nwt_substrate.condensate import (
    bfs_shortest_walks, walk_sigma_sequence, walk_to_3d_curve,
)
from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_k4_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


# Reference Hamilton cycles (skip-d Hamilton)
def hamilton_walk_skip(d: int) -> list[int]:
    return [(k * d) % 7 for k in range(7)] + [0]

HAM_D1 = hamilton_walk_skip(1)  # QR matter scaffold
HAM_D2 = hamilton_walk_skip(2)  # QR secondary
HAM_D3 = hamilton_walk_skip(3)  # NR CP-channel


# ---------------------------------------------------------------------------
# Approach 1: Gauss linking via pyknotid Link()
# ---------------------------------------------------------------------------

def gauss_linking(walk: list[int], ref: list[int]) -> int | None:
    """Compute Gauss linking number between two K_7 walks (lifted to 3D)."""
    if not PYKNOTID:
        return None
    try:
        curve_a = walk_to_3d_curve(walk, n_per_edge=20)
        curve_b = walk_to_3d_curve(ref, n_per_edge=20)
        # Shift the reference curve slightly so it's not coincident
        curve_b = curve_b + np.array([0.05, 0.05, 0.05])
        link = Link([curve_a, curve_b])
        # Linking number from Gauss code
        try:
            lk = link.linking_number()
            return int(lk) if lk is not None else 0
        except Exception:
            return None
    except Exception as e:
        return None


# ---------------------------------------------------------------------------
# Approach 2: spectral invariants of σ-transition matrix
# ---------------------------------------------------------------------------

def sigma_transition_matrix(walk: list[int]) -> np.ndarray:
    """7x7 matrix M[i,j] = # transitions σ_i → σ_j in the walk."""
    sig_seq = walk_sigma_sequence(walk)
    n = len(sig_seq)
    M = np.zeros((7, 7), dtype=int)
    for i in range(n):
        cur = sig_seq[i]
        nxt = sig_seq[(i + 1) % n]
        M[cur, nxt] += 1
    return M


def spectral_invariants(walk: list[int]) -> dict:
    """Various spectral / matrix invariants of the σ-transition matrix."""
    M = sigma_transition_matrix(walk)
    # Symmetrize for undirected counts
    M_sym = (M + M.T) / 2
    # Eigenvalues of symmetric part
    eigvals = np.linalg.eigvalsh(M_sym)
    # Non-zero eigenvalues
    n_nz_eigs = int(np.sum(np.abs(eigvals) > 1e-9))
    # Trace, determinant, frobenius norm
    trace = float(np.trace(M))      # = total self-loops
    trace_sym = float(np.trace(M_sym))
    spectral_radius = float(np.max(np.abs(eigvals)))
    rank = int(np.linalg.matrix_rank(M_sym, tol=1e-9))
    # # of distinct nonzero entries
    return {
        "trace": trace,
        "trace_sym": trace_sym,
        "rank": rank,
        "n_nz_eigs": n_nz_eigs,
        "spectral_radius": round(spectral_radius, 3),
        "max_eig_int": int(round(spectral_radius)),
    }


# ---------------------------------------------------------------------------
# Approach 3: Heffter face traversal
# ---------------------------------------------------------------------------

def heffter_faces() -> list[tuple[int, int, int]]:
    """14 triangular faces of the K_7 Heffter embedding."""
    from nwt_substrate.topology.K7 import heffter_rotation, trace_K7_faces
    rot = heffter_rotation()
    faces = trace_K7_faces(rot)
    return [tuple(sorted(f[:3])) for f in faces if len(f) == 3]


HEFFTER_FACES = heffter_faces()


def edge_to_faces(a: int, b: int) -> list[tuple[int, int, int]]:
    """Two Heffter faces that share edge (a, b)."""
    e = (min(a, b), max(a, b))
    out = []
    for f in HEFFTER_FACES:
        edges = [(min(f[i], f[(i+1)%3]), max(f[i], f[(i+1)%3])) for i in range(3)]
        if e in edges:
            out.append(f)
    return out


def face_traversal_invariants(walk: list[int]) -> dict:
    """Heffter face traversal for the walk."""
    faces_used = set()
    face_transitions = []
    for i in range(len(walk) - 1):
        faces = edge_to_faces(walk[i], walk[i + 1])
        for f in faces:
            faces_used.add(f)
        face_transitions.append(faces)
    return {
        "n_faces_used": len(faces_used),
        "n_faces_total": len(HEFFTER_FACES),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE K-4 — three new approaches to σ-orbit → carrier")
    print("=" * 78)
    print()
    print(f"  pyknotid available: {PYKNOTID}")
    print(f"  HAM_D1 (matter): {HAM_D1}")
    print(f"  HAM_D2 (QR sec): {HAM_D2}")
    print(f"  HAM_D3 (NR CP):  {HAM_D3}")
    print()
    print(f"  Heffter faces: {len(HEFFTER_FACES)} (expect 14)")
    print(f"  First few: {HEFFTER_FACES[:3]}")
    print()

    walks = bfs_shortest_walks(max_length=25)
    print(f"Loaded {len(walks)} (|p|, |q|) shortest walks")
    print()

    print(f"{'particle':<12} {'n_q':<4} {'lk_d1':<7} {'lk_d2':<7} {'lk_d3':<7} "
          f"{'rank':<5} {'sp_rad':<7} {'n_face':<7}")
    print("-" * 70)
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        # Gauss linking (may be slow)
        lk1 = gauss_linking(walk, HAM_D1) if PYKNOTID else None
        lk2 = gauss_linking(walk, HAM_D2) if PYKNOTID else None
        lk3 = gauss_linking(walk, HAM_D3) if PYKNOTID else None
        spec = spectral_invariants(walk)
        face = face_traversal_invariants(walk)
        rows.append({
            **entry, "walk": walk,
            "lk_d1": lk1, "lk_d2": lk2, "lk_d3": lk3,
            **spec, **face,
        })
        lk1_s = str(lk1) if lk1 is not None else "?"
        lk2_s = str(lk2) if lk2 is not None else "?"
        lk3_s = str(lk3) if lk3 is not None else "?"
        print(f"{entry['name']:<12} {entry['n_q']:<4} {lk1_s:<7} {lk2_s:<7} "
              f"{lk3_s:<7} {spec['rank']:<5} "
              f"{spec['spectral_radius']:<7} {face['n_faces_used']:<7}")
    print()

    # ---- Match counts -------------------------------------------------
    print("=" * 78)
    print("MATCH COUNTS — does any single invariant give n_q?")
    print("=" * 78)
    print()
    candidates = {
        "lk_d1 (linking with d=1 Ham)": "lk_d1",
        "lk_d2 (linking with d=2 Ham)": "lk_d2",
        "lk_d3 (linking with d=3 Ham)": "lk_d3",
        "|lk_d1|": lambda r: abs(r["lk_d1"]) if r["lk_d1"] is not None else None,
        "lk_d1 + lk_d2 + lk_d3": lambda r: sum(r[f"lk_d{i}"] for i in [1,2,3]) if all(r[f"lk_d{i}"] is not None for i in [1,2,3]) else None,
        "|lk_d1| + |lk_d2| + |lk_d3|": lambda r: sum(abs(r[f"lk_d{i}"]) for i in [1,2,3]) if all(r[f"lk_d{i}"] is not None for i in [1,2,3]) else None,
        "spectral_radius (rounded)": "max_eig_int",
        "rank": "rank",
        "n_faces_used": "n_faces_used",
        "n_faces_used / 2": lambda r: r["n_faces_used"] // 2,
        "n_faces_used − 14": lambda r: abs(r["n_faces_used"] - 14),
    }
    print(f"  {'invariant':<35} {'exact':<8} {'Pearson r'}")
    print("  " + "-" * 55)
    for name, key in candidates.items():
        if callable(key):
            vals = [key(r) for r in rows]
        else:
            vals = [r[key] for r in rows]
        # Skip if any None
        valid = [(r["n_q"], v) for r, v in zip(rows, vals) if v is not None]
        if not valid:
            print(f"  {name:<35} (all None)")
            continue
        n_q_arr = np.array([x[0] for x in valid])
        v_arr = np.array([x[1] for x in valid])
        exact = int(sum(1 for nq, v in valid if nq == v))
        try:
            r_p = np.corrcoef(n_q_arr, v_arr)[0, 1] if len(set(v_arr)) > 1 else float('nan')
        except Exception:
            r_p = float('nan')
        print(f"  {name:<35} {exact}/{len(valid)}    {r_p:+.3f}")
    print()

    # ---- By-sector breakdown for promising candidates -----------------
    from collections import defaultdict
    print("=" * 78)
    print("BY-SECTOR breakdown for spectral invariants")
    print("=" * 78)
    print()
    for attr in ["rank", "max_eig_int", "n_faces_used"]:
        print(f"  {attr}:")
        by_nq = defaultdict(list)
        for r in rows:
            by_nq[r["n_q"]].append(r[attr])
        for nq in sorted(by_nq.keys()):
            vals = sorted(set(by_nq[nq]))
            sector = {0:'leptons', 2:'mesons', 3:'hyperons', 5:'nucleons'}.get(nq, '?')
            print(f"    n_q={nq} ({sector}): values = {vals}")
        print()

    # Save data
    np.savez(OUT_DIR / "phase_k4_three_approaches.npz",
             names=np.array([r["name"] for r in rows]),
             n_q=np.array([r["n_q"] for r in rows]),
             lk_d1=np.array([r["lk_d1"] if r["lk_d1"] is not None else 999 for r in rows]),
             lk_d2=np.array([r["lk_d2"] if r["lk_d2"] is not None else 999 for r in rows]),
             lk_d3=np.array([r["lk_d3"] if r["lk_d3"] is not None else 999 for r in rows]),
             rank=np.array([r["rank"] for r in rows]),
             spectral_radius=np.array([r["spectral_radius"] for r in rows]),
             n_faces_used=np.array([r["n_faces_used"] for r in rows]))
    print(f"  data saved {OUT_DIR / 'phase_k4_three_approaches.npz'}")


if __name__ == "__main__":
    main()
