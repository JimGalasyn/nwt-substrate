"""
K_7 graph state and toroidal embedding (Heffter 1891).

K_7 (the complete graph on 7 vertices) admits a triangular embedding on
the torus with V=7, E=21, F=14, satisfying V - E + F = 0 (genus 1).

This embedding is the natural "background geometry" for both:
  - matter sector: (p, q) torus knots host vortex solitons
  - gravity sector: Wilson amplitude on 21 edges -> alpha^(21/2)

The Heffter rotation system at vertex i has cyclic order
  (i+1, i+3, i+2, i+6, i+4, i+5) mod 7.
"""

from __future__ import annotations


def heffter_rotation():
    """
    Heffter's rotation system for K_7's toroidal triangular embedding.

    Returns dict {v: [n1, ..., n6]} where each list is the cyclic order
    of v's neighbors in the embedding.
    """
    base = [1, 3, 2, 6, 4, 5]
    return {v: [(v + b) % 7 for b in base] for v in range(7)}


def trace_K7_faces(rotation):
    """
    Trace the face structure of an embedding from a rotation system.

    Algorithm: each directed edge (u, v) belongs to one face.  At v, find
    u's index in v's rotation and back-step to get the next edge.

    Returns list of faces; each face is a tuple of vertices.
    """
    visited = set()
    faces = []
    for u in rotation:
        for v in rotation[u]:
            if (u, v) in visited:
                continue
            face = []
            cur_u, cur_v = u, v
            while (cur_u, cur_v) not in visited:
                visited.add((cur_u, cur_v))
                face.append(cur_u)
                rot = rotation[cur_v]
                idx = rot.index(cur_u)
                next_v = rot[(idx - 1) % len(rot)]
                cur_u, cur_v = cur_v, next_v
            faces.append(tuple(face))
    return faces


def is_genus_one_embedding(rotation) -> tuple[bool, dict]:
    """
    Check whether an embedding has Euler characteristic 0 (torus, genus 1)
    and all triangular faces.

    Returns (passes, info_dict).
    """
    faces = trace_K7_faces(rotation)
    V = len(rotation)
    E = sum(len(rotation[v]) for v in rotation) // 2
    F = len(faces)
    chi = V - E + F
    sizes = sorted(len(f) for f in faces)
    info = {
        "V": V, "E": E, "F": F, "chi": chi,
        "face_sizes": sizes,
        "all_triangles": all(s == 3 for s in sizes),
        "genus": (2 - chi) // 2,
    }
    return (chi == 0 and info["all_triangles"]), info
