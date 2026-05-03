"""Tests for nwt_substrate.topology.K7."""

from nwt_substrate.topology.K7 import (
    heffter_rotation,
    trace_K7_faces,
    is_genus_one_embedding,
)


def test_heffter_has_seven_vertices():
    rot = heffter_rotation()
    assert len(rot) == 7
    for v in range(7):
        assert v in rot
        assert len(rot[v]) == 6   # K_7 -> each vertex has 6 neighbors


def test_K7_torus_embedding():
    """Heffter 1891: K_7 embeds on torus with V=7, E=21, F=14, chi=0."""
    rot = heffter_rotation()
    passes, info = is_genus_one_embedding(rot)
    assert passes, f"Embedding failed: {info}"
    assert info["V"] == 7
    assert info["E"] == 21
    assert info["F"] == 14
    assert info["chi"] == 0
    assert info["genus"] == 1
    assert info["all_triangles"]


def test_K7_face_count_matches_traced_faces():
    rot = heffter_rotation()
    faces = trace_K7_faces(rot)
    assert len(faces) == 14
    assert all(len(f) == 3 for f in faces)
