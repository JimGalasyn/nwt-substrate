"""Tests for nwt_substrate.condensate.walks."""
import pytest

from nwt_substrate.condensate.walks import bfs_shortest_walks


def test_bfs_rejects_invalid_start():
    """bfs_shortest_walks() guards the K_7 vertex range (0 .. N_VERTICES_K7-1)."""
    with pytest.raises(ValueError, match="start vertex"):
        bfs_shortest_walks(start=7)
    with pytest.raises(ValueError, match="start vertex"):
        bfs_shortest_walks(start=-1)
