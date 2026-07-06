"""
Regression test for `mesh.metrics.build_face_adjacency`.

Not one of the G0.1-G0.4 gates (design.md/roadmap.md don't require face
adjacency until P4's upwind element search), but it's part of P0's stated
metrics.py deliverable and had never actually been exercised: it crashed
under @njit (numba typed Dict can't hold growable-list values in nopython
mode) until fixed. This locks the fix in.
"""

import numpy as np
import pytest

from pyfp3d.mesh.metrics import build_face_adjacency

from .mesh_utils import generate_structured_cube_mesh

_FACE_LOCAL_NODES = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


def _reference_face_adjacency(elements: np.ndarray) -> np.ndarray:
    """Brute-force O(n_tets) dict-of-lists reference (plain Python, no numba)."""
    face_to_tets = {}
    for e, tet in enumerate(elements):
        for f, local in enumerate(_FACE_LOCAL_NODES):
            key = tuple(sorted(tet[local]))
            face_to_tets.setdefault(key, []).append((e, f))

    ref = np.full((len(elements), 4), -1, dtype=np.int32)
    for key, occurrences in face_to_tets.items():
        assert len(occurrences) <= 2, f"non-manifold face {key}: {occurrences}"
        if len(occurrences) == 2:
            (e1, f1), (e2, f2) = occurrences
            ref[e1, f1] = e2
            ref[e2, f2] = e1
    return ref


class TestFaceAdjacency:
    def test_runs_under_njit(self):
        """Smoke test: build_face_adjacency must not crash (see docstring)."""
        elements = np.array([
            [0, 1, 3, 4],
            [1, 2, 3, 6],
            [1, 3, 4, 6],
            [3, 4, 6, 7],
            [1, 4, 5, 6],
        ], dtype=np.int32)
        face_neighbors, _ = build_face_adjacency(elements)
        assert face_neighbors.shape == (5, 4)

    def test_matches_reference_on_structured_cube(self):
        _, elements = generate_structured_cube_mesh(n=3, L=1.0)
        face_neighbors, _ = build_face_adjacency(elements)
        reference = _reference_face_adjacency(elements)
        assert np.array_equal(face_neighbors, reference)

    def test_neighbor_relation_is_symmetric(self):
        _, elements = generate_structured_cube_mesh(n=3, L=1.0)
        face_neighbors, _ = build_face_adjacency(elements)

        for e in range(len(elements)):
            for f in range(4):
                nb = face_neighbors[e, f]
                if nb == -1:
                    continue
                assert e in face_neighbors[nb], f"element {e} face {f} -> {nb} not reciprocated"

    def test_boundary_face_count_matches_surface_area_expectation(self):
        # A structured n x n x n cube has 6 faces, each split into 2n^2
        # boundary triangles -> 12 n^2 boundary faces total.
        n = 4
        _, elements = generate_structured_cube_mesh(n=n, L=1.0)
        face_neighbors, _ = build_face_adjacency(elements)
        n_boundary = int(np.sum(face_neighbors == -1))
        assert n_boundary == 12 * n * n


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
