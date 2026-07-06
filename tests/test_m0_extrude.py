"""
M0 unit tests: single-layer extrusion + globally consistent prism -> 3-tet
split (pyfp3d/meshgen/extrude.py). Pure numpy, no Gmsh dependency.

Covers the M0 spec items (docs/roadmap.md Track M):
  - exactly one cell layer between two parallel planes, both tagged "symmetry"
  - 3 tets per prism, all positive volume, exact total volume
  - min-global-index diagonal rule -> every interior quad face split
    identically from both sides (assert_quad_split_consistency)
  - tagged interior sheets ("wake") coincide with conforming tet faces
  - the consistency assert actually *fires* on a broken split
"""

import numpy as np
import pytest

from pyfp3d.mesh.metrics import compute_tet_volumes
from pyfp3d.meshgen.extrude import (
    assert_quad_split_consistency,
    extrude_single_layer,
)


def unit_square_two_triangles():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    triangles = np.array([[0, 1, 2], [0, 2, 3]])
    boundary = {
        "wall": np.array([[0, 1], [1, 2]]),
        "farfield": np.array([[2, 3], [3, 0]]),
    }
    return points, triangles, boundary


def random_disk_mesh(n_points=120, seed=7):
    """Delaunay triangulation of random points in a disk (convex domain, so
    every Delaunay triangle is inside)."""
    from scipy.spatial import Delaunay

    rng = np.random.default_rng(seed)
    r = np.sqrt(rng.uniform(0.05, 1.0, n_points))
    t = rng.uniform(0, 2 * np.pi, n_points)
    interior = np.stack([r * np.cos(t), r * np.sin(t)], axis=1)
    ring = np.stack(
        [np.cos(np.linspace(0, 2 * np.pi, 48, endpoint=False)),
         np.sin(np.linspace(0, 2 * np.pi, 48, endpoint=False))], axis=1
    )
    points = np.vstack([ring, interior])
    dela = Delaunay(points)
    triangles = dela.simplices
    hull = dela.convex_hull  # (m, 2) boundary edges
    return points, triangles, {"farfield": hull}


class TestSingleLayerExtrusion:
    def test_two_triangle_square(self):
        points, triangles, boundary = unit_square_two_triangles()
        dz = 0.25
        mesh = extrude_single_layer(points, triangles, boundary, dz=dz, z0=-0.5)

        assert len(mesh.nodes) == 8
        assert len(mesh.elements) == 6  # 2 prisms x 3 tets

        vols = compute_tet_volumes(mesh.nodes, mesh.elements)
        assert np.all(vols > 0), "prism split produced inverted tets"
        assert abs(vols.sum() - 1.0 * dz) < 1e-14

        z = mesh.nodes[:, 2]
        assert set(np.round(z, 12)) == {-0.5, -0.25}
        sym = mesh.boundary_faces["symmetry"]
        assert len(sym) == 4  # 2 triangles per plane
        for tri in sym:
            assert len(set(np.round(z[tri], 12))) == 1, "symmetry face not planar"

        assert len(mesh.boundary_faces["wall"]) == 4  # 2 edges x 2 triangles
        assert len(mesh.boundary_faces["farfield"]) == 4

    def test_random_disk_positive_volumes_and_consistency(self):
        points, triangles, boundary = random_disk_mesh()
        dz = 0.1
        mesh = extrude_single_layer(points, triangles, boundary, dz=dz)

        vols = compute_tet_volumes(mesh.nodes, mesh.elements)
        assert np.all(vols > 0)

        # Exact prism volume: sum of 2D triangle areas x dz
        p = points
        t = triangles
        area = 0.5 * np.abs(
            (p[t[:, 1], 0] - p[t[:, 0], 0]) * (p[t[:, 2], 1] - p[t[:, 0], 1])
            - (p[t[:, 1], 1] - p[t[:, 0], 1]) * (p[t[:, 2], 0] - p[t[:, 0], 0])
        )
        assert abs(vols.sum() - area.sum() * dz) < 1e-12 * area.sum()

        # The constructor already runs the assert; run it again explicitly.
        assert_quad_split_consistency(mesh)

    def test_interior_wake_sheet_conforms(self):
        """Tag an interior edge as 'wake': its extruded faces must be interior
        tet faces shared by exactly two tets."""
        points, triangles, boundary = unit_square_two_triangles()
        wake = {"wake": np.array([[0, 2]])}  # the shared diagonal edge
        mesh = extrude_single_layer(points, triangles, boundary,
                                    interior_edge_groups=wake, dz=0.25)
        assert len(mesh.boundary_faces["wake"]) == 2
        assert_quad_split_consistency(mesh, interior_groups=("wake",))

    def test_consistency_assert_fires_on_broken_split(self):
        """Re-split one prism with the OPPOSITE diagonal rule: the shared
        interior quad face is now split differently from its two sides and
        the M0 assert must fire."""
        points, triangles, boundary = unit_square_two_triangles()
        mesh = extrude_single_layer(points, triangles, boundary, dz=0.25)
        n2 = 4

        # Prism over triangle (0, 1, 2): rebuild its 3 tets with the diagonal
        # on quad (1, 2, 2', 1') flipped from 1-2' to 2-1'.
        broken = mesh.elements.copy()
        v0, v1, v2 = 0, 1, 2
        v3, v4, v5 = v0 + n2, v1 + n2, v2 + n2
        broken[0] = [v0, v1, v2, v4]
        broken[1] = [v0, v4, v2, v5]
        broken[2] = [v0, v4, v5, v3]
        # (this is the *other* valid single-prism split, so each tet is fine
        # in isolation -- only the shared-face consistency is violated)
        mesh.elements = broken
        with pytest.raises(
            AssertionError,
            match="split inconsistently|not owned|not faces of any tet",
        ):
            assert_quad_split_consistency(mesh)

    def test_rejects_bad_input(self):
        points, triangles, boundary = unit_square_two_triangles()
        with pytest.raises(ValueError, match="dz"):
            extrude_single_layer(points, triangles, boundary, dz=0.0)
        with pytest.raises(ValueError, match="referenced by no"):
            extrude_single_layer(
                np.vstack([points, [[5.0, 5.0]]]), triangles, boundary, dz=0.1
            )


class TestFreestreamSpanwiseZero:
    """G2.5(a) preview: the nodal interpolant of phi = x on a single-layer
    mesh has machine-zero spanwise gradient in every element."""

    def test_phi_x_gradient_z_component(self):
        points, triangles, boundary = random_disk_mesh(n_points=80, seed=3)
        mesh = extrude_single_layer(points, triangles, boundary, dz=0.07)
        phi = mesh.nodes[:, 0].copy()

        el = mesh.elements
        p = mesh.nodes
        e = np.stack(
            [p[el[:, 1]] - p[el[:, 0]],
             p[el[:, 2]] - p[el[:, 0]],
             p[el[:, 3]] - p[el[:, 0]]], axis=1
        )  # (n_tets, 3 edges, 3 xyz)
        d = np.stack(
            [phi[el[:, 1]] - phi[el[:, 0]],
             phi[el[:, 2]] - phi[el[:, 0]],
             phi[el[:, 3]] - phi[el[:, 0]]], axis=1
        )
        grad = np.linalg.solve(e, d[:, :, None])[:, :, 0]  # (n_tets, 3)
        assert np.max(np.abs(grad[:, 2])) < 1e-12
        assert np.max(np.abs(grad[:, 0] - 1.0)) < 1e-12


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
