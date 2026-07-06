"""
Regression test for `mesh.reader.read_mesh` / `write_mesh` round-tripping
tagged boundary surfaces.

Gate G0.4 (tests/test_io_vtk.py) only exercises `post/vtk_out.py`'s
write_vtu/read_vtu pair on an untagged synthetic mesh -- a different code
path from `mesh/reader.py`'s read_mesh/write_mesh, which is what actually
loads/saves the named ("wall", "farfield", ...) physical-surface tags used
throughout the solver (Dirichlet BCs, wake topology asserts, etc). That gap
let write_mesh silently drop every named boundary group (it only recognized
a legacy "all_triangles" block) and let ".msh" resolve to meshio's "ansys"
writer instead of "gmsh" (dropping all tag data outright) go unnoticed.
This test locks in the fix using the real committed sphere-shell case.
"""

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh, write_mesh


@pytest.fixture
def sphere_coarse_path():
    from pathlib import Path
    path = Path(__file__).parent.parent / "cases" / "meshes" / "sphere_shell" / "coarse.msh"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return path


def _sorted_triangle_set(faces: np.ndarray) -> np.ndarray:
    rows = np.sort(faces, axis=1)
    return rows[np.lexsort(rows.T[::-1])]


class TestMeshReaderRoundTrip:
    def test_tagged_boundary_groups_survive_roundtrip(self, sphere_coarse_path, tmp_path):
        mesh = read_mesh(sphere_coarse_path)
        assert set(mesh.boundary_faces) == {"wall", "farfield"}, (
            "fixture assumption changed: sphere_shell/coarse.msh should have "
            "exactly the 'wall' and 'farfield' physical surfaces"
        )

        out_path = tmp_path / "roundtrip.msh"
        write_mesh(mesh, out_path)
        mesh2 = read_mesh(out_path)

        assert set(mesh2.boundary_faces) == set(mesh.boundary_faces)
        assert np.allclose(mesh.nodes, mesh2.nodes)
        assert np.array_equal(mesh.elements, mesh2.elements)
        assert np.array_equal(mesh.element_tags, mesh2.element_tags)

        for name, faces in mesh.boundary_faces.items():
            assert np.array_equal(
                _sorted_triangle_set(faces),
                _sorted_triangle_set(mesh2.boundary_faces[name]),
            ), f"boundary group '{name}' did not round-trip"

    def test_msh_extension_resolves_to_gmsh_not_ansys(self, sphere_coarse_path, tmp_path):
        """".msh" is ambiguous in meshio (matches both "ansys" and "gmsh");
        write_mesh must not silently fall back to a writer that can't
        represent our cell/physical-tag data."""
        mesh = read_mesh(sphere_coarse_path)
        out_path = tmp_path / "format_check.msh"
        write_mesh(mesh, out_path)

        # An ansys-format write would not round-trip as tetra/triangle cells
        # with gmsh:physical tags at all -- read_mesh would either raise or
        # come back with no boundary_faces.
        mesh2 = read_mesh(out_path)
        assert len(mesh2.elements) == len(mesh.elements)
        assert mesh2.boundary_faces  # non-empty


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
