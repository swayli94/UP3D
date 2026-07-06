"""
M0 gate checks for the single-layer extruded NACA0012 quasi-2D mesh family
(cases/meshes/naca0012_2.5d, generate_naca0012.py; docs/roadmap.md Track M).

Mesh-side M0 gate items covered here:
  - solver preprocessor (mesh/reader.py) ingests the family; tags complete
    (wall / farfield / symmetry / wake)
  - globally consistent prism->3-tet split: quad-split consistency assert
  - wake is ONE continuous planar interior sheet from the TE (x = 1) to the
    far field, conforming to tet faces, nodes NOT duplicated
  - symmetry planes planar, non-overlapping with the wall
  - single spanwise station (exactly two z-planes) -- Gamma will be a single
    scalar for the P2 Kutta loop
  - a non-lifting (alpha = 0) Laplace solve runs end to end on the coarse
    mesh (the lifting G2.3/G2.5(b) checks need P2's wake cut + Kutta loop)

The P2 wake-cut topology asserts (one +/- side element per wake face, TE
nodes not duplicated, ...) live in mesh/wake_cut.py when P2 opens; they are
NOT duplicated here.
"""

from collections import defaultdict

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.extrude import assert_quad_split_consistency
from pyfp3d.solve.picard import solve_laplace

R_FAR = 15.0
LEVELS = ["coarse", "medium"]


@pytest.fixture(scope="module")
def meshes(request):
    from pathlib import Path

    mesh_dir = Path(__file__).parent.parent / "cases" / "meshes" / "naca0012_2.5d"
    return {level: read_mesh(mesh_dir / f"{level}.msh") for level in LEVELS}


class TestNacaMeshFamily:
    @pytest.mark.parametrize("level", LEVELS)
    def test_ingestion_tags_and_consistency(self, meshes, level):
        mesh = meshes[level]
        assert set(mesh.boundary_faces) == {"wall", "farfield", "symmetry", "wake"}
        assert_quad_split_consistency(mesh, interior_groups=("wake",))

    @pytest.mark.parametrize("level", LEVELS)
    def test_single_layer_two_zplanes(self, meshes, level):
        mesh = meshes[level]
        z = np.round(mesh.nodes[:, 2], 12)
        planes = sorted(set(z))
        assert len(planes) == 2, "must be exactly one cell layer in z"
        for tri in mesh.boundary_faces["symmetry"]:
            assert len(set(z[tri])) == 1, "symmetry face not planar"

    @pytest.mark.parametrize("level", LEVELS)
    def test_symmetry_wall_disjoint(self, meshes, level):
        mesh = meshes[level]
        sym = set(map(tuple, np.sort(mesh.boundary_faces["symmetry"], axis=1).tolist()))
        wall = set(map(tuple, np.sort(mesh.boundary_faces["wall"], axis=1).tolist()))
        assert not (sym & wall), "symmetry and wall face sets overlap"

    @pytest.mark.parametrize("level", LEVELS)
    def test_wake_sheet_topology(self, meshes, level):
        """Wake: one continuous planar sheet, TE -> farfield, on both z-planes,
        no duplicated nodes (duplication is solver-side, P2)."""
        mesh = meshes[level]
        nodes = mesh.nodes
        wake = mesh.boundary_faces["wake"]
        wake_nodes = np.unique(wake)

        # planar (y = 0 exactly: the wake line is a straight Gmsh line and
        # extrusion does not move nodes)
        assert np.max(np.abs(nodes[wake_nodes, 1])) == 0.0

        # spans TE (x = 1) to the far-field circle (x = 0.5 + R_FAR)
        assert abs(nodes[wake_nodes, 0].min() - 1.0) < 1e-12
        assert abs(nodes[wake_nodes, 0].max() - (0.5 + R_FAR)) < 1e-9

        # covers both z-planes
        assert len(set(np.round(nodes[wake_nodes, 2], 12))) == 2

        # single edge-connected component
        edge_owner = defaultdict(list)
        for i, tri in enumerate(np.sort(wake, axis=1)):
            a, b, c = map(int, tri)
            for e in ((a, b), (a, c), (b, c)):
                edge_owner[e].append(i)
        parent = list(range(len(wake)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for owners in edge_owner.values():
            for other in owners[1:]:
                ra, rb = find(owners[0]), find(other)
                if ra != rb:
                    parent[ra] = rb
        n_components = len({find(i) for i in range(len(wake))})
        assert n_components == 1, f"wake sheet has {n_components} components"

        # nodes not duplicated: every wake node is also used by the volume mesh
        vol_nodes = np.unique(mesh.elements)
        assert np.all(np.isin(wake_nodes, vol_nodes))

    @pytest.mark.parametrize("level", LEVELS)
    def test_wall_is_unit_chord_airfoil(self, meshes, level):
        mesh = meshes[level]
        wall_nodes = np.unique(mesh.boundary_faces["wall"])
        x = mesh.nodes[wall_nodes, 0]
        y = mesh.nodes[wall_nodes, 1]
        assert abs(x.min()) < 1e-6 and abs(x.max() - 1.0) < 1e-12
        assert np.max(np.abs(y)) < 0.0605  # NACA0012 max half-thickness ~ 0.06

    def test_refinement_is_monotone(self, meshes):
        n = [len(meshes[level].elements) for level in LEVELS]
        assert n[0] < n[1], "medium must refine coarse"


class TestNacaNonliftingSolve:
    """End-to-end ingestion into the P1 solver: alpha = 0, no wake cut yet
    (non-lifting), freestream Dirichlet at the far field."""

    def test_alpha0_solve_and_spanwise(self, meshes):
        mesh = meshes["coarse"]
        nodes, elements = mesh.nodes, mesh.elements
        farfield_nodes = np.unique(mesh.boundary_faces["farfield"])

        result = solve_laplace(
            nodes, elements, farfield_nodes, nodes[farfield_nodes, 0],
            rtol=1e-11, maxiter=5000,
        )
        assert result["residual_norm"] < 1e-8

        # G2.5(a): interpolated freestream is exactly spanwise-zero
        from tests.test_m0_cylinder import element_gradients_all

        grad_x = element_gradients_all(nodes, elements, nodes[:, 0].copy())
        assert np.max(np.abs(grad_x[:, 2])) < 1e-12

        # solved field: O(h) spanwise noise, inherent to the 3-tet prism
        # split (see tests/test_m0_cylinder.py docstring; measured 5.4e-2)
        grad = element_gradients_all(nodes, elements, result["phi"])
        assert np.max(np.abs(grad[:, 2])) < 0.10


class TestNacaArtifacts:
    """Headless evidence per roadmap Sec 0.1 / M0 visual tests: edge-length
    histogram across levels on one figure (+ per-level layer PNGs are written
    by generate_naca0012.py at generation time)."""

    def test_export_m0_naca_edge_histogram(self, gate_artifacts_dir, meshes):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from pyfp3d.mesh.metrics import compute_edge_lengths

        fig, ax = plt.subplots(figsize=(8, 6))
        summary = []
        for level in LEVELS:
            mesh = meshes[level]
            lengths = compute_edge_lengths(mesh.nodes, mesh.elements)
            ax.hist(np.log10(lengths), bins=60, histtype="step",
                    label=f"{level} (n_tets={len(mesh.elements)})")
            summary.append((level, len(mesh.elements),
                            lengths.min(), lengths.max(), lengths.mean()))
        ax.set_xlabel("log10(edge length)")
        ax.set_ylabel("count")
        ax.set_title("M0 naca0012_2.5d: edge-length distribution by level")
        ax.legend()
        ax.grid(True, alpha=0.3)
        out_png = gate_artifacts_dir / "edge_length_histogram.png"
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert out_png.exists()

        with open(gate_artifacts_dir / "summary.csv", "w") as f:
            f.write("level,n_tets,min_edge,max_edge,mean_edge\n")
            for row in summary:
                f.write(f"{row[0]},{row[1]},{row[2]:.6e},{row[3]:.6e},{row[4]:.6e}\n")


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
