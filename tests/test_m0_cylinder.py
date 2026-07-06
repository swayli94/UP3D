"""
M0 cylinder-flow test case (quasi-2D pipeline validation).

Incompressible potential flow past a circular cylinder on the single-layer
extruded annulus mesh (cases/meshes/cylinder_2.5d, generate_cylinder.py):

    phi = U x (1 + a^2 / r^2),    surface Cp = 1 - 4 sin^2(theta)

Dirichlet exact potential at the far field, natural (zero-flux) BC on the
wall and both symmetry planes. This is the simplest end-to-end check that
the M0 meshing pipeline (vanilla-Gmsh 2D + consistent prism->3-tet split)
feeds the existing P1 solver correctly, complementing the lifting NACA0012
mesh which is the actual M0 deliverable (tests/test_m0_naca0012.py).

Accuracy expectations are calibrated, not aspirational: the wall is a
flat-faceted approximation of the true circle, i.e. the same curved-wall
variational crime already root-caused on the G1.6 sphere (see
PROJECT_STRUCTURE.md "Known gaps"), so wall Cp converges at ~O(h):
measured max |Cp err| ~ 0.091 (coarse) -> 0.045 (medium) with quadratic
surface recovery. Thresholds below assert those measured levels with
margin, plus monotone improvement under refinement.

Spanwise-velocity finding (G2.5 evidence, documented in roadmap ledger):
the *interpolated* freestream phi = x has machine-zero spanwise gradient
(G2.5(a) passes), but the *solved* field does NOT -- measured
max |w|/U ~ 2.9e-2 (coarse) -> 1.5e-2 (medium), converging at ~O(h). This
is inherent to any 3-tet prism subdivision: the split is necessarily
asymmetric under the z-mirror (on a lateral quad, integral of a vertex hat
function is S/3 or S/6 depending on the diagonal), so a z-invariant field
is not a discrete solution and the discrete minimizer picks up O(h)
spanwise noise. G2.5(b)'s < 1e-12 criterion is therefore unachievable as
written with single-layer 3-tet prisms; the tests assert the honest
behavior (small and decreasing under refinement) instead.
"""

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.extrude import assert_quad_split_consistency

from .mesh_utils import (
    CYLINDER_RADIUS as A,
    element_gradients_all,
    run_cylinder_case,
)


class TestCylinderMeshIngestion:
    """M0 gate items: solver preprocessor ingests the mesh, tags are complete,
    the quad-split consistency assert passes."""

    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_tags_and_quad_split_consistency(self, mesh_dir, level):
        mesh = read_mesh(mesh_dir / "cylinder_2.5d" / f"{level}.msh")
        assert set(mesh.boundary_faces) == {"wall", "farfield", "symmetry"}
        assert_quad_split_consistency(mesh, interior_groups=())

        z = np.round(mesh.nodes[:, 2], 12)
        assert len(set(z)) == 2, "not a single-layer mesh"
        sym_nodes = np.unique(mesh.boundary_faces["symmetry"])
        for tri in mesh.boundary_faces["symmetry"]:
            assert len(set(z[tri])) == 1, "symmetry face not planar"
        assert len(sym_nodes) == len(mesh.nodes), (
            "single layer: every node lies on one of the two symmetry planes"
        )

    def test_volume_matches_annulus(self, mesh_dir):
        from pyfp3d.mesh.metrics import compute_tet_volumes

        mesh = read_mesh(mesh_dir / "cylinder_2.5d" / "coarse.msh")
        dz = np.ptp(mesh.nodes[:, 2])
        vols = compute_tet_volumes(mesh.nodes, mesh.elements)
        assert np.all(vols > 0)
        analytic = np.pi * (20.0**2 - A**2) * dz
        # polygonal circles undershoot the analytic area by O(h^2)
        assert 0.99 * analytic < vols.sum() < analytic


class TestCylinderFlow:
    """End-to-end: solve, spanwise behavior, Cp vs 1 - 4 sin^2(theta)."""

    def test_freestream_interpolant_spanwise_zero(self, mesh_dir):
        """G2.5(a): nodal phi = x has machine-zero spanwise gradient."""
        mesh = read_mesh(mesh_dir / "cylinder_2.5d" / "coarse.msh")
        grad = element_gradients_all(mesh.nodes, mesh.elements,
                                     mesh.nodes[:, 0].copy())
        assert np.max(np.abs(grad[:, 2])) < 1e-12
        assert np.max(np.abs(grad[:, 0] - 1.0)) < 1e-12

    def test_solved_flow_coarse(self, mesh_dir):
        case = run_cylinder_case(mesh_dir / "cylinder_2.5d" / "coarse.msh")
        assert case["residual_norm"] < 1e-8

        # solved-field spanwise noise: O(h), inherent to the 3-tet prism
        # split (see module docstring); assert the measured level w/ margin
        grad = element_gradients_all(case["mesh"].nodes,
                                     case["mesh"].elements, case["phi"])
        assert np.max(np.abs(grad[:, 2])) < 0.05  # measured 2.9e-2

        # Cp accuracy: measured 0.091 max / 0.045 mean (curved-wall crime)
        assert case["error"].max() < 0.12
        assert case["error"].mean() < 0.06

        # physics sanity: stagnation Cp ~ +1 at theta = 0, pi;
        # suction peak Cp ~ -3 at theta = +-pi/2
        near_stag = case["sin2"] < 0.02
        assert np.all(case["cp_numeric"][near_stag] > 0.8)
        near_peak = case["sin2"] > 0.98
        assert np.all(case["cp_numeric"][near_peak] < -2.5)

    def test_refinement_improves_medium(self, mesh_dir):
        coarse = run_cylinder_case(mesh_dir / "cylinder_2.5d" / "coarse.msh")
        medium = run_cylinder_case(mesh_dir / "cylinder_2.5d" / "medium.msh")
        assert medium["residual_norm"] < 1e-8

        # Cp error and spanwise noise both drop under refinement (~O(h))
        assert medium["error"].max() < 0.06  # measured 0.045
        assert medium["error"].max() < 0.75 * coarse["error"].max()

        gz_c = np.max(np.abs(element_gradients_all(
            coarse["mesh"].nodes, coarse["mesh"].elements, coarse["phi"])[:, 2]))
        gz_m = np.max(np.abs(element_gradients_all(
            medium["mesh"].nodes, medium["mesh"].elements, medium["phi"])[:, 2]))
        assert gz_m < 0.75 * gz_c

    def test_top_bottom_mirror_pairing(self, mesh_dir):
        """Top-layer solution tracks the bottom layer to the same O(h) level
        as the spanwise noise (they would be exactly equal only for a
        z-mirror-symmetric subdivision)."""
        case = run_cylinder_case(mesh_dir / "cylinder_2.5d" / "coarse.msh")
        n2 = len(case["mesh"].nodes) // 2
        # extruder layout: node i and i + n2 are the same 2D point
        assert np.allclose(case["mesh"].nodes[:n2, :2],
                           case["mesh"].nodes[n2:, :2])
        dphi = np.max(np.abs(case["phi"][n2:] - case["phi"][:n2]))
        assert dphi < 0.02  # measured 5.8e-3 on coarse


class TestCylinderArtifacts:
    """Headless evidence per roadmap Sec 0.1: Cp(theta) plot + CSV."""

    def test_export_m0_cylinder_cp(self, gate_artifacts_dir, mesh_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        cases = {
            level: run_cylinder_case(mesh_dir / "cylinder_2.5d" / f"{level}.msh")
            for level in ["coarse", "medium"]
        }

        theta_line = np.linspace(0, 180, 200)
        cp_line = 1.0 - 4.0 * np.sin(np.radians(theta_line)) ** 2

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(theta_line, cp_line, "k-", linewidth=2,
                label="exact: 1 - 4 sin^2(theta)")
        for level, case in cases.items():
            nodes = case["mesh"].nodes[case["wall_nodes"]]
            theta = np.degrees(np.arctan2(np.abs(nodes[:, 1]), nodes[:, 0]))
            order = np.argsort(theta)
            ax.plot(theta[order], case["cp_numeric"][order], ".",
                    markersize=3, alpha=0.5, label=f"numeric ({level})")
        ax.set_xlabel("theta (deg from +x stagnation point)")
        ax.set_ylabel("Cp")
        ax.set_title("M0 cylinder (quasi-2D single layer): surface Cp")
        ax.legend()
        ax.grid(True, alpha=0.3)
        out_png = gate_artifacts_dir / "cylinder_cp.png"
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert out_png.exists()

        with open(gate_artifacts_dir / "summary.csv", "w") as f:
            f.write("level,max_cp_error,mean_cp_error,n_cg_iterations,residual_norm\n")
            for level, case in cases.items():
                f.write(
                    f"{level},{case['error'].max():.6e},{case['error'].mean():.6e},"
                    f"{case['n_cg_iterations']},{case['residual_norm']:.3e}\n"
                )


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
