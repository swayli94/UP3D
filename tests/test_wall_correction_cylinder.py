"""
Gate G1.3 (formerly G1.2-a0): cylinder oracle pre-study for the Option A
true-normal weak-flux correction (design.md §5.1/§5.1.1, roadmap G1.3).

OUTCOME (2026-07-06) -- the pre-study is complete and its acceptance
criterion is NOT met, for a reason the fix plan did not anticipate: on a
body-fitted mesh (wall vertices on the true surface) there is essentially
nothing for a boundary-DATA correction to correct.

1. The exact solution's net flux through EACH flat wall facet is exactly
   zero: the sliver between a chord facet and the true arc is a closed
   region bounded only by the facet and the (zero-flux) true surface, so
   the divergence theorem forces the facet flux to vanish for a harmonic
   potential. The weak consistency term <N_i, grad(phi_exact).n_facet> is
   therefore only a first-moment residue -- measured max |b_i| ~ 2e-5
   (coarse) scaling ~O(h^4).
2. Option A's t-decomposition form <N_i, grad(phi_exact).t> assembles to
   machine zero (~1e-17) on this mesh: contributions from the facets
   sharing each boundary edge cancel exactly on the uniformly-spaced
   circle (the assembly itself is verified correct on a hand-computed
   single-facet case below).
3. Consequently the oracle-corrected solve is bit-close to the uncorrected
   one and the Cp(theta) error is unchanged -- no "significant reduction",
   no order recovery. Uncorrected Cp error already converges at ~1.0 order.
4. Error anatomy overturns the cylinder's "same variational crime as the
   sphere" designation: feeding the EXACT potential into the surface
   recovery reproduces ~3/4 of the total Cp error at every level (recovery
   on the quasi-2D single-layer sliver strip is the dominant term), and
   the wall nodal phi error converges at ~1.2 order (healthy, unlike the
   sphere's decreasing sub-first order).

DP1 decision (recorded in roadmap.md): G1.3 fails -> ">5%" branch. A
direct sphere-oracle run (G1.4, absorbed into cases/demo/p1_laplace/run_demo.py)
confirms it on the gate geometry: full-flux oracle ceiling ~11.3% vs 11.6%
uncorrected on the medium sphere mesh.

The tests below lock in the measured facts as regressions; the gate's
acceptance criterion itself is a strict xfail, mirroring the G1.6 idiom.
"""

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.wall_correction import (
    assemble_wall_flux_correction_rhs,
    cylinder_closest_point_normal,
    facet_normal_deviation_angles,
    wall_correction_geometry,
    wall_face_adjacent_tets,
)

from .mesh_utils import (
    cylinder_grad_exact,
    cylinder_phi_exact,
    element_gradients_all,
    run_cylinder_case,
)

LEVELS = {"coarse": 0.10, "medium": 0.05, "fine": 0.025}  # h_wall per level


@pytest.fixture(scope="module")
def cylinder_oracle_cases():
    """Solve every level once, uncorrected and oracle-corrected."""
    from pathlib import Path

    mesh_dir = Path(__file__).parent.parent / "cases" / "meshes"
    cases = {}
    for level in LEVELS:
        path = mesh_dir / "cylinder_2.5d" / f"{level}.msh"
        mesh = read_mesh(path)
        wall_faces = mesh.boundary_faces["wall"]

        geo = wall_correction_geometry(
            mesh.nodes, mesh.elements, wall_faces, cylinder_closest_point_normal
        )
        grad_qp = cylinder_grad_exact(geo["qp"].reshape(-1, 3)).reshape(geo["qp"].shape)
        rhs = assemble_wall_flux_correction_rhs(
            len(mesh.nodes), wall_faces, geo, grad_qp
        )

        uncorrected = run_cylinder_case(path, mesh=mesh)
        corrected = run_cylinder_case(path, body_source_rhs=rhs, mesh=mesh)

        def spanwise_floor(phi, m=mesh):
            return np.max(np.abs(element_gradients_all(m.nodes, m.elements, phi)[:, 2]))

        cases[level] = {
            "mesh": mesh,
            "geometry": geo,
            "rhs": rhs,
            "uncorrected": uncorrected,
            "corrected": corrected,
            "w_uncorrected": spanwise_floor(uncorrected["phi"]),
            "w_corrected": spanwise_floor(corrected["phi"]),
        }
    return cases


class TestCorrectionAssembly:
    """The zero cylinder result is real, not an assembly bug: the assembly
    reproduces a hand-computed value on a fabricated single-facet case."""

    def _single_tet(self):
        nodes = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        elements = np.array([[0, 1, 2, 3]])
        wall_faces = np.array([[0, 1, 2]])
        return nodes, elements, wall_faces

    def test_hand_computed_single_facet(self):
        nodes, elements, wall_faces = self._single_tet()
        alpha = 0.1  # tilt the "true" normal by a known angle in the xz-plane
        n_true = np.array([np.sin(alpha), 0.0, -np.cos(alpha)])

        def cpn(pts):
            return pts.copy(), np.tile(n_true, (len(pts), 1))

        geo = wall_correction_geometry(nodes, elements, wall_faces, cpn)
        # facet normal must be domain-outward (-z, away from the tet)
        assert np.allclose(geo["n_facet"], [0.0, 0.0, -1.0])

        g = np.array([2.0, 0.0, 0.5])
        rhs = assemble_wall_flux_correction_rhs(
            4, wall_faces, geo, np.tile(g, (1, 3, 1)).reshape(1, 3, 3)
        )
        # constant integrand: b_i = (g . t) * area / 3 at each facet node
        t = geo["n_facet"][0] - np.dot(geo["n_facet"][0], n_true) * n_true
        expected = np.dot(g, t) * (0.5 / 3.0)
        assert np.allclose(rhs[:3], expected)
        assert rhs[3] == 0.0

    def test_flat_wall_gives_exactly_zero(self):
        """If the true surface IS the facet plane, t = 0 and the correction
        vanishes identically (the V0/no-curvature no-op property)."""
        nodes, elements, wall_faces = self._single_tet()

        def cpn(pts):
            n = np.tile(np.array([0.0, 0.0, -1.0]), (len(pts), 1))
            return pts.copy(), n

        geo = wall_correction_geometry(nodes, elements, wall_faces, cpn)
        rhs = assemble_wall_flux_correction_rhs(
            4, wall_faces, geo,
            np.random.default_rng(0).normal(size=(1, 3, 3)),
        )
        assert np.all(rhs == 0.0)

    def test_adjacent_tet_lookup_rejects_interior_face(self):
        nodes = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], float
        )
        elements = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
        with pytest.raises(ValueError, match="owned by two tets"):
            wall_face_adjacent_tets(elements, np.array([[1, 2, 3]]))


class TestCylinderOracle:
    """The measured pre-study facts, locked in as regressions."""

    def test_correction_rhs_is_machine_zero(self, cylinder_oracle_cases):
        """Fact 2 of the module docstring: adjacent-facet cancellation on the
        uniformly-spaced circle assembles the t-form RHS to machine zero,
        even though single quadrature points carry O(h) flux defects."""
        for level, case in cylinder_oracle_cases.items():
            assert np.abs(case["rhs"]).max() < 1e-12, level
            flux_defect = np.einsum(
                "fqk,fqk->fq",
                cylinder_grad_exact(
                    case["geometry"]["qp"].reshape(-1, 3)
                ).reshape(case["geometry"]["qp"].shape),
                case["geometry"]["t"],
            )
            assert np.abs(flux_defect).max() > 1e-3, (
                f"{level}: pointwise defects should be O(h), "
                "otherwise the zero RHS would be trivial"
            )

    def test_correction_leaves_solution_unchanged(self, cylinder_oracle_cases):
        """Fact 3: the oracle correction has no effect on the solution."""
        for level, case in cylinder_oracle_cases.items():
            dphi = np.abs(case["corrected"]["phi"] - case["uncorrected"]["phi"]).max()
            assert dphi < 1e-8, level
            assert case["corrected"]["residual_norm"] < 1e-8, level

    def test_uncorrected_error_first_order(self, cylinder_oracle_cases):
        """The uncorrected Cp error already converges at ~first order
        (measured slope 1.02) -- there is no sub-first-order pathology to
        recover from on this geometry (unlike the G1.6 sphere)."""
        h = np.array([LEVELS[lv] for lv in LEVELS])
        e = np.array(
            [cylinder_oracle_cases[lv]["uncorrected"]["error"].max() for lv in LEVELS]
        )
        slope = np.polyfit(np.log(h), np.log(e), 1)[0]
        assert 0.85 < slope < 1.25, slope
        # measured levels with margin: 0.091 / 0.045 / 0.022
        assert e[0] < 0.12 and e[1] < 0.06 and e[2] < 0.03

    def test_recovery_dominates_cylinder_error(self, cylinder_oracle_cases):
        """Fact 4: the exact potential fed through the surface recovery
        reproduces most of the total Cp error (measured ~76% at every
        level) -- the quasi-2D sliver-strip recovery, not the volume solve,
        dominates here. This de-designates the cylinder as a testbed for
        the sphere's variational-crime pathology."""
        from pyfp3d.post.surface import wall_tangential_gradient_quadratic

        for level, case in cylinder_oracle_cases.items():
            mesh = case["mesh"]
            wall_nodes = case["uncorrected"]["wall_nodes"]
            grad = wall_tangential_gradient_quadratic(
                mesh.nodes, mesh.boundary_faces["wall"], cylinder_phi_exact(mesh.nodes)
            )
            cp_oracle = 1.0 - np.sum(grad[wall_nodes] ** 2, axis=1)
            oracle_err = np.abs(cp_oracle - case["uncorrected"]["cp_exact"]).max()
            full_err = case["uncorrected"]["error"].max()
            assert oracle_err > 0.6 * full_err, (level, oracle_err, full_err)

    def test_spanwise_noise_floor(self, cylinder_oracle_cases):
        """G2.5(b) magnitude reference (by-product required by the gate
        spec): the corrected solve has the same O(h) spanwise-noise floor,
        decreasing under refinement."""
        w = {lv: cylinder_oracle_cases[lv]["w_corrected"] for lv in LEVELS}
        for lv in LEVELS:
            assert np.isclose(
                w[lv], cylinder_oracle_cases[lv]["w_uncorrected"], rtol=1e-6
            )
        # measured: 2.9e-2 / 1.5e-2 / 7.4e-3
        assert w["coarse"] < 0.05
        assert w["medium"] < 0.75 * w["coarse"]
        assert w["fine"] < 0.75 * w["medium"]

    @pytest.mark.xfail(
        strict=True,
        reason="G1.3 acceptance NOT met: the oracle correction does not "
        "reduce the Cp error (nothing to correct on a body-fitted mesh; "
        "see module docstring). DP1 '>5%' branch taken.",
    )
    def test_gate_acceptance_criterion(self, cylinder_oracle_cases):
        """The gate's acceptance as specified: corrected error significantly
        below uncorrected."""
        medium = cylinder_oracle_cases["medium"]
        assert (
            medium["corrected"]["error"].max()
            < 0.5 * medium["uncorrected"]["error"].max()
        )


class TestG13Artifacts:
    """Headless artifacts per the G1.3 gate spec (roadmap §0.1)."""

    def test_export_g13_artifacts(self, cylinder_oracle_cases, artifacts_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        out = artifacts_dir / "G1.3"
        out.mkdir(parents=True, exist_ok=True)

        # --- cp_theta_overlay.png: exact / uncorrected / oracle-corrected ---
        fig, ax = plt.subplots(figsize=(8, 6))
        theta_line = np.linspace(0, 180, 200)
        ax.plot(theta_line, 1.0 - 4.0 * np.sin(np.radians(theta_line)) ** 2,
                "k-", linewidth=2, label="exact: 1 - 4 sin^2(theta)")
        for level in ("coarse", "medium"):
            case = cylinder_oracle_cases[level]
            nodes = case["mesh"].nodes[case["uncorrected"]["wall_nodes"]]
            theta = np.degrees(np.arctan2(np.abs(nodes[:, 1]), nodes[:, 0]))
            order = np.argsort(theta)
            ax.plot(theta[order], case["uncorrected"]["cp_numeric"][order],
                    ".", markersize=3, alpha=0.5, label=f"uncorrected ({level})")
            ax.plot(theta[order], case["corrected"]["cp_numeric"][order],
                    "x", markersize=3, alpha=0.5,
                    label=f"oracle-corrected ({level})")
        ax.set_xlabel("theta (deg from +x stagnation point)")
        ax.set_ylabel("Cp")
        ax.set_title("G1.3: cylinder Cp -- correction has no effect "
                     "(curves coincide; that IS the finding)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.savefig(out / "cp_theta_overlay.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # --- error_vs_h.png ---
        h = np.array([LEVELS[lv] for lv in LEVELS])
        e_un = np.array([cylinder_oracle_cases[lv]["uncorrected"]["error"].max()
                         for lv in LEVELS])
        e_co = np.array([cylinder_oracle_cases[lv]["corrected"]["error"].max()
                         for lv in LEVELS])
        s_un = np.polyfit(np.log(h), np.log(e_un), 1)[0]
        s_co = np.polyfit(np.log(h), np.log(e_co), 1)[0]
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.loglog(h, e_un, "o-", label=f"uncorrected (slope {s_un:.2f})")
        ax.loglog(h, e_co, "x--", label=f"oracle-corrected (slope {s_co:.2f})")
        ax.set_xlabel("h_wall")
        ax.set_ylabel("max |Cp error|")
        ax.set_title("G1.3: Cp error vs h (lines coincide)")
        ax.legend()
        ax.grid(True, which="both", alpha=0.3)
        fig.savefig(out / "error_vs_h.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # --- normal_deviation.png: geometric-error self-check ---
        fig, ax = plt.subplots(figsize=(8, 6))
        for level in ("coarse", "medium"):
            geo = cylinder_oracle_cases[level]["geometry"]
            qp = geo["qp"].reshape(-1, 3)
            theta = np.degrees(np.arctan2(np.abs(qp[:, 1]), qp[:, 0]))
            dev = np.degrees(facet_normal_deviation_angles(geo)).ravel()
            ax.plot(theta, dev, ".", markersize=3, alpha=0.5, label=level)
        ax.set_xlabel("theta (deg)")
        ax.set_ylabel("angle(n_facet, n_true) (deg)")
        ax.set_title("G1.3: facet-normal deviation (O(h) as expected)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.savefig(out / "normal_deviation.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # --- section_symmetry_plane.png: section_cut interface pre-warm ---
        from pyfp3d.post.section_cut import plot_section_field, section_cut

        case = cylinder_oracle_cases["coarse"]
        mesh = case["mesh"]
        w_nodal = np.zeros(len(mesh.nodes))
        grad_z = element_gradients_all(mesh.nodes, mesh.elements,
                                       case["corrected"]["phi"])[:, 2]
        counts = np.zeros(len(mesh.nodes))
        for k in range(4):
            np.add.at(w_nodal, mesh.elements[:, k], grad_z)
            np.add.at(counts, mesh.elements[:, k], 1.0)
        w_nodal /= np.maximum(counts, 1.0)
        z0 = float(mesh.nodes[:, 2].min())
        section = section_cut(mesh, {"w": w_nodal}, z=z0)
        plot_section_field(
            section, "w", out / "section_symmetry_plane.png",
            title="G1.3: spanwise velocity w on z = z0 (O(h) prism-split noise)",
            cmap="RdBu_r", xlim=(-4, 4), ylim=(-4, 4),
        )

        # --- summary.csv ---
        with open(out / "summary.csv", "w") as f:
            f.write("level,h_wall,max_cp_err_uncorrected,mean_cp_err_uncorrected,"
                    "max_cp_err_corrected,mean_cp_err_corrected,"
                    "max_rhs,spanwise_floor_uncorrected,spanwise_floor_corrected,"
                    "n_cg_iterations,residual_norm\n")
            for lv in LEVELS:
                c = cylinder_oracle_cases[lv]
                f.write(
                    f"{lv},{LEVELS[lv]},"
                    f"{c['uncorrected']['error'].max():.6e},"
                    f"{c['uncorrected']['error'].mean():.6e},"
                    f"{c['corrected']['error'].max():.6e},"
                    f"{c['corrected']['error'].mean():.6e},"
                    f"{np.abs(c['rhs']).max():.3e},"
                    f"{c['w_uncorrected']:.6e},{c['w_corrected']:.6e},"
                    f"{c['corrected']['n_cg_iterations']},"
                    f"{c['corrected']['residual_norm']:.3e}\n"
                )
            f.write(f"slope_uncorrected,,{s_un:.3f},,,,,,,,\n")
            f.write(f"slope_corrected,,{s_co:.3f},,,,,,,,\n")

        for name in ("cp_theta_overlay.png", "error_vs_h.png",
                     "normal_deviation.png", "section_symmetry_plane.png",
                     "summary.csv"):
            assert (out / name).exists()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
