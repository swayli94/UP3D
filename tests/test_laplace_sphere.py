"""
Gate G1.6 (formerly G1.2): Incompressible sphere Cp vs 1 - 9/4 sin^2(theta).

Setup: uniform flow past a unit sphere (cases/meshes/sphere_shell), Dirichlet
BC from the exact potential phi = x(1 + a^3/(2r^3)) at the far field, natural
(zero-flux) BC at the wall. Surface Cp = 1 - q^2 is evaluated from the wall's
own tangential gradient (post/surface.py::wall_tangential_gradient_quadratic)
rather than by extrapolating the volume-element gradient -- see that
function's docstring for why the volume-based approach was rejected
(systematic underestimate at the boundary, and an ill-conditioned failure
mode for a naive least-squares fix).

Status (see docs/roadmap.md P1 ledger and PROJECT_STRUCTURE.md "Known gaps"):
this does NOT yet meet the <2% max-error target from design.md Sec 10 (V2)
on the medium mesh. Measured max error is ~11.6% on medium (down from ~12.0%
with the older linear-only recovery) and only reaches ~4% even at h_min=0.02
(~540k nodes) on a clean, single-variable h_min refinement sweep -- so the
remaining gap is not simply "refine the mesh more," and is left open rather
than papered over with a loosened threshold.

Root cause -- ERRATUM 2026-07-19 (phase P11; see PROJECT_STRUCTURE.md "Known
gaps" P11 block and roadmap track_p.md Sec P11 for the measurements): the
paragraph below was this file's original attribution; its "variational
crime dominates" mechanism is OVERTURNED by measurement. What P11 found:
a verified curved wall-adjacent layer (solve/curved_wall.py) moves the
medium error only 11.56% -> 11.33% (= the G1.4 boundary-data oracle
ceiling -- the geometric crime is worth ~0.2 pp); the sub-first-order
"signature" (orders 0.88/0.56/0.42) is the fixed-bulk-mesh pollution floor
of a sweep that shrinks h_min while h_max stays 3.0 (refining ONLY the far
mesh at h_min=0.03 drops the wall phi error 3.17x and restores order 1.89);
and a structured icosphere shell with the SAME flat facets converges at
~2nd order. The 11.6% is essentially the intrinsic P1-field max-norm
capability at h=0.08. The gate stays open; the recorded route fork (Option
C re-spec / isoparametric P2 wall layer / accept) is the user's call.

Original (superseded) attribution, kept for the record: NOT the surface
recovery scheme -- an oracle test that feeds
the *exact* analytic potential through the recovery step with no FEM solve at
all shows the recovery operator's own bias is a small fraction (well under a
tenth) of the total error at every mesh level tested, for both the linear and
quadratic recovery schemes [STANDS]. NOT under-resolved bulk mesh either --
tightening
the far-field mesh at fixed h_min helps a little then plateaus [SUPERSEDED:
that measurement stopped at ~4% max Cp because the MAX-norm is resolution-
limited at the wall; the wall-phi floor it left behind IS the bulk mesh,
per P11/E8]. The dominant
error source is the volume PDE solve's own accuracy next to the wall [STANDS]:
the
natural (zero-flux) BC is satisfied on the flat polyhedral wall-facet
approximation, not the true curved sphere, a geometric/variational-crime
inconsistency that pollutes the whole domain through ellipticity and shows up
as sub-first-order convergence (order ~0.4-0.9, decreasing as h shrinks) of
the raw nodal potential itself, not just its recovered gradient [OVERTURNED,
see the erratum above]. A first
attempt at a direct fix (a Nitsche/penalty term weakly forcing each
wall-adjacent tet's own volumetric gradient toward zero along the true
surface normal) was tried and rejected: it made the error *worse* with
increasing penalty strength, because a P1 tet spanning from the wall inward
necessarily has a nonzero radial gradient component representing the
interior falloff of tangential velocity -- that's correct FEM behavior, not a
BC violation, so penalizing it fights the physically-correct solution
[STANDS]. Closing
this gate for real looks like it needs genuine curved/isoparametric boundary
elements [EXECUTED by P11 2026-07-19: measured NEGATIVE -- superparametric
mapped-P1 curving swaps one O(h) error for another; an isoparametric P2
layer (field order raised) is the only remaining route to the literal
criterion].

The xfail below tracks the *actual* gate criterion so it turns into a hard
failure (strict=True) the day someone fixes the underlying accuracy issue and
forgets to remove the marker.
"""

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.surface import wall_tangential_gradient_quadratic
from pyfp3d.solve.picard import solve_laplace

A = 1.0  # sphere radius


def run_sphere_case(mesh_path):
    mesh = read_mesh(mesh_path)
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])

    r = np.linalg.norm(nodes, axis=1)
    phi_exact = nodes[:, 0] * (1.0 + 0.5 * A**3 / r**3)

    result = solve_laplace(
        nodes, elements, farfield_nodes, phi_exact[farfield_nodes], rtol=1e-11, maxiter=3000,
    )
    phi = result["phi"]

    grad_wall = wall_tangential_gradient_quadratic(nodes, wall_faces, phi)
    q_squared = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
    cp_numeric = 1.0 - q_squared

    cos_theta = nodes[wall_nodes, 0] / r[wall_nodes]
    cp_exact = 1.0 - 2.25 * (1.0 - cos_theta**2)

    return {
        "wall_nodes": wall_nodes,
        "cos_theta": cos_theta,
        "cp_numeric": cp_numeric,
        "cp_exact": cp_exact,
        "error": np.abs(cp_numeric - cp_exact),
        "n_cg_iterations": result["n_cg_iterations"],
        "residual_norm": result["residual_norm"],
    }


class TestSphereCp:
    """Gate G1.6: incompressible sphere Cp, max error < 2% on medium mesh."""

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Known open gap: max Cp error on the medium mesh is ~11.6% (quadratic "
            "surface recovery), not <2%. P11 (2026-07-19) re-attribution: this is "
            "the intrinsic P1-field max-norm capability at h=0.08 (curved wall "
            "elements measured a ~0.2 pp geometric share; the 2% bar needs an "
            "isoparametric P2 wall layer or the Option C gate re-spec) -- see "
            "module docstring erratum and PROJECT_STRUCTURE.md 'Known gaps'."
        ),
    )
    def test_sphere_cp_medium_mesh(self, mesh_dir):
        case = run_sphere_case(mesh_dir / "sphere_shell" / "medium.msh")
        assert case["residual_norm"] < 1e-6

        max_error = case["error"].max()
        assert max_error < 0.02, f"max |Cp error| = {max_error:.4f} >= 2% target"

    def test_sphere_stagnation_and_symmetry_sanity(self, mesh_dir):
        """Coarser sanity checks that should hold regardless of the open
        accuracy gap: front/rear stagnation Cp ~ 1, and Cp is symmetric
        about the equator (fore-aft symmetry of potential flow past a sphere)."""
        case = run_sphere_case(mesh_dir / "sphere_shell" / "coarse.msh")

        near_stagnation = np.abs(case["cos_theta"]) > 0.98
        assert np.all(case["cp_numeric"][near_stagnation] > 0.8), (
            "Cp near the fore/aft stagnation points should be close to +1"
        )

        near_equator = np.abs(case["cos_theta"]) < 0.05
        assert np.all(case["cp_numeric"][near_equator] < -1.0), (
            "Cp near the equator should be close to the -1.25 suction peak"
        )


class TestSphereCpArtifacts:
    """Generate visual artifacts for G1.6."""

    def test_export_sphere_cp_meridian(self, gate_artifacts_dir, mesh_dir):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        case = run_sphere_case(mesh_dir / "sphere_shell" / "medium.msh")
        theta_deg = np.degrees(np.arccos(np.clip(case["cos_theta"], -1.0, 1.0)))
        order = np.argsort(theta_deg)

        theta_line = np.linspace(0, 180, 200)
        cp_exact_line = 1.0 - 2.25 * np.sin(np.radians(theta_line)) ** 2

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(theta_line, cp_exact_line, "k-", linewidth=2, label="exact: 1 - 9/4 sin^2(theta)")
        ax.plot(theta_deg[order], case["cp_numeric"][order], ".", markersize=3, alpha=0.5,
                label="numeric (medium mesh)")
        ax.set_xlabel("theta (deg from +x stagnation point)")
        ax.set_ylabel("Cp")
        ax.set_title("G1.6: incompressible sphere surface Cp")
        ax.legend()
        ax.grid(True, alpha=0.3)

        output_file = gate_artifacts_dir / "sphere_cp_meridian.png"
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)
        assert output_file.exists()

        csv_file = gate_artifacts_dir / "summary.csv"
        with open(csv_file, "w") as f:
            f.write("metric,value\n")
            f.write(f"max_error,{case['error'].max():.6e}\n")
            f.write(f"mean_error,{case['error'].mean():.6e}\n")
            f.write(f"n_wall_nodes,{len(case['wall_nodes'])}\n")
            f.write(f"gate_target,0.02\n")
            f.write(f"gate_status,OPEN (see test_sphere_cp_medium_mesh xfail)\n")
        assert csv_file.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
