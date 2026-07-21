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
capability at h=0.08. Route fork RESOLVED 2026-07-22 (user-directed): (a)
Option C re-spec ADOPTED -- the active G1.6 gate is now the achievable,
measured criterion (all-scales-refined order >= 1.8 + mean Cp < 1% at
h_min 0.03; see class TestG16Respec below, which PASSES on P11's committed
sweep). The literal 2%-max-at-medium xfail below STAYS as the recorded P1
limitation (it is beyond any P1-field method; the isoparametric-P2 route (b)
was not taken). See PROJECT_STRUCTURE 'Known gaps' + roadmap track_p G1.6.

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
            "Recorded P1 limitation (NOT the active G1.6 gate since the 2026-07-22 "
            "Option C re-spec -- see TestG16Respec, which PASSES). The literal "
            "'max Cp < 2% at the medium mesh' is ~11.6% (quadratic recovery); P11 "
            "(2026-07-19) proved it is the intrinsic P1-field max-norm capability "
            "at h=0.08 (geometric share ~0.2 pp), beyond any P1-field method (the "
            "literal 2% needs an isoparametric P2 wall layer, route (b), not taken)"
            " -- see module docstring + PROJECT_STRUCTURE.md 'Known gaps'."
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


def _read_p11_csv(name):
    """Load a committed P11 sweep CSV (rows as dicts). The .msh are gitignored
    but these CSVs ARE the committed evidence (CLAUDE.md rule 3)."""
    import csv as _csv
    from pathlib import Path

    path = (Path(__file__).resolve().parent.parent
            / "cases" / "demo" / "p11_curved_walls" / "results" / name)
    if not path.exists():
        pytest.skip(f"P11 evidence {name} absent; regenerate via "
                    f"cases/demo/p11_curved_walls/run_demo.py")
    with open(path) as f:
        return list(_csv.DictReader(f))


class TestG16Respec:
    """G1.6 gate RE-SPEC (Option C, P11 route fork adopted 2026-07-22,
    user-directed).

    The literal '2% max Cp at the medium mesh' criterion is beyond ANY
    P1-field method (P11: it needs O(h^2) wall velocity at h=0.08) -- that
    remains recorded as a P1 limitation by the strict xfail above. The
    ACTIVE, ACHIEVABLE G1.6 gate this class defines, measured PASSING in
    P11's committed sweep, is: on an ALL-SCALES-refined family (h_min AND
    h_far scaling together, not the single-variable h_min sweep whose order
    collapse is the fixed-bulk-mesh pollution floor) the wall recovery
    converges at order >= 1.8 and reaches mean Cp < 1% at h_min = 0.03.

    Evidence (committed, no re-solve): cases/demo/p11_curved_walls/results/
    {e8_bulk_floor.csv, e6_ico_control.csv}.
    """

    def test_allscales_mean_cp_below_one_percent(self):
        """Option C criterion 1: all-scales-refined mean Cp < 1% at h_min 0.03
        (measured 0.60% -- the E8 h03_far10 mesh, h_min 0.03 AND h_max 1.0)."""
        rows = _read_p11_csv("e8_bulk_floor.csv")
        far = next(r for r in rows if r["mesh"].startswith("h03_far10"))
        cp_mean = float(far["cp_mean"])
        assert cp_mean < 0.01, (
            f"all-scales mean Cp {cp_mean:.4f} >= 1% (re-spec criterion)")
        # corroboration: refining ONLY the far mesh at fixed h_min=0.03 drops
        # the wall-phi error >= 2x -- the single-variable sweep's order
        # collapse was bulk-mesh pollution, not a wall/geometry error.
        h03 = next(r for r in rows if r["mesh"].startswith("h03 "))
        drop = float(h03["phi_err_wall_max"]) / float(far["phi_err_wall_max"])
        assert drop >= 2.0, f"far-only refinement drop {drop:.2f}x < 2x"

    def test_allscales_order_at_least_1p8(self):
        """Option C criterion 2: the all-scales-refined recovery order -> >= 1.8
        (asymptotic pair). Measured on the structured icosphere control with the
        SAME flat facets (s4->s5, clean 2x ratio): 1.98."""
        rows = _read_p11_csv("e6_ico_control.csv")
        by_level = {r["level"]: r for r in rows}
        s4, s5 = by_level["s4"], by_level["s5"]
        h_ratio = float(s4["h_wall_median"]) / float(s5["h_wall_median"])
        err_ratio = (float(s4["phi_err_wall_max"])
                     / float(s5["phi_err_wall_max"]))
        order = np.log(err_ratio) / np.log(h_ratio)
        assert order >= 1.8, (
            f"asymptotic wall-phi order {order:.2f} < 1.8 (re-spec criterion)")

    def test_intrinsic_p1_capability_recorded(self):
        """Recorded: the structured control reaches ~2% max Cp at h~0.036
        (s5) -- the intrinsic P1-field max-norm capability the literal 2%
        criterion runs into. This is a floor, not a bug (documents WHY the
        literal gate is xfail)."""
        rows = _read_p11_csv("e6_ico_control.csv")
        s5 = next(r for r in rows if r["level"] == "s5")
        assert float(s5["cp_max"]) < 0.03  # ~2.14% at h~0.036


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
