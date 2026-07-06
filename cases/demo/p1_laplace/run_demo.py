"""
P1 demo -- Laplace solver evidence (G1.1 + G1.2 closed; G1.3/G1.4 oracle
negative results; G1.6 OPEN, shown honestly as the documented gap).

What this shows, per docs/roadmap.md P1 and docs/demo_report.md:
  1. V0 freestream preservation: phi = x gives machine-zero interior
     residual on every mesh type -- the assembly is consistent.
  2. G1.1 MMS convergence: L2 error is O(h^2) (slope >= 1.9) -- the
     discretization has its design order.
  3. G1.2 CG+AMG mesh independence: iteration count nearly flat over a
     64x node-count increase -- the linear solver scales.
  4. Physical plausibility: the 3D sphere flow field (speed + Cp slice)
     has the right stagnation points, suction band, and fore-aft symmetry.
  5. G1.6 (OPEN, strict xfail in tests): wall Cp vs the analytic
     1 - (9/4)sin^2(theta); the ~11-12% max error is the root-caused
     curved-wall variational crime, NOT an unknown bug.
  6. G1.4 oracle negative result (evidence for DP1): even with the EXACT
     analytic gradient fed into the Option A boundary-data correction,
     the Cp error barely moves -- boundary-data corrections are ruled
     out; the sanctioned route is Option C + curved elements.

Standalone + self-checking:  python cases/demo/p1_laplace/run_demo.py
Outputs: cases/demo/p1_laplace/results/{*.png, summary.csv, checks.csv}
Exit code 0 iff every check passes (the G1.6 gate check is XFAIL by
design and does not fail the run). Runtime ~2-3 min.
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, DIV_BLUE_RED, GRID, INK, INK_2, MESH_DIR, MUTED,
    REPO_ROOT, S1_BLUE, S2_AQUA, S3_YELLOW, S5_VIOLET, SEQ_BLUE, SURFACE,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402
import pyamg  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from pyfp3d.kernels.residual import assemble_residual, assemble_stiffness_matrix  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.post.surface import (  # noqa: E402
    nodal_gradient_recovery, wall_tangential_gradient_quadratic,
)
from pyfp3d.solve.linear import apply_dirichlet  # noqa: E402
from pyfp3d.solve.picard import solve_laplace  # noqa: E402
from pyfp3d.solve.wall_correction import (  # noqa: E402
    assemble_wall_flux_correction_rhs, sphere_closest_point_normal,
    wall_correction_geometry,
)
from tests.mesh_utils import (  # noqa: E402
    cube_boundary_mask, generate_structured_cube_mesh, generate_sphere_shell_mesh,
)
from tests.test_laplace_mms import phi_exact_fn, run_mms_case  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
SPHERE_DIR = MESH_DIR / "sphere_shell"

# 7-point degree-5 triangle quadrature for the full-flux oracle variant
_A1, _B1 = 0.0597158717, 0.4701420641
_A2, _B2 = 0.7974269853, 0.1012865073
BARY7 = np.array([
    [1 / 3, 1 / 3, 1 / 3],
    [_A1, _B1, _B1], [_B1, _A1, _B1], [_B1, _B1, _A1],
    [_A2, _B2, _B2], [_B2, _A2, _B2], [_B2, _B2, _A2],
])
W7 = np.array([0.225] + [0.1323941527] * 3 + [0.1259391805] * 3)


def sphere_phi_exact(nodes, a=1.0):
    r = np.linalg.norm(nodes, axis=1)
    return nodes[:, 0] * (1.0 + 0.5 * a**3 / r**3)


def sphere_grad_exact(p, a=1.0):
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    r = np.linalg.norm(p, axis=1)
    c = 0.5 * a**3
    g = np.zeros_like(p)
    g[:, 0] = 1.0 + c / r**3 - 3 * c * x * x / r**5
    g[:, 1] = -3 * c * x * y / r**5
    g[:, 2] = -3 * c * x * z / r**5
    return g


def sphere_cp_exact(nodes, wall_nodes):
    r = np.linalg.norm(nodes[wall_nodes], axis=1)
    cos_t = nodes[wall_nodes, 0] / r
    return 1.0 - 2.25 * (1.0 - cos_t**2), cos_t


def solve_sphere(mesh, body_source_rhs=None):
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])
    phi_ex = sphere_phi_exact(mesh.nodes)
    return solve_laplace(mesh.nodes, mesh.elements, farfield_nodes,
                         phi_ex[farfield_nodes], body_source_rhs=body_source_rhs,
                         rtol=1e-11, maxiter=3000)


def sphere_cp_numeric(mesh, phi, wall_faces, wall_nodes):
    grad = wall_tangential_gradient_quadratic(mesh.nodes, wall_faces, phi)
    return 1.0 - np.sum(grad[wall_nodes] ** 2, axis=1)


# ---------------------------------------------------------------------------
# 1. V0 freestream preservation
# ---------------------------------------------------------------------------
def demo_freestream(checks):
    print("[1/6] V0 freestream preservation")
    cases = []

    nodes, elements = generate_structured_cube_mesh(n=8)
    interior = ~cube_boundary_mask(nodes)
    R = assemble_residual(nodes, elements, nodes[:, 0].copy())
    cases.append((f"structured cube (Kuhn, {len(elements)/1e3:.1f}k tets)",
                  np.abs(R[interior]).max()))

    nodes, elements, wall, far = generate_sphere_shell_mesh(subdivisions=2, n_layers=12)
    interior = np.ones(len(nodes), dtype=bool)
    interior[wall] = False
    interior[far] = False
    R = assemble_residual(nodes, elements, nodes[:, 0].copy())
    cases.append((f"icosphere shell ({len(elements)/1e3:.1f}k tets)",
                  np.abs(R[interior]).max()))

    # gmsh meshes: exclude only wall/farfield rows (they carry the physical
    # through-boundary flux of phi = x). Symmetry-plane rows stay in the
    # check: their normal is +-z and grad(x).z = 0, so they must be machine
    # zero too -- on the quasi-2D NACA mesh every node is on some boundary,
    # making this the stronger statement.
    for family, level in (("sphere_shell", "medium"), ("naca0012_2.5d", "medium")):
        mesh = read_mesh(MESH_DIR / family / f"{level}.msh")
        flux_nodes = np.unique(np.concatenate(
            [mesh.boundary_faces[tag].ravel() for tag in ("wall", "farfield")]))
        interior = np.ones(len(mesh.nodes), dtype=bool)
        interior[flux_nodes] = False
        R = assemble_residual(mesh.nodes, mesh.elements, mesh.nodes[:, 0].copy())
        cases.append((f"gmsh {family} {level} ({len(mesh.elements)/1e3:.0f}k tets)",
                      np.abs(R[interior]).max()))

    vals = np.array([c[1] for c in cases])
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    y = np.arange(len(cases))[::-1]
    x_lo = vals.min() * 0.35
    ax.hlines(y, x_lo, vals, color=GRID, linewidth=1.0, zorder=2)
    ax.plot(vals, y, "o", color=S1_BLUE, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
    ax.axvline(1e-12, color=CRITICAL, linewidth=1.2, linestyle="--", zorder=2)
    ax.set_ylim(-0.75, len(cases) - 0.4)
    ax.text(1.4e-12, -0.55, "gate: 1e-12", color=CRITICAL, fontsize=9, va="center")
    for yi, v in zip(y, vals):
        ax.text(v * 2.2, yi, f"{v:.1e}", va="center", color=INK_2, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(x_lo, 3e-11)
    ax.set_yticks(y)
    ax.set_yticklabels([c[0] for c in cases], color=INK_2)
    ax.set_xlabel("max interior |R(phi = x)|")
    ax.set_title("V0: freestream preserved to machine zero on every mesh type")
    ax.grid(axis="y", visible=False)
    finish(fig, OUT, "v0_freestream.png")
    write_csv(OUT, "v0_freestream.csv", "mesh,max_interior_residual",
              [(c[0].replace(",", ";"), f"{c[1]:.3e}") for c in cases])

    checks.add("V0", "max interior residual, phi=x", f"{vals.max():.2e}",
               "< 1e-12", bool(vals.max() < 1e-12))


# ---------------------------------------------------------------------------
# 2. G1.1 MMS convergence
# ---------------------------------------------------------------------------
def demo_mms(checks):
    print("[2/6] G1.1 MMS convergence")
    levels = [run_mms_case(n) for n in (4, 8, 16, 32)]
    h = np.array([lvl["h"] for lvl in levels])
    err = np.array([lvl["l2_error"] for lvl in levels])
    slope = np.polyfit(np.log(h), np.log(err), 1)[0]

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ref = err[0] * (h / h[0]) ** 2
    ax.loglog(h, ref, "--", color=MUTED, linewidth=1.4, zorder=2)
    ax.text(h[2] * 1.07, ref[2] * 0.62, "O(h$^2$) reference", color=MUTED, fontsize=9.5)
    ax.loglog(h, err, "o-", color=S1_BLUE, markeredgecolor=SURFACE,
              markeredgewidth=2, zorder=3)
    ax.text(0.05, 0.9, f"measured slope = {slope:.2f}  (gate: >= 1.9)",
            transform=ax.transAxes, color=INK_2, fontsize=10)
    ax.set_xlabel("h (cube edge / n)")
    ax.set_ylabel("L2 error of phi")
    ax.set_title("G1.1: 2nd-order MMS convergence "
                 "(phi = sin pix sin piy sin piz, 4 levels)")
    ax.grid(True, which="both")
    finish(fig, OUT, "g11_mms_convergence.png")
    write_csv(OUT, "g11_mms.csv", "n,h,l2_error,cg_iterations",
              [(l["n"], l["h"], f"{l['l2_error']:.6e}", l["n_cg_iterations"])
               for l in levels])

    checks.add("G1.1", "MMS L2 slope (4 levels)", f"{slope:.3f}",
               ">= 1.9", bool(slope >= 1.9))
    return slope


# ---------------------------------------------------------------------------
# 3. G1.2 CG+AMG mesh independence
# ---------------------------------------------------------------------------
def demo_cg_amg(checks):
    print("[3/6] G1.2 CG+AMG mesh independence")
    series = [(8, S1_BLUE), (16, S2_AQUA), (32, S3_YELLOW)]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.6))
    rows = []
    for n, color in series:
        nodes, elements = generate_structured_cube_mesh(n=n)
        dirichlet = np.where(cube_boundary_mask(nodes))[0]
        phi_d = phi_exact_fn(nodes[:, 0], nodes[:, 1], nodes[:, 2])[dirichlet]
        A = assemble_stiffness_matrix(nodes, elements)
        A_free, b_free, _, _ = apply_dirichlet(A, np.zeros(len(nodes)), dirichlet, phi_d)
        M = pyamg.smoothed_aggregation_solver(A_free).aspreconditioner()
        hist = []
        spla.cg(A_free, b_free, M=M, rtol=1e-10, maxiter=2000,
                callback=lambda xk: hist.append(np.linalg.norm(b_free - A_free @ xk)))
        hist = np.array(hist) / np.linalg.norm(b_free)
        ax1.semilogy(np.arange(1, len(hist) + 1), np.maximum(hist, 1e-16), "o-",
                     color=color, markeredgecolor=SURFACE, markeredgewidth=1.6,
                     markersize=6.5, label=f"n={n} ({len(nodes):,} nodes)")
        rows.append((n, len(nodes), len(hist)))

    ax1.set_xlabel("CG iteration")
    ax1.set_ylabel("relative residual")
    ax1.set_title("Convergence histories overlap across mesh levels")
    ax1.legend(loc="upper right")
    ax1.grid(True, which="both")

    n_nodes = np.array([r[1] for r in rows])
    iters = np.array([r[2] for r in rows])
    ax2.semilogx(n_nodes, iters, "o-", color=S1_BLUE,
                 markeredgecolor=SURFACE, markeredgewidth=2)
    for x, it in zip(n_nodes, iters):
        ax2.annotate(str(it), (x, it), xytext=(0, 9), textcoords="offset points",
                     ha="center", color=INK_2, fontsize=10)
    ax2.axhline(iters.min() * 2, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax2.text(n_nodes[0], iters.min() * 2 * 1.05, "gate: 2x growth bound",
             color=CRITICAL, fontsize=9)
    ax2.set_ylim(0, max(iters.max(), iters.min() * 2) * 1.35)
    ax2.set_xlabel("number of nodes")
    ax2.set_ylabel("CG iterations to 1e-10")
    ax2.set_title("Iteration count nearly flat over 64x more nodes")
    fig.suptitle("G1.2: AMG-preconditioned CG is mesh-independent",
                 fontsize=13, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    finish(fig, OUT, "g12_amg_cg_scaling.png")
    write_csv(OUT, "g12_cg_iterations.csv", "n,n_nodes,cg_iterations", rows)

    growth = iters.max() / iters.min()
    checks.add("G1.2", "CG iteration growth over 64x nodes", f"{growth:.2f}x",
               "< 2.0x", bool(growth < 2.0))


# ---------------------------------------------------------------------------
# 4. sphere flow-field slice (physical plausibility)
# ---------------------------------------------------------------------------
def demo_flowfield(mesh, result, checks):
    print("[4/6] sphere flow-field slice")
    import matplotlib.tri as mtri
    import pyvista as pv

    nodes, elements = mesh.nodes, mesh.elements
    grad = nodal_gradient_recovery(nodes, elements, result["phi"])
    q = np.linalg.norm(grad, axis=1)

    n_tets = len(elements)
    cells = np.hstack([np.full((n_tets, 1), 4, dtype=np.int64),
                       elements.astype(np.int64)]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(n_tets, pv.CellType.TETRA), nodes)
    grid.point_data["q"] = q
    grid.point_data["cp"] = 1.0 - q**2
    sl = grid.slice(normal=(0, 0, 1), origin=(0, 0, 0)).triangulate()
    pts = sl.points
    tris = sl.faces.reshape(-1, 4)[:, 1:]
    triang = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.4))
    for ax, field, title, cmap, norm, levels in (
        (axes[0], sl.point_data["q"], "speed |grad phi| / U_inf", SEQ_BLUE,
         None, np.linspace(0.0, 1.5, 16)),
        (axes[1], sl.point_data["cp"], "Cp (red = compression, blue = suction)",
         DIV_BLUE_RED, TwoSlopeNorm(vmin=-1.3, vcenter=0.0, vmax=1.0),
         np.linspace(-1.3, 1.0, 24)),
    ):
        cs = ax.tricontourf(triang, field, levels=levels, cmap=cmap, norm=norm,
                            extend="both")
        circle = plt.Circle((0, 0), 1.0, facecolor=SURFACE, edgecolor=INK_2,
                            linewidth=1.4, zorder=3)
        ax.add_patch(circle)
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-2.6, 2.6)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(False)
        cb = fig.colorbar(cs, ax=ax, shrink=0.85, pad=0.02)
        cb.ax.tick_params(labelsize=9, colors=MUTED)
        cb.outline.set_edgecolor(BASELINE)
    axes[0].annotate("flow", xy=(-2.15, 2.15), xytext=(-3.0, 2.15), color=INK_2,
                     fontsize=10, va="center",
                     arrowprops=dict(arrowstyle="-|>", color=INK_2, linewidth=1.4))
    fig.suptitle("Physical plausibility: 3D sphere potential field "
                 "(z = 0 slice, medium gmsh mesh, 95k tets)",
                 fontsize=13, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    finish(fig, OUT, "sphere_flowfield.png")

    # physics sanity on the volume solution: stagnation + fore-aft symmetry
    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    r = np.linalg.norm(nodes[wall_nodes], axis=1)
    cos_t = nodes[wall_nodes, 0] / r
    q_wall = q[wall_nodes]
    q_stag = q_wall[np.abs(cos_t) > 0.995].max()
    q_eq = q_wall[np.abs(cos_t) < 0.05]
    checks.add("V2-sanity", "residual norm (medium sphere solve)",
               f"{result['residual_norm']:.1e}", "< 1e-8",
               bool(result["residual_norm"] < 1e-8))
    checks.add("V2-sanity", "wall speed at stagnation poles", f"{q_stag:.3f}",
               "< 0.2 (exact 0)", bool(q_stag < 0.2))
    checks.add("V2-sanity", "wall speed at equator (exact 1.5)",
               f"{q_eq.mean():.3f}", "in [1.3, 1.6]",
               bool(1.3 < q_eq.mean() < 1.6))


# ---------------------------------------------------------------------------
# 5. G1.6 sphere Cp vs analytic (OPEN gate, XFAIL)
# ---------------------------------------------------------------------------
def demo_sphere_cp(checks):
    print("[5/6] G1.6 sphere Cp (open gate)")
    mesh = read_mesh(SPHERE_DIR / "medium.msh")
    result = solve_sphere(mesh)
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    cp_num = sphere_cp_numeric(mesh, result["phi"], wall_faces, wall_nodes)
    cp_ex, cos_t = sphere_cp_exact(mesh.nodes, wall_nodes)
    err = np.abs(cp_num - cp_ex)
    theta = np.degrees(np.arccos(np.clip(cos_t, -1, 1)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.8))
    tl = np.linspace(0, 180, 300)
    ax1.plot(tl, 1.0 - 2.25 * np.sin(np.radians(tl)) ** 2, "-", color=INK,
             linewidth=2, label="exact: 1 - (9/4) sin$^2$theta", zorder=3)
    ax1.plot(theta, cp_num, ".", color=S1_BLUE, markersize=3.5, alpha=0.35,
             label=f"FEM + quadratic recovery ({len(wall_nodes):,} wall nodes)",
             zorder=2)
    ax1.set_xlabel("theta from +x stagnation point (deg)")
    ax1.set_ylabel("Cp")
    ax1.set_title("Shape is right; an error layer remains")
    ax1.legend(loc="upper center")

    ax2.semilogy(theta, np.maximum(err, 1e-6), ".", color=S1_BLUE,
                 markersize=3.5, alpha=0.35)
    ax2.axhline(0.02, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax2.text(183, 0.02, "gate: 2%", color=CRITICAL, fontsize=9.5, va="center",
             clip_on=False)
    ax2.axhline(err.max(), color=MUTED, linewidth=1.0, linestyle=":")
    ax2.text(183, err.max(), f"max = {err.max()*100:.1f}%", color=INK_2,
             fontsize=9.5, va="center", clip_on=False)
    ax2.set_xlabel("theta (deg)")
    ax2.set_ylabel("|Cp error|")
    ax2.set_title("Flat-facet wall limit, not a solver bug")
    fig.suptitle("G1.6 OPEN gate: incompressible sphere Cp, medium mesh "
                 "-- documented gap, route = Option C + curved elements",
                 fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    finish(fig, OUT, "g16_sphere_cp_open_gate.png")

    checks.add("G1.6", "max |Cp err| sphere medium", f"{err.max()*100:.1f}%",
               "< 2% (gate)", bool(err.max() < 0.02), xfail=True,
               note="open gate; strict xfail in tests/test_laplace_sphere.py")
    return err


# ---------------------------------------------------------------------------
# 6. G1.4 oracle negative result (Option A ceiling)
# ---------------------------------------------------------------------------
def full_flux_rhs(mesh, wall_faces, geometry):
    tri = mesh.nodes[wall_faces]
    qp = np.einsum("qi,fik->fqk", BARY7, tri)
    g = sphere_grad_exact(qp.reshape(-1, 3)).reshape(qp.shape)
    flux = np.einsum("fqk,fk->fq", g, geometry["n_facet"])
    rhs = np.zeros(len(mesh.nodes))
    for v in range(3):
        contrib = geometry["area"][:, None] * W7[None, :] * BARY7[None, :, v] * flux
        np.add.at(rhs, wall_faces[:, v], contrib.sum(axis=1))
    return rhs


def demo_oracle(checks):
    print("[6/6] G1.4 oracle experiment (negative result)")
    rows, results = [], {}
    for level in ("coarse", "medium"):
        mesh = read_mesh(SPHERE_DIR / f"{level}.msh")
        wall_faces = mesh.boundary_faces["wall"]
        wall_nodes = np.unique(wall_faces)
        cp_ex, _ = sphere_cp_exact(mesh.nodes, wall_nodes)

        geometry = wall_correction_geometry(
            mesh.nodes, mesh.elements, wall_faces, sphere_closest_point_normal)
        grad_qp = sphere_grad_exact(geometry["qp"].reshape(-1, 3)).reshape(
            geometry["qp"].shape)
        rhs_t = assemble_wall_flux_correction_rhs(
            len(mesh.nodes), wall_faces, geometry, grad_qp)
        rhs_full = full_flux_rhs(mesh, wall_faces, geometry)

        errs = {}
        for name, rhs in (("uncorrected", None), ("t-form", rhs_t),
                          ("full-flux", rhs_full)):
            res = solve_sphere(mesh, body_source_rhs=rhs)
            cp = sphere_cp_numeric(mesh, res["phi"], wall_faces, wall_nodes)
            errs[name] = np.abs(cp - cp_ex).max()
            rows.append((level, name, f"{errs[name]:.4f}",
                         "" if rhs is None else f"{np.abs(rhs).max():.2e}"))
            print(f"    {level:7s} {name:12s} max|Cp err| = {errs[name]:.4f}")
        results[level] = {"errs": errs, "rhs_t_max": np.abs(rhs_t).max()}

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    variants = ["uncorrected", "t-form", "full-flux"]
    colors = [S1_BLUE, S2_AQUA, S5_VIOLET]
    width = 0.24
    x = np.arange(2)
    for i, (v, c) in enumerate(zip(variants, colors)):
        vals = [results[lv]["errs"][v] * 100 for lv in ("coarse", "medium")]
        bars = ax.bar(x + (i - 1) * width, vals, width * 0.92, color=c, label=v)
        for b, val in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, val + 0.25, f"{val:.1f}%",
                    ha="center", color=INK_2, fontsize=9.5)
    ax.axhline(2.0, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax.text(1.42, 2.25, "gate: 2%", color=CRITICAL, fontsize=9.5)
    ax.axhline(5.0, color=MUTED, linewidth=1.2, linestyle=":")
    ax.text(1.42, 5.25, "DP1 branch point: 5%", color=INK_2, fontsize=9.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["coarse (17.5k tets)", "medium (95k tets)"], color=INK_2)
    ax.set_ylabel("max |Cp error| at wall (%)")
    ax.set_ylim(0, max(results["coarse"]["errs"].values()) * 100 * 1.28)
    ax.set_title("G1.4 oracle: even the EXACT-gradient boundary correction "
                 "cannot close the gap\n(correction RHS is tiny on body-fitted "
                 "walls: max ~7e-5 medium -- almost no boundary-data defect to fix)")
    ax.legend(loc="upper right", title="correction variant")
    ax.grid(axis="x", visible=False)
    finish(fig, OUT, "g14_oracle_negative_result.png")
    write_csv(OUT, "g14_oracle.csv", "level,variant,max_cp_err,max_rhs", rows)

    med = results["medium"]
    best = min(med["errs"]["t-form"], med["errs"]["full-flux"])
    delta = abs(med["errs"]["uncorrected"] - best)
    checks.add("G1.4", "Option A t-form RHS magnitude (medium)",
               f"{med['rhs_t_max']:.1e}", "< 1e-3 (near-zero data defect)",
               bool(med["rhs_t_max"] < 1e-3),
               note="body-fitted wall leaves ~no flux defect; G1.3 cylinder "
                    "analogue is exactly machine zero (1.5e-17)")
    checks.add("G1.4", "correction moves max |Cp err| by",
               f"{delta*100:.2f} pp", "< 0.5 pp (correction ineffective)",
               bool(delta < 0.005))
    checks.add("DP1", "best oracle-corrected max |Cp err|", f"{best*100:.1f}%",
               "> 5% (confirms '>5%' branch)", bool(best > 0.05),
               note="boundary-data corrections ruled out; route = Option C")


def main():
    apply_style()
    t0 = time.time()
    checks = CheckList("P1 Laplace solver (V0, G1.1, G1.2 closed; G1.4/DP1 "
                       "negative-result evidence; G1.6 open)")

    demo_freestream(checks)
    demo_mms(checks)
    demo_cg_amg(checks)

    mesh_med = read_mesh(SPHERE_DIR / "medium.msh")
    result_med = solve_sphere(mesh_med)
    demo_flowfield(mesh_med, result_med, checks)
    demo_sphere_cp(checks)
    demo_oracle(checks)

    write_csv(OUT, "summary.csv", "metric,value",
              [(c["gate"] + " " + c["name"].replace(",", ";"),
                str(c["value"]) + " [" + c["status"] + "]")
               for c in checks.checks] + [("runtime_seconds", f"{time.time()-t0:.1f}")])
    code = checks.report(OUT)
    print(f"done in {time.time() - t0:.1f}s -> {OUT.relative_to(REPO_ROOT)}/")
    sys.exit(code)


if __name__ == "__main__":
    main()
