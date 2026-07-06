"""
P1 capability demonstration: headless PNG/CSV artifacts showing what the
solver can do today (P0+P1: Laplace on unstructured tets, CG+AMG, surface
recovery) and where the known gaps are (G1.6 curved-wall accuracy, no wake/
compressibility yet).

Not a gate test -- gates live in tests/. This script exists so a human can
*see* the current state in one pass. Everything runs headless (matplotlib
Agg; pyvista used only as a slice filter, no rendering).

Usage:
    python cases/demo/p1_capability_demo.py            # full run, ~2-4 min
    python cases/demo/p1_capability_demo.py --fast     # skip the finest family level

Outputs: artifacts/demo_P1/*.png + *.csv
"""

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyamg
import scipy.sparse.linalg as spla
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.kernels.residual import assemble_residual, assemble_stiffness_matrix
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.surface import (
    nodal_gradient_recovery,
    wall_tangential_gradient,
    wall_tangential_gradient_quadratic,
)
from pyfp3d.solve.linear import apply_dirichlet
from pyfp3d.solve.picard import solve_laplace

from tests.mesh_utils import (
    cube_boundary_mask,
    generate_structured_cube_mesh,
    generate_sphere_shell_mesh,
    icosphere,
)

OUT_DIR = REPO_ROOT / "artifacts" / "demo_P1"
MESH_DIR = REPO_ROOT / "cases" / "meshes" / "sphere_shell"

# ---------------------------------------------------------------------------
# Chart style (single light-mode palette for static PNG artifacts)
# ---------------------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
S1_BLUE = "#2a78d6"
S2_AQUA = "#1baf7a"
S3_YELLOW = "#eda100"
S5_VIOLET = "#4a3aa7"
CRITICAL = "#d03b3b"

SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue",
    ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
)
DIV_BLUE_RED = LinearSegmentedColormap.from_list(
    "div_blue_red",
    ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a9a8", "#e34948", "#8f1f1f"],
)

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": BASELINE,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "semibold",
        "axes.labelsize": 10.5,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "grid.linestyle": "-",
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.frameon": False,
        "legend.fontsize": 9.5,
        "legend.labelcolor": INK_2,
        "font.family": "sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 8,
    }
)


def _finish(fig, name):
    path = OUT_DIR / name
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.relative_to(REPO_ROOT)}")


def _write_csv(name, header, rows):
    path = OUT_DIR / name
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(v) for v in row) + "\n")
    print(f"  wrote {path.relative_to(REPO_ROOT)}")


# ---------------------------------------------------------------------------
# Shared sphere-case helpers (mirrors tests/test_laplace_sphere.py)
# ---------------------------------------------------------------------------
def sphere_phi_exact(nodes, a=1.0):
    r = np.linalg.norm(nodes, axis=1)
    return nodes[:, 0] * (1.0 + 0.5 * a**3 / r**3)


def sphere_cp_exact(nodes, wall_nodes):
    r = np.linalg.norm(nodes[wall_nodes], axis=1)
    cos_t = nodes[wall_nodes, 0] / r
    return 1.0 - 2.25 * (1.0 - cos_t**2), cos_t


def solve_gmsh_sphere(mesh_path):
    mesh = read_mesh(mesh_path)
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])
    phi_exact = sphere_phi_exact(nodes)
    result = solve_laplace(
        nodes, elements, farfield_nodes, phi_exact[farfield_nodes], rtol=1e-11, maxiter=3000
    )
    return mesh, wall_faces, wall_nodes, phi_exact, result


def cp_from_gradient(grad_wall, wall_nodes):
    q2 = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
    return 1.0 - q2


# ---------------------------------------------------------------------------
# 1. V0 freestream preservation across mesh types
# ---------------------------------------------------------------------------
def demo_freestream():
    print("[1/6] V0 freestream preservation")
    cases = []

    nodes, elements = generate_structured_cube_mesh(n=8)
    interior = ~cube_boundary_mask(nodes)
    R = assemble_residual(nodes, elements, nodes[:, 0].copy())
    cases.append(("structured cube (Kuhn, 3.1k tets)", len(elements), np.abs(R[interior]).max()))

    nodes, elements, wall, far = generate_sphere_shell_mesh(subdivisions=2, n_layers=12)
    interior = np.ones(len(nodes), dtype=bool)
    interior[wall] = False
    interior[far] = False
    R = assemble_residual(nodes, elements, nodes[:, 0].copy())
    cases.append(("icosphere shell (9.7k tets)", len(elements), np.abs(R[interior]).max()))

    for level in ("coarse", "medium"):
        mesh = read_mesh(MESH_DIR / f"{level}.msh")
        boundary_nodes = np.unique(
            np.concatenate([f.ravel() for f in mesh.boundary_faces.values()])
        )
        interior = np.ones(len(mesh.nodes), dtype=bool)
        interior[boundary_nodes] = False
        R = assemble_residual(mesh.nodes, mesh.elements, mesh.nodes[:, 0].copy())
        cases.append(
            (f"gmsh sphere shell, {level} ({len(mesh.elements)/1e3:.0f}k tets)",
             len(mesh.elements), np.abs(R[interior]).max())
        )

    labels = [c[0] for c in cases]
    vals = np.array([c[2] for c in cases])
    x_lo = vals.min() * 0.35

    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    y = np.arange(len(cases))[::-1]
    ax.hlines(y, x_lo, vals, color=GRID, linewidth=1.0, zorder=2)
    ax.plot(vals, y, "o", color=S1_BLUE, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
    ax.axvline(1e-12, color=CRITICAL, linewidth=1.2, linestyle="--", zorder=2)
    ax.text(1.4e-12, 3.42, "gate: 1e-12", color=CRITICAL, fontsize=9, va="center")
    for yi, v in zip(y, vals):
        ax.text(v * 2.2, yi, f"{v:.1e}", va="center", color=INK_2, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(x_lo, 3e-11)
    ax.set_ylim(-0.5, 3.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=INK_2)
    ax.set_xlabel("max interior |R(phi = x)|")
    ax.set_title("V0 capability: freestream preserved to machine zero on every mesh type")
    ax.grid(axis="y", visible=False)
    _finish(fig, "01_v0_freestream.png")
    _write_csv("freestream_residuals.csv", "mesh,n_tets,max_interior_residual",
               [(l.replace(",", ";"), n, f"{v:.3e}") for (l, n, v) in cases])


# ---------------------------------------------------------------------------
# 2. G1.1 MMS convergence
# ---------------------------------------------------------------------------
def demo_mms():
    print("[2/6] G1.1 MMS convergence")
    from tests.test_laplace_mms import run_mms_case

    levels = [run_mms_case(n) for n in (4, 8, 16, 32)]
    h = np.array([lvl["h"] for lvl in levels])
    err = np.array([lvl["l2_error"] for lvl in levels])
    slope = np.polyfit(np.log(h), np.log(err), 1)[0]

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ref = err[0] * (h / h[0]) ** 2
    ax.loglog(h, ref, "--", color=MUTED, linewidth=1.4, zorder=2)
    ax.text(h[2] * 1.07, ref[2] * 0.62, "O(h$^2$) reference", color=MUTED, fontsize=9.5)
    ax.loglog(h, err, "o-", color=S1_BLUE, markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
    ax.text(0.05, 0.9, f"measured, fitted slope = {slope:.2f}",
            transform=ax.transAxes, color=INK_2, fontsize=10)
    ax.set_xlabel("h (cube edge / n)")
    ax.set_ylabel("L2 error of phi")
    ax.set_title("G1.1 capability: 2nd-order MMS convergence (sin-product solution)")
    ax.grid(True, which="both")
    _finish(fig, "02_g11_mms_convergence.png")
    _write_csv("mms_convergence.csv", "n,h,l2_error,cg_iterations",
               [(l["n"], l["h"], f"{l['l2_error']:.6e}", l["n_cg_iterations"]) for l in levels])
    return slope


# ---------------------------------------------------------------------------
# 3. G1.2 CG+AMG mesh independence
# ---------------------------------------------------------------------------
def demo_cg_amg():
    print("[3/6] G1.2 CG+AMG mesh independence")
    from tests.test_laplace_mms import phi_exact_fn

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
        ax1.annotate(f"n={n}", (len(hist), hist[-1]), xytext=(6, 0),
                     textcoords="offset points", color=INK_2, fontsize=9, va="center")
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
    ax2.set_ylim(0, max(iters) * 1.6)
    ax2.set_xlabel("number of nodes")
    ax2.set_ylabel("CG iterations to 1e-10")
    ax2.set_title("Iteration count nearly flat over 64x more nodes")
    fig.suptitle("G1.2 capability: AMG-preconditioned CG is mesh-independent",
                 fontsize=13, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _finish(fig, "03_g13_cg_amg.png")
    _write_csv("cg_iterations.csv", "n,n_nodes,cg_iterations", rows)


# ---------------------------------------------------------------------------
# 4. G1.6 sphere Cp: current accuracy vs the analytic solution
# ---------------------------------------------------------------------------
def demo_sphere_cp():
    print("[4/6] G1.6 sphere Cp (medium gmsh mesh)")
    mesh, wall_faces, wall_nodes, phi_exact, result = solve_gmsh_sphere(MESH_DIR / "medium.msh")
    nodes = mesh.nodes

    grad_wall = wall_tangential_gradient_quadratic(nodes, wall_faces, result["phi"])
    cp_num = cp_from_gradient(grad_wall, wall_nodes)
    cp_ex, cos_t = sphere_cp_exact(nodes, wall_nodes)
    err = np.abs(cp_num - cp_ex)
    theta = np.degrees(np.arccos(np.clip(cos_t, -1, 1)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.8))
    tl = np.linspace(0, 180, 300)
    ax1.plot(tl, 1.0 - 2.25 * np.sin(np.radians(tl)) ** 2, "-", color=INK,
             linewidth=2, label="exact: 1 - (9/4) sin$^2$theta", zorder=3)
    ax1.plot(theta, cp_num, ".", color=S1_BLUE, markersize=3.5, alpha=0.35,
             label=f"FEM + quadratic recovery ({len(wall_nodes):,} wall nodes)", zorder=2)
    ax1.set_xlabel("theta from +x stagnation point (deg)")
    ax1.set_ylabel("Cp")
    ax1.set_title("Shape is right; boundary layer of error remains")
    ax1.legend(loc="upper center")

    ax2.semilogy(theta, np.maximum(err, 1e-6), ".", color=S1_BLUE, markersize=3.5, alpha=0.35)
    ax2.axhline(0.02, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax2.text(183, 0.02, "gate: 2%", color=CRITICAL, fontsize=9.5, va="center", clip_on=False)
    ax2.axhline(err.max(), color=MUTED, linewidth=1.0, linestyle=":")
    ax2.text(183, err.max(), f"max = {err.max()*100:.1f}%", color=INK_2, fontsize=9.5,
             va="center", clip_on=False)
    ax2.set_xlabel("theta (deg)")
    ax2.set_ylabel("|Cp error|")
    ax2.set_title(f"Gap: max {err.max()*100:.1f}%, mean {err.mean()*100:.2f}% (target < 2%)")
    fig.suptitle("G1.6 open gate: incompressible sphere Cp on the medium mesh",
                 fontsize=13, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _finish(fig, "04_g12_sphere_cp.png")

    return mesh, wall_faces, wall_nodes, phi_exact, result, err


# ---------------------------------------------------------------------------
# 5. Flow-field slice (visual sanity of the 3D volume solution)
# ---------------------------------------------------------------------------
def demo_flowfield(mesh, result):
    print("[5/6] flow-field slice through the volume solution")
    import pyvista as pv

    nodes, elements = mesh.nodes, mesh.elements
    grad = nodal_gradient_recovery(nodes, elements, result["phi"])
    q = np.linalg.norm(grad, axis=1)
    cp = 1.0 - q**2

    n_tets = len(elements)
    cells = np.hstack([np.full((n_tets, 1), 4, dtype=np.int64),
                       elements.astype(np.int64)]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(n_tets, pv.CellType.TETRA), nodes)
    grid.point_data["q"] = q
    grid.point_data["cp"] = cp
    sl = grid.slice(normal=(0, 0, 1), origin=(0, 0, 0)).triangulate()

    pts = sl.points
    tris = sl.faces.reshape(-1, 4)[:, 1:]
    import matplotlib.tri as mtri

    triang = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.4))
    for ax, field, title, cmap, norm, levels in (
        (axes[0], sl.point_data["q"], "speed |grad phi| / U_inf", SEQ_BLUE,
         None, np.linspace(0.0, 1.5, 16)),
        (axes[1], sl.point_data["cp"], "Cp (red = compression, blue = suction)",
         DIV_BLUE_RED, TwoSlopeNorm(vmin=-1.3, vcenter=0.0, vmax=1.0),
         np.linspace(-1.3, 1.0, 24)),
    ):
        cs = ax.tricontourf(triang, field, levels=levels, cmap=cmap, norm=norm, extend="both")
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
    fig.suptitle("Capability: full 3D potential field (z = 0 slice, medium gmsh mesh)",
                 fontsize=13, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _finish(fig, "05_flowfield_slice.png")


# ---------------------------------------------------------------------------
# 6. G1.6 error anatomy: family convergence + oracle decomposition
# ---------------------------------------------------------------------------
def demo_error_anatomy(gmsh_results, fast=False):
    print("[6/6] G1.6 error anatomy (this is the slow section)")
    subdivs = (1, 2, 3) if fast else (1, 2, 3, 4)
    n_layers = {1: 12, 2: 20, 3: 31, 4: 49}

    rows = []
    for s in subdivs:
        t0 = time.time()
        nodes, elements, wall_nodes, far_nodes = generate_sphere_shell_mesh(
            subdivisions=s, n_layers=n_layers[s], r_inner=1.0, r_outer=25.0
        )
        _, wall_faces = icosphere(s)
        phi_exact = sphere_phi_exact(nodes)
        result = solve_laplace(nodes, elements, far_nodes, phi_exact[far_nodes],
                               rtol=1e-11, maxiter=3000)

        grad_wall = wall_tangential_gradient_quadratic(nodes, wall_faces, result["phi"])
        cp_num = cp_from_gradient(grad_wall, wall_nodes)
        cp_ex, _ = sphere_cp_exact(nodes, wall_nodes)
        cp_err = np.abs(cp_num - cp_ex).max()
        phi_err = np.abs(result["phi"][wall_nodes] - phi_exact[wall_nodes]).max()

        edges = np.concatenate([
            np.linalg.norm(nodes[wall_faces[:, i]] - nodes[wall_faces[:, j]], axis=1)
            for i, j in ((0, 1), (1, 2), (2, 0))
        ])
        h_wall = edges.mean()
        rows.append((s, len(nodes), len(elements), h_wall, cp_err, phi_err))
        print(f"    subdiv={s}: {len(nodes):,} nodes, {len(elements):,} tets, "
              f"h_wall={h_wall:.3f}, max|Cp err|={cp_err:.3f}, "
              f"max|phi err|={phi_err:.2e}  ({time.time()-t0:.1f}s)")

    h = np.array([r[3] for r in rows])
    cp_e = np.array([r[4] for r in rows])
    phi_e = np.array([r[5] for r in rows])
    orders_cp = np.log(cp_e[:-1] / cp_e[1:]) / np.log(h[:-1] / h[1:])
    orders_phi = np.log(phi_e[:-1] / phi_e[1:]) / np.log(h[:-1] / h[1:])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.0))

    ax1.set_xlim(h[-1] * 0.72, h[0] * 1.55)
    ax1.loglog(h, cp_e, "o-", color=S1_BLUE, markeredgecolor=SURFACE,
               markeredgewidth=2, label="max |Cp error| at wall")
    ax1.loglog(h, phi_e, "s-", color=S5_VIOLET, markeredgecolor=SURFACE,
               markeredgewidth=2, label="max |phi error| at wall")
    for ref_p, anchor in ((1, cp_e[0] * 1.7), (2, cp_e[0] * 0.55)):
        ref = anchor * (h / h[0]) ** ref_p
        ax1.loglog(h, ref, "--", color=MUTED, linewidth=1.2)
        ax1.text(h[0] * 1.12, ref[0], f"O(h$^{ref_p}$)", color=MUTED, fontsize=9, va="center")
    ax1.axhline(0.02, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax1.text(h[0] * 1.12, 0.02, "gate: 2%", color=CRITICAL, fontsize=9, va="center")
    ax1.annotate(f"observed order {orders_cp[-1]:.2f}", (h[-1], cp_e[-1]),
                 xytext=(-4, 12), textcoords="offset points", color=INK_2, fontsize=9.5)
    ax1.annotate(f"observed order {orders_phi[-1]:.2f}", (h[-1], phi_e[-1]),
                 xytext=(-4, -18), textcoords="offset points", color=INK_2, fontsize=9.5)
    ax1.set_xlabel("mean wall-edge length h")
    ax1.set_ylabel("max error at wall")
    ax1.set_title("Refinement saturates: sub-first-order at the curved wall\n"
                  "(icosphere shell family, single-variable sweep)")
    ax1.legend(loc="lower right")
    ax1.grid(True, which="both")

    # Oracle decomposition on the committed gmsh meshes: total pipeline error
    # vs the recovery operator's own error when fed the exact potential.
    labels, full_err, oracle_err = [], [], []
    for level, (mesh, wall_faces, wall_nodes, phi_exact, result, err) in gmsh_results.items():
        grad_oracle = wall_tangential_gradient_quadratic(mesh.nodes, wall_faces, phi_exact)
        cp_oracle = cp_from_gradient(grad_oracle, wall_nodes)
        cp_ex, _ = sphere_cp_exact(mesh.nodes, wall_nodes)
        labels.append(f"gmsh {level}")
        full_err.append(err.max() * 100)
        oracle_err.append(np.abs(cp_oracle - cp_ex).max() * 100)

    y = np.arange(len(labels))[::-1]
    for yi, fe, oe in zip(y, full_err, oracle_err):
        ax2.plot([oe, fe], [yi, yi], "-", color=GRID, linewidth=2.5, zorder=2)
        ax2.annotate("", xy=(fe, yi - 0.22), xytext=(oe, yi - 0.22),
                     arrowprops=dict(arrowstyle="<->", color=MUTED, linewidth=1.0))
        ax2.text((fe + oe) / 2, yi - 0.34, "volume-solve error", ha="center",
                 va="top", color=INK_2, fontsize=9)
    ax2.plot(full_err, y, "o", color=S1_BLUE, markersize=11,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=3,
             label="full pipeline (FEM solve + recovery)")
    ax2.plot(oracle_err, y, "o", color=S2_AQUA, markersize=11,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=3,
             label="recovery only (exact phi fed in)")
    for yi, fe, oe in zip(y, full_err, oracle_err):
        ax2.text(fe, yi + 0.16, f"{fe:.1f}%", ha="center", color=INK_2, fontsize=9.5)
        ax2.text(oe, yi + 0.16, f"{oe:.2f}%", ha="center", color=INK_2, fontsize=9.5)
    ax2.set_yticks(y)
    ax2.set_yticklabels(labels, color=INK_2)
    ax2.set_xlabel("max |Cp error| at wall (%)")
    ax2.set_xlim(-0.5, max(full_err) * 1.25)
    ax2.set_ylim(-0.8, len(labels) - 0.2 + 0.6)
    ax2.axvline(2.0, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax2.text(2.0, len(labels) - 0.25, " gate: 2%", color=CRITICAL, fontsize=9)
    ax2.set_title("Error split: the PDE solve dominates, not the recovery\n"
                  "(oracle test on committed gmsh meshes)")
    ax2.legend(loc="lower right")
    ax2.grid(axis="y", visible=False)

    fig.suptitle("G1.6 root cause, visualized: flat-facet wall (variational crime) "
                 "limits accuracy -- needs isoparametric wall elements",
                 fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _finish(fig, "06_g12_error_anatomy.png")

    _write_csv("sphere_family_convergence.csv",
               "subdiv,n_nodes,n_tets,h_wall,max_cp_err,max_phi_err_wall",
               [(s, nn, nt, f"{hw:.4f}", f"{ce:.4e}", f"{pe:.4e}")
                for (s, nn, nt, hw, ce, pe) in rows])
    _write_csv("oracle_decomposition.csv",
               "mesh,full_pipeline_max_cp_err_pct,recovery_only_max_cp_err_pct",
               [(l, f"{fe:.3f}", f"{oe:.4f}") for l, fe, oe in
                zip(labels, full_err, oracle_err)])
    return rows, orders_cp, orders_phi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast", action="store_true",
                        help="skip the finest icosphere family level")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    demo_freestream()
    mms_slope = demo_mms()
    demo_cg_amg()

    gmsh_results = {}
    mesh, wall_faces, wall_nodes, phi_exact, result, err_med = demo_sphere_cp()
    gmsh_results["medium"] = (mesh, wall_faces, wall_nodes, phi_exact, result, err_med)
    demo_flowfield(mesh, result)

    print("    solving gmsh coarse for the oracle comparison...")
    mesh_c, wf_c, wn_c, pe_c, res_c = solve_gmsh_sphere(MESH_DIR / "coarse.msh")
    grad_c = wall_tangential_gradient_quadratic(mesh_c.nodes, wf_c, res_c["phi"])
    cp_c = cp_from_gradient(grad_c, wn_c)
    cp_ex_c, _ = sphere_cp_exact(mesh_c.nodes, wn_c)
    gmsh_results["coarse"] = (mesh_c, wf_c, wn_c, pe_c, res_c, np.abs(cp_c - cp_ex_c))
    gmsh_results = {k: gmsh_results[k] for k in ("coarse", "medium")}

    rows, orders_cp, orders_phi = demo_error_anatomy(gmsh_results, fast=args.fast)

    _write_csv(
        "summary.csv", "metric,value",
        [
            ("mms_l2_slope", f"{mms_slope:.3f}"),
            ("sphere_medium_max_cp_err", f"{err_med.max():.4f}"),
            ("sphere_medium_mean_cp_err", f"{err_med.mean():.4f}"),
            ("family_observed_cp_order_finest", f"{orders_cp[-1]:.3f}"),
            ("family_observed_phi_order_finest", f"{orders_phi[-1]:.3f}"),
            ("g12_gate_target", "0.02"),
            ("g12_gate_status", "OPEN"),
            ("runtime_seconds", f"{time.time() - t0:.1f}"),
        ],
    )
    print(f"done in {time.time() - t0:.1f}s -> {OUT_DIR.relative_to(REPO_ROOT)}/")


if __name__ == "__main__":
    main()
