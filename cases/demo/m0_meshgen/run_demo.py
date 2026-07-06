"""
M0 demo -- quasi-2D meshing pipeline evidence (Track M, closed 2026-07-06).

M0 delivers pyfp3d/meshgen/ (vanilla-Gmsh planar mesh + single-layer
extrusion + globally consistent prism->3-tet split) and the committed
mesh families under cases/meshes/. Its formal acceptance gate is G2.5
(the Track M <-> Track P link, shown in cases/demo/p2_kutta_lifting/);
this demo shows the mesh-side evidence:

  1. What the meshes look like: airfoil + embedded wake line + cylinder,
     with LE/TE resolution where the physics needs it.
  2. Topology soundness: quad-split consistency and wake-cut topology
     asserts pass on EVERY committed mesh (agent-rules hard rule 7).
  3. Mesh quality: edge-length grading and bounded aspect ratios across
     refinement levels.
  4. The meshes solve and converge: cylinder wall Cp vs the analytic
     1 - 4 sin^2(theta) over coarse/medium/fine, max error dropping at
     the expected ~O(h) (curved-wall recovery limit, same mechanism as
     G1.6 -- documented, not a pipeline defect).

Standalone + self-checking:  python cases/demo/m0_meshgen/run_demo.py
Outputs: cases/demo/m0_meshgen/results/{*.png, summary.csv, checks.csv}
Exit code 0 iff every acceptance check passes. Runtime ~1 min.
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    CRITICAL, GRID, INK, INK_2, MESH_DIR, MUTED, REPO_ROOT,
    S1_BLUE, S2_AQUA, S3_YELLOW, S5_VIOLET, SURFACE,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.metrics import (  # noqa: E402
    compute_aspect_ratios, compute_edge_lengths, compute_tet_volumes,
)
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import assert_wake_topology, cut_wake  # noqa: E402
from pyfp3d.meshgen.extrude import assert_quad_split_consistency  # noqa: E402
from tests.mesh_utils import run_cylinder_case  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"


def bottom_plane_triangulation(mesh):
    """2D triangulation of the bottom symmetry plane of an extruded mesh."""
    z = np.round(mesh.nodes[:, 2], 12)
    z_bot = min(set(z))
    tris = [tri for tri in mesh.boundary_faces["symmetry"]
            if np.all(z[tri] == z_bot)]
    return np.asarray(tris), z_bot


def boundary_polyline(mesh, tag, z_bot):
    """Nodes of a tagged boundary lying on the bottom plane, as 2D points."""
    z = np.round(mesh.nodes[:, 2], 12)
    nodes = np.unique(mesh.boundary_faces[tag])
    nodes = nodes[z[nodes] == z_bot]
    return mesh.nodes[nodes][:, :2]


# ---------------------------------------------------------------------------
# 1. mesh gallery
# ---------------------------------------------------------------------------
def demo_gallery():
    print("[1/4] mesh gallery")
    naca = read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh")
    cyl = read_mesh(MESH_DIR / "cylinder_2.5d" / "coarse.msh")

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.8),
                             gridspec_kw={"width_ratios": [1.15, 1.15, 1.0]})

    for ax, (x0, x1, y0, y1), title in (
        (axes[0], (-2.0, 6.0, -3.0, 3.0), "NACA0012 + embedded wake line"),
        (axes[1], (-0.25, 1.45, -0.55, 0.55), "airfoil zoom: LE/TE clustering"),
    ):
        tris, z_bot = bottom_plane_triangulation(naca)
        ax.triplot(naca.nodes[:, 0], naca.nodes[:, 1], tris,
                   color=GRID, linewidth=0.45)
        wall = boundary_polyline(naca, "wall", z_bot)
        order = np.argsort(np.arctan2(wall[:, 1] - 0.0, wall[:, 0] - 0.5))
        ax.plot(np.append(wall[order, 0], wall[order[0], 0]),
                np.append(wall[order, 1], wall[order[0], 1]),
                "-", color=INK, linewidth=1.3)
        wake = boundary_polyline(naca, "wake", z_bot)
        worder = np.argsort(wake[:, 0])
        ax.plot(wake[worder, 0], wake[worder, 1], "-", color=S5_VIOLET,
                linewidth=1.8)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.grid(False)
    axes[0].text(3.2, 0.18, "wake sheet (embedded,\nconforming tet faces)",
                 color=S5_VIOLET, fontsize=9)

    ax = axes[2]
    tris, z_bot = bottom_plane_triangulation(cyl)
    ax.triplot(cyl.nodes[:, 0], cyl.nodes[:, 1], tris, color=GRID, linewidth=0.45)
    wall = boundary_polyline(cyl, "wall", z_bot)
    order = np.argsort(np.arctan2(wall[:, 1], wall[:, 0]))
    ax.plot(np.append(wall[order, 0], wall[order[0], 0]),
            np.append(wall[order, 1], wall[order[0], 1]),
            "-", color=INK, linewidth=1.3)
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)
    ax.set_aspect("equal")
    ax.set_title("cylinder annulus (G1.3 testbed)")
    ax.grid(False)
    fig.suptitle("M0 mesh family: vanilla-Gmsh planar mesh, single-layer "
                 "extrusion, consistent prism->3-tet split "
                 "(bottom-plane view, coarse level)",
                 fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    finish(fig, OUT, "mesh_gallery.png")


# ---------------------------------------------------------------------------
# 2. topology asserts across the whole committed family (hard rule 7)
# ---------------------------------------------------------------------------
def demo_topology(checks):
    print("[2/4] topology asserts on every committed mesh")
    rows = []
    all_ok = True
    for msh in sorted(MESH_DIR.glob("*/*.msh")):
        name = f"{msh.parent.name}/{msh.stem}"
        mesh = read_mesh(msh)
        extruded = msh.parent.name.endswith("_2.5d")
        has_wake = "wake" in mesh.boundary_faces

        def run(fn):
            try:
                fn()
                return "pass"
            except AssertionError as exc:  # pragma: no cover - evidence path
                nonlocal_ok[0] = False
                return f"FAIL: {exc}"

        nonlocal_ok = [True]
        tags = "pass" if {"wall", "farfield"} <= set(mesh.boundary_faces) \
            else "FAIL"
        quad = run(lambda: assert_quad_split_consistency(
            mesh, interior_groups=("wake",) if has_wake else ())) \
            if extruded else "n/a"

        if has_wake:
            def wake_check():
                mesh_cut, wc = cut_wake(mesh)
                assert_wake_topology(mesh_cut, wc)
            wake = run(wake_check)
        else:
            wake = "n/a"
        ok = nonlocal_ok[0] and tags == "pass"
        all_ok &= ok
        rows.append((name, len(mesh.elements), tags, quad, wake))
        print(f"    {name}: tags={tags} quad-split={quad} wake-cut={wake}")

    fig, ax = plt.subplots(figsize=(9.2, 0.55 * len(rows) + 1.6))
    ax.axis("off")
    cols = ["mesh", "n_tets", "tags complete", "quad-split\nconsistency",
            "wake-cut topology\n(cut + asserts)"]
    cell_text = [[r[0], f"{r[1]:,}", r[2], r[3], r[4]] for r in rows]
    table = ax.table(cellText=cell_text, colLabels=cols, loc="center",
                     cellLoc="center", colWidths=[0.3, 0.14, 0.17, 0.17, 0.22])
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.55)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.get_text().set_color(INK_2)
        if r == 0:
            cell.get_text().set_color(INK)
            cell.get_text().set_weight("semibold")
            cell.set_facecolor(SURFACE)
        elif c >= 2:
            txt = cell.get_text().get_text()
            if txt == "pass":
                cell.get_text().set_color(S2_AQUA)
                cell.get_text().set_weight("semibold")
            elif txt.startswith("FAIL"):
                cell.get_text().set_color(CRITICAL)
                cell.get_text().set_weight("semibold")
            else:
                cell.get_text().set_color(MUTED)
    ax.set_title("Hard rule 7: topology asserts on every mesh in cases/meshes/ "
                 "(n/a = check does not apply to that family)",
                 color=INK, fontsize=12, fontweight="semibold", pad=18)
    finish(fig, OUT, "topology_asserts.png")
    write_csv(OUT, "topology_asserts.csv",
              "mesh,n_tets,tags_complete,quad_split,wake_cut",
              [(r[0].replace(",", ";"), r[1], r[2], r[3], r[4]) for r in rows])

    checks.add("M0", f"topology asserts on {len(rows)} committed meshes",
               "all pass", "no assert fires (hard rule 7)", bool(all_ok))


# ---------------------------------------------------------------------------
# 3. mesh quality across refinement levels
# ---------------------------------------------------------------------------
def demo_quality(checks):
    print("[3/4] mesh quality")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.4))
    rows = []

    for level, color in (("coarse", S1_BLUE), ("medium", S2_AQUA)):
        mesh = read_mesh(MESH_DIR / "naca0012_2.5d" / f"{level}.msh")
        lengths = compute_edge_lengths(mesh.nodes, mesh.elements)
        ax1.hist(np.log10(lengths), bins=60, histtype="step", linewidth=2.0,
                 color=color,
                 label=f"{level} ({len(mesh.elements)/1e3:.0f}k tets)")
        rows.append(("naca0012 " + level, len(mesh.elements),
                     f"{lengths.min():.4f}", f"{lengths.max():.3f}", ""))
    ax1.set_xlabel("log10(edge length / chord)")
    ax1.set_ylabel("edge count")
    ax1.set_title("NACA0012: graded, wall-clustered edges")
    ax1.legend(loc="upper right")

    min_volume = np.inf
    sphere_aspect_max = 0.0
    for (family, level), color in ((("naca0012_2.5d", "coarse"), S1_BLUE),
                                   (("cylinder_2.5d", "coarse"), S3_YELLOW),
                                   (("sphere_shell", "coarse"), S5_VIOLET)):
        mesh = read_mesh(MESH_DIR / family / f"{level}.msh")
        aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
        volumes = compute_tet_volumes(mesh.nodes, mesh.elements)
        min_volume = min(min_volume, float(volumes.min()))
        if family == "sphere_shell":
            sphere_aspect_max = float(aspect.max())
        ax2.hist(np.log10(aspect), bins=60, histtype="step", linewidth=2.0,
                 color=color, label=f"{family} {level} (max {aspect.max():.1f})")
        rows.append((f"{family} {level}", len(mesh.elements), "", "",
                     f"{aspect.max():.2f}"))
    ax2.set_xlabel("log10(tet aspect ratio)")
    ax2.set_ylabel("element count")
    ax2.set_yscale("log")
    ax2.set_title("Anisotropy sits where the design puts it\n"
                  "(quasi-2D far field: fixed dz by design; sphere < 4)")
    ax2.legend(loc="upper right", fontsize=9)
    fig.suptitle("M0 mesh quality", fontsize=12.5, fontweight="semibold",
                 color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    finish(fig, OUT, "mesh_quality.png")
    write_csv(OUT, "mesh_quality.csv",
              "mesh,n_tets,min_edge,max_edge,max_aspect", rows)

    checks.add("M0", "min tet volume (3 families)", f"{min_volume:.1e}",
               "> 0 (no degenerate/inverted tets)", bool(min_volume > 0.0))
    checks.add("M0", "isotropic sphere family aspect max", f"{sphere_aspect_max:.1f}",
               "< 5 (matches sphere_stats 3.5)", bool(sphere_aspect_max < 5.0),
               note="quasi-2D families are anisotropic in the far field by "
                    "design (fixed dz)")


# ---------------------------------------------------------------------------
# 4. cylinder Cp convergence: the meshes solve, and errors shrink
# ---------------------------------------------------------------------------
def demo_cylinder_convergence(checks):
    print("[4/4] cylinder Cp convergence (coarse/medium/fine)")
    levels = ("coarse", "medium", "fine")
    h = {"coarse": 0.1, "medium": 0.05, "fine": 0.025}  # h_wall from stats
    cases, rows = {}, []
    for level in levels:
        case = run_cylinder_case(MESH_DIR / "cylinder_2.5d" / f"{level}.msh")
        cases[level] = case
        rows.append((level, len(case["mesh"].elements), h[level],
                     f"{case['error'].max():.4e}", f"{case['error'].mean():.4e}",
                     f"{case['residual_norm']:.2e}"))
        print(f"    {level}: max|Cp err| = {case['error'].max():.4f}, "
              f"residual = {case['residual_norm']:.1e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8))
    tl = np.linspace(0, 180, 300)
    ax1.plot(tl, 1.0 - 4.0 * np.sin(np.radians(tl)) ** 2, "-", color=INK,
             linewidth=2, label="exact: 1 - 4 sin$^2$theta", zorder=4)
    for level, color in (("coarse", S1_BLUE), ("medium", S2_AQUA),
                         ("fine", S3_YELLOW)):
        case = cases[level]
        nodes = case["mesh"].nodes[case["wall_nodes"]]
        theta = np.degrees(np.arctan2(np.abs(nodes[:, 1]), nodes[:, 0]))
        order = np.argsort(theta)
        ax1.plot(theta[order], case["cp_numeric"][order], ".", color=color,
                 markersize=3, alpha=0.5, label=f"{level}")
    ax1.set_xlabel("theta (deg from +x stagnation point)")
    ax1.set_ylabel("Cp")
    ax1.set_title("Wall Cp closes on the analytic curve under refinement")
    ax1.legend(loc="upper center", fontsize=9)

    hv = np.array([h[lv] for lv in levels])
    ev = np.array([cases[lv]["error"].max() for lv in levels])
    slope = np.polyfit(np.log(hv), np.log(ev), 1)[0]
    ax2.loglog(hv, ev, "o-", color=S1_BLUE, markeredgecolor=SURFACE,
               markeredgewidth=2, label="max |Cp error| at wall")
    ref = ev[0] * hv / hv[0]
    ax2.loglog(hv, ref, "--", color=MUTED, linewidth=1.2)
    ax2.text(hv[1], ref[1] * 1.18, "O(h) reference", color=MUTED, fontsize=9)
    ax2.text(0.97, 0.06, f"measured slope = {slope:.2f}\n"
             "(curved-wall recovery limit, same\nmechanism as G1.6 -- documented)",
             transform=ax2.transAxes, color=INK_2, fontsize=9.5,
             ha="right", va="bottom")
    ax2.set_xlabel("wall spacing h")
    ax2.set_ylabel("max |Cp error|")
    ax2.set_title("Errors shrink monotonically at the expected order")
    ax2.legend(loc="upper left", fontsize=9)
    ax2.grid(True, which="both")
    fig.suptitle("M0 meshes are solvable and refinable: cylinder "
                 "quasi-2D end-to-end", fontsize=12.5,
                 fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    finish(fig, OUT, "cylinder_cp_convergence.png")
    write_csv(OUT, "cylinder_convergence.csv",
              "level,n_tets,h_wall,max_cp_err,mean_cp_err,residual_norm", rows)

    checks.add("M0", "cylinder solve residual (all levels)",
               f"{max(cases[lv]['residual_norm'] for lv in levels):.1e}",
               "< 1e-8", bool(max(cases[lv]["residual_norm"]
                                  for lv in levels) < 1e-8))
    mono = all(cases[b]["error"].max() < 0.75 * cases[a]["error"].max()
               for a, b in (("coarse", "medium"), ("medium", "fine")))
    checks.add("M0", "max |Cp err| drops >= 25% per level", mono,
               "monotone refinement", bool(mono))
    checks.add("M0", "Cp error slope vs h", f"{slope:.2f}",
               "in [0.85, 1.25] (~O(h), documented limit)",
               bool(0.85 <= slope <= 1.25))


def main():
    apply_style()
    t0 = time.time()
    checks = CheckList("M0 quasi-2D meshing pipeline (closed; acceptance "
                       "link = G2.5 in the P2 demo)")

    demo_gallery()
    demo_topology(checks)
    demo_quality(checks)
    demo_cylinder_convergence(checks)

    write_csv(OUT, "summary.csv", "metric,value",
              [(c["gate"] + " " + c["name"].replace(",", ";"),
                str(c["value"]) + " [" + c["status"] + "]")
               for c in checks.checks] + [("runtime_seconds", f"{time.time()-t0:.1f}")])
    code = checks.report(OUT)
    print(f"done in {time.time() - t0:.1f}s -> {OUT.relative_to(REPO_ROOT)}/")
    sys.exit(code)


if __name__ == "__main__":
    main()
