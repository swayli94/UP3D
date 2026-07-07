"""
M1 demo -- swept/tapered wing meshing evidence (Track M).

M1 delivers pyfp3d/meshgen/wing3d.py (OCC ruled-loft ONERA M6 half wing,
spherical far field ~15 MAC, planar wake sheet swept from the sharp TE
and embedded conformally via occ.fragment + mesh.embed) plus the
committed mesh family cases/meshes/onera_m6 (coarse/medium; fine is
regenerable on demand, stats committed) and the wake_cut.py extensions
that make the solver preprocessor ingest a 3D swept wake: per-node TE
stations, single-valued sheet FREE-edge nodes at the wing tip
(Gamma(tip) = 0 discretely) and the off-plane Kutta-probe fallback.

This demo shows the mesh-side evidence on the COMMITTED meshes:

  1. What the mesh looks like: wake sheet + wing wireframe, tip zoom
     (wake-tip closure) and cut planes just in/outboard of the tip.
  2. Ingestion: cut_wake + the P2 topology asserts pass on both levels;
     station/free-edge/Kutta-probe semantics are the M1-specific ones.
  3. Element quality (M1 gate: min dihedral, max aspect within bounds)
     across the refinement ladder, fine level from its committed stats.
  4. G2.1-style freestream preservation on the CUT coarse mesh.

Standalone + self-checking:  python cases/demo/m1_wing_mesh/run_demo.py
Outputs: cases/demo/m1_wing_mesh/results/{*.png, summary.csv, checks.csv}
Exit code 0 iff every acceptance check passes. Runtime ~1 min.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    GRID, INK, MESH_DIR, REPO_ROOT, S1_BLUE, S3_YELLOW, S5_VIOLET,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402

from pyfp3d.kernels.residual import assemble_stiffness_matrix  # noqa: E402
from pyfp3d.constraints.wake import WakeConstraint  # noqa: E402
from pyfp3d.mesh.metrics import (  # noqa: E402
    compute_aspect_ratios, compute_min_dihedral_angles,
)
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, x_te  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
M6_DIR = MESH_DIR / "onera_m6"
LEVELS = ["coarse", "medium"]
BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}


# ---------------------------------------------------------------------------
# 1. mesh gallery: wake + tip wireframe, cut planes at the tip
# ---------------------------------------------------------------------------
def _tri_edges(tris):
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    return np.unique(np.sort(e, axis=1), axis=0)


def _slice_segments(nodes, faces, z0):
    p = nodes[faces]
    d = p[:, :, 2] - z0
    cross = (d.min(axis=1) < 0) & (d.max(axis=1) > 0)
    p, d = p[cross], d[cross]
    segs = []
    for pk, dk in zip(p, d):
        pts = [pk[a] + dk[a] / (dk[a] - dk[b]) * (pk[b] - pk[a])
               for a, b in ((0, 1), (1, 2), (2, 0)) if (dk[a] < 0) != (dk[b] < 0)]
        if len(pts) == 2:
            segs.append([pts[0][:2], pts[1][:2]])
    return np.asarray(segs) if segs else np.empty((0, 2, 2))


def demo_gallery(mesh):
    print("[1/4] wing + wake gallery (coarse)")
    wall = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)

    fig = plt.figure(figsize=(14.4, 5.4))
    for i, (title, win) in enumerate((
        ("wake sheet swept TE -> far field", None),
        ("tip zoom: wake-tip closure", (x_te(B_SEMI) - 0.35, x_te(B_SEMI) + 0.55,
                                        B_SEMI - 0.45, B_SEMI + 0.15)),
    ), start=1):
        ax = fig.add_subplot(1, 2, i, projection="3d")
        for tris, color, lw in ((wall, INK, 0.3), (wake, S5_VIOLET, 0.25)):
            e = _tri_edges(tris)
            if win is not None:
                x0, x1, z0, z1 = win
                keep = ((mesh.nodes[e, 0].min(axis=1) > x0 - 0.1)
                        & (mesh.nodes[e, 0].max(axis=1) < x1 + 0.1)
                        & (mesh.nodes[e, 2].min(axis=1) > z0 - 0.1)
                        & (mesh.nodes[e, 2].max(axis=1) < z1 + 0.1))
                e = e[keep]
            ax.add_collection3d(
                Line3DCollection(mesh.nodes[e], colors=color, linewidths=lw))
        if win is None:
            ax.set_xlim(-0.5, 3.5); ax.set_ylim(-2, 2); ax.set_zlim(-0.5, 3.5)
        else:
            x0, x1, z0, z1 = win
            ax.set_xlim(x0, x1)
            ax.set_ylim(-0.5 * (x1 - x0), 0.5 * (x1 - x0))
            ax.set_zlim(z0, z1)
        ax.view_init(elev=28, azim=-125)
        ax.set_title(title)
    fig.suptitle("M1 ONERA M6 half wing: ruled loft, chord-plane wake sheet "
                 "(wall=dark, wake=violet; coarse level)",
                 fontsize=12.5, fontweight="semibold", color=INK)
    finish(fig, OUT, "wing_wake_gallery.png")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    el = np.asarray(mesh.elements, dtype=np.int64)
    for ax, z0 in zip(axes, (B_SEMI - 0.05, B_SEMI + 0.05)):
        crossing = el[(mesh.nodes[el, 2].min(axis=1) < z0)
                      & (mesh.nodes[el, 2].max(axis=1) > z0)]
        faces = np.concatenate([crossing[:, [0, 1, 2]], crossing[:, [0, 1, 3]],
                                crossing[:, [0, 2, 3]], crossing[:, [1, 2, 3]]])
        faces = np.unique(np.sort(faces, axis=1), axis=0)
        from matplotlib.collections import LineCollection
        ax.add_collection(LineCollection(
            _slice_segments(mesh.nodes, faces, z0), colors=GRID, linewidths=0.35))
        ws = _slice_segments(mesh.nodes, wake, z0)
        if len(ws):
            ax.add_collection(LineCollection(ws, colors=S5_VIOLET, linewidths=1.8))
        ax.set_xlim(x_te(B_SEMI) - 0.9, x_te(B_SEMI) + 1.1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_aspect("equal")
        ax.grid(False)
        side = "inboard" if z0 < B_SEMI else "outboard"
        ax.set_title(f"z = {z0:.3f} ({side} of tip)")
    fig.suptitle("M1 visual gate: valid volume fill at the tip; the wake "
                 "slice (violet) exists inboard and ends AT the tip",
                 fontsize=12.5, fontweight="semibold", color=INK)
    finish(fig, OUT, "tip_cut_planes.png")


# ---------------------------------------------------------------------------
# 2. ingestion: cut_wake topology + M1 station/free-edge semantics
# ---------------------------------------------------------------------------
def demo_ingestion(checks, meshes):
    print("[2/4] cut_wake ingestion on the committed family")
    rows = []
    cut_coarse = None
    for level in LEVELS:
        mesh = meshes[level]
        t0 = time.time()
        mesh_cut, wc = cut_wake(mesh)  # topology asserts run inside
        dt = time.time() - t0
        rows.append((level, len(mesh.elements), len(wc.slave_nodes),
                     len(wc.free_nodes), len(wc.te_nodes), wc.n_stations,
                     f"{dt:.1f}"))
        if level == "coarse":
            cut_coarse = (mesh_cut, wc)

        ok_tags = set(mesh.boundary_faces) == {"wall", "farfield",
                                               "symmetry", "wake"}
        checks.add("M1", f"{level}: tags wall/farfield/symmetry/wake",
                   sorted(mesh.boundary_faces), "complete", ok_tags)
        checks.add("M1", f"{level}: P2 topology asserts through cut_wake",
                   "pass", "no assert fires", True)
        checks.add("M1", f"{level}: per-node TE stations (swept TE)",
                   wc.n_stations, "== n_te_nodes",
                   wc.n_stations == len(wc.te_nodes))
        tip_ok = (len(wc.free_nodes) > 0
                  and np.abs(mesh.nodes[wc.free_nodes, 2] - B_SEMI).max() < 1e-6
                  and not np.isin(wc.free_nodes, wc.master_nodes).any())
        checks.add("M1", f"{level}: tip free edge single-valued",
                   len(wc.free_nodes), "all at z=B_SEMI, none duplicated",
                   bool(tip_ok))
        probes_ok = (np.all(mesh.nodes[wc.kutta_upper, 1] > 0)
                     and np.all(mesh.nodes[wc.kutta_lower, 1] < 0))
        checks.add("M1", f"{level}: Kutta probes found both sides",
                   f"{len(wc.kutta_upper)} pairs", "y>0 upper / y<0 lower",
                   bool(probes_ok))
    write_csv(OUT, "ingestion.csv",
              "level,n_tets,n_dup,n_free,n_te,n_stations,cut_seconds", rows)
    return cut_coarse


# ---------------------------------------------------------------------------
# 3. element quality across the refinement ladder
# ---------------------------------------------------------------------------
def demo_quality(checks, meshes):
    print("[3/4] element quality (M1 gate: min dihedral / max aspect)")
    rows = []
    for level in LEVELS:
        mesh = meshes[level]
        dih = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
        asp = compute_aspect_ratios(mesh.nodes, mesh.elements)
        rows.append((level, len(mesh.elements), f"{dih.min():.2f}",
                     f"{asp.max():.2f}", "recomputed"))

    # Fine level: 2.5M tets, .msh regenerable on demand; quality numbers
    # come from its committed generation-time stats.
    fine_stats = {}
    with open(M6_DIR / "fine_stats.csv") as f:
        for row in csv.DictReader(f):
            fine_stats[row["metric"]] = row["value"]
    rows.append(("fine", int(fine_stats["n_elements"]),
                 f"{float(fine_stats['min_dihedral_deg']):.2f}",
                 f"{float(fine_stats['max_aspect_ratio']):.2f}",
                 "committed fine_stats.csv"))
    write_csv(OUT, "mesh_quality.csv",
              "level,n_tets,min_dihedral_deg,max_aspect,source", rows)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    x = np.arange(len(rows))
    names = [r[0] for r in rows]
    for ax, vals, bound, label, better in (
        (axes[0], [float(r[2]) for r in rows], BOUNDS["min_dihedral_deg"],
         "min dihedral angle [deg]", "higher is better"),
        (axes[1], [float(r[3]) for r in rows], BOUNDS["max_aspect_ratio"],
         "max aspect ratio", "lower is better"),
    ):
        ax.bar(x, vals, 0.55, color=S1_BLUE)
        ax.axhline(bound, color=S3_YELLOW, linewidth=1.6,
                   label=f"gate bound {bound}")
        ax.set_xticks(x, names)
        ax.set_title(f"{label} ({better})")
        ax.legend()
    axes[1].set_yscale("log")
    fig.suptitle("M1 element quality across the refinement ladder "
                 "(55k / 351k / 2513k tets)",
                 fontsize=12.5, fontweight="semibold", color=INK)
    finish(fig, OUT, "mesh_quality.png")

    ok = all(float(r[2]) >= BOUNDS["min_dihedral_deg"]
             and float(r[3]) <= BOUNDS["max_aspect_ratio"] for r in rows)
    worst = min(float(r[2]) for r in rows)
    checks.add("M1", "quality within bounds on all 3 levels",
               f"worst min-dihedral {worst:.2f} deg",
               f">= {BOUNDS['min_dihedral_deg']} deg, aspect <= "
               f"{BOUNDS['max_aspect_ratio']}", ok)
    n = [r[1] for r in rows]
    checks.add("M1", "refinement ladder monotone (one h parameter)",
               f"{n[0]}/{n[1]}/{n[2]}", "n_tets grows ~2^3 per level",
               n[0] * 3 < n[1] < n[2] and n[1] * 3 < n[2])


# ---------------------------------------------------------------------------
# 4. freestream preservation on the CUT mesh (G2.1 analogue)
# ---------------------------------------------------------------------------
def demo_freestream(checks, cut_coarse):
    print("[4/4] freestream preservation on the cut coarse mesh")
    mesh_cut, wc = cut_coarse
    A = assemble_stiffness_matrix(mesh_cut.nodes, mesh_cut.elements)
    con = WakeConstraint(A, wc)
    R = con.A_reduced @ mesh_cut.nodes[: con.n_reduced, 0]
    to_master = np.arange(len(mesh_cut.nodes))
    to_master[wc.slave_nodes] = wc.master_nodes
    check = np.ones(con.n_reduced, dtype=bool)
    for tag in ("wall", "farfield", "symmetry"):
        check[np.unique(to_master[np.unique(mesh_cut.boundary_faces[tag])])] = False
    r_max = float(np.max(np.abs(R[check])))
    checks.add("M1", "freestream residual on cut mesh (interior rows)",
               f"{r_max:.1e}", "< 1e-10 (G2.1 analogue)", r_max < 1e-10)


def main():
    apply_style()
    t0 = time.time()
    checks = CheckList("M1 swept-wing meshing pipeline (ONERA M6; consumed "
                       "by P5)")

    meshes = {level: read_mesh(M6_DIR / f"{level}.msh") for level in LEVELS}
    demo_gallery(meshes["coarse"])
    cut_coarse = demo_ingestion(checks, meshes)
    demo_quality(checks, meshes)
    demo_freestream(checks, cut_coarse)

    write_csv(OUT, "summary.csv", "metric,value",
              [(c["gate"] + " " + c["name"].replace(",", ";"),
                str(c["value"]) + " [" + c["status"] + "]")
               for c in checks.checks] + [("runtime_seconds", f"{time.time()-t0:.1f}")])
    code = checks.report(OUT)
    print(f"done in {time.time() - t0:.1f}s -> {OUT.relative_to(REPO_ROOT)}/")
    sys.exit(code)


if __name__ == "__main__":
    main()
