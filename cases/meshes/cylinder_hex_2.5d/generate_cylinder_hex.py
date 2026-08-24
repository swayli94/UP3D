"""Structured body-fitted O-grid cylinder, 2.5-D (phase 3 task 2, route A).

The structured counterpart of `cases/meshes/cylinder_2.5d/` (gmsh unstructured annulus).
SAME geometry, SAME extrusion machinery, SAME output format -- so the two families give a
clean same-geometry A/B, which is why the cylinder is the first leg: it has an EXACT
potential solution, and the unstructured baseline is already measured (uncorrected wall Cp
error converging at ~1.0 order, tests/C/test_C02_cylinder_wall_correction.py).

Pipeline: `structured.cylinder_o_grid_2d` (analytic, no gmsh size field) -> each quad cut on
a fixed diagonal into 2 triangles -> `extrude_single_layer` (prism -> 3 tets, globally
consistent orientation, already locked by tests) => 6 tets per hex.

★ The ladder refines the WALL-NORMAL and WALL-TANGENTIAL resolutions together by 2x per
level, with the grading rate and the domain radius HELD FIXED. That is only meaningful
because the construction makes those four knobs independent -- see the module docstring of
pyfp3d/meshgen/structured.py and the G0 evidence in bench/gate_results/hex_g0_single_var.csv.

    python cases/meshes/cylinder_hex_2.5d/generate_cylinder_hex.py --all

Outputs: <level>.msh (gitignored, like every other mesh family) + <level>_stats.csv (tracked).
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import mesh_stats, write_mesh          # noqa: E402
from pyfp3d.meshgen.extrude import extrude_single_layer        # noqa: E402
from pyfp3d.meshgen.structured import cylinder_o_grid_2d       # noqa: E402

OUT_DIR = Path(__file__).resolve().parent
RADIUS = 1.0
#: r_far and growth are FIXED across the ladder on purpose: the refinement axis is wall
#: resolution, and holding the other two lets the convergence order be attributed to it.
R_FAR, GROWTH = 20.0, 1.15
LEVELS = {
    "coarse": dict(n_theta=48, h_wall_normal=0.04, dz=0.20),
    "medium": dict(n_theta=96, h_wall_normal=0.02, dz=0.10),
    "fine": dict(n_theta=192, h_wall_normal=0.01, dz=0.05),
}


def build(level: str):
    p = LEVELS[level]
    pts, tris, edges, info = cylinder_o_grid_2d(
        radius=RADIUS, r_far=R_FAR, n_theta=p["n_theta"],
        h_wall_normal=p["h_wall_normal"], growth=GROWTH)
    mesh = extrude_single_layer(pts, tris, edges, dz=p["dz"],
                               name=f"cylinder_hex_{level}")
    return mesh, info, (pts, tris, edges)


def wall_facet_normal_error(mesh) -> dict:
    """max/mean angle (deg) between each wall facet normal and the EXACT cylinder normal.

    ★ RECORDED with no criterion, per the pre-registration's risk 2: an O-grid wall is still
    FLAT FACETS, so this O(h) error is NOT expected to improve over the unstructured family.
    It is measured so that nobody later counts it as a gain of route (A)."""
    wall = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    p = mesh.nodes[wall]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    c = p.mean(axis=1)
    exact = c.copy(); exact[:, 2] = 0.0
    exact /= np.linalg.norm(exact, axis=1, keepdims=True)
    cos = np.abs(np.einsum("ij,ij->i", n, exact))
    ang = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    return dict(normal_err_max_deg=float(ang.max()),
                normal_err_mean_deg=float(ang.mean()))


def write_inspection_png(path, points2d, triangles, edge_groups, level, r_far):
    """Headless artifact (hard rule: every visual gate needs one; matplotlib Agg).

    ★ Added 2026-08-11 after the user pointed out this family shipped WITHOUT one while every
    sibling family has its `*_layer.png`. Two panels: the full domain, and a near-wall zoom --
    the zoom is the one that matters here, because the whole claim of route (A) is about the
    ORTHOGONAL graded near-wall layer, and that is invisible at domain scale.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    colors = {"wall": "tab:red", "farfield": "tab:blue"}
    for ax, lim in zip(axes, [r_far * 1.05, 1.35]):
        ax.triplot(points2d[:, 0], points2d[:, 1], triangles,
                   linewidth=0.25, color="0.6")
        for tag, edges in edge_groups.items():
            for seg in points2d[np.asarray(edges)]:
                ax.plot(seg[:, 0], seg[:, 1], color=colors.get(tag, "k"), linewidth=1.2)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    axes[0].set_title(f"cylinder_hex_2.5d {level}: structured O-grid + tags")
    axes[1].set_title("near-wall zoom -- the orthogonal graded layer (wall=red)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", action="append", choices=sorted(LEVELS))
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    levels = sorted(LEVELS) if a.all else (a.level or ["coarse", "medium"])
    for level in levels:
        mesh, info, planar = build(level)
        path = OUT_DIR / f"{level}.msh"
        write_mesh(mesh, path)
        st = mesh_stats(mesh)
        #: built by successive update() rather than one dict(**a, **b, ...) call: `info` and
        #: `mesh_stats` both carry keys this row also sets (growth, h_wall_normal, n_nodes,
        #: n_tets), and dict() raises on a duplicate keyword -- which is how the first two runs
        #: failed. mesh_stats keys are PREFIXED so a collision cannot silently overwrite.
        row = {"level": level}
        row.update({k: info[k] for k in sorted(info)})
        row.update(dz=LEVELS[level]["dz"], r_far_requested=R_FAR,
                   n_nodes=len(mesh.nodes), n_tets=len(mesh.elements),
                   n_tris_wall=len(mesh.boundary_faces["wall"]),
                   n_tris_farfield=len(mesh.boundary_faces["farfield"]))
        row.update(wall_facet_normal_error(mesh))
        row.update({f"stats_{k}": v for k, v in st.items()
                    if isinstance(v, (int, float))})
        write_inspection_png(OUT_DIR / f"{level}_layer.png", planar[0], planar[1],
                             planar[2], level, R_FAR)
        with open(OUT_DIR / f"{level}_stats.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row)); w.writeheader(); w.writerow(row)
        print(f"  {level:7} n_theta {info['n_theta']:>4} n_radial {info['n_radial']:>3} "
              f"nodes {len(mesh.nodes):>7} tets {len(mesh.elements):>8} "
              f"wall tris {row['n_tris_wall']:>5}  h_t {info['h_wall_tangential']:.6f} "
              f"h_n {info['h_wall_normal']:.4f}  normal err max "
              f"{row['normal_err_max_deg']:.4f} deg  -> {path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
