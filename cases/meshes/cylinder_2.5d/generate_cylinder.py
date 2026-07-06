"""
Quasi-2D circular-cylinder flow test case (M0 pipeline demonstrator).

A single-layer extruded annulus around a unit-radius cylinder: the 2D
circle-in-circle domain is meshed with vanilla Gmsh (graded from the wall),
then extruded exactly one cell layer in z with the globally consistent
prism -> 3-tet split (pyfp3d/meshgen/extrude.py). No wake: incompressible
cylinder flow is non-lifting, and the analytic solution

    phi = U x (1 + a^2 / r^2),   surface Cp = 1 - 4 sin^2(theta)

makes this the simplest end-to-end validation of the quasi-2D meshing
pipeline + Laplace solver (tests/test_m0_cylinder.py), complementing the
lifting NACA0012 case which is the actual M0 deliverable.

Tags: wall (cylinder), farfield (outer circle), symmetry (both z-planes).

Usage:
    python generate_cylinder.py --level coarse --level medium
    python generate_cylinder.py --all
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.meshgen.extrude import extrude_single_layer
from pyfp3d.meshgen.planar import cylinder_annulus_2d

RADIUS = 1.0
R_FAR = 20.0

LEVELS = {
    # h_wall: wall-adjacent 2D size; h_far: bulk size; dz: layer thickness
    "coarse": dict(h_wall=0.10, h_far=2.5, dist_max=10.0, dz=0.20),
    "medium": dict(h_wall=0.05, h_far=2.0, dist_max=10.0, dz=0.10),
    "fine": dict(h_wall=0.025, h_far=1.5, dist_max=10.0, dz=0.05),
}


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    points2d, triangles, edge_groups = cylinder_annulus_2d(
        radius=RADIUS, r_far=R_FAR,
        h_wall=p["h_wall"], h_far=p["h_far"], dist_max=p["dist_max"],
    )
    mesh = extrude_single_layer(
        points2d, triangles, edge_groups, dz=p["dz"], z0=0.0,
        name=f"cylinder_2.5d_{level}",
    )

    out_path = out_dir / f"{level}.msh"
    write_mesh(mesh, out_path)

    stats = mesh_stats(mesh)
    stats["n_2d_triangles"] = len(triangles)
    stats["dz"] = p["dz"]
    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    if inspect:
        _write_inspection_png(out_dir / f"{level}_layer.png", points2d,
                              triangles, edge_groups, level)
    return out_path


def _write_inspection_png(path: Path, points2d, triangles, edge_groups, level):
    """Headless artifact (roadmap Sec 0.1): 2D layer triangulation with tagged
    boundaries, full domain + near-wall zoom."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    colors = {"wall": "tab:red", "farfield": "tab:blue"}
    for ax, lim in zip(axes, [R_FAR * 1.05, 3.0]):
        ax.triplot(points2d[:, 0], points2d[:, 1], triangles,
                   linewidth=0.25, color="0.6")
        for tag, edges in edge_groups.items():
            seg = points2d[np.asarray(edges)]
            for s in seg:
                ax.plot(s[:, 0], s[:, 1], color=colors.get(tag, "k"),
                        linewidth=1.2)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
    axes[0].set_title(f"cylinder_2.5d {level}: 2D layer + tags")
    axes[1].set_title("near-wall zoom (wall=red, farfield=blue)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", action="append", choices=sorted(LEVELS))
    parser.add_argument("--all", action="store_true", help="generate all levels")
    args = parser.parse_args()
    levels = sorted(LEVELS) if args.all else (args.level or ["coarse", "medium"])

    out_dir = Path(__file__).parent
    for level in levels:
        out_path = generate_level(out_dir, level)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
