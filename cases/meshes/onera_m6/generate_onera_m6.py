"""
M1 deliverable: ONERA M6 half-wing tet mesh family with embedded wake
sheet (docs/roadmap.md Track M, gate M1; consumed by solver phase P5).

Geometry, axis convention and wake-sheet construction live in
pyfp3d/meshgen/wing3d.py (vanilla Gmsh/OCC; node duplication stays in the
solver preprocessor, hard rule 8). This script owns the refinement-level
family, the element-quality report and the headless inspection artifacts:

  <level>.msh            solver-ready mesh (tags: wall, farfield,
                         symmetry, wake)
  <level>_stats.csv      size + quality report (M1 gate: min dihedral,
                         max aspect within QUALITY_BOUNDS)
  <level>_wake_tip.png   wireframe of wake sheet + wing tip (M1 visual
                         check: wake-tip closure, no cracks)
  <level>_cutplane.png   z = const cut planes just inboard/outboard of
                         the tip (M1 visual check: valid volume fill)

One parameter per level (h = h_wall on the wing surface, 2x ladder
0.030 / 0.015 / 0.0075 m); everything else scales with it. Sizes are
runtime-driven: the 2.5D P4 gates showed solver wall time (not memory)
is the binding constraint, so `coarse` is the day-to-day P5 development
mesh and `medium` the gate mesh. Measured (2026-07-07):

    level   tets     nodes   gen time   cut_wake ingest
    coarse  55.5 k   11.0 k     ~7 s        ~2 s
    medium  350.7 k  63.2 k    ~23 s       ~12 s
    fine    2.513 M  428 k    ~264 s       ~94 s

No .msh file is committed (cases/meshes/onera_m6/*.msh is gitignored --
they are large; coarse+medium regenerate in ~30 s, fine in ~4.5 min).
The per-level stats CSVs and inspection PNGs ARE committed as the M1
evidence record. Note the default-suite implications: the 13 tests in
tests/test_m1_onera_m6.py skip when the meshes are absent, and the
hard-rule-7 glob sweep (test_p2_wake_cut.py) ingests whatever M6 meshes
exist locally -- fine.msh alone adds ~94 s there, another reason to
regenerate it only for mesh-convergence studies.

Usage:
    python generate_onera_m6.py --level coarse
    python generate_onera_m6.py --all          # includes fine
(default: coarse + medium, same as the M0 family)
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.meshgen.wing3d import (
    B_SEMI, MAC, onera_m6_wing_mesh, x_te,
)

R_FAR = 15.0 * MAC

# M1 gate: element quality "within bounds". Bounds chosen from what
# gmsh Delaunay + optimization reliably delivers on this geometry; a
# violation means the sizing fields or the geometry stitching regressed.
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}


def _level_params(h_wall: float, clamp_h_far: bool = True) -> dict:
    """Everything derived from the single wall-size parameter h.

    ★ THE h_far CLAMP BREAKS SELF-SIMILARITY (found 2026-07-13, P13/G13.3).
    Every length scales with h_wall EXCEPT the far field, which is capped:

        h_far = min(2.5, 120 * h_wall)

    and that cap BITES ONLY AT `coarse` (120*0.030 = 3.6 > 2.5), never at
    medium (1.8) or fine (0.9). So coarse->medium refines the far field by
    just 1.39x while the wall refines by 2x, whereas medium->fine is a clean
    2x throughout: **`coarse` is NOT on the same refinement ray as the other
    two.** Consequences, measured:

      * Any THREE-POINT RICHARDSON over (coarse, medium, fine) is INVALID --
        including the one P9/G9.1 attempted. Two of the points are not
        members of the same family.
      * It shows up as a spanwise loading "slosh" coarse->medium (-6.3% at the
        root against +5.4% mid-span, nearly cancelling in the total) that
        looks like convergence noise but is really a comparison between two
        different meshes families.

    The clamp is KEPT BY DEFAULT so that `coarse`/`medium`/`fine` stay
    BIT-IDENTICAL -- they are the shipped dev/gate meshes and carry the
    regression locks of P5, P8/G8.2, B7 and the M1 asserts. For grid
    convergence use RICHARDSON_LADDER below, which swaps in an UNCLAMPED
    coarse ("coarse_ss") and reuses medium/fine unchanged (they already
    satisfy h_far = 120*h_wall exactly), giving a family that is self-similar
    by a factor 2 in EVERY length scale.
    """
    h_far = 120.0 * h_wall
    return dict(
        h_wall=h_wall,
        h_wake=3.0 * h_wall,
        h_edge=0.5 * h_wall,
        h_far=min(2.5, h_far) if clamp_h_far else h_far,
    )


LEVELS = {
    # shipped dev/gate family (clamp retained => bit-identical, locks intact)
    "coarse": _level_params(0.030),
    "medium": _level_params(0.015),
    "fine": _level_params(0.0075),
    # self-similar coarse: the ONLY level the clamp ever touched, regenerated
    # without it (h_far 2.5 -> 3.6) so that {coarse_ss, medium, fine} is a
    # legitimate 2x refinement ladder in every length scale.
    "coarse_ss": _level_params(0.030, clamp_h_far=False),
}

#: The only member set valid for a three-point grid-convergence / Richardson
#: study (P13/G13.3). medium and fine are the SAME meshes as in the shipped
#: family -- only the coarse end had to be re-cut.
RICHARDSON_LADDER = ("coarse_ss", "medium", "fine")


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    t0 = time.perf_counter()
    mesh = onera_m6_wing_mesh(
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        h_edge=p["h_edge"], r_far=R_FAR, name=f"onera_m6_{level}",
    )
    gen_seconds = time.perf_counter() - t0

    out_path = out_dir / f"{level}.msh"
    write_mesh(mesh, out_path)

    stats = mesh_stats(mesh)
    dihedral = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
    aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
    stats["min_dihedral_deg"] = float(dihedral.min())
    stats["p01_dihedral_deg"] = float(np.percentile(dihedral, 1))
    stats["max_aspect_ratio"] = float(aspect.max())
    stats["p99_aspect_ratio"] = float(np.percentile(aspect, 99))
    for tag, faces in mesh.boundary_faces.items():
        stats[f"n_{tag}_faces"] = len(faces)
    stats["h_wall"] = p["h_wall"]
    stats["gen_seconds"] = round(gen_seconds, 1)

    ok = (stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"]
          and stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"])
    stats["quality_within_bounds"] = ok

    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    _check_wake_tip_closure(mesh)

    if inspect:
        _write_wake_tip_png(out_dir / f"{level}_wake_tip.png", mesh, level)
        _write_cutplane_png(out_dir / f"{level}_cutplane.png", mesh, level)

    if not ok:
        raise AssertionError(
            f"{level}: element quality outside bounds "
            f"(min dihedral {stats['min_dihedral_deg']:.2f} deg, "
            f"max aspect {stats['max_aspect_ratio']:.1f}; "
            f"bounds {QUALITY_BOUNDS})"
        )
    return out_path


def _check_wake_tip_closure(mesh) -> None:
    """Programmatic M1 wake-tip closure check (the PNG is for eyeballing).

    The sheet's tip edge (wake-sheet boundary edges at z = B_SEMI) must
    form one connected open chain that starts at the wing tip TE corner --
    a node shared EXACTLY between wall and wake -- i.e. no crack and no
    dangling pieces between the sheet and the tip.
    """
    wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)
    wall_nodes = np.unique(mesh.boundary_faces["wall"])

    # Wake-sheet boundary edges = edges used by exactly one wake triangle.
    edges = np.sort(np.concatenate([wake[:, [0, 1]], wake[:, [1, 2]],
                                    wake[:, [2, 0]]]), axis=1)
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = uniq[counts == 1]

    z = mesh.nodes[:, 2]
    at_tip = np.abs(z - B_SEMI) < 1e-6
    tip_edges = boundary_edges[at_tip[boundary_edges].all(axis=1)]
    assert len(tip_edges) > 0, "wake sheet has no tip edge at z = B_SEMI"

    # Chain check: every node in <= 2 tip edges, exactly two endpoints.
    nodes, deg = np.unique(tip_edges, return_counts=True)
    assert deg.max() <= 2, "wake tip edge self-intersects (node degree > 2)"
    endpoints = nodes[deg == 1]
    assert len(endpoints) == 2, (
        f"wake tip edge is not one open chain ({len(endpoints)} endpoints)"
    )

    # One endpoint is the tip TE corner (on the wall), and it sits at the
    # exact planform TE tip -- shared geometry, not merely nearby.
    on_wall = np.isin(endpoints, wall_nodes)
    assert on_wall.sum() == 1, (
        "wake tip edge does not attach to the wall at exactly one end"
    )
    corner = endpoints[on_wall][0]
    xy = mesh.nodes[corner]
    assert abs(xy[0] - x_te(B_SEMI)) < 1e-9 and abs(xy[1]) < 1e-9, (
        f"tip TE corner at {xy}, expected ({x_te(B_SEMI):.6f}, 0, {B_SEMI})"
    )


def _write_wake_tip_png(path, mesh, level):
    """Headless artifact (roadmap Sec 0.1 / M1 visual check): wireframe of
    the wake sheet + wing wall near the tip, plus a full-sheet overview."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    def tri_edges(tris):
        e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
        return np.unique(np.sort(e, axis=1), axis=0)

    wall = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)

    fig = plt.figure(figsize=(15, 6))
    views = [
        ("full wake sheet + wing", None),
        ("tip zoom (wake-tip closure)",
         (x_te(B_SEMI) - 0.35, x_te(B_SEMI) + 0.55,
          B_SEMI - 0.45, B_SEMI + 0.15)),
    ]
    for i, (title, win) in enumerate(views, start=1):
        ax = fig.add_subplot(1, 2, i, projection="3d")
        for tris, color, lw in ((wall, "tab:red", 0.3),
                                (wake, "tab:green", 0.25)):
            e = tri_edges(tris)
            if win is not None:
                x0, x1, z0, z1 = win
                keep = ((mesh.nodes[e, 0].min(axis=1) > x0 - 0.1)
                        & (mesh.nodes[e, 0].max(axis=1) < x1 + 0.1)
                        & (mesh.nodes[e, 2].min(axis=1) > z0 - 0.1)
                        & (mesh.nodes[e, 2].max(axis=1) < z1 + 0.1))
                e = e[keep]
            segs = mesh.nodes[e]
            ax.add_collection3d(Line3DCollection(segs, colors=color,
                                                 linewidths=lw))
        if win is None:
            ax.set_xlim(-0.5, 3.5)
            ax.set_ylim(-2.0, 2.0)
            ax.set_zlim(-0.5, 3.5)
        else:
            x0, x1, z0, z1 = win
            ax.set_xlim(x0, x1)
            ax.set_ylim(-0.5 * (x1 - x0), 0.5 * (x1 - x0))
            ax.set_zlim(z0, z1)
        ax.view_init(elev=28, azim=-125)
        ax.set_title(f"onera_m6 {level}: {title}\n(wall=red, wake=green)")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _write_cutplane_png(path, mesh, level):
    """M1 visual check: tet-slice wireframe on z = const planes just
    inboard and just outboard of the tip -- verifies valid volume fill
    around the tip and that the wake sheet ends AT the tip."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    planes = [B_SEMI - 0.05, B_SEMI + 0.05]
    for ax, z0 in zip(axes, planes):
        segs = _slice_segments(mesh.nodes, np.asarray(mesh.elements), z0)
        ax.add_collection(LineCollection(segs, colors="0.4", linewidths=0.3))
        wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)
        wsegs = _slice_segments(mesh.nodes, None, z0, tris=wake)
        if len(wsegs):
            ax.add_collection(LineCollection(wsegs, colors="tab:green",
                                             linewidths=1.5))
        ax.set_xlim(x_te(B_SEMI) - 0.9, x_te(B_SEMI) + 1.1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_aspect("equal")
        side = "inboard" if z0 < B_SEMI else "outboard"
        ax.set_title(f"z = {z0:.3f} ({side} of tip); wake slice = green")
    fig.suptitle(f"onera_m6 {level}: cut planes near the tip")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _slice_segments(nodes, elements, z0, tris=None):
    """Intersection segments of tet faces (or of given triangles) with the
    plane z = z0, projected to (x, y). Returns (n, 2, 2) array."""
    if tris is None:
        el = elements[
            (nodes[elements, 2].min(axis=1) < z0)
            & (nodes[elements, 2].max(axis=1) > z0)
        ]
        faces = np.concatenate([el[:, [0, 1, 2]], el[:, [0, 1, 3]],
                                el[:, [0, 2, 3]], el[:, [1, 2, 3]]])
        faces = np.unique(np.sort(faces, axis=1), axis=0)
    else:
        faces = tris

    p = nodes[faces]                      # (n, 3, 3)
    d = p[:, :, 2] - z0                   # signed distance of each vertex
    cross = (d.min(axis=1) < 0) & (d.max(axis=1) > 0)
    p, d = p[cross], d[cross]

    segs = []
    edge_pairs = [(0, 1), (1, 2), (2, 0)]
    for pk, dk in zip(p, d):
        pts = []
        for a, b in edge_pairs:
            if (dk[a] < 0) != (dk[b] < 0):
                t = dk[a] / (dk[a] - dk[b])
                q = pk[a] + t * (pk[b] - pk[a])
                pts.append(q[:2])
        if len(pts) == 2:
            segs.append(pts)
    return np.asarray(segs) if segs else np.empty((0, 2, 2))


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
