"""
M5 deliverable: ONERA M6 half-wing tet mesh family with a ROUNDED TIP CAP
(docs/roadmap.md Track M gate M5; consumed by solver phase P13/G13.3).

WHY THIS FAMILY EXISTS. The M1 family closes the wing with a FLAT tip cap --
a documented deliberate simplification, standard for FP validation meshes.
P13/G13.3 measured what it costs: the flat cap meets the upper and lower
surfaces at a sharp convex edge, and in potential flow that edge is a
singularity. On the self-similar M1 ladder its box peak Mach DIVERGES under
uniform refinement (exponent p = +0.321) while the wake free edge (+0.045,
fixed by G13.2's spanwise loading taper) and the wing interior (-0.014) stay
bounded. A diverging singularity keeps the lift sequence out of the asymptotic
range, which is why no three-point Richardson has ever been earned on ONERA M6.

The cap is not an element-order problem -- isoparametric (P11) elements cannot
regularize a genuinely sharp geometric EDGE. The geometry is what is wrong, so
the fix is here: `wing3d.py`'s `tip_cap="round"` closes the wing with the half
body of revolution swept by the tip section about its own chord line. See that
module for why the TE line, the wake sheet, the Kutta stations and B_SEMI are
all UNCHANGED by it -- this family differs from M1 in the tip WALL and nothing
else, which is what makes the A/B against M1 a controlled one.

  <level>.msh            solver-ready mesh (tags: wall, farfield, symmetry,
                         wake -- same as M1)
  <level>_stats.csv      size + quality report, plus the M5 gate metric:
                         the crease angle on the tip-section seam (below)
  <level>_tipcap.png     the cap itself: wall wireframe at the tip, M1's flat
                         cap beside this family's rounded one when the M1 mesh
                         exists locally
  <level>_cutplane.png   z = const cut planes inboard of / at / outboard of
                         z = B_SEMI -- the third one is empty for M1 and shows
                         the cap for M5

★ THE GATE METRIC: SEAM CREASE ANGLE. On the tip-section profile at z = B_SEMI
-- the locus that used to BE the sharp edge -- take the turning angle between
the outward normals of adjacent wall triangles (post/surface.py::
wall_crease_angles), away from the LE and TE (both are sharp BY DESIGN in either
family; the TE carries the Kutta condition). Measured:

    h_wall      flat (M1)      round (M5)
    0.030        91.9 deg       46.8 deg
    0.015        92.1 deg       25.0 deg

The flat cap's angle does not move: it is a real edge, and refinement resolves
it rather than removing it. The rounded cap's HALVES when h halves -- the O(h)
signature of facets approximating a smooth surface, i.e. no edge in the limit.
That is the discrete statement of what G13.3 needs.

REFINEMENT LADDER. Unlike the M1 family this one is self-similar from the
start: h_far is NOT clamped (the M1b defect -- `h_far = min(2.5, 120 h_wall)`
bit only at coarse and knocked it off the refinement ray), so
RICHARDSON_LADDER = (coarse, medium, fine) refines by exactly 2 in every length
scale, h_tip included. There is no `coarse_ss` here because there is nothing to
repair.

No .msh committed (gitignored, M1 policy). The stats CSVs and PNGs ARE the
committed evidence.

Usage:
    python generate_onera_m6_roundtip.py --level coarse
    python generate_onera_m6_roundtip.py --all          # includes fine
(default: coarse + medium)
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.mesh.reader import mesh_stats, write_mesh
from pyfp3d.meshgen.wing3d import (
    B_SEMI, C_TIP, MAC, TIP_CAP_RADIUS, onera_m6_wing_mesh, x_le, x_te,
)
from pyfp3d.post.surface import wall_crease_angles

R_FAR = 15.0 * MAC

#: Same bounds as M1 (cases/meshes/onera_m6/generate_onera_m6.py). The rounded
#: cap is a small-radius feature and does cost some quality -- these bounds are
#: what says "some" is still inside what the solver has always accepted.
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}

#: The seam is only measured away from the LE and the TE: the section has zero
#: thickness at both, so the cap radius vanishes there and the geometry is
#: sharp in BOTH families by design (that is the airfoil, not the cap).
SEAM_XI = (0.05, 0.95)

#: M5 gate: the seam crease must be an O(h) artifact of faceting, not an edge.
#: Bounds are set well inside the measured flat/round separation (92 deg vs
#: 47/25/12), so a regression that reintroduced an edge would fire immediately.
SEAM_MAX_DEG = {"coarse": 60.0, "medium": 35.0, "fine": 20.0}


def _level_params(h_wall: float) -> dict:
    """Everything from the one wall-size parameter. No h_far clamp: see the
    module docstring and the M1b ledger row."""
    return dict(
        h_wall=h_wall,
        h_wake=3.0 * h_wall,
        h_edge=0.5 * h_wall,
        h_far=120.0 * h_wall,
        h_tip=0.25 * h_wall,
    )


LEVELS = {
    "coarse": _level_params(0.030),
    "medium": _level_params(0.015),
    "fine": _level_params(0.0075),
}

#: Self-similar by construction -- every length scale halves per level.
RICHARDSON_LADDER = ("coarse", "medium", "fine")


def seam_crease_angles(mesh) -> np.ndarray:
    """Turning angle of the wall across the tip-section seam at z = B_SEMI --
    the locus that IS the sharp edge in the flat-cap family."""
    ang, mid = wall_crease_angles(mesh.nodes, mesh.elements,
                                  mesh.boundary_faces["wall"])
    xi = (mid[:, 0] - x_le(B_SEMI)) / C_TIP
    on_seam = ((np.abs(mid[:, 2] - B_SEMI) < 1e-9)
               & (xi > SEAM_XI[0]) & (xi < SEAM_XI[1]))
    return ang[on_seam]


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    t0 = time.perf_counter()
    mesh = onera_m6_wing_mesh(
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        h_edge=p["h_edge"], h_tip=p["h_tip"], r_far=R_FAR,
        tip_cap="round", name=f"onera_m6_roundtip_{level}",
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
    for k, v in p.items():
        stats[k] = v
    stats["gen_seconds"] = round(gen_seconds, 1)

    seam = seam_crease_angles(mesh)
    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    stats["n_seam_edges"] = len(seam)
    stats["seam_crease_max_deg"] = float(seam.max())
    stats["seam_crease_p99_deg"] = float(np.percentile(seam, 99))
    stats["wall_z_max"] = float(mesh.nodes[wall_nodes, 2].max())
    stats["tip_cap_apex_z"] = B_SEMI + TIP_CAP_RADIUS

    quality_ok = (stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"]
                  and stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"])
    seam_ok = stats["seam_crease_max_deg"] <= SEAM_MAX_DEG[level]
    stats["quality_within_bounds"] = quality_ok
    stats["seam_within_bounds"] = seam_ok

    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    _check_wake_tip_closure(mesh)

    if inspect:
        _write_tipcap_png(out_dir / f"{level}_tipcap.png", mesh, level)
        _write_cutplane_png(out_dir / f"{level}_cutplane.png", mesh, level)

    if not quality_ok:
        raise AssertionError(
            f"{level}: element quality outside bounds "
            f"(min dihedral {stats['min_dihedral_deg']:.2f} deg, "
            f"max aspect {stats['max_aspect_ratio']:.1f}; "
            f"bounds {QUALITY_BOUNDS})"
        )
    if not seam_ok:
        raise AssertionError(
            f"{level}: tip-section seam crease {stats['seam_crease_max_deg']:.1f} "
            f"deg exceeds {SEAM_MAX_DEG[level]:.1f} -- the cap is not rounded, "
            "or h_tip no longer resolves it (the flat cap measures ~92 deg)"
        )
    return out_path


def _check_wake_tip_closure(mesh) -> None:
    """The M1 check, unchanged and deliberately re-run here: the rounded cap
    must not have disturbed the wake attachment. The sheet's tip edge must be
    one open chain starting at the wing tip TE corner -- which is exactly where
    the cap radius vanishes, so it should still be the same point it always was.
    """
    wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)
    wall_nodes = np.unique(mesh.boundary_faces["wall"])

    edges = np.sort(np.concatenate([wake[:, [0, 1]], wake[:, [1, 2]],
                                    wake[:, [2, 0]]]), axis=1)
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = uniq[counts == 1]

    at_tip = np.abs(mesh.nodes[:, 2] - B_SEMI) < 1e-6
    tip_edges = boundary_edges[at_tip[boundary_edges].all(axis=1)]
    assert len(tip_edges) > 0, "wake sheet has no tip edge at z = B_SEMI"

    nodes, deg = np.unique(tip_edges, return_counts=True)
    assert deg.max() <= 2, "wake tip edge self-intersects (node degree > 2)"
    endpoints = nodes[deg == 1]
    assert len(endpoints) == 2, (
        f"wake tip edge is not one open chain ({len(endpoints)} endpoints)"
    )

    on_wall = np.isin(endpoints, wall_nodes)
    assert on_wall.sum() == 1, (
        "wake tip edge does not attach to the wall at exactly one end"
    )
    xyz = mesh.nodes[endpoints[on_wall][0]]
    assert abs(xyz[0] - x_te(B_SEMI)) < 1e-9 and abs(xyz[1]) < 1e-9, (
        f"tip TE corner at {xyz}, expected ({x_te(B_SEMI):.6f}, 0, {B_SEMI})"
    )


def _write_tipcap_png(path, mesh, level):
    """The M5 visual check: the cap. Wall wireframe at the tip, coloured by the
    crease angle across each edge -- the flat cap draws its own sharp edge as a
    bright ~90 deg line, the rounded one has nothing to draw."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    from pyfp3d.mesh.reader import read_mesh

    ang, mid = wall_crease_angles(mesh.nodes, mesh.elements,
                                  mesh.boundary_faces["wall"])
    panels = [("M5 rounded cap (this family)", mesh, ang, mid)]

    m1 = path.parent.parent / "onera_m6" / f"{level}.msh"
    if m1.exists():
        flat = read_mesh(m1)
        fang, fmid = wall_crease_angles(flat.nodes, flat.elements,
                                        flat.boundary_faces["wall"])
        panels.insert(0, ("M1 flat cap (sharp edge)", flat, fang, fmid))

    norm = Normalize(0, 90)
    fig = plt.figure(figsize=(7.2 * len(panels), 6.0))
    for i, (title, m, a, c) in enumerate(panels, start=1):
        ax = fig.add_subplot(1, len(panels), i, projection="3d")
        wall = np.asarray(m.boundary_faces["wall"], dtype=np.int64)
        e = np.unique(np.sort(np.concatenate(
            [wall[:, [0, 1]], wall[:, [1, 2]], wall[:, [2, 0]]]), axis=1), axis=0)
        near = ((m.nodes[e, 2].min(axis=1) > B_SEMI - 0.16)
                & (m.nodes[e, 0].max(axis=1) < x_te(B_SEMI) + 0.05))
        ax.add_collection3d(Line3DCollection(m.nodes[e[near]], colors="0.72",
                                             linewidths=0.35))
        # the crease edges, coloured -- only the ones worth looking at
        hot = a > 8.0
        if hot.any():
            ax.scatter(c[hot, 0], c[hot, 1], c[hot, 2], c=a[hot], cmap="inferno",
                       norm=norm, s=5, depthshade=False)
        ax.set_xlim(x_le(B_SEMI) - 0.05, x_te(B_SEMI) + 0.05)
        ax.set_ylim(-0.25, 0.25)
        ax.set_zlim(B_SEMI - 0.16, B_SEMI + 0.34)
        ax.view_init(elev=22, azim=-132)
        ax.set_title(f"{title}\n{level}: wall wireframe + crease angle > 8 deg")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
    fig.colorbar(ScalarMappable(norm=norm, cmap="inferno"),
                 ax=fig.axes, shrink=0.6, label="wall crease angle [deg]")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _write_cutplane_png(path, mesh, level):
    """z = const tet slices inboard of / at / outboard of the tip. The third
    plane is the one that matters: it is EMPTY for the flat cap and shows the
    rounded cap's section for this family."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    planes = [B_SEMI - 0.05, B_SEMI - 0.005, B_SEMI + 0.5 * TIP_CAP_RADIUS]
    labels = ["inboard of the tip", "just inboard of z = B_SEMI",
              "OUTBOARD of z = B_SEMI (the rounded cap)"]
    for ax, z0, lab in zip(axes, planes, labels):
        segs = _slice_segments(mesh.nodes, np.asarray(mesh.elements), z0)
        if len(segs):
            ax.add_collection(LineCollection(segs, colors="0.4", linewidths=0.3))
        wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)
        wsegs = _slice_segments(mesh.nodes, None, z0, tris=wake)
        if len(wsegs):
            ax.add_collection(LineCollection(wsegs, colors="tab:green",
                                             linewidths=1.5))
        ax.set_xlim(x_te(B_SEMI) - 0.9, x_te(B_SEMI) + 0.6)
        ax.set_ylim(-0.55, 0.55)
        ax.set_aspect("equal")
        ax.set_title(f"z = {z0:.4f}\n{lab}")
    fig.suptitle(f"onera_m6_roundtip {level}: cut planes through the tip "
                 f"(wake slice = green)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _slice_segments(nodes, elements, z0, tris=None):
    """Intersection segments of tet faces (or of given triangles) with the
    plane z = z0, projected to (x, y). Returns (n, 2, 2)."""
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

    p = nodes[faces]
    d = p[:, :, 2] - z0
    cross = (d.min(axis=1) < 0) & (d.max(axis=1) > 0)
    p, d = p[cross], d[cross]

    segs = []
    for pk, dk in zip(p, d):
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if (dk[a] < 0) != (dk[b] < 0):
                t = dk[a] / (dk[a] - dk[b])
                pts.append((pk[a] + t * (pk[b] - pk[a]))[:2])
        if len(pts) == 2:
            segs.append(pts)
    return np.asarray(segs) if segs else np.empty((0, 2, 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", action="append", choices=sorted(LEVELS))
    parser.add_argument("--all", action="store_true", help="generate all levels")
    args = parser.parse_args()
    levels = (["coarse", "medium", "fine"] if args.all
              else (args.level or ["coarse", "medium"]))

    out_dir = Path(__file__).parent
    for level in levels:
        out_path = generate_level(out_dir, level)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
