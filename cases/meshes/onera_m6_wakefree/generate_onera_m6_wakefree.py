"""
M4 deliverable: wake-FREE ONERA M6 half-wing tet mesh family (roadmap
Track M M4; design_track_b.md section 5.7 dual-mesh rule -- the 3D half).

Identical geometry, sizing policy and refinement ladder to the M1 family
(cases/meshes/onera_m6) except `embed_wake=False`: the chord-plane sheet is
built ONLY as the source of the wake-corridor Distance size field and is
neither fragmented into the fluid BRep nor embedded, so the tet mesh does
not conform to it and no "wake" tag exists. This is Track B's deliverable
form in 3D -- generic cuts through generic elements, with the swept TE and
the wing tip that the 2.5D M3 family cannot exercise (design_track_b.md D9
TE-polyline ruled surface; the spanwise clip that keeps the level set from
cutting beyond the tip).

Because the corridor size field is the SAME field the M1 family uses
(same h_wake, same Distance thresholds, sheet at the same chord-plane
location), element counts land close to M1's at equal h_wall -- which is
what makes the B7 A/B against the P5/P8 baselines a controlled
comparison. A wide alpha-sweep corridor (the 2.5D M3 wedge) is deliberately
NOT used here: in 3D the wedge volume scales with span and the element cost
is prohibitive (recorded in the roadmap M4 row) -- 3D alpha sweeps re-aim
the level set within the near-nominal band the corridor already resolves.

Tags: wall, farfield, symmetry. NO wake tag.

No .msh committed (gitignored like the M1 family; coarse+medium regenerate
in ~30 s). The per-level stats CSVs and the inspection PNGs ARE committed:

  <level>_stats.csv        size + quality + corridor-sizing report
  <level>_tags.png         boundary tags in 3D + the level-set wake plane
                           drawn as a REFERENCE overlay (it is not a mesh
                           entity) -- the counterpart of the M1 family's
                           <level>_wake_tip.png, which shows a real meshed
                           sheet there
  <level>_cutplane.png     z = const slices just inboard/outboard of the
                           tip, same framing as the M1 family's file --
                           the wake line is dashed (no mesh edges follow it)
  <level>_wake_corridor.png  mid-span slice with the level-set CUT elements
                           highlighted, side by side with the same slice of
                           the wake-embedded M1 mesh when it exists locally
                           -- the direct "wake in the topology vs wake in
                           the level set" comparison

Usage:
    python generate_onera_m6_wakefree.py --level coarse
    python generate_onera_m6_wakefree.py --all          # includes fine
(default: coarse + medium)
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import write_mesh, mesh_stats
from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, onera_m6_wing_mesh, x_te

R_FAR = 15.0 * MAC
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}


def _level_params(h_wall: float) -> dict:
    """Same policy as the M1 family (one parameter per level)."""
    return dict(
        h_wall=h_wall,
        h_wake=3.0 * h_wall,
        h_edge=0.5 * h_wall,
        h_far=min(2.5, 120.0 * h_wall),
    )


LEVELS = {
    "coarse": _level_params(0.030),
    "medium": _level_params(0.015),
    "fine": _level_params(0.0075),
}


def _corridor_stats(mesh, h_wake: float) -> dict:
    """Edge lengths of elements whose centroid lies in the wake corridor
    (behind the TE, within the sheet's span, |y| < 0.15 m) -- the M4 gate's
    sizing evidence against the M1 family's sheet-adjacent h_wake."""
    cen = mesh.nodes[mesh.elements].mean(axis=1)
    z = cen[:, 2]
    in_span = (z > 0.0) & (z < B_SEMI)
    behind = cen[:, 0] > np.array([x_te(zi) for zi in z])
    near_plane = np.abs(cen[:, 1]) < 0.15
    near_wing = cen[:, 0] < x_te(B_SEMI) + 2.0
    sel = in_span & behind & near_plane & near_wing
    el = mesh.elements[sel].astype(np.int64)
    if len(el) == 0:
        return {"corridor_n_elements": 0}
    p = mesh.nodes
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    e = np.concatenate([
        np.linalg.norm(p[el[:, i]] - p[el[:, j]], axis=1) for i, j in pairs
    ])
    return {
        "corridor_n_elements": int(len(el)),
        "corridor_median_edge": float(np.median(e)),
        "corridor_p90_edge": float(np.percentile(e, 90)),
        "h_wake_target": h_wake,
    }


def generate_level(out_dir: Path, level: str, inspect: bool = True) -> Path:
    p = LEVELS[level]
    t0 = time.perf_counter()
    mesh = onera_m6_wing_mesh(
        h_wall=p["h_wall"], h_far=p["h_far"], h_wake=p["h_wake"],
        h_edge=p["h_edge"], r_far=R_FAR,
        name=f"onera_m6_wakefree_{level}", embed_wake=False,
    )
    gen_seconds = time.perf_counter() - t0
    assert "wake" not in mesh.boundary_faces, "M4 mesh must carry no wake tag"

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
    stats.update(_corridor_stats(mesh, p["h_wake"]))

    ok = (stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"]
          and stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"])
    stats["quality_within_bounds"] = ok

    with open(out_dir / f"{level}_stats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in stats.items():
            writer.writerow([k, v])

    if inspect:
        _write_tags_png(out_dir / f"{level}_tags.png", mesh, level)
        _write_cutplane_png(out_dir / f"{level}_cutplane.png", mesh, level)
        _write_corridor_png(out_dir / f"{level}_wake_corridor.png",
                            mesh, level)

    if not ok:
        raise AssertionError(
            f"{level}: element quality outside bounds "
            f"(min dihedral {stats['min_dihedral_deg']:.2f} deg, "
            f"max aspect {stats['max_aspect_ratio']:.1f}; "
            f"bounds {QUALITY_BOUNDS})"
        )
    return out_path


# ---------------------------------------------------------------------------
# Headless inspection artifacts (roadmap Sec 0.1). The framing deliberately
# mirrors the M1 family's PNGs so the two can be put side by side.
# ---------------------------------------------------------------------------

WAKE_X_END = 6.0          # how far downstream to draw the wake reference


def _wake_levelset():
    """The chord-plane ruled wake the solver's Track B path uses on this
    geometry: swept TE line, freestream direction (alpha = 0 -> +x)."""
    from pyfp3d.wake import WakeLevelSet
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    return WakeLevelSet(te, direction=(1.0, 0.0, 0.0))


def _tri_edges(tris):
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    return np.unique(np.sort(e, axis=1), axis=0)


def _write_tags_png(path, mesh, level):
    """Boundary tags in 3D (wall / symmetry / far field) plus the level-set
    wake plane as a REFERENCE overlay. Compare with the M1 family's
    <level>_wake_tip.png, where the green sheet is an actual meshed surface
    with its own triangles and a conforming tip edge; here the green quad is
    drawn by matplotlib and exists nowhere in the mesh."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

    wall = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    sym = np.asarray(mesh.boundary_faces["symmetry"], dtype=np.int64)
    far = np.asarray(mesh.boundary_faces["farfield"], dtype=np.int64)

    # reference wake quad: TE line swept downstream in the chord plane
    quad = np.array([
        [x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI],
        [WAKE_X_END, 0.0, B_SEMI], [WAKE_X_END, 0.0, 0.0],
    ])

    fig = plt.figure(figsize=(15, 6))
    views = [
        ("boundary tags + the level-set wake plane", None),
        ("tip zoom: no sheet, no tip edge in the mesh",
         (x_te(B_SEMI) - 0.35, x_te(B_SEMI) + 0.55,
          B_SEMI - 0.45, B_SEMI + 0.15)),
    ]
    for i, (title, win) in enumerate(views, start=1):
        ax = fig.add_subplot(1, 2, i, projection="3d")
        layers = [(wall, "tab:red", 0.3)]
        if win is None:
            layers += [(sym, "tab:blue", 0.15), (far, "0.75", 0.12)]
        for tris, color, lw in layers:
            e = _tri_edges(tris)
            if win is not None:
                x0, x1, z0, z1 = win
                keep = ((mesh.nodes[e, 0].min(axis=1) > x0 - 0.1)
                        & (mesh.nodes[e, 0].max(axis=1) < x1 + 0.1)
                        & (mesh.nodes[e, 2].min(axis=1) > z0 - 0.1)
                        & (mesh.nodes[e, 2].max(axis=1) < z1 + 0.1))
                e = e[keep]
            ax.add_collection3d(Line3DCollection(mesh.nodes[e], colors=color,
                                                 linewidths=lw))
        if win is None:
            ax.add_collection3d(Poly3DCollection(
                [quad], facecolors="tab:green", alpha=0.18,
                edgecolors="tab:green", linewidths=1.0, linestyles="--"))
            ax.set_xlim(-0.5, 3.5)
            ax.set_ylim(-2.0, 2.0)
            ax.set_zlim(-0.5, 3.5)
        else:
            # In the zoom draw the sheet's OUTLINE only (a filled plane
            # would hide the mesh): the TE line and, crucially, the sheet's
            # tip edge -- both exist in the level set alone. On the M1 mesh
            # the same two lines are chains of real mesh edges.
            x0, x1, z0, z1 = win
            ax.plot([x_te(B_SEMI), x1], [0, 0], [B_SEMI, B_SEMI],
                    "--", color="tab:green", linewidth=1.6,
                    label="sheet tip edge (level set only)")
            ax.plot([x_te(B_SEMI - 0.4), x_te(B_SEMI)], [0, 0],
                    [B_SEMI - 0.4, B_SEMI], "-", color="tab:green",
                    linewidth=1.2, label="TE line (wake origin)")
            ax.legend(loc="upper left", fontsize=8)
            ax.set_xlim(x0, x1)
            ax.set_ylim(-0.5 * (x1 - x0), 0.5 * (x1 - x0))
            ax.set_zlim(z0, z1)
        ax.view_init(elev=28, azim=-125)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
    fig.suptitle(
        f"onera_m6_wakefree {level}: wall=red, symmetry=blue, farfield=grey; "
        "green = level-set wake, NOT a mesh entity\n"
        "(compare cases/meshes/onera_m6/"
        f"{level}_wake_tip.png, where the sheet is a real meshed surface "
        "with a conforming tip edge)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _write_cutplane_png(path, mesh, level):
    """Same framing as the M1 family's <level>_cutplane.png: z = const cut
    planes just inboard and just outboard of the tip. In the M1 figure the
    wake slice is a green POLYLINE of conforming mesh edges inboard of the
    tip and disappears outboard; here the wake is drawn dashed because no
    mesh edge follows it -- and outboard of the tip the level set does not
    cut at all (the spanwise clip, Gamma(tip) = 0)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, z0 in zip(axes, [B_SEMI - 0.05, B_SEMI + 0.05]):
        segs = _slice_segments(mesh.nodes, np.asarray(mesh.elements), z0)
        ax.add_collection(LineCollection(segs, colors="0.4", linewidths=0.3))
        inboard = z0 < B_SEMI
        if inboard:
            ax.plot([x_te(z0), x_te(B_SEMI) + 1.1], [0.0, 0.0],
                    "--", color="tab:green", linewidth=1.6,
                    label="level-set wake (no mesh edges on it)")
        else:
            ax.plot([], [], " ", label="outboard of tip: no wake, no cut")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_xlim(x_te(B_SEMI) - 0.9, x_te(B_SEMI) + 1.1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_aspect("equal")
        side = "inboard" if inboard else "outboard"
        ax.set_title(f"z = {z0:.3f} ({side} of tip)")
    fig.suptitle(f"onera_m6_wakefree {level}: cut planes near the tip "
                 f"(compare cases/meshes/onera_m6/{level}_cutplane.png, "
                 "where the wake slice is a chain of real mesh edges)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _write_corridor_png(path, mesh, level):
    """The direct comparison: a mid-span z = const slice through the wake
    corridor. Left = this wake-free mesh, with the elements the level set
    CUTS highlighted (they straddle the wake plane -- the jump is carried by
    extra DOFs, not by split nodes). Right = the same slice of the
    wake-embedded M1 mesh if it exists locally, where mesh edges lie ON the
    wake plane and the sheet is a real surface in the topology."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.wake import CutElementMap

    z0 = 0.5 * B_SEMI
    wall_nodes = np.unique(mesh.boundary_faces["wall"])
    cm = CutElementMap(mesh.nodes, mesh.elements, _wake_levelset(),
                       wall_nodes=wall_nodes)
    cut = np.asarray(mesh.elements, dtype=np.int64)[cm.cut_elems]

    m1_path = path.parent.parent / "onera_m6" / f"{level}.msh"
    m1 = read_mesh(m1_path) if m1_path.exists() else None

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharex=True, sharey=True)

    ax = axes[0]
    segs = _slice_segments(mesh.nodes, np.asarray(mesh.elements), z0)
    ax.add_collection(LineCollection(segs, colors="0.75", linewidths=0.3))
    csegs = _slice_segments(mesh.nodes, cut, z0)
    ax.add_collection(LineCollection(csegs, colors="tab:orange", linewidths=0.6,
                                     label="level-set cut elements"))
    ax.plot([x_te(z0), WAKE_X_END], [0.0, 0.0], "--", color="tab:green",
            linewidth=1.6, label="level-set wake (implicit)")
    ax.set_title(f"M4 wake-free: wake lives in the level set\n"
                 f"{len(cm.cut_elems)} cut elements, "
                 f"{cm.n_ext_dofs} extra DOFs (no node splitting)")
    ax.legend(loc="upper right", fontsize=8)

    ax = axes[1]
    if m1 is not None:
        segs = _slice_segments(m1.nodes, np.asarray(m1.elements), z0)
        ax.add_collection(LineCollection(segs, colors="0.75", linewidths=0.3))
        wsegs = _slice_segments(m1.nodes, None, z0,
                                tris=np.asarray(m1.boundary_faces["wake"],
                                                dtype=np.int64))
        ax.add_collection(LineCollection(wsegs, colors="tab:green",
                                         linewidths=1.8,
                                         label="conforming wake sheet (mesh faces)"))
        ax.set_title("M1 wake-embedded: wake is a surface in the topology\n"
                     "(mesh edges lie ON the sheet; nodes get duplicated)")
        ax.legend(loc="upper right", fontsize=8)
    else:
        ax.text(0.5, 0.5, "onera_m6/%s.msh not generated locally\n"
                          "(regenerate for the side-by-side)" % level,
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("M1 wake-embedded (mesh absent)")

    for ax in axes:
        ax.set_xlim(x_te(z0) - 0.5, x_te(z0) + 2.2)
        ax.set_ylim(-0.5, 0.5)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.suptitle(f"onera_m6 {level}: wake corridor at z = {z0:.3f} "
                 "(mid-span) -- same geometry, same sizing, different topology")
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
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    levels = sorted(LEVELS) if args.all else (args.level or ["coarse", "medium"])

    out_dir = Path(__file__).parent
    for level in levels:
        out_path = generate_level(out_dir, level)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
