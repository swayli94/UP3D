"""
Track B / B9 (re-spec 2026-07-17): the CONFORMING wing-body mesh family --
ONERA M6 wing + fuselage with an EMBEDDED wake sheet, for the conforming
solver path (`cut_wake` + coupled Newton + P14 pressure Kutta).

Run:  python cases/meshes/onera_m6_wingbody_conforming/generate_onera_m6_wingbody_conforming.py
      [--levels coarse medium]

Artifacts per level:
  <level>.msh                 volume mesh (GITIGNORED -- M1/M5 policy)
  <level>_stats.csv           counts, quality, junction/fuselage/wake metrics
  <level>_wingbody.png        overview + junction zoom + wake topology panel

Same geometry, same h_wall ladder and same size fields as the wake-free M2
family (`../onera_m6_wingbody/`) -- the ONLY difference is
`embed_wake=True`: the sheet spans z_lo..tip, pierces the far-field sphere,
and is fragment+embed'ed. The fragment trims its inboard boundary to the
fuselage WATERLINE (the y=0 top meridian) alongside the body and to the
symmetry plane aft of it, which is what keeps the cut simply connected
around body+wake; downstream it reaches the sphere so the branch-cut
far-field Dirichlet data has duplicated wake nodes to land on.

GENERATION-TIME cut_wake INGEST GATE (the crack detector): a sheet-body
stitch failure turns waterline edges into interior FREE edges, so we run
cut_wake right here and assert every free node sits at the tip
(z ~ B_SEMI). That converts the riskiest failure mode of this family --
an OCC fragment that silently fails to stitch -- into a generator error
with the geometry still in scope, instead of a wrong solve much later.
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.mesh.reader import write_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.fuselage import FuselageParams, radius_at
from pyfp3d.meshgen.wing3d import B_SEMI, C_ROOT, x_te
from pyfp3d.meshgen.wingbody import (
    H_FAR_IN_H_WALL,
    R_FAR,
    junction_z,
    onera_m6_wingbody_mesh,
)
from pyfp3d.post.surface import wall_crease_angles

OUT_DIR = Path(__file__).resolve().parent

# Same wall sizes as the M1/M5/M2 ladders (controlled A/B at equal h_wall).
# coarse + medium ONLY -- the B9 scope (no fine; the conforming fine is the
# recorded splu trap anyway).
LEVELS = {"coarse": 0.030, "medium": 0.015}

QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}
FUSELAGE_CREASE_MAX_DEG = {"coarse": 25.0, "medium": 15.0}

FUSELAGE = FuselageParams()


def _level_params(h_wall: float) -> dict:
    """Identical to the wake-free family's _level_params (self-similar)."""
    return {
        "h_wall": h_wall,
        "h_far": H_FAR_IN_H_WALL * h_wall,
        "h_wake": 3.0 * h_wall,
        "h_edge": 0.5 * h_wall,
        "h_tip": 0.25 * h_wall,
        "h_junction": 0.5 * h_wall,
        "h_body": 2.0 * h_wall,
        "h_body_tip": 0.25 * h_wall,
    }


def _fuselage_crease(mesh) -> tuple:
    ang, _ = wall_crease_angles(mesh.nodes, mesh.elements,
                                mesh.boundary_faces["fuselage"])
    return float(np.percentile(ang, 99)), float(ang.max()), float(np.median(ang))


def _bad_tet_localization(mesh, dihedral, deg: float = 5.0) -> dict:
    """Where the worst tets live (x-centroid stats), so a dihedral-gate
    failure is diagnosable from the CSV alone (downstream pole vs junction
    vs waterline)."""
    bad = dihedral < deg
    if not bad.any():
        return {"n_tets_below_5deg": 0, "bad_tet_x_min": float("nan"),
                "bad_tet_x_med": float("nan"), "bad_tet_x_max": float("nan")}
    xc = mesh.nodes[mesh.elements[bad]].mean(axis=1)[:, 0]
    return {
        "n_tets_below_5deg": int(bad.sum()),
        "bad_tet_x_min": float(xc.min()),
        "bad_tet_x_med": float(np.median(xc)),
        "bad_tet_x_max": float(xc.max()),
    }


def _ingest_gate(mesh, p: FuselageParams) -> dict:
    """Run cut_wake NOW and assert the B9 wake topology (the crack detector).

    - free nodes only at the tip edge: a stitch failure anywhere (sheet LE
      vs wing TE, waterline vs body, aft edge vs symmetry) turns those edges
      into interior free edges at z < B_SEMI, and this assert names it;
    - the innermost Kutta station is the junction (z = r_f) and stations are
      one per TE node (3D swept TE);
    - every waterline (wake∩fuselage) node is duplicated (carries the jump);
    - TE nodes stay WING-only: `wall_tag="wall"` -- widening it to include
      the fuselage would mint spurious Kutta stations along the waterline.
    """
    t0 = time.perf_counter()
    mesh_cut, wc = cut_wake(mesh)          # internal assert_wake_topology
    dt = time.perf_counter() - t0

    z_junc = junction_z(p)
    if len(wc.free_nodes):
        zf = mesh.nodes[wc.free_nodes, 2]
        assert zf.min() > B_SEMI - 1e-6, (
            f"free (single-valued) sheet nodes at z_min={zf.min():.4f} < "
            f"B_SEMI -- the sheet is NOT stitched (crack at the TE, the "
            "waterline, or the aft symmetry edge)"
        )
    assert len(wc.te_nodes) > 0, "no TE (wall∩wake) nodes"
    st_z = np.sort(wc.station_z)
    assert abs(st_z[0] - z_junc) < 0.05, (
        f"innermost Kutta station z={st_z[0]:.4f}, expected the junction "
        f"{z_junc:.4f}"
    )
    assert wc.n_stations == len(wc.te_nodes), (
        f"{wc.n_stations} stations != {len(wc.te_nodes)} TE nodes -- the 3D "
        "swept TE should give one station per TE node"
    )
    # Waterline nodes (wake∩fuselage) must all be duplicated -- they carry
    # the potential jump onto the body (same rule as the symmetry root edge
    # on the wing-alone family).
    wake_nodes = np.unique(mesh.boundary_faces["wake"])
    fus_nodes = set(np.unique(mesh.boundary_faces["fuselage"]).tolist())
    waterline = np.array([n for n in wake_nodes.tolist() if n in fus_nodes],
                         dtype=np.int64)
    master_set = set(wc.master_nodes.tolist())
    n_wl_dup = sum(1 for n in waterline.tolist() if n in master_set)
    assert n_wl_dup == len(waterline), (
        f"only {n_wl_dup}/{len(waterline)} waterline nodes duplicated -- "
        "the jump does not reach the fuselage"
    )
    # TE stations are wing-only by construction; the innermost TE node is
    # the junction vertex itself.
    te_z = mesh.nodes[wc.te_nodes, 2]
    assert te_z.min() > z_junc - 1e-6, "TE station inboard of the junction"

    return {
        "cut_n_stations": int(wc.n_stations),
        "cut_n_te_nodes": int(len(wc.te_nodes)),
        "cut_n_free": int(len(wc.free_nodes)),
        "cut_n_dup": int(len(wc.master_nodes)),
        "n_waterline_nodes": int(len(waterline)),
        "cut_ingest_s": dt,
    }


def generate_level(level: str, render: bool = True) -> dict:
    prm = _level_params(LEVELS[level])
    print(f"[{level}] h_wall={prm['h_wall']:.4f} (embed_wake=True) ...",
          flush=True)
    mesh = onera_m6_wingbody_mesh(fuselage=FUSELAGE, tip_cap="round",
                                  embed_wake=True, **prm)

    dihedral = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
    aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
    cre_p99, cre_max, cre_med = _fuselage_crease(mesh)

    wake_nodes = np.unique(mesh.boundary_faces["wake"])
    stats = {
        "level": level,
        **prm,
        "n_nodes": len(mesh.nodes),
        "n_tets": len(mesh.elements),
        "min_dihedral_deg": float(dihedral.min()),
        "max_aspect_ratio": float(aspect.max()),
        "fuselage_crease_p99_deg": cre_p99,
        "fuselage_crease_max_deg": cre_max,
        "fuselage_crease_median_deg": cre_med,
        **{k: len(v) for k, v in
           (("n_tris_" + g, f) for g, f in mesh.boundary_faces.items())},
        "n_wake_nodes": int(len(wake_nodes)),
        "wake_z_min": float(mesh.nodes[wake_nodes, 2].min()),
        "wake_z_max": float(mesh.nodes[wake_nodes, 2].max()),
        "wake_x_max": float(mesh.nodes[wake_nodes, 0].max()),
        **_bad_tet_localization(mesh, dihedral),
        **_ingest_gate(mesh, FUSELAGE),
    }

    # --- gates -----------------------------------------------------------
    assert "wake" in mesh.boundary_faces, "conforming family must embed"
    assert stats["min_dihedral_deg"] >= QUALITY_BOUNDS["min_dihedral_deg"], (
        f"min dihedral {stats['min_dihedral_deg']:.2f} deg below bound "
        f"{QUALITY_BOUNDS['min_dihedral_deg']} (bad tets at x in "
        f"[{stats['bad_tet_x_min']:.2f}, {stats['bad_tet_x_max']:.2f}])"
    )
    assert stats["max_aspect_ratio"] <= QUALITY_BOUNDS["max_aspect_ratio"], (
        f"max aspect {stats['max_aspect_ratio']:.1f} above bound "
        f"{QUALITY_BOUNDS['max_aspect_ratio']}"
    )
    assert cre_p99 <= FUSELAGE_CREASE_MAX_DEG[level], (
        f"fuselage crease p99 {cre_p99:.1f} deg above bound "
        f"{FUSELAGE_CREASE_MAX_DEG[level]}"
    )

    write_mesh(mesh, str(OUT_DIR / f"{level}.msh"))
    with open(OUT_DIR / f"{level}_stats.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "value"])
        for k, v in stats.items():
            w.writerow([k, v])
    if render:
        _write_png(OUT_DIR / f"{level}_wingbody.png", mesh, level)

    print(f"[{level}] {stats['n_tets']} tets | dihedral "
          f"{stats['min_dihedral_deg']:.2f} deg | aspect "
          f"{stats['max_aspect_ratio']:.1f} | crease p99 {cre_p99:.1f} deg | "
          f"wake tris {stats['n_tris_wake']} | TE/stations "
          f"{stats['cut_n_te_nodes']}/{stats['cut_n_stations']} | "
          f"free {stats['cut_n_free']} | dup {stats['cut_n_dup']} | "
          f"waterline {stats['n_waterline_nodes']}", flush=True)
    return stats


def _write_png(path: Path, mesh, level: str) -> None:
    """Overview + junction zoom + wake topology panel (waterline + TE +
    tip edge + aft symmetry edge -- the objects the ingest gate asserts)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    p = FUSELAGE
    z_junc = junction_z(p)

    def edges(tris):
        t = np.asarray(tris, dtype=np.int64)
        return np.unique(np.sort(np.concatenate(
            [t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1), axis=0)

    def true_aspect(ax):
        ax.set_box_aspect([hi - lo for lo, hi in
                           (ax.get_xlim(), ax.get_ylim(), ax.get_zlim())])

    fig = plt.figure(figsize=(20, 6.6))

    ax = fig.add_subplot(1, 3, 1, projection="3d")
    for g, col, lw in (("fuselage", "0.62", 0.2), ("wall", "tab:blue", 0.2),
                       ("wake", "tab:red", 0.15)):
        e = edges(mesh.boundary_faces[g])
        near = mesh.nodes[e, 0].max(axis=1) < p.x_tail_tip + 1.5
        ax.add_collection3d(Line3DCollection(mesh.nodes[e[near]], colors=col,
                                             linewidths=lw))
    ax.set_xlim(p.x_nose_tip - 0.1, p.x_tail_tip + 1.5)
    ax.set_ylim(-0.6, 0.6)
    ax.set_zlim(0.0, 1.3)
    ax.view_init(elev=25, azim=-60)
    ax.set_title(f"{level}: wing (blue) + fuselage (grey) + EMBEDDED wake "
                 "(red)\nsheet rides the waterline, reaches the symmetry "
                 "plane aft of the body")
    true_aspect(ax)

    ax = fig.add_subplot(1, 3, 2, projection="3d")
    for g, col, lw in (("fuselage", "0.62", 0.35), ("wall", "tab:blue", 0.35),
                       ("wake", "tab:red", 0.3)):
        e = edges(mesh.boundary_faces[g])
        near = ((mesh.nodes[e, 2].max(axis=1) < 3 * z_junc)
                & (mesh.nodes[e, 0].max(axis=1) < x_te(z_junc) + 0.6)
                & (mesh.nodes[e, 0].min(axis=1) > 0.3))
        ax.add_collection3d(Line3DCollection(mesh.nodes[e[near]], colors=col,
                                             linewidths=lw))
    ax.set_xlim(0.3, x_te(z_junc) + 0.6)
    ax.set_ylim(-0.2, 0.2)
    ax.set_zlim(0.0, 3 * z_junc)
    ax.view_init(elev=28, azim=-70)
    ax.set_title("JUNCTION zoom: the sheet must leave the wing TE and ride\n"
                 "the fuselage waterline with no crack")
    true_aspect(ax)

    # Wake sheet in plan view (x-z), with the boundary classes marked.
    ax = fig.add_subplot(1, 3, 3)
    wake_nodes = np.unique(mesh.boundary_faces["wake"])
    fus_nodes = set(np.unique(mesh.boundary_faces["fuselage"]).tolist())
    wall_nodes = set(np.unique(mesh.boundary_faces["wall"]).tolist())
    sym_nodes = set(np.unique(mesh.boundary_faces["symmetry"]).tolist())
    wn = mesh.nodes[wake_nodes]
    ax.plot(wn[:, 0], wn[:, 2], ".", ms=1.0, color="0.8", label="wake nodes")
    for sel_set, col, lab in ((wall_nodes, "tab:blue", "TE (wall∩wake)"),
                              (fus_nodes, "tab:red", "waterline (fus∩wake)"),
                              (sym_nodes, "tab:green", "symmetry edge")):
        sel = np.array([n in sel_set for n in wake_nodes.tolist()])
        if sel.any():
            ax.plot(wn[sel, 0], wn[sel, 2], ".", ms=3.0, color=col, label=lab)
    xs = np.linspace(p.x_nose_tip, p.x_tail_tip, 200)
    ax.plot(xs, [radius_at(p, float(x)) for x in xs], "k-", lw=0.8,
            label="body waterline R(x)")
    ax.axhline(B_SEMI, color="0.5", ls=":", lw=0.8)
    ax.set_xlim(0.0, p.x_tail_tip + 2.0)
    ax.set_ylim(-0.02, B_SEMI + 0.1)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_title("wake sheet plan view: boundary classes\n"
                 "(zoomed to the body; sheet continues to the sphere)")

    fig.suptitle(f"onera_m6_wingbody_conforming {level} (embedded wake)")
    fig.subplots_adjust(wspace=0.3)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse"],
                    choices=list(LEVELS))
    ap.add_argument("--no-render", action="store_true")
    args = ap.parse_args()

    print(f"fuselage: r_f={FUSELAGE.r_f} length={FUSELAGE.length:.3f} "
          f"junction_z={junction_z(FUSELAGE):.3f}  r_far={R_FAR:.2f}")
    rows = [generate_level(lv, render=not args.no_render) for lv in args.levels]

    with open(OUT_DIR / "ladder_stats.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT_DIR}/ladder_stats.csv")


if __name__ == "__main__":
    main()
