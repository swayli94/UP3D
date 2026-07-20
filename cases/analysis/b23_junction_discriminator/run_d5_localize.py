"""D5 (post-hoc localization, NOT pre-registered): WHERE is the pocket?

Runs no new solves -- reloads the cached D1/D2 states and classifies the
highest-M^2 elements against the wake geometry:

  * |s|  signed distance to the wake sheet (WakeLevelSet.evaluate)
  * q    spanwise arclength along the TE polyline; q ~ 0 = the sheet's
         INBOARD FREE EDGE (the junction-TE end), q ~ span_length = tip edge
  * d    downstream coordinate from the TE
  * cut  whether the element is a cut (wake) element, a TE-fan element,
         or a plain element

Pre-registered branch-2 evidence (lift/wake-coupled): a pocket sitting at
s ~ 0, q ~ 0, d > 0 lives on the wake's inboard free edge AFT of the body --
NOT on the wing-fuselage wall crease (x in [0, 0.806]).

Run:  python cases/analysis/b23_junction_discriminator/run_d5_localize.py
Artifacts: results/d5_localize.csv, results/d5_pocket_map.png
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb_common import (ALPHA_REF, FUS, LS_MESH_DIR, OUT,
                       build_ls, load_mesh)
from pyfp3d.meshgen.fuselage import radius_at

TOPN = 200

#: (case_tag, level, mesh_path, cache_tag, alpha)
CASES = [
    ("d1_coarse_a3.06", "coarse", LS_MESH_DIR / "coarse.msh", "d1_coarse", ALPHA_REF),
    ("d1_medium_a3.06", "medium", LS_MESH_DIR / "medium.msh", "d1_medium", ALPHA_REF),
    ("d1_medium_a2.0", "medium", LS_MESH_DIR / "medium.msh", "d1_medium", 2.0),
    ("d1_medium_a1.0", "medium", LS_MESH_DIR / "medium.msh", "d1_medium", 1.0),
    ("d1_medium_a0.0", "medium", LS_MESH_DIR / "medium.msh", "d1_medium", 0.0),
]


def analyze(case, mesh):
    tag, level, mp, cache_tag, alpha = case
    cache = OUT / f"{cache_tag}_a{alpha:.2f}_m0.50.npz"
    if not cache.exists():
        print(f"[skip] {tag}: no cache {cache.name}", flush=True)
        return None, None
    d = np.load(cache, allow_pickle=True)
    phi_ext = d["phi_ext"]
    wls, cm, mvop = build_ls(mesh, alpha)
    m2 = np.asarray(mvop.element_mach2(phi_ext, 0.5, 1.4, 1.0))
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    top = np.argsort(m2)[::-1][:TOPN]
    tc = cents[top]
    s, dd, q = wls.evaluate(tc)

    cut_set = set(cm.cut_elems.tolist())
    fan_set = set(cm.te_lower_elems.tolist())
    is_cut = np.array([(int(e) in cut_set) for e in top])
    is_fan = np.array([(int(e) in fan_set) for e in top])

    rr = np.array([radius_at(FUS, float(x)) for x in tc[:, 0]])
    dist_fus = np.abs(np.hypot(tc[:, 1], tc[:, 2]) - rr)

    sup = m2[top] > 1.0
    n_sup = int(sup.sum())
    agg = dict(
        case=tag, alpha=alpha, n_tets=len(mesh.elements), topn=TOPN,
        n_sup_topn=n_sup,
        mmax=float(np.sqrt(m2[top[0]])),
        frac_cut_topn=float(is_cut.mean()),
        frac_fan_topn=float(is_fan.mean()),
        med_abs_s_topn=float(np.median(np.abs(s))),
        med_q_topn=float(np.median(q)),
        med_d_topn=float(np.median(dd)),
        med_dist_fus_topn=float(np.median(dist_fus)),
        med_x_topn=float(np.median(tc[:, 0])),
        med_z_topn=float(np.median(tc[:, 2])),
        frac_x_in_crease=float(((tc[:, 0] >= 0.0) & (tc[:, 0] <= 0.806)
                                & (np.abs(tc[:, 2]) < 0.30)).mean()),
    )
    if n_sup:
        agg.update(dict(
            med_abs_s_sup=float(np.median(np.abs(s[sup]))),
            med_q_sup=float(np.median(q[sup])),
            med_dist_fus_sup=float(np.median(dist_fus[sup])),
            med_x_sup=float(np.median(tc[sup, 0])),
            frac_cut_sup=float(is_cut[sup].mean()),
        ))
    else:
        agg.update(dict(med_abs_s_sup=np.nan, med_q_sup=np.nan,
                        med_dist_fus_sup=np.nan, med_x_sup=np.nan,
                        frac_cut_sup=np.nan))
    cloud = dict(x=tc[:, 0], z=tc[:, 2], m2=m2[top], s=s, q=q,
                 cut=is_cut, alpha=alpha, tag=tag)
    print(f"  [{tag}] Mmax={agg['mmax']:.2f} n_sup(top{TOPN})={n_sup} "
          f"cut={agg['frac_cut_topn']:.2f} med|s|={agg['med_abs_s_topn']:.4f} "
          f"med q={agg['med_q_topn']:.3f} med x={agg['med_x_topn']:.2f} "
          f"crease-frac={agg['frac_x_in_crease']:.2f}", flush=True)
    return agg, cloud


def main():
    meshes = {}
    rows, clouds = [], []
    for case in CASES:
        level = case[1]
        if level not in meshes:
            if not case[2].exists():
                print(f"[skip] {level}: mesh absent")
                continue
            meshes[level] = load_mesh(case[2])
        agg, cloud = analyze(case, meshes[level])
        if agg:
            rows.append(agg)
            clouds.append(cloud)

    if not rows:
        print("no cached solves found; run D1 first")
        sys.exit(0)

    keys = ["case", "alpha", "n_tets", "topn", "n_sup_topn", "mmax",
            "frac_cut_topn", "frac_fan_topn", "med_abs_s_topn", "med_q_topn",
            "med_d_topn", "med_dist_fus_topn", "med_x_topn", "med_z_topn",
            "frac_x_in_crease", "med_abs_s_sup", "med_q_sup",
            "med_dist_fus_sup", "med_x_sup", "frac_cut_sup"]
    with open(OUT / "d5_localize.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'d5_localize.csv'}")

    # pocket map: top-N M^2 elements in the (x, z) plane, per case
    n = len(clouds)
    fig, ax = plt.subplots(1, n, figsize=(5.2 * n, 4.8), squeeze=False)
    for j, cl in enumerate(clouds):
        a = ax[0][j]
        sc = a.scatter(cl["x"], cl["z"], c=np.log10(cl["m2"] + 1e-12),
                       cmap="inferno", s=14)
        a.scatter(cl["x"][cl["cut"]], cl["z"][cl["cut"]],
                  facecolors="none", edgecolors="tab:cyan", s=26, lw=0.7,
                  label="cut (wake) elem")
        a.axvspan(0.0, 0.806, color="tab:green", alpha=0.08)
        a.axhline(0.15, color="0.5", ls="--", lw=0.8)
        a.axhline(1.1963, color="0.5", ls=":", lw=0.8)
        a.set_xlim(-0.1, 3.2); a.set_ylim(-0.05, 1.35)
        a.set_xlabel("x"); a.set_ylabel("z")
        a.set_title(f"{cl['tag']}", fontsize=9)
        a.grid(alpha=0.3)
        if j == 0:
            a.legend(fontsize=7, loc="upper right")
    fig.suptitle("D5 pocket map: top-200 M^2 elements (green band = wall "
                 "crease chord; dashed = z_junc; dotted = tip)")
    fig.tight_layout()
    fig.savefig(OUT / "d5_pocket_map.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'d5_pocket_map.png'}")


if __name__ == "__main__":
    main()
