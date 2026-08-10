"""D2b -- region-split Mmax for the D2 junction-refinement ladder.

The D2 global-Mmax trend (4.71 / 14.66 / 23.73 for h_junction 0.015 /
0.0075 / 0.00375) is POLLUTED by the wingtip free-edge singularity: on
j100 and j025 the argmax element sits at the tip TE (z ~ 1.197, x ~ 1.135,
|r-R| ~ 1.05), not at the junction. The pre-registered D2 question -- "is
the junction pocket a function of h_junction, and which way does it move?"
-- must be answered on a JUNCTION-REGION metric, with the tip region
reported separately as the control.

Region split on element centroids (reuses the cached D2 solves; no re-solve):
  * junction corridor: z < 0.5  AND  x > 0.8   (the W1 pocket's home)
  * tip region:        z > 0.8                  (P13 tip free edge)
  * sheet corridor:    |s| < 0.03 (level-set) AND z < 0.5 AND x > 0.8
    (the pocket proper, on the wake sheet)

Run:  python cases/analysis/b23_junction_discriminator/run_d2b_regionmax.py
Artifacts: results/d2b_regionmax.csv
"""

import csv

import numpy as np

from wb_common import (ALPHA_REF, FUS, LS_MESH_DIR, OUT, build_ls, load_mesh)

HERE = OUT.parent
MESH_DIR = HERE / "meshes"

CASES = {"j100": (MESH_DIR / "j100.msh", 0.015),
         "j050": (LS_MESH_DIR / "medium.msh", 0.0075),
         "j025": (MESH_DIR / "j025.msh", 0.00375)}


def region_max(mesh, mvop, phi_ext, mask):
    m2 = np.asarray(mvop.element_mach2(phi_ext, 0.5, 1.4, 1.0))
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    if not mask.any():
        return np.nan, np.nan, np.nan, 0
    m2m = np.where(mask, m2, -1.0)
    i0 = int(np.argmax(m2m))
    c = cents[i0]
    return float(np.sqrt(m2[i0])), float(c[0]), float(c[2]), int(mask.sum())


def main():
    rows = []
    for tag, (mp, h_j) in CASES.items():
        cache = OUT / f"d2_{tag}_a{ALPHA_REF:.2f}_m0.50.npz"
        if not (mp.exists() and cache.exists()):
            print(f"[skip] {tag}: mesh or cache absent", flush=True)
            continue
        mesh = load_mesh(mp)
        wls, _, mvop = build_ls(mesh, ALPHA_REF)
        d = np.load(cache, allow_pickle=True)
        phi_ext = d["phi_ext"]
        cents = mesh.nodes[mesh.elements].mean(axis=1)
        z, x = cents[:, 2], cents[:, 0]
        s, _, _ = wls.evaluate(cents)

        junc = (z < 0.5) & (x > 0.8)
        tip = z > 0.8
        sheet = (np.abs(s) < 0.03) & junc
        mj = region_max(mesh, mvop, phi_ext, junc)
        mt = region_max(mesh, mvop, phi_ext, tip)
        ms = region_max(mesh, mvop, phi_ext, sheet)
        rows.append(dict(tag=tag, h_junction=h_j,
                         junc_mmax=mj[0], junc_x=mj[1], junc_z=mj[2],
                         tip_mmax=mt[0], tip_x=mt[1], tip_z=mt[2],
                         sheet_mmax=ms[0], sheet_x=ms[1], sheet_z=ms[2],
                         n_sheet=ms[3]))
        print(f"  [{tag} h_j={h_j:g}] junction={mj[0]:.2f}@({mj[1]:.2f},{mj[2]:.2f}) "
              f"tip={mt[0]:.2f}@({mt[1]:.2f},{mt[2]:.2f}) "
              f"sheet={ms[0]:.2f}@({ms[1]:.2f},{ms[2]:.2f})", flush=True)

    keys = ["tag", "h_junction", "junc_mmax", "junc_x", "junc_z",
            "tip_mmax", "tip_x", "tip_z",
            "sheet_mmax", "sheet_x", "sheet_z", "n_sheet"]
    with open(OUT / "d2b_regionmax.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'d2b_regionmax.csv'}")


if __name__ == "__main__":
    main()
