"""D1 -- alpha sweep on the EXISTING wing-body LS meshes (no new meshes).

Question (PRE_REGISTRATION.md): does the junction pocket need lift/wake, or
is it already there at alpha = 0 (pure geometry)?

  * alpha = 0 pocket (Mmax >> 1)          -> geometry-driven
  * pocket appears/grows with alpha       -> lift/wake/Kutta-coupled
  * self-check: cl_fus(alpha=0) ~ 0 (symmetry)

Run:  python cases/analysis/wb_junction_discriminator/run_d1.py
Artifacts: results/d1_summary.csv, results/d1_topk_*.csv, results/d1_alpha.png
"""

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb_common import (ALPHA_REF, LS_MESH_DIR, OUT, load_mesh, measure,
                       solve_ls)

#: (level, alphas). Medium carries the full sweep (the pocket's home mesh);
#: coarse gets the two endpoints as the refinement-trend anchor.
SWEEP = {"medium": [0.0, 1.0, 2.0, ALPHA_REF],
         "coarse": [0.0, ALPHA_REF]}


def main():
    levels = sys.argv[1:] or list(SWEEP)
    rows = []
    for level, alphas in ((lv, SWEEP[lv]) for lv in levels):
        mp = LS_MESH_DIR / f"{level}.msh"
        if not mp.exists():
            print(f"[skip] {level}: mesh absent ({mp})")
            continue
        mesh = load_mesh(mp)
        print(f"=== D1 {level}: {len(mesh.elements)} tets, "
              f"alphas={alphas} ===", flush=True)
        for a in alphas:
            rec, mvop = solve_ls(mesh, f"d1_{level}", a)
            m, top = measure(mesh, mvop, rec["phi_ext"], a)
            rows.append(dict(level=level, alpha=a, n_tets=len(mesh.elements),
                             converged=rec["conv"], n_outer=rec["n"],
                             wall_s=rec["wall_s"], **m))
            print(f"  [{level} a={a:>4}] Mmax={m['mmax']:.2f} at "
                  f"(x={m['argmax_x']:.3f}, y={m['argmax_y']:.3f}, "
                  f"z={m['argmax_z']:.3f}; |r-R|={m['argmax_dist_fus']:.4f}, "
                  f"z-z_j={m['argmax_z_minus_zjunc']:+.4f}, "
                  f"d_te={m['argmax_dist_te_junc']:.3f}) "
                  f"n_sup={m['n_sup']} cl_fus/wing={m['cl_fus_over_wing']:.3f}",
                  flush=True)
            _write_topk(level, a, top)

    if not rows:
        print("no meshes found; nothing ran")
        sys.exit(0)
    _write_summary(rows)
    _write_fig(rows)


def _write_topk(level, alpha, top):
    import csv
    p = OUT / f"d1_topk_{level}_a{alpha:.2f}.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(top[0]))
        w.writeheader()
        w.writerows(top)


def _write_summary(rows):
    import csv
    keys = ["level", "alpha", "n_tets", "converged", "n_outer", "wall_s",
            "mmax", "m2_max", "argmax_elem", "argmax_x", "argmax_y", "argmax_z",
            "argmax_dist_fus", "argmax_z_minus_zjunc", "argmax_dist_te_junc",
            "n_sup", "sup_x_min", "sup_x_max", "sup_y_min", "sup_y_max",
            "sup_z_min", "sup_z_max", "cl_p", "cl_fus", "cl_fus_over_wing"]
    with open(OUT / "d1_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'d1_summary.csv'}")


def _write_fig(rows):
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.4))
    for level, col in (("coarse", "tab:blue"), ("medium", "tab:red")):
        rr = sorted((r for r in rows if r["level"] == level),
                    key=lambda r: r["alpha"])
        if not rr:
            continue
        a = [r["alpha"] for r in rr]
        ax[0].plot(a, [r["mmax"] for r in rr], "o-", color=col, label=level)
        ax[1].plot(a, [r["n_sup"] for r in rr], "o-", color=col, label=level)
        ax[2].plot(a, [r["cl_fus_over_wing"] for r in rr], "o-", color=col,
                   label=level)
    ax[0].axhline(1.0, color="0.5", ls="--", lw=0.8)
    ax[0].set_ylabel("Mmax (element M^2 max, sqrt)")
    ax[0].set_title("pocket strength vs alpha")
    ax[1].set_ylabel("# supersonic elements (M^2 > 1)")
    ax[1].set_title("pocket extent vs alpha")
    ax[2].set_ylabel("|cl_fus| / cl_p_wing")
    ax[2].set_title("fuselage spurious lift vs alpha")
    for a in ax[:3]:
        a.set_xlabel("alpha (deg)")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    # pocket-location map: argmax (x, z) per case against the planform
    for level, col, mk in (("coarse", "tab:blue", "o"), ("medium", "tab:red", "s")):
        for r in rows:
            if r["level"] != level:
                continue
            ax[3].plot(r["argmax_x"], r["argmax_z"], mk, color=col, ms=9,
                       fillstyle="none", mew=1.6,
                       label=f"{level} a={r['alpha']}")
    ax[3].axvspan(0.0, 0.806, color="tab:green", alpha=0.10,
                  label="wing root chord (the CREASE)")
    ax[3].axhline(0.15, color="0.5", ls="--", lw=0.8, label="z_junc = r_f")
    ax[3].axhline(1.1963, color="0.5", ls=":", lw=0.8, label="tip")
    ax[3].set_xlabel("x of Mmax element"); ax[3].set_ylabel("z of Mmax element")
    ax[3].set_title("WHERE is the pocket? (crease = x<0.81)")
    ax[3].grid(alpha=0.3); ax[3].legend(fontsize=7)
    fig.suptitle("D1: does the junction pocket need lift? "
                 "(LS Picard M0.5, freestream+pin_gamma)")
    fig.tight_layout()
    fig.savefig(OUT / "d1_alpha.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'d1_alpha.png'}")


if __name__ == "__main__":
    main()
