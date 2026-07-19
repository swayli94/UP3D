"""E1 -- waterline-extension A/B: does killing the near-field inboard free
edge kill the junction pocket? (PRE_REGISTRATION.md, 2026-07-19)

Legs (LS Picard M0.5, freestream+pin_gamma, the b23 recipe):
  A  (control):   current 2-point TE polyline; solves = committed b23 D1
                  caches, loaded read-only, measured with the SAME code.
  B1 (treatment): te_polyline(extend="waterline") -- the sheet's inboard
                  free edge pushed to the far field.

Pre-registered criteria:
  E1: B corridor Mmax <= 1.3 AND corridor n_sup = 0 (pocket GONE), or the
      corridor peak moves past x_tail (pocket MOVED to the far field).
  E2: |cl_p| within 2% of A, dgamma <= 5%, tip Mmax within [0.5, 2]x,
      alpha=0 clean, n_outer <= 1.5x A, n_te_nodes unchanged (asserted in
      wb24), nlim+nflr within +10 of the A reference (B20 committed anchor
      3/3 at medium M0.5; the b23 caches carry no limiter counters).
  E3: band cl_fus collapses (<= 0.3x A), outside within +/-20%.

Run:  python cases/analysis/b24_wake_inboard_end/run_e1.py
Artifacts: results/e1_summary.csv, results/e1_w2.csv, results/e1_pocket.png
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb24 import (ALPHA_REF, BW_LADDER, LS_MESH_DIR, M_INF, OUT, X_TAIL,
                  load_mesh, measure_e1, solve_side)

#: (level, alphas) per side. alpha=0 legs are the self-check (no new
#: contamination source from the extension); alpha=1 skipped (D1 showed the
#: transition already).
LEGS = {"A": {"medium": [0.0, 2.0, ALPHA_REF], "coarse": [0.0, ALPHA_REF]},
        "B1": {"medium": [0.0, 2.0, ALPHA_REF], "coarse": [0.0, ALPHA_REF]}}


def main():
    only = sys.argv[1:] or None
    rows = []
    for side in ("A", "B1"):
        for level, alphas in LEGS[side].items():
            if only and f"{side}_{level}" not in only and side not in only:
                continue
            mesh = load_mesh(LS_MESH_DIR / f"{level}.msh")
            print(f"=== E1 {side} {level}: {len(mesh.elements)} tets, "
                  f"alphas={alphas} ===", flush=True)
            for a in alphas:
                rec, mvop, wls = solve_side(mesh, level, side, a)
                m = measure_e1(mesh, mvop, wls, rec["phi_ext"], a, level)
                rows.append(dict(side=side, level=level, alpha=a,
                                 converged=rec["conv"], n_outer=rec["n"],
                                 wall_s=rec["wall_s"], gamma=rec["gamma"],
                                 nlim=rec["nlim"], nflr=rec["nflr"],
                                 res_final=rec["res"][-1] if rec["res"] else np.nan,
                                 **m))
                print(f"  [{side} {level} a={a:>4}] corrM={m['corr_mmax']:.2f}"
                      f"@x={m['corr_x']:.2f} sheetM={m['sheet_mmax']:.2f} "
                      f"tipM={m['tip_mmax']:.2f} n_sup={m['n_sup']} "
                      f"cl_p={m['cl_p']:.4f} cl_fus={m['cl_fus']:+.5f} "
                      f"(band {m['cl_fus_band']:+.5f}) gamma={rec['gamma']:.4f}",
                      flush=True)

    if not rows:
        print("nothing ran")
        sys.exit(0)
    _write_summary(rows)
    _write_w2(rows)
    _write_fig(rows)


def _write_summary(rows):
    keys = ["side", "level", "alpha", "converged", "n_outer", "wall_s",
            "gamma", "nlim", "nflr", "res_final",
            "n_te_nodes", "n_aux_farfield", "n_cut",
            "corr_mmax", "corr_x", "corr_z", "sheet_mmax", "sheet_x",
            "sheet_z", "tip_mmax", "tip_x", "tip_z", "n_sup", "n_sup_corr",
            "top_med_abs_s", "top_med_q", "top_med_x",
            "pocket_peak_x", "pocket_past_tail",
            "cl_p", "cl_fus", "cl_fus_band", "cl_fus_out", "cl_fus_poles"]
    with open(OUT / "e1_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'e1_summary.csv'}")


def _write_w2(rows):
    keys = ["side", "level", "alpha", "cl_fus", "cl_fus_band", "cl_fus_out",
            "cl_fus_poles"]
    with open(OUT / "e1_w2.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: r[k] for k in keys} for r in rows])
    print(f"wrote {OUT/'e1_w2.csv'}")


def _write_fig(rows):
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    style = {"A": ("o", "--"), "B1": ("s", "-")}
    for level, col in (("coarse", "tab:blue"), ("medium", "tab:red")):
        for side in ("A", "B1"):
            rr = sorted((r for r in rows
                         if r["level"] == level and r["side"] == side),
                        key=lambda r: r["alpha"])
            if not rr:
                continue
            a = [r["alpha"] for r in rr]
            mk, ls = style[side]
            lab = f"{level} {side}"
            ax[0][0].plot(a, [r["corr_mmax"] for r in rr], marker=mk, ls=ls,
                          color=col, label=lab)
            ax[0][1].plot(a, [r["corr_x"] for r in rr], marker=mk, ls=ls,
                          color=col, label=lab)
            ax[1][0].plot(a, [r["cl_p"] for r in rr], marker=mk, ls=ls,
                          color=col, label=lab)
            ax[1][1].plot(a, [r["cl_fus_band"] for r in rr], marker=mk,
                          ls=ls, color=col, label=lab + " band")
            ax[1][1].plot(a, [r["cl_fus_out"] for r in rr], marker="^",
                          ls=":", color=col, alpha=0.6, label=lab + " out")
    ax[0][0].axhline(1.3, color="0.5", ls="--", lw=0.8, label="E1 threshold")
    ax[0][0].set_ylabel("junction-corridor Mmax (z<0.5, x>0.8)")
    ax[0][0].set_title("E1 primary: pocket strength")
    ax[0][1].axvline(X_TAIL, color="0.5", ls="--", lw=0.8,
                     label=f"x_tail_start={X_TAIL:.2f}")
    ax[0][1].set_ylabel("x of corridor Mmax element")
    ax[0][1].set_title("E1 primary: pocket location (right of line = far field)")
    ax[1][0].set_ylabel("cl_p (wing)")
    ax[1][0].set_title("E2 guardrail: wing load unchanged")
    ax[1][1].axhline(0.0, color="0.5", lw=0.8)
    ax[1][1].set_ylabel("cl_fus contribution")
    ax[1][1].set_title("E3: band collapses, outside (carryover) unchanged")
    for a in ax.flat:
        a.set_xlabel("alpha (deg)")
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    fig.suptitle("E1 A/B: waterline extension vs current free-edge sheet "
                 "(LS Picard M0.5, freestream+pin_gamma)")
    fig.tight_layout()
    fig.savefig(OUT / "e1_pocket.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'e1_pocket.png'}")


if __name__ == "__main__":
    main()
