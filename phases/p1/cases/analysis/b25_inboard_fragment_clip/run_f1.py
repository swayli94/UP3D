"""F1 -- inboard fragment clip A/C: does moving the sheet's inboard end onto
the fuselage surface / symmetry plane (the conforming topology) kill the
junction pocket? (PRE_REGISTRATION.md, 2026-07-19)

Legs (LS Picard M0.5, freestream+pin_gamma, the b23 recipe; no knob -- a
clean A/C pair):
  A (control):   current q >= 0 clip; solves = committed b23 D1 caches,
                 loaded read-only, measured with the SAME code. The b23 set
                 has no coarse alpha=2 leg: that one control is solved with
                 the unchanged default path (bit-identical code, S6) into
                 b25's own cache.
  C (treatment): inboard_clip = make_inboard_clip(FUS) -- the sheet runs
                 inboard to the fuselage surface (trace on the wall) and,
                 aft of the body, to the z = 0 symmetry plane.
  alpha=0 legs are the self-check (new cut strips beside the body / on the
  symmetry plane, but Gamma ~ 0 -- must be inert: |d cl_p| <= 1%, no new
  pocket, normal convergence).

Pre-registered criteria (section 3):
  C-A (healed):  medium a=3.06 corridor Mmax <= 1.3 AND corridor n_sup = 0;
                 cl_p within 2%, |dgamma| <= 5%, root profile distortion
                 <= 5%; alpha=0 inert; n_outer <= 1.5x A, res <= 1e-7;
                 outside-band carryover |d| <= 20%.
  C-B:           pocket gone but a NEW peak at the body trace / sliver
                 metrics degrade together ("topology right, wall-grazing
                 implementation sick").
  C-C:           pocket unchanged -> free-edge hypothesis challenged a
                 third time (very unlikely; B23/B24 double-confirmed).
  C-D:           alpha=0 leg NOT inert -> the clip introduces a new
                 contamination source; roll back and analyse.

Run:  python cases/analysis/b25_inboard_fragment_clip/run_f1.py [A_coarse C_coarse ...]
Artifacts: results/f1_summary.csv, results/f1_pocket.png
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb25 import (ALPHA_REF, LS_MESH_DIR, M_INF, OUT, X_TAIL_TIP, load_mesh,
                  measure_f1, root_distortion, solve_side, te_profile)

ALPHAS = [0.0, 2.0, ALPHA_REF]
SIDES = ("A", "C")
LEVELS = ("coarse", "medium")


def main():
    only = sys.argv[1:] or None
    rows = []
    profs = {}
    for level in LEVELS:
        mesh = load_mesh(LS_MESH_DIR / f"{level}.msh")
        for side in SIDES:
            if only and f"{side}_{level}" not in only and side not in only:
                continue
            print(f"=== F1 {side} {level}: {len(mesh.elements)} tets, "
                  f"alphas={ALPHAS} ===", flush=True)
            for a in ALPHAS:
                rec, mvop, wls = solve_side(mesh, level, side, a)
                m = measure_f1(mesh, mvop, wls, rec["phi_ext"], a, level)
                profs[(side, level, a)] = te_profile(mesh, mvop,
                                                     rec["phi_ext"])
                rows.append(dict(side=side, level=level, alpha=a,
                                 converged=rec["conv"], n_outer=rec["n"],
                                 wall_s=rec["wall_s"], gamma=rec["gamma"],
                                 nlim=rec["nlim"], nflr=rec["nflr"],
                                 res_final=rec["res"][-1] if rec["res"]
                                 else np.nan,
                                 **m))
                print(f"  [{side} {level} a={a:>4}] corrM={m['corr_mmax']:.2f}"
                      f"@x={m['corr_x']:.2f} q={m['corr_q']:+.3f} "
                      f"n_sup={m['n_sup']} cl_p={m['cl_p']:.4f} "
                      f"gamma={rec['gamma']:.4f} n_cut={m['n_cut']} "
                      f"(wing {m['n_cut_wing']} body {m['n_cut_inboard_body']}"
                      f" sym {m['n_cut_inboard_sym']}) "
                      f"dih_min={m['sliver_dih_min_deg']:.2f}deg "
                      f"strip_jmp={m['strip_aux_jump_max']:.2e}",
                      flush=True)

    if not rows:
        print("nothing ran")
        sys.exit(0)
    _add_deltas(rows, profs)
    _write_summary(rows)
    _write_fig(rows)


def _add_deltas(rows, profs):
    """Pairwise C - A deltas on the C rows (pre-reg 2.3.3 + section 3)."""
    for row in rows:
        if row["side"] != "C":
            continue
        key_a = ("A", row["level"], row["alpha"])
        row_a = next((r for r in rows
                      if (r["side"], r["level"], r["alpha"]) == key_a), None)
        if row_a is None or key_a not in profs:
            continue
        d_cl = row["cl_p"] - row_a["cl_p"]
        row["d_cl_p"] = d_cl
        row["d_cl_p_rel"] = (abs(d_cl) / abs(row_a["cl_p"])
                             if abs(row_a["cl_p"]) > 1e-8 else np.nan)
        g_a, g_c = row_a["gamma"], row["gamma"]
        row["d_gamma_rel"] = (abs(g_c - g_a) / abs(g_a)
                              if abs(g_a) > 1e-12 else np.nan)
        row["root_jump_distortion"] = root_distortion(
            *profs[key_a], *profs[("C", row["level"], row["alpha"])])
        row["d_cl_fus_out_rel"] = (
            abs(row["cl_fus_out"] - row_a["cl_fus_out"])
            / max(abs(row_a["cl_fus_out"]), 1e-12))
        print(f"  [delta {row['level']} a={row['alpha']:>4}] "
              f"d_cl_p={d_cl:+.5f} ({row['d_cl_p_rel']:+.2%}) "
              f"d_gamma={row['d_gamma_rel']:+.2%} "
              f"root_distort={row['root_jump_distortion']:.2%}",
              flush=True)


def _write_summary(rows):
    keys = ["side", "level", "alpha", "converged", "n_outer", "wall_s",
            "gamma", "nlim", "nflr", "res_final",
            "d_cl_p", "d_cl_p_rel", "d_gamma_rel", "root_jump_distortion",
            "d_cl_fus_out_rel",
            "n_te_nodes", "n_aux_farfield", "n_aux_symmetry", "n_cut",
            "n_cut_wing", "n_cut_inboard_body", "n_cut_inboard_sym",
            "sliver_dih_min_deg", "sliver_dih_p05_deg",
            "sliver_vol_min", "sliver_vol_p05",
            "strip_aux_jump_max", "strip_aux_jump_p95",
            "strip_aux_jump_over_gamma",
            "corr_mmax", "corr_x", "corr_z", "corr_q", "sheet_mmax",
            "sheet_x", "sheet_z", "tip_mmax", "tip_x", "tip_z",
            "n_sup", "n_sup_corr", "top_med_abs_s", "top_med_q",
            "top_med_x", "pocket_peak_x", "pocket_past_tail",
            "cl_p", "cl_fus", "cl_fus_band", "cl_fus_out", "cl_fus_poles"]
    with open(OUT / "f1_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'f1_summary.csv'}")


def _write_fig(rows):
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    style = {"A": ("o", "--"), "C": ("s", "-")}
    for level, col in (("coarse", "tab:blue"), ("medium", "tab:red")):
        for side in SIDES:
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
            ax[0][2].plot(a, [r["n_sup"] for r in rr], marker=mk, ls=ls,
                          color=col, label=lab)
            ax[1][0].plot(a, [r["cl_p"] for r in rr], marker=mk, ls=ls,
                          color=col, label=lab)
            ax[1][1].plot(a, [r["cl_fus_band"] for r in rr], marker=mk,
                          ls=ls, color=col, label=lab + " band")
            ax[1][1].plot(a, [r["cl_fus_out"] for r in rr], marker="^",
                          ls=":", color=col, alpha=0.6, label=lab + " out")
            ax[1][2].plot(a, [r["n_cut"] for r in rr], marker=mk, ls=ls,
                          color=col, label=lab + " total")
            if side == "C":
                ax[1][2].plot(a, [r["n_cut_inboard_body"] for r in rr],
                              marker="^", ls=":", color=col,
                              label=lab + " inboard(body)")
                ax[1][2].plot(a, [r["n_cut_inboard_sym"] for r in rr],
                              marker="v", ls=":", color=col,
                              label=lab + " inboard(sym)")
    ax[0][0].axhline(1.3, color="0.5", ls="--", lw=0.8, label="C-A threshold")
    ax[0][0].set_ylabel("junction-corridor Mmax (z<0.5, x>0.8)")
    ax[0][0].set_title("F1 primary: pocket strength")
    ax[0][1].axvline(X_TAIL_TIP, color="0.5", ls="--", lw=0.8,
                     label=f"x_tail_tip={X_TAIL_TIP:.2f}")
    ax[0][1].set_ylabel("x of corridor Mmax element")
    ax[0][1].set_title("F1 primary: pocket location")
    ax[0][2].set_ylabel("n supersonic elements")
    ax[0][2].set_title("F1 primary: supersonic census")
    ax[1][0].set_ylabel("cl_p (wing)")
    ax[1][0].set_title("guardrail: wing load unchanged (C-A <= 2%)")
    ax[1][1].axhline(0.0, color="0.5", lw=0.8)
    ax[1][1].set_ylabel("cl_fus contribution")
    ax[1][1].set_title("band vs outside (carryover) decomposition")
    ax[1][2].set_ylabel("cut-element census")
    ax[1][2].set_title("the new inboard strips (body / symmetry)")
    for a in ax.flat:
        a.set_xlabel("alpha (deg)")
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    fig.suptitle("F1 A/C: inboard fragment clip (conforming topology) vs "
                 "q>=0 control (LS Picard M0.5, freestream+pin_gamma)")
    fig.tight_layout()
    fig.savefig(OUT / "f1_pocket.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'f1_pocket.png'}")


if __name__ == "__main__":
    main()
