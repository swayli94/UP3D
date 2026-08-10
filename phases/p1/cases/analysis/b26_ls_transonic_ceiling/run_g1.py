"""G1 -- the pocket-healed LS transonic ceiling re-measurement (B26
PRE_REGISTRATION.md, 2026-07-19).

Legs (pre-registered execution order):
  A_coarse  same-code control + the B18 reproduction gate (checks.csv row 3:
            died at m=0.55, Mmax=1.31)
  C_coarse  fragment clip (smoke + trend)
  A_medium  same-code control + the GB20.5 anchor (checks.csv row 5 /
            legb_b18_hypothesis.csv main: died at m=0.50, Mmax=5.22,
            nlim/nflr 3/3, res 1.1e-13)
  C_medium  the MAIN judgement leg

Both sides share the nominal ladder m_start=0.50 -> m_target=0.84, dm=0.05
(0.84 = the conforming coarse ceiling, the physical upper reference); the
ramp stops honestly (a failed level halves dm down to dm_min=0.01). Note the
B21/B22 freeze-capture fixes moved the code since B18 was committed: the A
re-run is the PRIMARY control (P14 same-code discipline), the committed CSVs
are the historical anchor -- a significant A-vs-anchor drift is recorded as
an independent B21/B22 finding, not a B26 failure (risk T1).

Verdicts (pre-reg section 3):
  B26-A (ceiling raised):  C medium m_last_converged >= 0.60 (past B18's
      0.50 death point by >= 2 rungs, warm-start chain intact), converged
      levels free of class-(a) failures, cl trend consistent with the
      conforming anchors -> the pocket WAS the ceiling limiter.
  B26-B (ceiling unmoved): C medium dies at the same level as A medium with
      class-(a) features at the junction site -> re-attribute (G1.6/P11,
      wake-LS conditioning, upwind-strip interaction).
  B26-C (new failure mode): C climbs past A's death point but dies of
      class (b)/(c) or at a new site.

Run:  python cases/analysis/b26_ls_transonic_ceiling/run_g1.py [A_coarse C_coarse ...]
Artifacts: results/g1_summary.csv, results/g1_levels.csv,
           results/g1_peaks.csv, results/g1_ceiling.png
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb26 import (ALPHA, B18_DEATH, B18_DEATH_MMAX, CONF_CL, DM, LS_MESH_DIR,
                  M_START, M_TARGET, NOMINAL, OUT, classify, corridor_peaks,
                  load_mesh, measure_leg, run_ramp)

LEGS = [("A", "coarse"), ("C", "coarse"), ("A", "medium"), ("C", "medium")]


def main():
    only = sys.argv[1:] or None
    rows, lv_rows, pk_rows = [], [], []
    for side, level in LEGS:
        if only and f"{side}_{level}" not in only:
            continue
        mesh = load_mesh(LS_MESH_DIR / f"{level}.msh")
        print(f"=== G1 {side} {level}: {len(mesh.elements)} tets, "
              f"alpha={ALPHA}, {M_START}->{M_TARGET} dm={DM} ===", flush=True)
        rec, wls, mvop = run_ramp(mesh, level, side)
        m = measure_leg(mesh, wls, mvop, rec["phi_ext"], rec["m_final"], level)
        cls = classify(rec)
        die = rec["levels"][-1]
        mlc = rec["m_last"]
        rows.append(dict(
            side=side, level=level, cls=cls,
            target_reached=rec["target_reached"],
            m_last_converged=(mlc if mlc is not None else ""),
            m_final=rec["m_final"],
            n_levels_run=len(rec["levels"]),
            n_halvings=len(rec["mach_schedule"]) - len(NOMINAL),
            wall_s=round(rec["wall_s"], 1),
            die_m=die["m_inf"], die_tag=die["tag"],
            die_converged=die["converged"], die_res=die["residual_norm"],
            die_nlim=die["n_limited"], die_nflr=die["n_floored"],
            die_n_newton=die["n_newton"], die_mmax=die["mach_max"],
            die_gamma=die["gamma"], die_accept=die["accept_reason"],
            **m))
        for lv in rec["levels"]:
            lv_rows.append(dict(side=side, level=level, **lv))
        for p in corridor_peaks(mesh, wls, mvop, rec["phi_ext"],
                                rec["m_final"]):
            pk_rows.append(dict(side=side, level=level, **p))
        print(f"  [{side} {level}] cls={cls} reached={rec['target_reached']} "
              f"m_last={mlc} m_final={rec['m_final']} "
              f"die(res={die['residual_norm']:.1e} nlim={die['n_limited']} "
              f"nflr={die['n_floored']} Mmax={die['mach_max']:.2f}) "
              f"corrM={m['corr_mmax']:.2f}@x={m['corr_x']:.2f},"
              f"q={m['corr_q']:+.3f} n_sup={m['n_sup']} "
              f"pkM={m['pk_mmax']:.2f} dTEj={m['pk_dist_te_junc']:.3f} "
              f"cl_p={m['cl_p']:.4f} strip_jmp={m['strip_aux_jump_max']:.2e}",
              flush=True)

    if not rows:
        print("nothing ran")
        sys.exit(0)
    _write_summary(rows)
    _write_levels(lv_rows)
    _write_peaks(pk_rows)
    _write_fig(rows)
    _verdict_hint(rows)


def _write_summary(rows):
    keys = ["side", "level", "cls", "target_reached", "m_last_converged",
            "m_final", "n_levels_run", "n_halvings", "wall_s",
            "die_m", "die_tag", "die_converged", "die_res", "die_nlim",
            "die_nflr", "die_n_newton", "die_mmax", "die_gamma", "die_accept",
            "pk_mmax", "pk_x", "pk_y", "pk_z", "pk_dist_fus",
            "pk_z_minus_zjunc", "pk_dist_te_junc",
            "corr_mmax", "corr_x", "corr_z", "corr_q", "n_sup", "n_sup_corr",
            "sheet_mmax", "sheet_x", "sheet_z",
            "tip_mmax", "tip_x", "tip_z",
            "top_med_abs_s", "top_med_q", "top_med_x",
            "pocket_peak_x", "pocket_past_tail",
            "cl_p", "cl_fus", "cl_fus_band", "cl_fus_out", "cl_fus_poles",
            "n_te_nodes", "n_aux_farfield", "n_aux_symmetry", "n_cut",
            "n_cut_wing", "n_cut_inboard_body", "n_cut_inboard_sym",
            "sliver_dih_min_deg", "sliver_dih_p05_deg",
            "sliver_vol_min", "sliver_vol_p05",
            "strip_aux_jump_max", "strip_aux_jump_p95",
            "strip_aux_jump_over_gamma"]
    with open(OUT / "g1_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'g1_summary.csv'}")


def _write_levels(lv_rows):
    keys = ["side", "level", "m_inf", "tag", "converged", "accept_reason",
            "n_newton", "residual_norm", "n_limited", "n_floored",
            "mach_max", "gamma", "cl_kj", "froze", "n_freeze_refresh",
            "n_freeze_reverts", "n_refactor", "n_schur_fallback",
            "n_lin_iters", "n_lin_solves", "wall_s"]
    with open(OUT / "g1_levels.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(lv_rows)
    print(f"wrote {OUT/'g1_levels.csv'}")


def _write_peaks(pk_rows):
    keys = ["side", "level", "rank", "mach", "x", "y", "z", "q",
            "dist_fus_surface", "z_minus_zjunc", "dist_te_junc"]
    with open(OUT / "g1_peaks.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(pk_rows)
    print(f"wrote {OUT/'g1_peaks.csv'}")


def _write_fig(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, level in zip(axes, ("coarse", "medium")):
        rr = {r["side"]: r for r in rows if r["level"] == level}
        xs = np.arange(2)
        ml = [rr[s]["m_last_converged"] if s in rr and
              rr[s]["m_last_converged"] != "" else 0.0 for s in "AC"]
        bars = ax.bar(xs - 0.05, ml, width=0.55,
                      color=["0.55", "tab:blue"], alpha=0.85)
        for x, v, s in zip(xs, ml, "AC"):
            if s not in rr:
                continue
            ax.text(x - 0.05, max(v, 0.405) + 0.008,
                    f"{v:.2f}" if v else "none", ha="center", fontsize=9)
            ax.plot(x - 0.05, float(rr[s]["m_final"]), "x",
                    color="tab:red", ms=10, mew=2.2)
        ax.axhline(B18_DEATH[level], color="tab:red", ls=":", lw=1.2,
                   label=f"B18 death M{B18_DEATH[level]:.2f} "
                         f"(Mmax {B18_DEATH_MMAX[level]})")
        ax.axhline(0.84 if level == "coarse" else 0.79, color="0.4",
                   ls="--", lw=1.2,
                   label="conforming ceiling "
                         f"({0.84 if level == 'coarse' else 0.79})")
        ax.axhline(M_TARGET, color="0.75", ls="-", lw=0.8)
        ax.set_xticks(xs)
        ax.set_xticklabels(["A (clip=None)", "C (fragment clip)"])
        ax.set_title(f"{level}: LS ceiling A vs C (x = died-at Mach)")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8, loc="lower right")
    axes[0].set_ylabel("m_last_converged")
    axes[0].set_ylim(0.4, 0.9)
    fig.suptitle(f"B26 G1: pocket-healed LS transonic ceiling, "
                 f"alpha={ALPHA}, {M_START}->{M_TARGET} dm={DM}")
    fig.tight_layout()
    fig.savefig(OUT / "g1_ceiling.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT/'g1_ceiling.png'}")


def _verdict_hint(rows):
    """Print the pre-reg section-3 verdict ingredients (the VERDICT.md writeup
    is a separate analysis step)."""
    print("\n--- verdict hint (pre-reg section 3) ---")
    by = {(r["side"], r["level"]): r for r in rows}
    for level in ("coarse", "medium"):
        a, c = by.get(("A", level)), by.get(("C", level))
        if a is None or c is None:
            continue
        ml_c = c["m_last_converged"]
        print(f"  {level}: A m_last={a['m_last_converged'] or 'none'} "
              f"cls={a['cls']} die@{a['die_m']} | C m_last={ml_c or 'none'} "
              f"cls={c['cls']} die@{c['die_m']} "
              f"(B18 anchor death M{B18_DEATH[level]:.2f})")
    med = by.get(("C", "medium"))
    if med is not None and med["m_last_converged"] != "":
        ml = float(med["m_last_converged"])
        if ml >= 0.60:
            print(f"  C medium m_last={ml:.2f} >= 0.60 -> B26-A candidate "
                  f"(check cl_kj trend vs conforming anchors {CONF_CL} in "
                  f"g1_levels.csv)")
        elif med["die_m"] == by.get(("A", "medium"), {}).get("die_m") \
                and med["cls"].startswith("a"):
            print(f"  C medium died at the same level as A with class-(a) "
                  f"features -> B26-B candidate (re-attribute)")
        else:
            print(f"  C medium: neither A nor B pattern -> B26-C candidate")
    elif med is not None:
        a = by.get(("A", "medium"))
        same = a is not None and med["die_m"] == a["die_m"]
        site = "same level as A" if same else "NOT the A death point"
        cand = ("B26-B candidate" if same and med["cls"].startswith("a")
                else "B26-C candidate")
        print(f"  C medium m_last=none, died@{med['die_m']} "
              f"cls={med['cls']} ({site}) -> {cand}")


if __name__ == "__main__":
    main()
