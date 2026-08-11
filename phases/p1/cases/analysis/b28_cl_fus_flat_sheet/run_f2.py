"""B28 F2: the flat-vs-tilted decoupling legs + the pre-registered verdicts.

Legs (PRE_REGISTRATION.md section 3.1): F medium a=3.06 (decisive),
F medium a=2.0 (trend), F coarse a=3.06 (level endpoint). Comparators are
the committed b25 f1_summary.csv (A/C) and w2_conf.csv (conforming oracle),
read-only.

Decision tree (section 3.3, TOL = 15%):
  F1 (position sole factor): |F_out - oracle| <= .15|oracle| AND
                             |F_out - C_out| > .15|C_out|
  F2 (position exonerated):  |F_out - C_out| <= .15|C_out| AND
                             |F_out - oracle| > .15|oracle|
  F3 (mixed):                otherwise; report r = |C_out - F_out| /
                             |C_out - oracle|

Artifacts: results/f2_summary.csv, results/f2_decomposition.png
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))

from wb28 import (B25_RESULTS, OUT, load_c_te_profile, load_mesh, measure_f1,
                  solve_flat)
from wb25 import root_distortion, te_profile

LEGS = (("medium", 3.06), ("medium", 2.0), ("coarse", 3.06))
TOL = 0.15
MESH = {lv: REPO_ROOT / "cases/meshes/onera_m6_wingbody" / f"{lv}.msh"
        for lv in ("coarse", "medium")}


def committed_comparators():
    """b25 committed A/C rows keyed by (side, level, alpha) + the oracle."""
    comp = {}
    with open(B25_RESULTS / "f1_summary.csv") as fh:
        for row in csv.DictReader(fh):
            key = (row["side"], row["level"], float(row["alpha"]))
            comp[key] = {k: float(v) for k, v in row.items()
                         if k not in ("side", "level", "converged",
                                      "pocket_past_tail") and v != ""}
    with open(B25_RESULTS / "w2_conf.csv") as fh:
        oracle = next(r for r in csv.DictReader(fh)
                      if r["path"] == "conforming_pressure"
                      and r["level"] == "medium")
    oracle = {"all": float(oracle["all"]), "band": float(oracle["band_bw0.06"]),
              "out": float(oracle["out_bw0.06"]),
              "poles": float(oracle["poles"])}
    return comp, oracle


def main():
    comp, oracle = committed_comparators()
    rows = []
    for level, alpha in LEGS:
        mesh = load_mesh(MESH[level])
        rec, mvop, wls = solve_flat(mesh, level, alpha)
        m = measure_f1(mesh, mvop, wls, rec["phi_ext"], alpha, level)
        row = dict(side="F", level=level, alpha=alpha,
                   converged=rec["conv"], n_outer=rec["n"],
                   wall_s=rec["wall_s"], gamma=rec["gamma"],
                   nlim=rec["nlim"], nflr=rec["nflr"],
                   res_final=rec["res"][-1], **m)
        c = comp[("C", level, alpha)]
        row["d_cl_p_rel_vs_c"] = abs(row["cl_p"] - c["cl_p"]) / abs(c["cl_p"])
        row["d_gamma_rel_vs_c"] = (abs(row["gamma"] - c["gamma"])
                                   / abs(c["gamma"]))
        row["d_cl_fus_out_rel_vs_c"] = (abs(row["cl_fus_out"]
                                            - c["cl_fus_out"])
                                        / abs(c["cl_fus_out"]))
        prof_f = te_profile(mesh, mvop, rec["phi_ext"])
        prof_c = load_c_te_profile(mesh, level, alpha)
        row["root_jump_distortion_vs_c"] = root_distortion(*prof_c, *prof_f)
        if (level, alpha) == ("medium", 3.06):
            row["d_cl_fus_out_rel_vs_oracle"] = (
                abs(row["cl_fus_out"] - oracle["out"]) / abs(oracle["out"]))
        rows.append(row)
        print(f"  [F {level} a={alpha:>4}] cl_p={row['cl_p']:.5f} "
              f"cl_fus={row['cl_fus']:.5f} band={row['cl_fus_band']:.5f} "
              f"out={row['cl_fus_out']:.5f} poles={row['cl_fus_poles']:.6f} "
              f"| d_cl_p_vs_c={row['d_cl_p_rel_vs_c']:.2%} "
              f"d_out_vs_c={row['d_cl_fus_out_rel_vs_c']:.2%}",
              flush=True)

    _write_summary(rows)
    _write_fig(rows, comp, oracle)
    _verdict(rows, comp, oracle)


def _write_summary(rows):
    keys = ["side", "level", "alpha", "converged", "n_outer", "wall_s",
            "gamma", "nlim", "nflr", "res_final",
            "d_cl_p_rel_vs_c", "d_gamma_rel_vs_c", "d_cl_fus_out_rel_vs_c",
            "root_jump_distortion_vs_c", "d_cl_fus_out_rel_vs_oracle",
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
    with open(OUT / "f2_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'f2_summary.csv'}")


def _write_fig(rows, comp, oracle):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    width = 0.2
    for j, (level, alpha) in enumerate((("medium", 3.06), ("medium", 2.0))):
        a = ax[j]
        groups = ["A", "C", "F"]
        series = {g: comp.get((g, level, alpha)) for g in groups}
        series["F"] = next(r for r in rows
                           if r["level"] == level and r["alpha"] == alpha)
        if j == 0:
            series["oracle"] = {"cl_fus_band": oracle["band"],
                                "cl_fus_out": oracle["out"],
                                "cl_fus": oracle["all"]}
        labels = list(series)
        for i, lab in enumerate(labels):
            s = series[lab]
            vals = [s["cl_fus_band"], s["cl_fus_out"], s["cl_fus"]]
            a.bar(np.arange(3) + (i - (len(labels) - 1) / 2) * width, vals,
                  width, label=lab)
        a.set_xticks(range(3))
        a.set_xticklabels(["band", "out-band", "total"])
        a.axhline(0.0, color="0.5", lw=0.8)
        a.grid(alpha=0.3, axis="y")
        a.set_title(f"{level} a={alpha}")
        if j == 0:
            a.set_ylabel("cl_fus contribution")
            a.legend(fontsize=8)
    fig.suptitle("B28 F2: fuselage-lift decomposition -- A (tilt, q>=0) / "
                 "C (tilt, fragment) / F (FLAT, fragment) / oracle (conf)")
    fig.tight_layout()
    fig.savefig(OUT / "f2_decomposition.png", dpi=150)
    print(f"wrote {OUT/'f2_decomposition.png'}")


def _verdict(rows, comp, oracle):
    dec = next(r for r in rows if r["level"] == "medium"
               and r["alpha"] == 3.06)
    c = comp[("C", "medium", 3.06)]
    f_out, c_out, o_out = dec["cl_fus_out"], c["cl_fus_out"], oracle["out"]
    vs_oracle = abs(f_out - o_out) / abs(o_out)
    vs_c = abs(f_out - c_out) / abs(c_out)
    r_factor = abs(c_out - f_out) / abs(c_out - o_out)
    print("\n=== B28 F2 pre-registered verdict (decisive leg medium a=3.06) "
          "===")
    print(f"  F_out={f_out:.5f}  C_out={c_out:.5f}  oracle_out={o_out:.5f}")
    print(f"  |F-oracle|/|oracle| = {vs_oracle:.2%} (F1 needs <= {TOL:.0%})")
    print(f"  |F-C|/|C|           = {vs_c:.2%} (F2 needs <= {TOL:.0%})")
    print(f"  position factor r   = {r_factor:.2f}")
    if vs_oracle <= TOL and vs_c > TOL:
        print("  BRANCH F1: the out-band gap is the sheet POSITION alone.")
    elif vs_c <= TOL and vs_oracle > TOL:
        print("  BRANCH F2: position EXONERATED -- topology/root-vortex "
              "strength (M2 wiring) is the next suspect.")
    else:
        print(f"  BRANCH F3: mixed/ambiguous (r = {r_factor:.2f}).")


if __name__ == "__main__":
    main()
