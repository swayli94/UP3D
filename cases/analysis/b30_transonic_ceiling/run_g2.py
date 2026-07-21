"""GB30.2 -- (b)-class ceiling attribution census (B30 PRE_REGISTRATION.md).

Two legs, each a warm-started single-level STRICT solve at the production
recipe's dying level, with the pre-registered continuation clause (if the
single level converges, climb by dm=0.01 to the first FAILURE level,
<= 0.84 -- path-dependence vs the cascade is itself an independent finding):

  LS    m0=0.7875 from ls_flat_medium_084.npz (phi_ext @ 0.775)
  CONF  m0=0.80   from conf_medium_079.npz     (phi @ 0.79)

Census at the first failing level per leg: limited/floored element masks
(the solvers' own accounting, replicated and VALIDATED against the
recorded counts), region histogram (tip_box / corridor / near_fus /
field), count-level oscillation (step_records / clamp_history), peak-M
attribution, corridor corrM, cl_p / gamma vs the committed anchors.

Verdict tree (pre-registered):
  B30-SAME   both paths' clamp sets >= 2/3 in tip_box and total < 50
  B30-SPLIT  CONF clamps junction/shock-located while LS tip-located
  B30-AMBIG  in between -> full histogram to the user

Run:  python cases/analysis/b30_transonic_ceiling/run_g2.py [LS CONF]
Artifacts: results/g2_census.csv, results/g2_regions.csv,
           results/g2_oscillation.png
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wb30 import (ALPHA, DEMO_OUT, LS_MESH_DIR, OUT, UPWIND_DEFAULT,
                  build_conf, build_ls_flat, census_rows, clamp_masks_conf,
                  clamp_masks_ls, cl_p_conf, cl_p_ls, corridor_mmax,
                  load_mesh, mach2_conf, mach2_ls, oscillation_stats,
                  region_codes, run_conf_level, run_ls_level)
from pyfp3d.solve.newton import NewtonWorkspace
from pyfp3d.constraints.wake import kutta_targets

M0 = {"LS": 0.7875, "CONF": 0.80}
SEED = {"LS": ("ls_flat_medium_084.npz", "phi_ext"),
        "CONF": ("conf_medium_079.npz", "phi")}
# committed anchors for the trend check (checks.csv / B26-B29)
CL_ANCHOR = {"LS": 0.2475, "CONF": 0.2579}     # LS @0.775-dying, CONF @0.79
TIP_FRAC_SAME = 2.0 / 3.0
CLAMPS_LOCALIZED = 50


def leg_ls(mesh):
    _wls, _cm, mvop = build_ls_flat(mesh)
    d = np.load(DEMO_OUT / SEED["LS"][0])
    phi0 = d[SEED["LS"][1]]
    m, phi = M0["LS"], phi0
    chain = []
    while True:
        rec = run_ls_level(mesh, mvop, m, phi, tag="g2")
        chain.append(rec)
        if not rec["converged"] or m >= 0.84 - 1e-12:
            break
        phi = rec["phi"]
        m = round(m + 0.01, 4)
    top = chain[-1]
    m_inf = top["m_inf"]
    phi_ext = top["phi"]
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    code, _df = region_codes(cents)
    m2 = mach2_ls(mvop, phi_ext, m_inf)
    lim, flr = clamp_masks_ls(mesh, mvop, phi_ext, m_inf,
                              upwind_c=top["upwind_c"],
                              m_crit=UPWIND_DEFAULT["m_crit"],
                              m_cap=UPWIND_DEFAULT["m_cap"],
                              rho_floor=UPWIND_DEFAULT["rho_floor"])
    # pre-reg self-check: the replicated counts must equal the solver's
    check_ok = (int(lim.sum()) == top["n_limited"]
                and int(flr.sum()) == top["n_floored"])
    rows = (census_rows("LS", "limited", lim, m2, cents, code)
            + census_rows("LS", "floored", flr, m2, cents, code))
    corr_m, corr_xyz = corridor_mmax(m2, cents)
    i = int(np.argmax(m2))
    summary = dict(leg="LS", m_census=m_inf, converged=top["converged"],
                   accept_reason=top["accept_reason"],
                   res=top["residual_norm"], n_newton=top["n_newton"],
                   n_limited=top["n_limited"], n_floored=top["n_floored"],
                   mask_lim=int(lim.sum()), mask_flr=int(flr.sum()),
                   self_check=check_ok,
                   n_climb=len(chain),
                   mmax=float(np.sqrt(m2.max())), pk_x=float(cents[i, 0]),
                   pk_y=float(cents[i, 1]), pk_z=float(cents[i, 2]),
                   corr_mmax=corr_m, corr_x=corr_xyz[0], corr_z=corr_xyz[2],
                   cl_p=cl_p_ls(mesh, mvop, phi_ext, m_inf),
                   gamma=float(np.mean(top["gamma"])),
                   wall_s=sum(c["wall_s"] for c in chain))
    summary.update(oscillation_stats(top))
    return summary, rows, top


def leg_conf():
    mc, wc = build_conf("medium")
    d = np.load(DEMO_OUT / SEED["CONF"][0])
    phi0 = d[SEED["CONF"][1]]
    # T6 (measured): a zeros-Gamma seed kills the wake jump + vortex far
    # field and diverges -- reconstruct Gamma from the cached cut phi.
    gamma0 = kutta_targets(phi0, wc)
    m, phi, gamma = M0["CONF"], phi0, gamma0
    chain = []
    while True:
        rec = run_conf_level(mc, wc, m, phi, gamma_init=gamma, tag="g2")
        chain.append(rec)
        if not rec["converged"] or m >= 0.84 - 1e-12:
            break
        phi, gamma = rec["phi"], rec["gamma"]
        m = round(m + 0.01, 4)
    top = chain[-1]
    m_inf = top["m_inf"]
    ws = NewtonWorkspace(mc, wc, ALPHA, 1.0, 1.4, (0.25, 0.0), True,
                         kutta_estimator="pressure")
    ws.set_mach(m_inf)
    cents = mc.nodes[mc.elements].mean(axis=1)
    code, _df = region_codes(cents)
    m2 = mach2_conf(ws, top["phi"], top["gamma"], m_inf)
    lim, flr = clamp_masks_conf(ws, top["phi"], top["gamma"], m_inf,
                                upwind_c=top["upwind_c"],
                                m_crit=UPWIND_DEFAULT["m_crit"],
                                m_cap=UPWIND_DEFAULT["m_cap"],
                                rho_floor=UPWIND_DEFAULT["rho_floor"])
    # pre-reg self-check: (i) the replication must read 0/0 on the last
    # CONVERGED chain state (the conforming converged gate refuses clamped
    # states -- probe-verified 2026-07-21 on the 0.82 cache) and (ii) it
    # matches the dying level within the final-step bookkeeping drift:
    # clamp_history is one update stale on this path (newton.py:672 vs
    # :992), while the LS solver re-evaluates at the final state
    # (newton_ls.py:928-932), so the LS leg checks EXACT.
    conv = [c for c in chain[:-1] if c["converged"]]
    clean_ok = True
    if conv:
        cc = conv[-1]
        ws.set_mach(cc["m_inf"])
        l0, f0 = clamp_masks_conf(ws, cc["phi"], cc["gamma"], cc["m_inf"],
                                  upwind_c=cc["upwind_c"],
                                  m_crit=UPWIND_DEFAULT["m_crit"],
                                  m_cap=UPWIND_DEFAULT["m_cap"],
                                  rho_floor=UPWIND_DEFAULT["rho_floor"])
        clean_ok = (int(l0.sum()) == 0 and int(f0.sum()) == 0)
    drift_ok = (abs(int(lim.sum()) - top["n_limited"]) <= 2
                and abs(int(flr.sum()) - top["n_floored"]) <= 2)
    check_ok = clean_ok and drift_ok
    rows = (census_rows("CONF", "limited", lim, m2, cents, code)
            + census_rows("CONF", "floored", flr, m2, cents, code))
    corr_m, corr_xyz = corridor_mmax(m2, cents)
    i = int(np.argmax(m2))
    summary = dict(leg="CONF", m_census=m_inf, converged=top["converged"],
                   accept_reason=top["accept_reason"],
                   res=top["residual_norm"], n_newton=top["n_newton"],
                   n_limited=top["n_limited"], n_floored=top["n_floored"],
                   mask_lim=int(lim.sum()), mask_flr=int(flr.sum()),
                   self_check=check_ok,
                   n_climb=len(chain),
                   mmax=float(np.sqrt(m2.max())), pk_x=float(cents[i, 0]),
                   pk_y=float(cents[i, 1]), pk_z=float(cents[i, 2]),
                   corr_mmax=corr_m, corr_x=corr_xyz[0], corr_z=corr_xyz[2],
                   cl_p=cl_p_conf(mc, top["phi"], m_inf),
                   gamma=float(np.mean(top["gamma"])),
                   wall_s=sum(c["wall_s"] for c in chain))
    summary.update(oscillation_stats(top))
    return summary, rows, top


def tip_fraction(rows, leg):
    tot = sum(r["count"] for r in rows
              if r["leg"] == leg and r["clamp"] == "limited") \
        + sum(r["count"] for r in rows
              if r["leg"] == leg and r["clamp"] == "floored")
    tip = sum(r["count"] for r in rows
              if r["leg"] == leg and r["region"] == "tip_box")
    return (tip / tot if tot else 0.0), tot


def main():
    only = sys.argv[1:] or ["LS", "CONF"]
    summaries, region_rows, tops = [], [], {}
    if "LS" in only:
        mesh = load_mesh(LS_MESH_DIR / "medium.msh")
        print(f"=== GB30.2 LS leg: {len(mesh.elements)} tets, "
              f"m0={M0['LS']} ===", flush=True)
        s, rr, top = leg_ls(mesh)
        summaries.append(s)
        region_rows += rr
        tops["LS"] = top
        print(f"  [LS] m_census={s['m_census']} conv={s['converged']} "
              f"clamps={s['n_limited']}+{s['n_floored']} "
              f"self_check={s['self_check']} mmax={s['mmax']:.2f} "
              f"@z={s['pk_z']:.3f} corrM={s['corr_mmax']:.2f} "
              f"cl_p={s['cl_p']:.4f}", flush=True)
    if "CONF" in only:
        print(f"=== GB30.2 CONF leg: m0={M0['CONF']} ===", flush=True)
        s, rr, top = leg_conf()
        summaries.append(s)
        region_rows += rr
        tops["CONF"] = top
        print(f"  [CONF] m_census={s['m_census']} conv={s['converged']} "
              f"clamps={s['n_limited']}+{s['n_floored']} "
              f"self_check={s['self_check']} mmax={s['mmax']:.2f} "
              f"@z={s['pk_z']:.3f} corrM={s['corr_mmax']:.2f} "
              f"cl_p={s['cl_p']:.4f}", flush=True)

    _write_census(summaries)
    _write_regions(region_rows)
    _write_fig(tops)
    _verdict_hint(summaries, region_rows)


def _write_census(rows):
    keys = ["leg", "m_census", "converged", "accept_reason", "res",
            "n_newton", "n_limited", "n_floored", "mask_lim", "mask_flr",
            "self_check", "n_climb", "mmax", "pk_x", "pk_y", "pk_z",
            "corr_mmax", "corr_x", "corr_z", "cl_p", "gamma",
            "step_nlim_min", "step_nlim_max", "step_nlim_last",
            "step_nlim_range", "step_nflr_min", "step_nflr_max",
            "step_nflr_last", "step_nflr_range", "wall_s"]
    with open(OUT / "g2_census.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g2_census.csv'}")


def _write_regions(rows):
    keys = ["leg", "clamp", "region", "count", "peak_mach",
            "pk_x", "pk_y", "pk_z"]
    with open(OUT / "g2_regions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / 'g2_regions.csv'}")


def _write_fig(tops):
    fig, axes = plt.subplots(1, len(tops) or 1, figsize=(6 * max(len(tops), 1), 4.2),
                             squeeze=False)
    for ax, (leg, top) in zip(axes[0], tops.items()):
        its = np.arange(1, len(top["step_res"]) + 1)
        ax2 = ax.twinx()
        ax.semilogy(its, top["step_res"], "k-", lw=1.0, label="|R|")
        ax2.plot(its, top["step_nlim"], "r-o", ms=3, label="n_limited")
        ax2.plot(its, top["step_nflr"], "b-s", ms=3, label="n_floored")
        ax.set_xlabel("Newton iteration")
        ax.set_ylabel("|R|_inf")
        ax2.set_ylabel("clamped cells")
        ax.set_title(f"{leg} dying level m={top['m_inf']} "
                     f"(conv={top['converged']})")
        ax.grid(alpha=0.3)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    fig.suptitle("GB30.2: clamp-count oscillation at the dying levels")
    fig.tight_layout()
    fig.savefig(OUT / "g2_oscillation.png", dpi=120)
    plt.close(fig)
    print(f"wrote {OUT / 'g2_oscillation.png'}")


def _verdict_hint(summaries, region_rows):
    print("\n--- verdict hint (pre-reg GB30.2) ---")
    fracs = {}
    for s in summaries:
        f, tot = tip_fraction(region_rows, s["leg"])
        fracs[s["leg"]] = (f, tot)
        print(f"  {s['leg']}: clamps total={tot}, tip_box fraction={f:.2f}, "
              f"self_check={s['self_check']}, m_census={s['m_census']}, "
              f"pk z={s['pk_z']:.3f}, corrM={s['corr_mmax']:.2f}")
    if len(fracs) == 2:
        same = all(f >= TIP_FRAC_SAME and tot < CLAMPS_LOCALIZED
                   for f, tot in fracs.values())
        split = fracs["LS"][0] >= TIP_FRAC_SAME and fracs["CONF"][0] < 1.0 / 3.0
        if same:
            print("  -> B30-SAME candidate: both clamp sets tip-localized; "
                  "GB30.3 runs BOTH legs")
        elif split:
            print("  -> B30-SPLIT candidate: LS tip-localized, CONF not; "
                  "GB30.3 runs the LS leg only")
        else:
            print("  -> B30-AMBIG: full histogram in g2_regions.csv; "
                  "user arbitration")
    bad = [s["leg"] for s in summaries if not s["self_check"]]
    if bad:
        print(f"  !! census self-check FAILED on {bad} -- fix the mask "
              f"replication before reading the histogram")


if __name__ == "__main__":
    main()
