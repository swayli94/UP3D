"""
Track A / A1 -- figure builders for the solver bottleneck study.

Shared by run_a1.py (2.5-D headline) and run_a1_m6.py (gated 3-D leg). Uses
the demo house style (`cases/demo/_common`): matplotlib Agg, the fixed light
palette, `finish` at dpi 130. One stable colour per method so every figure
reads the same.

Both legs write into the SAME results/ directory, so every figure name is
namespaced by `prefix`: "a1" for the 2.5-D leg (the default -- these are the
names demo_report/track_a.md cites), "a1_m6" for the gated 3-D leg. Without it
the 3-D run silently overwrites the committed 2.5-D figures with 3-D content
and the report's figure references quietly start lying.
"""

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

import matplotlib.pyplot as plt                                    # noqa: E402

from cases.demo._common import (                                  # noqa: E402
    BASELINE,
    CRITICAL,
    INK,
    INK_2,
    MUTED,
    S1_BLUE,
    S2_AQUA,
    S3_YELLOW,
    S4_ROSE,
    finish,
)
from pyfp3d.solve.timing import PHASES                            # noqa: E402

from _bench import METHOD_LABEL, METHODS                          # noqa: E402

METHOD_COLOR = {
    "conf_picard": S1_BLUE,
    "conf_newton": S2_AQUA,
    "ls_picard": S3_YELLOW,
    "ls_newton": S4_ROSE,
}
# one colour per timing phase for the stacked bars
PHASE_COLOR = {
    "seed": BASELINE,
    "assembly": S1_BLUE,
    "precond": S2_AQUA,
    "linsolve": S3_YELLOW,
    "residual": S4_ROSE,
    "kutta": INK_2,
    "other": MUTED,
}


def _iters_wall(rec):
    """Per-iteration (residual, gamma, cum-wall) arrays for a run, flattened
    across ramp levels when present."""
    steps = _all_steps(rec)
    res = np.array([s["residual"] for s in steps], dtype=float)
    gam = np.array([s.get("gamma_mean", np.nan) for s in steps], dtype=float)
    wall = np.array([s.get("wall_cum_s", np.nan) for s in steps], dtype=float)
    return res, gam, wall


def _all_steps(rec):
    if rec.get("_levels_steps"):
        return rec["_levels_steps"]
    return rec.get("_steps", [])


# --------------------------------------------------------------------------

def fig_convergence(recs, out_dir, title, prefix="a1"):
    """2x2: residual vs iteration | vs wall; gamma vs iteration | vs wall.

    The wall-clock column is the point of the study: it is the only view in
    which a Newton step and a Picard sweep are commensurable."""
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    for m in METHODS:
        rec = recs.get(m)
        if rec is None:
            continue
        res, gam, wall = _iters_wall(rec)
        if len(res) == 0:
            continue
        it = np.arange(1, len(res) + 1)
        res_pos = np.where(res > 0, res, np.nan)
        c = METHOD_COLOR[m]
        lbl = METHOD_LABEL[m]
        ax[0, 0].semilogy(it, res_pos, color=c, label=lbl, marker=".", ms=4)
        ax[0, 1].semilogy(wall, res_pos, color=c, label=lbl, marker=".", ms=4)
        ax[1, 0].plot(it, gam, color=c, label=lbl, marker=".", ms=4)
        ax[1, 1].plot(wall, gam, color=c, label=lbl, marker=".", ms=4)
    ax[0, 0].set(xlabel="outer / Newton iteration", ylabel="residual |R|")
    ax[0, 1].set(xlabel="wall-clock (s)", ylabel="residual |R|")
    ax[1, 0].set(xlabel="outer / Newton iteration", ylabel="circulation Γ")
    ax[1, 1].set(xlabel="wall-clock (s)", ylabel="circulation Γ")
    ax[0, 0].set_title("convergence per iteration")
    ax[0, 1].set_title("convergence per wall-second — the honest comparison")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle(title, fontweight="bold")
    finish(fig, out_dir, f"{prefix}_convergence")


def fig_time_breakdown(recs, out_dir, title, regime_label="", prefix="a1"):
    """Stacked per-phase bars per method, annotated with wall_s."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.set_axisbelow(True)              # gridlines under the bars, not striping
    labels = [METHOD_LABEL[m] for m in METHODS if m in recs]
    xs = np.arange(len(labels))
    bottoms = np.zeros(len(labels))
    used = [m for m in METHODS if m in recs]
    for p in list(PHASES) + ["other"]:
        vals = np.array([recs[m]["t_" + p] for m in used])
        ax.bar(xs, vals, bottom=bottoms, color=PHASE_COLOR[p], label=p,
               edgecolor="white", linewidth=0.5, zorder=3)
        bottoms += vals
    for x, m in zip(xs, used):
        ax.text(x, bottoms[list(xs).index(x)], f"{recs[m]['wall_s']:.1f}s",
                ha="center", va="bottom", fontsize=9, color=INK,
                fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=12, ha="right")
    ax.set_ylabel("wall-clock (s)")
    ax.legend(title="phase", fontsize=8, ncol=2)
    ax.set_title(f"{title}{(' — ' + regime_label) if regime_label else ''}",
                 fontweight="bold")
    finish(fig, out_dir, f"{prefix}_time_breakdown")


def fig_ramp(recs, out_dir, title, prefix="a1"):
    """Mach-ramp anatomy: outer iterations and wall seconds per level."""
    ramps = {m: recs[m] for m in METHODS
             if m in recs and recs[m].get("_level_rows")}
    if not ramps:
        return False
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for m, rec in ramps.items():
        lv = rec["_level_rows"]
        mach = [r["m"] for r in lv]
        niter = [r["n_iter"] for r in lv]
        wall = [r["wall_s"] for r in lv]
        c = METHOD_COLOR[m]
        ax[0].plot(mach, niter, color=c, marker="o", label=METHOD_LABEL[m])
        ax[1].plot(mach, wall, color=c, marker="o", label=METHOD_LABEL[m])
    ax[0].set(xlabel="Mach level", ylabel="iterations at level")
    ax[1].set(xlabel="Mach level", ylabel="wall-clock at level (s)")
    ax[0].set_title("iterations per Mach level")
    ax[1].set_title("cost per Mach level — where the ramp spends time")
    ax[0].legend(fontsize=8)
    fig.suptitle(title, fontweight="bold")
    finish(fig, out_dir, f"{prefix}_ramp")
    return True


def fig_linear_solver(recs, out_dir, title, prefix="a1"):
    """Linear-algebra anatomy: iters/solve, solves/outer, refactors, stalls."""
    used = [m for m in METHODS if m in recs]
    xs = np.arange(len(used))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for a in ax:
        a.set_axisbelow(True)
    ips = [recs[m]["n_lin_iters"] / max(recs[m]["n_lin_solves"], 1)
           for m in used]
    spo = [recs[m]["n_lin_solves"] / max(recs[m]["n_outer_total"], 1)
           for m in used]
    ax[0].bar(xs, ips, color=[METHOD_COLOR[m] for m in used], zorder=3)
    ax[1].bar(xs, spo, color=[METHOD_COLOR[m] for m in used], zorder=3)
    ax[1].axhline(1.0, color=CRITICAL, ls="--", lw=1,
                  label="1 solve / outer (Newton ideal)")
    for a in ax:
        a.set_xticks(xs)
        a.set_xticklabels([METHOD_LABEL[m] for m in used], rotation=12,
                          ha="right")
    ax[0].set_ylabel("linear iterations per solve")
    ax[1].set_ylabel("linear solves per outer")
    ax[0].set_title("Krylov/CG work per solve")
    ax[1].set_title("solves per outer — the Picard inner-Kutta amplifier")
    ax[1].legend(fontsize=8)
    fig.suptitle(title, fontweight="bold")
    finish(fig, out_dir, f"{prefix}_linear_solver")


def fig_cp(recs, out_dir, title, curves, ref=None, prefix="a1"):
    """Cp(x/c) overlay of the four methods (+ optional reference)."""
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for m in METHODS:
        cv = curves.get(m)
        if cv is None:
            continue
        c = METHOD_COLOR[m]
        ax.plot(cv["x_upper"], cv["cp_upper"], color=c, label=METHOD_LABEL[m])
        ax.plot(cv["x_lower"], cv["cp_lower"], color=c, ls="--")
    if ref is not None:
        ax.scatter(ref[0], ref[1], s=16, color=INK, zorder=5,
                   label="reference")
    ax.invert_yaxis()
    ax.set(xlabel="x/c", ylabel="Cp")
    ax.legend(fontsize=8)
    ax.set_title(title, fontweight="bold")
    finish(fig, out_dir, f"{prefix}_cp")


def fig_spanwise(recs, out_dir, title, prefix="a1"):
    """Gamma(z) / sectional loading vs z/b, four methods (3-D leg)."""
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    plotted = False
    for m in METHODS:
        rec = recs.get(m)
        if rec is None or "_span_z" not in rec:
            continue
        ax.plot(rec["_span_z"], rec["_span_gamma"], color=METHOD_COLOR[m],
                marker=".", ms=4, label=METHOD_LABEL[m])
        plotted = True
    ax.set(xlabel="z / b (span)", ylabel="circulation Γ(z)")
    ax.legend(fontsize=8)
    ax.set_title(title, fontweight="bold")
    finish(fig, out_dir, f"{prefix}_spanwise")
    return plotted


def fig_dashboard(rows, out_dir, title, prefix="a1"):
    """Speedup matrix: wall_s per (method, regime) with the dominant phase."""
    regimes = sorted({r["regime"] for r in rows},
                     key=lambda x: 0 if x == "subsonic" else 1)
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.set_axisbelow(True)
    used = [m for m in METHODS]
    xs = np.arange(len(used))
    w = 0.8 / max(len(regimes), 1)
    for j, reg in enumerate(regimes):
        vals, doms = [], []
        for m in used:
            rr = [r for r in rows if r["method"] == m and r["regime"] == reg
                  and r["wake_model_headline"]]
            vals.append(rr[0]["wall_s"] if rr else 0.0)
            doms.append(rr[0]["dominant_phase"] if rr else "")
        bars = ax.bar(xs + j * w, vals, w, label=reg,
                      color=[METHOD_COLOR[m] for m in used],
                      alpha=0.6 if j else 1.0, edgecolor="white")
        for b, d in zip(bars, doms):
            if d:
                ax.text(b.get_x() + b.get_width() / 2, b.get_height(), d,
                        ha="center", va="bottom", fontsize=7, rotation=90,
                        color=INK_2)
    ax.set_xticks(xs + w * (len(regimes) - 1) / 2)
    ax.set_xticklabels([METHOD_LABEL[m] for m in used], rotation=12,
                       ha="right")
    ax.set_ylabel("wall-clock (s)")
    ax.legend(title="regime", fontsize=8)
    ax.set_title(title, fontweight="bold")
    finish(fig, out_dir, f"{prefix}_dashboard")
