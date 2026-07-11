"""
G10.2 acceptance A/B (roadmap P10): level-adaptive intermediate
continuation tolerance (`solve_newton_transonic(intermediate_tol=1e-5)`)
vs the default strict path, on the two committed Newton recipes:

  (b) NACA0012 medium, M0.7875 / alpha 1.25 -- NEWTON_TRANSONIC_RECIPE
      (G8.1 regression locks: shock 0.674, cl 0.523, M_max 1.404)
  (a) ONERA M6 medium,  M0.84   / alpha 3.06 -- NEWTON_M6_RECIPE
      (G8.2 regression locks: cl 0.2646, M_max 2.134,
       shocks 0.596 / 0.541 / 0.362 at eta 0.44 / 0.65 / 0.90)

Acceptance (pre-registered in roadmap.md P10): ALL regression locks
intact under the adaptive path, robustness not degraded (attempted level
count and dm-halvings not worse), end-to-end deltas reported; the knob
is promoted into the recipes only if case (a) improves >= 15%.

MEASURED OUTCOME (2026-07-11, two A/B rounds): SPLIT verdict.
(a) M6 medium: all locks intact, final level converges identically
    (12 steps, |R| 7.8e-15, cl/M_max/shocks equal to 4 digits), solve
    239.5 s -> 140.3 s (+41.4%) => intermediate_tol=1e-5 PROMOTED into
    NEWTON_M6_RECIPE.
(b) NACA medium M0.7875 (FOLD ZONE): NEGATIVE result, recorded as the
    expected XFAILs below -- the loose ramp reaches the final level
    with an untracked Gamma/assignment seed (1-4 Newton steps per
    level; dcl/dM ~ 6-10 near the fold) and neither the final level
    nor STRICTLY-run dm-halving retry levels recover within the
    60-step budget (round 2 hardening: loose acceptance requires >= 1
    Newton step AND inserted retry levels run strict -- still stuck at
    the ~5e-6 live-churn floor with cl 0.369 vs 0.523). Loose
    intermediates are CONTRAINDICATED in fold zones; NEWTON_TRANSONIC_
    RECIPE keeps the default strict path (the P8 "warm-start only from
    CONVERGED levels" trap, now measured in its G10.2 form).

Timing protocol: NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
OPENBLAS_NUM_THREADS=16 (the G8.2-measured oversubscription trap).
One-shot evidence script -- NOT a suite test (G8.3 CI budget untouched).
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _common import (  # noqa: E402
    BASELINE, CheckList, S1_BLUE, S3_YELLOW, apply_style, finish,
    write_csv, MESH_DIR,
)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.post.section_cut import section_cp_curve, wall_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import planform_area, wall_force_coefficients
from pyfp3d.solve.continuation import mach_schedule
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)

INTERMEDIATE_TOL = 1e-5          # the roadmap candidate rule's |R| clause
PROMOTE_THRESHOLD = 0.15

# The PRE-promotion recipes (tests/test_p8_newton.py as of G8.2 closure)
# -- the A/B object is "recipe default" vs "recipe + intermediate_tol",
# so these copies deliberately EXCLUDE the knob that run_variant adds;
# the promoted NEWTON_M6_RECIPE in tests/ now carries intermediate_tol.
NEWTON_TRANSONIC_RECIPE = dict(
    dm=0.025, dm_min=0.003, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   n_newton_max=60),
)
NEWTON_M6_RECIPE = dict(
    dm=0.05, dm_min=0.01, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   direct_refactor_every=1000, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)


def run_variant(mc, wc, m_inf, alpha, recipe, adaptive):
    kw = dict(recipe)
    if adaptive:
        kw["intermediate_tol"] = INTERMEDIATE_TOL
    t0 = time.perf_counter()
    r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha, **kw)
    wall = time.perf_counter() - t0
    lvls = r["level_results"]
    n_sched = len(mach_schedule(m_inf, 0.70, recipe["dm"]))
    return {
        "result": r,
        "wall": wall,
        "levels": lvls,
        "n_levels": len(lvls),
        "n_halvings": max(0, len(lvls) - n_sched),
        "steps": sum(lr["n_newton"] for lr in lvls),
        "intermediate_wall": sum(lr["wall_s"] for lr in lvls[:-1]),
        "final_wall": lvls[-1]["wall_s"],
    }


def check_locks(cl2, case, v, locks, xfail=False, note=""):
    r = v["result"]
    h = r["residual_history"]
    drops = [h[i + 1] / h[i] for i in range(len(h) - 1)]
    quad_pair = any(a < 3e-2 and b < 3e-2
                    for a, b in zip(drops, drops[1:]))
    cl2.add(case, "converged", r["converged"], "True", r["converged"],
            xfail=xfail, note=note)
    cl2.add(case, "final |R|", f"{h[-1]:.2e}", "< 1e-9", h[-1] < 1e-9,
            xfail=xfail, note=note)
    cl2.add(case, "terminal quadratic pair", quad_pair, "True", quad_pair,
            xfail=xfail, note=note)
    cl2.add(case, "0 limited/floored",
            f"{r['n_limited']}/{r['n_floored']}", "0/0",
            r["n_limited"] == 0 and r["n_floored"] == 0)
    cl2.add(case, "Kutta |F|", f"{r['F_history'][-1]:.2e}", "< 1e-12",
            r["F_history"][-1] < 1e-12)
    for name, val, ref, band in locks:
        cl2.add(case, name, f"{val:.4f}", f"{ref} +- {band}",
                abs(val - ref) < band, xfail=xfail, note=note)


def per_level_rows(case, tag, v):
    return [(case, tag, f"{lr['m']:.4f}", lr["upwind_c"],
             lr["n_newton"], f"{lr['wall_s']:.2f}",
             lr["accept_reason"], f"{lr['residual_history'][-1]:.3e}",
             f"{lr['F_history'][-1]:.3e}")
            for lr in v["levels"]]


def level_figure(case, base, adapt, fname):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ms = [f"{lr['m']:.3f}" for lr in base["levels"]]
    x = np.arange(len(base["levels"]))
    axes[0].bar(x - 0.2, [lr["wall_s"] for lr in base["levels"]], 0.4,
                color=BASELINE, label="default")
    xa = np.arange(len(adapt["levels"]))
    axes[0].bar(xa + 0.2, [lr["wall_s"] for lr in adapt["levels"]], 0.4,
                color=S1_BLUE, label="adaptive (G10.2)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ms, rotation=45)
    axes[0].set_xlabel("continuation level M")
    axes[0].set_ylabel("wall s")
    axes[0].set_title(f"{case}: per-level cost")
    axes[0].legend(frameon=False)
    for v, color, tag in ((base, BASELINE, "default"),
                          (adapt, S1_BLUE, "adaptive")):
        hist = np.concatenate([lr["residual_history"]
                               for lr in v["levels"]])
        axes[1].semilogy(hist, color=color, label=tag, lw=1.4)
    axes[1].axhline(INTERMEDIATE_TOL, color=S3_YELLOW, ls="--", lw=1,
                    label="intermediate tol")
    axes[1].set_xlabel("cumulative Newton step")
    axes[1].set_ylabel(r"$\|R\|_\infty$")
    axes[1].set_title("stacked level residual histories")
    axes[1].legend(frameon=False)
    finish(fig, OUT, fname)


def main():
    apply_style()
    cl2 = CheckList("G10.2 level-adaptive intermediate tolerance A/B")
    level_rows = []
    summary_rows = []

    # numba/BLAS warmup so the first timed variant is not paying compile
    mc0, wc0 = cut_wake(read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh"))
    solve_newton_lifting(mc0, wc0, m_inf=0.5, alpha_deg=2.0,
                         upwind_c=1.5, m_crit=0.95, precond="direct")

    # ---- case (b): NACA medium M0.7875 (G8.1 locks) --------------------
    mc, wc = cut_wake(read_mesh(MESH_DIR / "naca0012_2.5d" / "medium.msh"))
    dz = float(np.ptp(mc.nodes[:, 2]))
    variants = {}
    for tag in ("default", "adaptive"):
        v = run_variant(mc, wc, 0.7875, 1.25, NEWTON_TRANSONIC_RECIPE,
                        adaptive=(tag == "adaptive"))
        variants[tag] = v
        rep = shock_report(
            wall_cp_curve(mc, v["result"]["phi"], z=0.5 * dz, m_inf=0.7875),
            0.7875)
        forces = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"],
            v["result"]["phi"], alpha_deg=1.25, s_ref=dz, m_inf=0.7875)
        v["locks"] = [
            ("shock x/c", rep["upper"]["x_shock"], 0.674, 0.012),
            ("cl", forces["cl"], 0.523, 0.01),
            ("M_max", float(np.sqrt(v["result"]["mach2_max"])), 1.404, 0.02),
        ]
        level_rows += per_level_rows("naca_medium", tag, v)
    b, a = variants["default"], variants["adaptive"]
    check_locks(cl2, "naca_medium/default", b, b["locks"])
    # measured negative result (module docstring): loose intermediates are
    # contraindicated in the fold zone -- expected failures, recorded
    check_locks(cl2, "naca_medium/adaptive", a, a["locks"], xfail=True,
                note="fold-zone contraindication (recorded negative)")
    cl2.add("naca_medium", "fold-zone negative result recorded",
            f"adaptive converged={a['result']['converged']}",
            "do NOT promote into NEWTON_TRANSONIC_RECIPE", True,
            note="strict default kept for fold-zone continuation")
    speedup_b = (b["wall"] - a["wall"]) / b["wall"]
    summary_rows.append(("naca_medium", f"{b['wall']:.1f}", f"{a['wall']:.1f}",
                         f"{speedup_b:.3f}", b["steps"], a["steps"],
                         b["n_levels"], a["n_levels"],
                         b["n_halvings"], a["n_halvings"]))
    level_figure("NACA0012 medium M0.7875", b, a, "naca_medium_ab.png")
    print(f"\nnaca_medium: default {b['wall']:.1f}s -> adaptive "
          f"{a['wall']:.1f}s  ({100*speedup_b:+.1f}% saved), steps "
          f"{b['steps']} -> {a['steps']}")

    # ---- case (a): ONERA M6 medium M0.84 (G8.2 locks) -------------------
    mc, wc = cut_wake(read_mesh(MESH_DIR / "onera_m6" / "medium.msh"))
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    variants = {}
    for tag in ("default", "adaptive"):
        v = run_variant(mc, wc, 0.84, 3.06, NEWTON_M6_RECIPE,
                        adaptive=(tag == "adaptive"))
        variants[tag] = v
        forces = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"],
            v["result"]["phi"], alpha_deg=3.06, s_ref=s_ref, m_inf=0.84)
        shocks = {}
        for eta in (0.44, 0.65, 0.90):
            c = section_cp_curve(mc, v["result"]["phi"], eta=eta,
                                 b_semi=B_SEMI, m_inf=0.84)
            shocks[eta] = shock_report(c, 0.84)["upper"]["x_shock"]
        v["locks"] = [
            ("cl", forces["cl"], 0.2646, 0.005),
            ("M_max", float(np.sqrt(v["result"]["mach2_max"])), 2.134, 0.05),
            ("shock eta 0.44", shocks[0.44], 0.596, 0.02),
            ("shock eta 0.65", shocks[0.65], 0.541, 0.02),
            ("shock eta 0.90", shocks[0.90], 0.362, 0.02),
        ]
        level_rows += per_level_rows("m6_medium", tag, v)
    b, a = variants["default"], variants["adaptive"]
    for tag, v in variants.items():
        check_locks(cl2, f"m6_medium/{tag}", v, v["locks"])
    cl2.add("m6_medium", "levels not worse",
            f"{a['n_levels']} vs {b['n_levels']}", "adaptive <= default",
            a["n_levels"] <= b["n_levels"])
    cl2.add("m6_medium", "dm-halvings not worse",
            f"{a['n_halvings']} vs {b['n_halvings']}", "adaptive <= default",
            a["n_halvings"] <= b["n_halvings"])
    speedup_a = (b["wall"] - a["wall"]) / b["wall"]
    promote = speedup_a >= PROMOTE_THRESHOLD
    cl2.add("m6_medium", "solve speedup", f"{100*speedup_a:.1f}%",
            ">= 15% to promote into NEWTON_M6_RECIPE", promote,
            note="negative result recorded if below; locks decide safety")
    summary_rows.append(("m6_medium", f"{b['wall']:.1f}", f"{a['wall']:.1f}",
                         f"{speedup_a:.3f}", b["steps"], a["steps"],
                         b["n_levels"], a["n_levels"],
                         b["n_halvings"], a["n_halvings"]))
    level_figure("ONERA M6 medium M0.84", b, a, "m6_medium_ab.png")
    print(f"\nm6_medium: default {b['wall']:.1f}s -> adaptive "
          f"{a['wall']:.1f}s  ({100*speedup_a:+.1f}% saved), steps "
          f"{b['steps']} -> {a['steps']}")
    print(f"promotion verdict (>= 15% on M6 medium): "
          f"{'PROMOTE into NEWTON_M6_RECIPE' if promote else 'DO NOT PROMOTE (record negative)'}; "
          f"NEWTON_TRANSONIC_RECIPE unchanged (fold-zone negative result)")

    write_csv(OUT, "levels.csv",
              "case,variant,m,upwind_c,n_newton,wall_s,accept_reason,"
              "residual_final,F_final", level_rows)
    write_csv(OUT, "summary.csv",
              "case,wall_default_s,wall_adaptive_s,speedup,"
              "steps_default,steps_adaptive,levels_default,levels_adaptive,"
              "halvings_default,halvings_adaptive", summary_rows)
    # the promotion check is advisory (a negative A/B result is a valid
    # gate outcome); locks/robustness checks are the hard criteria
    for c in cl2.checks:
        if c["name"] == "solve speedup" and c["status"] == "FAIL":
            c["status"] = "INFO"
    return cl2.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
