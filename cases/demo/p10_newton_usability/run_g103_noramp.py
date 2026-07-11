"""
G10.3 no-ramp direct-solve feasibility study (roadmap P10, user-directed
2026-07-11, scheduled BEFORE P9): can the coupled Newton converge AT the
target M_inf from a cold start -- no Mach continuation at all?

Post-G10.2 motivation: on M6 medium the (already loosened) intermediate
levels still cost ~89 s of the 140.3 s adaptive solve; the remaining
value of the ramp is (i) FP branch selection near the non-uniqueness
fold (design.md Sec 12 risk 2 -- why the ramp exists) and (ii) the
Gamma/active-set seed. This study measures both.

Protocol (pre-registered in roadmap.md G10.3): single-level
solve_newton_lifting at target M_inf with the current machinery, under
two seedings -- s1 = the standard short Picard seed (n_picard_seed=5),
s2 = a deeper Picard seed (n_picard_seed=40, the P5-era warm start) --
on four cases:

    m6_coarse    M0.84 /alpha3.06  (far from fold)   locks: ramp solution
    m6_medium    M0.84 /alpha3.06  (far from fold)   locks: G8.2
    naca_coarse  M0.80 /alpha1.25  (fold zone)       locks: G8.1
    naca_medium  M0.7875/alpha1.25 (fold zone)       locks: G8.1

Outcome classes per case (pre-registered):
    (A) SAME-solution: converged, terminal quadratic, 0 lim/flr at the
        end AND all ramp-solution locks met.  For PROMOTION the path
        must also be clamp-free (clamp_history all zero -- a transient
        clamp means the no-ramp transient left the isentropic model's
        validity on the way).
    (B) WRONG-BRANCH: converged cleanly but outside the locks -- hard
        evidence the ramp does branch selection; cl/shock recorded.
    (C) NO convergence (stall / divergence / clamped state).

Promotion rule: no-ramp goes into NEWTON_M6_RECIPE only if BOTH M6
cases land in (A) with end-to-end >= 20% below the current post-G10.2
recipe (ramp baselines re-measured in this same session for a fair
clock). The Mach ramp stays the shipped default for the 2.5D fold-zone
recipe REGARDLESS (G10.2's measured contraindication), and stays
available everywhere as the fallback.

Timing protocol: NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
OPENBLAS_NUM_THREADS=16. One-shot evidence script, not a suite test.
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
    BASELINE, CRITICAL, CheckList, MESH_DIR, S1_BLUE, S2_AQUA, S3_YELLOW,
    apply_style, finish, write_csv,
)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.post.section_cut import section_cp_curve, wall_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import planform_area, wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)

PROMOTE_GAIN = 0.20
SEEDINGS = {"s1_picard5": 5, "s2_picard40": 40}

# newton_kw shared by ramp + no-ramp (the shipped recipes' machinery)
NACA_KW = dict(freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
               n_newton_max=60)
M6_KW = dict(freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
             direct_refactor_every=1000, n_newton_max=60,
             farfield_spanwise_gamma=True)

CASES = {
    # case: (geom, level, m_inf, alpha, newton_kw, ramp_recipe, locks)
    # locks = [(name, ref, band)] around the RAMP solution (G8.1/G8.2 +
    # the capability-demo M6-coarse record)
    "m6_coarse": ("onera_m6", "coarse", 0.84, 3.06, M6_KW,
                  dict(dm=0.05, dm_min=0.01, freeze_tol=1e-6,
                       intermediate_tol=1e-5, newton_kw=M6_KW),
                  [("cl", 0.2560, 0.005), ("m_max", 1.397, 0.05),
                   ("shock_044", 0.600, 0.02), ("shock_065", 0.573, 0.02),
                   ("shock_090", 0.429, 0.02)]),
    "m6_medium": ("onera_m6", "medium", 0.84, 3.06, M6_KW,
                  dict(dm=0.05, dm_min=0.01, freeze_tol=1e-6,
                       intermediate_tol=1e-5, newton_kw=M6_KW),
                  [("cl", 0.2646, 0.005), ("m_max", 2.129, 0.05),
                   ("shock_044", 0.596, 0.02), ("shock_065", 0.541, 0.02),
                   ("shock_090", 0.362, 0.02)]),
    "naca_coarse": ("naca0012_2.5d", "coarse", 0.80, 1.25, NACA_KW, None,
                    [("cl", 0.459, 0.01), ("m_max", 1.408, 0.02),
                     ("shock", 0.658, 0.012)]),
    "naca_medium": ("naca0012_2.5d", "medium", 0.7875, 1.25, NACA_KW, None,
                    [("cl", 0.523, 0.01), ("m_max", 1.404, 0.02),
                     ("shock", 0.674, 0.012)]),
}


def measure(case, mc, wc, m_inf, alpha, r):
    """cl / M_max / shocks for the lock comparison."""
    out = {"m_max": float(np.sqrt(r["mach2_max"]))}
    if case.startswith("m6"):
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        forces = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
            alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
        out["cl"] = float(forces["cl"])
        for eta, key in ((0.44, "shock_044"), (0.65, "shock_065"),
                         (0.90, "shock_090")):
            c = section_cp_curve(mc, r["phi"], eta=eta, b_semi=B_SEMI,
                                 m_inf=m_inf)
            out[key] = shock_report(c, m_inf)["upper"]["x_shock"]
    else:
        dz = float(np.ptp(mc.nodes[:, 2]))
        forces = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
            alpha_deg=alpha, s_ref=dz, m_inf=m_inf)
        out["cl"] = float(forces["cl"])
        rep = shock_report(
            wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf), m_inf)
        out["shock"] = (rep["upper"]["x_shock"]
                        if rep["upper"]["has_shock"] else float("nan"))
    return out


def quad_pair(hist):
    drops = [hist[i + 1] / hist[i] for i in range(len(hist) - 1)]
    return any(a < 3e-2 and b < 3e-2 for a, b in zip(drops, drops[1:]))


def classify(r, vals, locks):
    """Pre-registered class (A)/(B)/(C) + the clamp-free-path flag."""
    clamps = r.get("clamp_history", [])
    clamp_free = all(l == 0 and f == 0 for l, f in clamps)
    locks_ok = all(np.isfinite(vals[k]) and abs(vals[k] - ref) < band
                   for k, ref, band in locks)
    converged_clean = (r["converged"] and r["n_limited"] == 0
                       and r["n_floored"] == 0
                       and quad_pair(r["residual_history"]))
    if converged_clean and locks_ok:
        return "A", clamp_free
    if converged_clean:
        return "B", clamp_free
    return "C", clamp_free


def main():
    apply_style()
    cl = CheckList("G10.3 no-ramp direct-solve feasibility")
    rows = []
    hist_data = {}
    ramp_walls = {}
    noramp_best = {}

    # numba warmup (cache-warm subsonic solve; not timed against anything)
    mc0, wc0 = cut_wake(read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh"))
    solve_newton_lifting(mc0, wc0, m_inf=0.5, alpha_deg=2.0,
                         precond="direct")

    for case, (geom, level, m_inf, alpha, nkw, ramp_recipe,
               locks) in CASES.items():
        mc, wc = cut_wake(read_mesh(MESH_DIR / geom / f"{level}.msh"))

        # fresh ramp baseline (fair same-session clock; M6 cases only --
        # the fold-zone NACA baselines are the committed G8.1/A-B numbers
        # and no promotion decision hangs on their clock)
        if ramp_recipe is not None:
            t0 = time.perf_counter()
            rb = solve_newton_transonic(mc, wc, m_inf=m_inf,
                                        alpha_deg=alpha, **ramp_recipe)
            ramp_walls[case] = time.perf_counter() - t0
            assert rb["converged"], f"ramp baseline diverged on {case}!"
            print(f"[{case}] ramp baseline: {ramp_walls[case]:.1f}s")

        for tag, n_seed in SEEDINGS.items():
            t0 = time.perf_counter()
            r = solve_newton_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                     n_picard_seed=n_seed, **nkw)
            wall = time.perf_counter() - t0
            r.pop("workspace", None)
            vals = measure(case, mc, wc, m_inf, alpha, r)
            klass, clamp_free = classify(r, vals, locks)
            max_clamp = (max((l + f) for l, f in r["clamp_history"])
                         if r["clamp_history"] else 0)
            hist_data[(case, tag)] = (r["residual_history"], klass)
            if klass == "A" and (case not in noramp_best
                                 or wall < noramp_best[case][0]):
                noramp_best[case] = (wall, tag, clamp_free)
            rows.append((case, tag, klass, clamp_free,
                         r["converged"], r["n_newton"],
                         f"{r['residual_history'][-1]:.3e}",
                         f"{r['F_history'][-1]:.1e}",
                         r["n_limited"], r["n_floored"], max_clamp,
                         r["n_freeze_reverts"],
                         f"{vals.get('cl', float('nan')):.4f}",
                         f"{vals.get('m_max', float('nan')):.3f}",
                         ";".join(f"{vals[k]:.3f}" for k in vals
                                  if k.startswith("shock")),
                         f"{wall:.1f}"))
            print(f"[{case}] {tag}: class {klass} "
                  f"(clamp_free={clamp_free}) wall {wall:.1f}s "
                  f"steps {r['n_newton']} cl {vals.get('cl'):.4f}")
            cl.add("G10.3", f"{case}/{tag} outcome class",
                   f"{klass} (clamp_free={clamp_free}, {wall:.0f}s)",
                   "pre-registered A/B/C -- any class is a valid result",
                   True, note="evidence gate: the class IS the result")

    # ---- promotion decision (pre-registered) ----------------------------
    both_a = all(c in noramp_best and noramp_best[c][2]
                 for c in ("m6_coarse", "m6_medium"))
    gains = {}
    if both_a:
        for c in ("m6_coarse", "m6_medium"):
            gains[c] = (ramp_walls[c] - noramp_best[c][0]) / ramp_walls[c]
        promote = all(g >= PROMOTE_GAIN for g in gains.values())
    else:
        promote = False
    verdict = ("PROMOTE no-ramp into NEWTON_M6_RECIPE" if promote else
               "KEEP the Mach ramp (record the evidence)")
    cl.add("G10.3", "promotion rule",
           f"both M6 class-A+clamp-free={both_a}, gains=" +
           (", ".join(f"{c}:{100*g:.0f}%" for c, g in gains.items())
            if gains else "n/a"),
           ">= 20% on BOTH M6 cases, class A, clamp-free path", True,
           note=verdict)
    print(f"\nG10.3 verdict: {verdict}")

    write_csv(OUT, "g103_noramp.csv",
              "case,seeding,class,clamp_free_path,converged,n_newton,"
              "residual_final,F_final,n_limited,n_floored,max_clamped,"
              "freeze_reverts,cl,m_max,shocks,wall_s", rows)
    write_csv(OUT, "g103_summary.csv",
              "case,ramp_wall_s,noramp_best_wall_s,noramp_best_seeding,"
              "gain",
              [(c, f"{ramp_walls.get(c, float('nan')):.1f}",
                f"{noramp_best[c][0]:.1f}" if c in noramp_best else "n/a",
                noramp_best[c][1] if c in noramp_best else "n/a",
                f"{gains[c]:.3f}" if c in gains else "n/a")
               for c in CASES])

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    colors = {"A": S2_AQUA, "B": S3_YELLOW, "C": CRITICAL}
    for (case, tag), (h, klass) in hist_data.items():
        style = "-" if tag.startswith("s1") else "--"
        axes[0 if case.startswith("m6") else 1].semilogy(
            h, style, lw=1.3, color=colors[klass],
            label=f"{case}/{tag} [{klass}]")
    for ax, title in zip(axes, ("ONERA M6 (far from fold)",
                                "NACA 2.5D (fold zone)")):
        ax.set_xlabel("Newton step")
        ax.set_ylabel(r"$\|R\|_\infty$")
        ax.set_title(f"no-ramp cold start: {title}")
        ax.legend(frameon=False, fontsize=7)
    finish(fig, OUT, "g103_noramp_convergence.png")
    # CheckList.report always writes checks.csv; that name belongs to the
    # committed G10.2 A/B evidence in this shared results/ dir, so stash
    # and restore it around the report and keep this gate's table under
    # its own name (first run learned this the hard way -- the rename
    # alone clobbered the A/B file, restored in commit 25b1bf4)
    ab_checks = (OUT / "checks.csv").read_bytes() \
        if (OUT / "checks.csv").exists() else None
    code = cl.report(OUT)
    (OUT / "checks.csv").rename(OUT / "g103_checks.csv")
    if ab_checks is not None:
        (OUT / "checks.csv").write_bytes(ab_checks)
    return code


if __name__ == "__main__":
    sys.exit(main())
