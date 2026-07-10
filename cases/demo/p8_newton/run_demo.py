"""
P8 demo -- fully-coupled (phi_red, Gamma) Newton (gates G8.1 + the N6
G8.2 ONERA M6 performance run; G8.3 is the suite-runtime CI budget,
recorded in roadmap.md).

The exact Jacobian (design.md Sec 6.3, frozen-selection walk flux from P7)
plus the coupled Gamma unknown (Sec 8.1) replace the Picard/secant with a
Newton iteration that converges to the ACTUAL discrete solution:

  * subsonic: a handful of quadratic steps to ||R||_inf ~ 1e-13, Kutta
    closed to machine precision (the P3 nested Picard/secant needs ~15
    density iterations to 1e-10 and closes Kutta to 1e-8);
  * transonic: Mach continuation with the N5 robustness chain -- direct
    exact steps (the shock-position soft mode stalls eta-accurate Krylov
    steps), stall-adaptive freeze of the upwind assignment with active-set
    refresh (the 2.5D prism-split mesh parks ~1e3 elements in the
    max(nu_e, nu_u) near-tie band; live Newton limit-cycles there), and
    divergence safety nets. Within a frozen assignment the residual is
    smooth and the terminal rate is quadratic (Lopez's frozen-selection
    architecture made persistent).

★ Baseline finding (2026-07-11, roadmap P4 erratum): the P4 Picard
"engineering-converged" states are NOT discrete solutions -- the coupled
Newton residual at the committed coarse M0.80 state is 2.2e-4, and Newton
started FROM it walks (quadratically) to the true solution: shock 0.658,
cl 0.459 -- aft-shock/high-lift, dissipation-robust, continuation-path
independent. Conservative FP overshoots Euler on strong-shock cases
(Holst PAS 2000), so the Euler-anchored G4.1 band does not bind these
solutions. On the medium mesh the solution family steepens toward the FP
non-uniqueness fold (cl 0.396 at M0.775 -> 0.523 at M0.7875; no reachable
isolated solution at M0.80) -- G8.1 was therefore re-specced
(user-approved) to coarse M0.80 + medium M0.7875.

  part 1 (always, ~2 min): subsonic Newton-vs-Picard agreement + coarse
      M0.80 continuation. V8.1a residual histories with terminal-quadratic
      annotation; N2 Jacobian-density measurement (nnz vs Picard).
  part 2 (PYFP3D_TRANSONIC_GATES=1, ~5 min): medium M0.7875 gate run --
      V8.1b convergence + per-stage runtime breakdown (V8.2-lite).
  part 3 (PYFP3D_TRANSONIC_GATES=1 + the gitignored ONERA M6 medium.msh,
      ~4 min): N6 / G8.2 -- M6 medium M0.84/alpha3.06 Newton end to end
      < 5 min (vs the P5 Picard 4539 s), V8.2 runtime breakdown, and the
      P5-caveat measurement (Newton residual at the P5 Picard state is
      ~8e-6 with Kutta |F| ~5e-4; the Newton true solution lifts cl_p
      +5.8% coarse / +7.9% medium with shock positions essentially
      unchanged -- the P4-erratum finding in kind, much milder in
      degree; cl_KJ moves 0.2499 -> 0.2692, narrowing the P9 gap to the
      Tranair/KRATOS 0.288 reference).

Headless; writes results/*.png + checks.csv; exits nonzero on FAIL.
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    CRITICAL, INK_2, MESH_DIR, S1_BLUE, S2_AQUA, S3_YELLOW, S4_ROSE,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve  # noqa: E402
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients  # noqa: E402
from pyfp3d.solve.newton import (  # noqa: E402
    solve_newton_lifting, solve_newton_transonic,
)
from pyfp3d.solve.picard import solve_subsonic_lifting  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
ALPHA = 1.25

#: the N5 transonic recipe (tests/test_p8_newton.py::NEWTON_TRANSONIC_RECIPE)
RECIPE = dict(
    dm=0.025, dm_min=0.003, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct", n_newton_max=60),
)


def _cl(mc, phi, m_inf, alpha):
    dz = float(np.ptp(mc.nodes[:, 2]))
    return float(wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
        alpha_deg=alpha, s_ref=dz, m_inf=m_inf)["cl"])


def _terminal_drops(h):
    """Best consecutive drop pair in the final-level history (the
    freeze-refresh honesty re-evaluations interleave live jumps with the
    quadratic frozen phases, so the last two raw entries can straddle a
    refresh -- same protocol as tests/test_p8_newton.py)."""
    drops = [h[i + 1] / h[i] for i in range(len(h) - 1)]
    best = (np.inf, np.inf)
    for i in range(len(drops) - 1):
        if max(drops[i], drops[i + 1]) < max(best):
            best = (drops[i], drops[i + 1])
    return best


def _plot_history(ax, h, color, label):
    ax.semilogy(range(len(h)), h, "o-", ms=3.5, lw=1.2, color=color,
                label=label)


def part1(cl: CheckList):
    print("\n[1/3] subsonic agreement + coarse M0.80 Newton continuation")
    mc, wc = cut_wake(read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh"))

    # -- subsonic: Newton vs the P3 nested Picard/secant ------------------
    rp = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0)
    rn = solve_newton_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0)
    cl_p = _cl(mc, rp["phi"], 0.5, 2.0)
    cl_n = _cl(mc, rn["phi"], 0.5, 2.0)
    dgamma = float(np.max(np.abs(np.asarray(rn["gamma"])
                                 - np.asarray(rp["gamma"]))))
    cl.add("G8.1", "subsonic cl Newton vs P3 Picard",
           f"{abs(cl_n / cl_p - 1):.2e}", "< 0.5%",
           abs(cl_n / cl_p - 1) < 5e-3,
           note=f"Newton {rn['n_newton']} steps to "
                f"{rn['residual_history'][-1]:.1e}; |dGamma| {dgamma:.1e}")

    # -- coarse transonic continuation ------------------------------------
    t0 = time.perf_counter()
    rt = solve_newton_transonic(mc, wc, m_inf=0.80, alpha_deg=ALPHA,
                                **RECIPE)
    wall = time.perf_counter() - t0
    h = rt["residual_history"]
    d2, d1 = _terminal_drops(h)
    dz = float(np.ptp(mc.nodes[:, 2]))
    rep = shock_report(wall_cp_curve(mc, rt["phi"], z=0.5 * dz, m_inf=0.80),
                       0.80)
    clv = _cl(mc, rt["phi"], 0.80, ALPHA)

    cl.add("G8.1", "coarse M0.80 converged", str(rt["converged"]), "True",
           bool(rt["converged"]), note=f"{wall:.0f} s end-to-end")
    cl.add("G8.1", "coarse terminal residual", f"{h[-1]:.1e}", "< 1e-9",
           h[-1] < 1e-9)
    cl.add("G8.1", "coarse terminal quadratic drops",
           f"{d2:.1e}, {d1:.1e}", "both < 3e-2", d2 < 3e-2 and d1 < 3e-2)
    cl.add("G8.1", "coarse Kutta closure |F|", f"{rt['F_history'][-1]:.1e}",
           "< 1e-12 (machine)", rt["F_history"][-1] < 1e-12)
    cl.add("G8.1", "coarse shock x/c (regression lock)",
           f"{rep['upper']['x_shock']:.4f}", "0.658 +/- 0.012",
           abs(rep["upper"]["x_shock"] - 0.658) < 0.012,
           note="TRUE discrete solution; P4 Picard stall state was 0.604 "
                "(Newton residual there: 2.2e-4)")
    cl.add("G8.1", "coarse cl (regression lock)", f"{clv:.4f}",
           "0.459 +/- 0.010", abs(clv - 0.459) < 0.010)
    cl.add("G8.1", "coarse clean pocket",
           f"lim {rt['n_limited']} flr {rt['n_floored']}", "0 / 0",
           rt["n_limited"] == 0 and rt["n_floored"] == 0)

    # N2 measurement: Jacobian density vs the Picard matrix
    jac_nnz = rt.get("jacobian_nnz")
    n_t3 = rt.get("n_term3_active")
    write_csv(OUT, "n2_jacobian_density.csv", "quantity,value",
              [("newton_nnz", jac_nnz), ("n_term3_active", n_t3),
               ("coarse_wall_s", f"{wall:.1f}"),
               ("n_gmres_stalled", rt.get("n_gmres_stalled"))])

    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    _plot_history(ax, h, S1_BLUE,
                  f"coarse M0.80 final level ({len(h)} its)")
    _plot_history(ax, rn["residual_history"], S2_AQUA,
                  f"subsonic M0.50 ({len(rn['residual_history'])} its)")
    ax.axhline(1e-10, color=CRITICAL, lw=1.0, ls="--", label="tol 1e-10")
    ax.annotate(f"terminal drops {d2:.0e}, {d1:.0e}\n(quadratic tail)",
                xy=(len(h) - 1, h[-1]), xytext=(-120, 30),
                textcoords="offset points", fontsize=9, color=INK_2,
                arrowprops=dict(arrowstyle="->", color=INK_2, lw=0.8))
    ax.set_xlabel("Newton iteration (final Mach level)")
    ax.set_ylabel(r"$\|R\|_\infty$ (free reduced dofs)")
    ax.set_title("V8.1a coupled-Newton convergence, NACA0012 coarse "
                 "(subsonic + M0.80 continuation final level)")
    ax.legend()
    finish(fig, OUT, "v81a_convergence_coarse.png")


def part2(cl: CheckList):
    print("\n[2/3] medium M0.7875 G8.1 gate run [heavy]")
    mc, wc = cut_wake(read_mesh(MESH_DIR / "naca0012_2.5d" / "medium.msh"))
    t0 = time.perf_counter()
    r = solve_newton_transonic(mc, wc, m_inf=0.7875, alpha_deg=ALPHA,
                               **RECIPE)
    wall = time.perf_counter() - t0
    h = r["residual_history"]
    d2, d1 = _terminal_drops(h)
    dz = float(np.ptp(mc.nodes[:, 2]))
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                     m_inf=0.7875), 0.7875)
    clv = _cl(mc, r["phi"], 0.7875, ALPHA)

    cl.add("G8.1", "medium M0.7875 converged", str(r["converged"]), "True",
           bool(r["converged"]), note=f"{wall:.0f} s end-to-end (Picard "
           "G4.1 medium was 16m39s and not a solution)")
    cl.add("G8.1", "medium terminal residual", f"{h[-1]:.1e}", "< 1e-9",
           h[-1] < 1e-9)
    cl.add("G8.1", "medium terminal quadratic drops",
           f"{d2:.1e}, {d1:.1e}", "both < 3e-2", d2 < 3e-2 and d1 < 3e-2)
    cl.add("G8.1", "medium shock x/c (regression lock)",
           f"{rep['upper']['x_shock']:.4f}", "0.674 +/- 0.012",
           abs(rep["upper"]["x_shock"] - 0.674) < 0.012)
    cl.add("G8.1", "medium cl (regression lock)", f"{clv:.4f}",
           "0.523 +/- 0.010", abs(clv - 0.523) < 0.010)
    cl.add("G8.1", "medium clean pocket",
           f"lim {r['n_limited']} flr {r['n_floored']}", "0 / 0",
           r["n_limited"] == 0 and r["n_floored"] == 0)
    ru = r["residual_unfrozen"]
    cl.add("G8.1", "assignment-discontinuity floor (reported)",
           "None" if ru is None else f"{ru:.1e}",
           "< 1e-5 (intrinsic C0-flux floor, honesty metric)",
           ru is None or ru < 1e-5,
           note=f"{r['n_assignment_stale']} stale assignments at accept")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.6),
                                   gridspec_kw={"width_ratios": [1.5, 1]})
    _plot_history(ax1, h, S4_ROSE, f"medium M0.7875 final level ({len(h)} its)")
    ax1.axhline(1e-10, color=CRITICAL, lw=1.0, ls="--", label="tol 1e-10")
    if ru is not None:
        ax1.axhline(ru, color=S3_YELLOW, lw=1.0, ls=":",
                    label=f"live assignment floor {ru:.0e}")
    ax1.set_xlabel("Newton iteration (final Mach level)")
    ax1.set_ylabel(r"$\|R\|_\infty$")
    ax1.set_title(f"V8.1b medium M0.7875 (drops {d2:.0e}, {d1:.0e})")
    ax1.legend()

    tm = r["timings"]
    stages = ["residual", "jacobian", "amg_setup", "gmres"]
    labels = ["residual sweeps", "Jacobian assembly",
              "LU / precond setup", "linear solves"]
    vals = [tm.get(s, 0.0) for s in stages]
    ax2.barh(labels[::-1], vals[::-1],
             color=[S1_BLUE, S2_AQUA, S3_YELLOW, S4_ROSE][::-1])
    for i, v in enumerate(vals[::-1]):
        ax2.text(v, i, f" {v:.1f}s", va="center", fontsize=9, color=INK_2)
    ax2.set_xlabel("wall-clock seconds (final level)")
    ax2.set_title("V8.2-lite per-stage runtime")
    finish(fig, OUT, "v81b_convergence_medium.png")

    write_csv(OUT, "g81_medium.csv", "quantity,value",
              [("m_inf", 0.7875), ("converged", r["converged"]),
               ("residual_final", f"{h[-1]:.3e}"),
               ("residual_unfrozen", ru),
               ("n_assignment_stale", r["n_assignment_stale"]),
               ("shock_x_c", f"{rep['upper']['x_shock']:.4f}"),
               ("cl_p", f"{clv:.4f}"),
               ("mach_max", f"{np.sqrt(r['mach2_max']):.4f}"),
               ("wall_s", f"{wall:.1f}"),
               ("levels", len(r["level_history"]))])


#: G8.2 physics regression locks = the measured M6 medium Newton solution
#: (2026-07-11 N6 run; same lock-band semantics as the G8.1 constants)
M6_LOCK = dict(cl_p=0.2646, shock44=0.596, shock65=0.541, shock90=0.362,
               m_max=2.134)

#: the N6 M6 recipe (tests/test_p8_newton.py::NEWTON_M6_RECIPE): N5 chain +
#: lagged-LU direct steps (true-3D LU fill: 18.6 s/refactor at 63k dofs) +
#: the P5 dm=0.05 schedule (far from the NACA-medium fold)
M6_RECIPE = dict(
    dm=0.05, dm_min=0.01, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   direct_refactor_every=1000, n_newton_max=60,
                   farfield_spanwise_gamma=True),
)


def part3(cl: CheckList):
    print("\n[3/3] N6/G8.2: ONERA M6 medium M0.84 Newton end to end [heavy]")
    from pyfp3d.meshgen.wing3d import B_SEMI
    from pyfp3d.post.section_cut import section_cp_curve
    from pyfp3d.post.surface import cl_kj_3d, planform_area
    from pyfp3d.solve.newton import NewtonWorkspace

    msh = MESH_DIR / "onera_m6" / "medium.msh"
    if not msh.exists():
        print("  onera_m6/medium.msh not generated "
              "(cases/meshes/onera_m6/generate_onera_m6.py); part 3 skipped")
        return
    m_inf, alpha = 0.84, 3.06

    # ---- G8.2 end-to-end clock: mesh -> cut -> solve -> forces + shocks --
    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(msh))
    r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                               **M6_RECIPE)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
    cl_kj = cl_kj_3d(r["gamma"], wc.station_z, s_ref=s_ref, b_semi=B_SEMI)
    shocks = {}
    for eta in (0.44, 0.65, 0.90):
        c = section_cp_curve(mc, r["phi"], eta=eta, b_semi=B_SEMI,
                             m_inf=m_inf)
        shocks[eta] = shock_report(c, m_inf)["upper"]["x_shock"]
    wall = time.perf_counter() - t0
    # ---- end-to-end clock stops -------------------------------------------

    h = r["residual_history"]
    d2, d1 = _terminal_drops(h)
    clv = float(forces["cl"])
    m_max = float(np.sqrt(r["mach2_max"]))

    cl.add("G8.2", "M6 medium end to end", f"{wall:.0f} s",
           "< 300 s single node", wall < 300.0,
           note=f"P5 Picard medium solve was 4539 s")
    cl.add("G8.2", "M6 medium converged", str(r["converged"]), "True",
           bool(r["converged"]),
           note=f"levels {len(r['level_history'])}, final |R| {h[-1]:.1e}")
    cl.add("G8.2", "M6 terminal quadratic drops", f"{d2:.1e}, {d1:.1e}",
           "both < 3e-2", d2 < 3e-2 and d1 < 3e-2)
    cl.add("G8.2", "M6 Kutta closure |F|", f"{r['F_history'][-1]:.1e}",
           "< 1e-12 (machine; P5 secant+polish was 5.8e-4)",
           r["F_history"][-1] < 1e-12)
    cl.add("G8.2", "M6 clean pocket",
           f"lim {r['n_limited']} flr {r['n_floored']}", "0 / 0",
           r["n_limited"] == 0 and r["n_floored"] == 0)
    cl.add("G8.2", "M6 cl_p (regression lock)", f"{clv:.4f}",
           f"{M6_LOCK['cl_p']} +/- 0.005",
           abs(clv - M6_LOCK["cl_p"]) < 0.005,
           note=f"cl_KJ {cl_kj:.4f}; P5 Picard cl_p 0.2453 (caveat, below)")
    cl.add("G8.2", "M6 shock x/c eta 0.44/0.65/0.90 (regression lock)",
           f"{shocks[0.44]:.3f}/{shocks[0.65]:.3f}/{shocks[0.90]:.3f}",
           f"{M6_LOCK['shock44']}/{M6_LOCK['shock65']}/{M6_LOCK['shock90']}"
           " +/- 0.02",
           all(abs(shocks[e] - M6_LOCK[k]) < 0.02 for e, k in
               [(0.44, "shock44"), (0.65, "shock65"), (0.90, "shock90")]))
    cl.add("G8.2", "M6 M_max (regression lock)", f"{m_max:.3f}",
           f"{M6_LOCK['m_max']} +/- 0.05",
           abs(m_max - M6_LOCK["m_max"]) < 0.05)

    # ---- P5-caveat measurement: Newton residual at the P5 Picard state --
    p5 = Path(__file__).resolve().parents[1] / "p5_onera_m6" / "results" \
        / "medium_solution.npz"
    res_p5 = None
    if p5.exists():
        d = np.load(p5)
        ws = NewtonWorkspace(mc, wc, alpha_deg=alpha,
                             farfield_spanwise_gamma=True)
        ws.set_mach(m_inf)
        R_free, F, _ = ws.eval_residual(
            d["phi"][:ws.n_red][ws.free], d["gamma"], 1.5, 0.95, 3.0, 0.05)
        res_p5 = float(np.max(np.abs(R_free)))
        cl.add("G8.2", "P5-caveat: Newton residual at P5 Picard state",
               f"|R| {res_p5:.1e}, |F| {float(np.max(np.abs(F))):.1e}",
               "reported (P4-erratum-in-kind measurement)", True,
               note=f"Newton true solution: cl_p {clv:.4f} vs P5 0.2453, "
                    f"M_max {m_max:.3f} vs 1.995")
    else:
        print("  p5 medium_solution.npz cache absent; caveat row skipped")

    # ---- V8.2: runtime breakdown chart (timings = FINAL level only; the
    # earlier continuation levels + mesh/cut/post are the remainder bar) --
    tm = r["timings"]
    stages = ["residual", "jacobian", "amg_setup", "gmres"]
    labels = ["residual sweeps", "Jacobian assembly",
              "LU factorization", "GMRES + LU solves"]
    vals = [tm.get(s, 0.0) for s in stages]
    other = wall - sum(vals)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.6),
                                   gridspec_kw={"width_ratios": [1.5, 1]})
    _plot_history(ax1, h, S1_BLUE, f"M6 medium M0.84 final level ({len(h)} its)")
    ax1.axhline(1e-10, color=CRITICAL, lw=1.0, ls="--", label="tol 1e-10")
    ax1.set_xlabel("Newton iteration (final Mach level)")
    ax1.set_ylabel(r"$\|R\|_\infty$")
    ax1.set_title(f"V8.1c ONERA M6 medium M0.84 (drops {d2:.0e}, {d1:.0e})")
    ax1.legend()
    bl = labels + ["earlier levels + mesh/cut/post"]
    bv = vals + [other]
    colors = [S1_BLUE, S2_AQUA, S3_YELLOW, S4_ROSE, INK_2]
    ax2.barh(bl[::-1], bv[::-1], color=colors[::-1])
    for i, v in enumerate(bv[::-1]):
        ax2.text(v, i, f" {v:.1f}s", va="center", fontsize=9, color=INK_2)
    ax2.set_xlabel("wall-clock seconds (final Mach level; remainder bar = rest)")
    ax2.set_title(f"V8.2 M6 medium runtime breakdown ({wall:.0f}s end to end)")
    finish(fig, OUT, "v82_m6_medium.png")

    write_csv(OUT, "g82_m6_medium.csv", "quantity,value",
              [("m_inf", m_inf), ("alpha_deg", alpha),
               ("wall_end_to_end_s", f"{wall:.1f}"),
               ("converged", r["converged"]),
               ("residual_final", f"{h[-1]:.3e}"),
               ("kutta_F", f"{r['F_history'][-1]:.3e}"),
               ("cl_p", f"{clv:.4f}"), ("cl_kj", f"{cl_kj:.4f}"),
               ("mach_max", f"{m_max:.4f}"),
               ("shock_eta44", f"{shocks[0.44]:.4f}"),
               ("shock_eta65", f"{shocks[0.65]:.4f}"),
               ("shock_eta90", f"{shocks[0.90]:.4f}"),
               ("newton_residual_at_p5_state",
                "n/a" if res_p5 is None else f"{res_p5:.3e}"),
               ("levels", len(r["level_history"])),
               ("p5_picard_solve_s", 4539)])


def main():
    apply_style()
    cl = CheckList("P8 fully-coupled Newton (G8.1 + N6 G8.2)")
    part1(cl)
    if os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1":
        part2(cl)
        part3(cl)
    else:
        print("\n[2/3] medium gate + [3/3] M6 runs skipped (~10 min); "
              "set PYFP3D_TRANSONIC_GATES=1 to run")
    sys.exit(cl.report(OUT))


if __name__ == "__main__":
    main()
