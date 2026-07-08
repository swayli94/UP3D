"""
P6 demo -- removing the non-physical surface-Cp sawtooth (gate G6.1/G6.2).

Root-cause finding (2026-07-08): the ~2-cell wall-Cp serration in the
supersonic run is a per-triangle gradient-RECOVERY artifact on the sliver
prism-split wall triangulation -- NOT the artificial-density flux. Decisive
evidence: nodal/edge-neighbour smoothing of the SAME solution's wall gradient
drops the G6.1 sawtooth metric ~330x, while a smoother artificial-density flux
(the streamline kernel) does NOT reduce it at all. So P6's G6.1 fix is a
normal-gated wall-gradient smoothing in post-processing (`smooth_passes`,
preserves the sharp TE), applied to BOTH the Cp curve and the force integral.

This re-runs the P4 and P5 validation CASES through the walk solver and shows,
on the exact gate cases:
  part 1 (always, ~6 min): NACA0012 M=0.80 alpha=1.25 coarse (P4 G4.1).
      raw per-triangle Cp vs smoothed Cp (G6.1 metric), physics preserved
      (G6.2: shock/M_max/cl unchanged -- smoothing is post-processing). Also
      runs the streamline-kernel flux to DOCUMENT that it does not fix the
      sawtooth (the negative result that redirected P6).
  part 2 (PYFP3D_TRANSONIC_GATES=1, heavy): ONERA M6 M=0.84 coarse (P5 G5.1),
      section Cp at eta=0.65 raw vs smoothed.

Headless; writes results/*.png + checks.csv; exits nonzero on unexpected FAIL.
"""

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, INK_2, MESH_DIR, S1_BLUE, S4_ROSE,
    CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.post.section_cut import (  # noqa: E402
    cp_oscillation_metric, section_cp_curve, wall_cp_curve,
)
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients  # noqa: E402
from pyfp3d.solve.continuation import solve_transonic_lifting  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
SMOOTH_PASSES = 1  # G6.1 recovery smoothing (normal-gated; preserves the TE)


def part1_naca(cl: CheckList):
    print("\n[1/2] NACA0012 M=0.80 alpha=1.25 coarse -- P4 G4.1 case")
    mesh = read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)
    dz = float(np.ptp(mc.nodes[:, 2]))

    # Same walk solution; only the Cp/force RECOVERY changes (post-processing).
    r = solve_transonic_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25,
                                max_gamma_evals=12, n_picard_eval=800)
    raw = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=0.80, smooth_passes=0)
    sm = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=0.80,
                       smooth_passes=SMOOTH_PASSES)
    rep_raw = shock_report(raw, 0.80)
    rep_sm = shock_report(sm, 0.80)
    cp_star = rep_raw["cp_critical"]
    m_raw = cp_oscillation_metric(raw["x_upper"], raw["cp_upper"], cp_star)
    m_sm = cp_oscillation_metric(sm["x_upper"], sm["cp_upper"], cp_star)
    f_raw = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], r["phi"],
                                    alpha_deg=1.25, s_ref=dz, m_inf=0.80)
    f_sm = wall_force_coefficients(mc.nodes, mc.elements,
                                   mc.boundary_faces["wall"], r["phi"],
                                   alpha_deg=1.25, s_ref=dz, m_inf=0.80,
                                   smooth_passes=SMOOTH_PASSES)

    # Negative control: the streamline-kernel flux (opt-in) -- does it reduce
    # the RAW sawtooth? (finding: no; the sawtooth is a recovery artifact.)
    rk = solve_transonic_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25,
                                 max_gamma_evals=12, n_picard_eval=800,
                                 upwind_weighted=True, upwind_mode="kernel",
                                 upwind_c=2.5, upwind_reach_frac=1.0)
    ck = wall_cp_curve(mc, rk["phi"], z=0.5 * dz, m_inf=0.80, smooth_passes=0)
    repk = shock_report(ck, 0.80)
    m_kernel = cp_oscillation_metric(ck["x_upper"], ck["cp_upper"],
                                     repk["cp_critical"])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(raw["x_upper"], raw["cp_upper"], ".-", ms=4, lw=0.8, color=BASELINE,
            label=f"raw per-triangle (metric {m_raw['metric']:.3f})")
    ax.plot(sm["x_upper"], sm["cp_upper"], ".-", ms=4, lw=1.2, color=S1_BLUE,
            label=f"G6.1 smoothed (metric {m_sm['metric']:.4f})")
    ax.axhline(cp_star, color=INK_2, lw=0.8, ls=":", label="Cp* (sonic)")
    ax.invert_yaxis(); ax.set_xlabel("x/c"); ax.set_ylabel("Cp")
    ax.set_title("V6.1 NACA0012 M0.80 wall Cp: recovery smoothing removes the sawtooth")
    ax.legend()
    finish(fig, OUT, "g61_cp_raw_vs_smoothed_coarse.png")

    clw_kj = 2.0 * float(r["gamma"][0])
    write_csv(OUT, "g61_naca_coarse.csv", "quantity,raw,smoothed",
              [("sawtooth_metric_upper", f"{m_raw['metric']:.4e}",
                f"{m_sm['metric']:.4e}"),
               ("n_reversals_upper", m_raw["n_reversals"], m_sm["n_reversals"]),
               ("upper_shock_x_c", f"{rep_raw['upper']['x_shock']:.4f}",
                f"{rep_sm['upper']['x_shock']:.4f}"),
               ("cl_pressure", f"{f_raw['cl']:.4f}", f"{f_sm['cl']:.4f}"),
               ("cd_pressure", f"{f_raw['cd_pressure']:.5f}",
                f"{f_sm['cd_pressure']:.5f}"),
               ("cl_kj (Gamma, unchanged)", f"{clw_kj:.4f}", f"{clw_kj:.4f}"),
               ("mach_max (solution, unchanged)",
                f"{np.sqrt(r['mach2_max']):.4f}", f"{np.sqrt(r['mach2_max']):.4f}")])
    write_csv(OUT, "g61_kernel_negative.csv", "quantity,walk_raw,kernel_raw",
              [("sawtooth_metric_upper", f"{m_raw['metric']:.4e}",
                f"{m_kernel['metric']:.4e}"),
               ("note", "flux-change", "does-not-reduce-sawtooth")])

    # G6.1: recovery smoothing drops the sawtooth metric far below the raw one.
    cl.add("G6.1", "sawtooth metric raw->smoothed",
           f"{m_raw['metric']:.3f} -> {m_sm['metric']:.4f}",
           "smoothed << raw (>=10x)", m_sm["metric"] <= 0.1 * m_raw["metric"])
    cl.add("G6.1", "reversal count raw->smoothed",
           f"{m_raw['n_reversals']} -> {m_sm['n_reversals']}",
           "few reversals left", m_sm["n_reversals"] <= 3)
    # G6.2: post-processing -- the SOLUTION (shock, M_max, cl_KJ) is unchanged;
    # smoothed shock stays in band and smoothed cl_p tracks the raw within tol.
    cl.add("G6.2", "smoothed shock x/c in band",
           f"{rep_sm['upper']['x_shock']:.3f}", "0.62 +/- 0.03",
           abs(rep_sm["upper"]["x_shock"] - 0.62) <= 0.03)
    cl.add("G6.2", "smoothed cl_p vs raw cl_p",
           f"{f_sm['cl']:.4f} vs {f_raw['cl']:.4f}", "within 3%",
           abs(f_sm["cl"] - f_raw["cl"]) <= 0.03 * abs(f_raw["cl"]))
    # Negative result (documented): the flux change does NOT fix the sawtooth.
    cl.add("finding", "kernel flux vs walk raw sawtooth",
           f"{m_kernel['metric']:.3f} vs {m_raw['metric']:.3f}",
           "kernel does NOT reduce it (recorded)",
           m_kernel["metric"] >= 0.5 * m_raw["metric"],
           note="sawtooth is a recovery artifact, not a flux artifact")


def part2_onera_m6(cl: CheckList):
    print("\n[2/2] ONERA M6 M=0.84 alpha=3.06 coarse -- P5 G5.1 case [heavy]")
    mesh_path = MESH_DIR / "onera_m6" / "coarse.msh"
    if not mesh_path.exists():
        print("  M6 coarse.msh missing; skipping (run generate_onera_m6.py)")
        return
    mc, wc = cut_wake(read_mesh(mesh_path))
    b_semi = float(mc.nodes[:, 2].max())
    r = solve_transonic_lifting(mc, wc, m_inf=0.84, alpha_deg=3.06,
                                n_picard_seed=40, n_picard_eval=300,
                                max_gamma_evals=10, rtol=1e-7, n_kutta_polish=4,
                                farfield_spanwise_gamma=True)
    eta = 0.65
    raw = section_cp_curve(mc, r["phi"], eta=eta, b_semi=b_semi, m_inf=0.84,
                           smooth_passes=0)
    sm = section_cp_curve(mc, r["phi"], eta=eta, b_semi=b_semi, m_inf=0.84,
                          smooth_passes=SMOOTH_PASSES)
    cps = shock_report(raw, 0.84)["cp_critical"]
    m_raw = cp_oscillation_metric(raw["x_upper"], raw["cp_upper"], cps)
    m_sm = cp_oscillation_metric(sm["x_upper"], sm["cp_upper"], cps)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(raw["x_upper"], raw["cp_upper"], ".-", ms=4, lw=0.8, color=BASELINE,
            label=f"raw (metric {m_raw['metric']:.3f})")
    ax.plot(sm["x_upper"], sm["cp_upper"], ".-", ms=4, lw=1.2, color=S1_BLUE,
            label=f"G6.1 smoothed (metric {m_sm['metric']:.4f})")
    ax.axhline(cps, color=INK_2, lw=0.8, ls=":", label="Cp*")
    ax.invert_yaxis(); ax.set_xlabel("x/c"); ax.set_ylabel("Cp")
    ax.set_title(f"V6.1 ONERA M6 coarse eta={eta}: recovery smoothing")
    ax.legend()
    finish(fig, OUT, "g61_m6_section_cp_coarse.png")
    write_csv(OUT, "g61_m6_coarse.csv", "quantity,raw,smoothed",
              [("sawtooth_metric_eta065", f"{m_raw['metric']:.4e}",
                f"{m_sm['metric']:.4e}"),
               ("mach_max (unchanged)", f"{np.sqrt(r['mach2_max']):.4f}",
                f"{np.sqrt(r['mach2_max']):.4f}")])
    cl.add("G6.1 M6", "section-Cp metric raw->smoothed",
           f"{m_raw['metric']:.3f} -> {m_sm['metric']:.4f}",
           "smoothed << raw", m_sm["metric"] < 0.5 * m_raw["metric"])


def main():
    apply_style()
    cl = CheckList("P6 differentiable-flux / sawtooth removal (G6.1-G6.2)")
    part1_naca(cl)
    if os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1":
        part2_onera_m6(cl)
    else:
        print("\n[2/2] ONERA M6 raw-vs-smoothed skipped (~10-20 min); set "
              "PYFP3D_TRANSONIC_GATES=1 to run")
    sys.exit(cl.report(OUT))


if __name__ == "__main__":
    main()
