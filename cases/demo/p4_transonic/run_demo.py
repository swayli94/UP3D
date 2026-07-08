"""
P4 demo -- transonic artificial-density evidence (gates G4.1-G4.3).

What this shows, per docs/roadmap.md P4 and docs/demo_report.md:
  1. G4.2: the upwind machinery is an EXACT no-op below critical Mach --
     M = 0.5 with the upstream walk + nu switch + rho_tilde sweep active
     is bit-identical to the P3 path (nu == 0.0 exactly).
  2. The dissipation-starvation root cause and its fix: on the M0
     prism-split meshes a single-hop upstream neighbor reaches only
     ~1/3 of an element's streamwise extent, starving the artificial
     density of exactly the dissipation the (M^2-1)/M^2 stability bound
     demands -- the multi-hop walk restores a genuine one-cell reach
     (histogram evidence).
  3. G4.1 (coarse evidence run): NACA0012 M = 0.80 alpha = 1.25 deg via
     Mach continuation + outer Gamma secant: monotone upper shock at
     x/c ~ 0.60 within the reference band (Euler anchor + documented
     conservative-FP shift, cases/reference_data/naca0012_m080/), weak
     lower shock ~ 0.36, no expansion shock, physical M_max ~ 1.36,
     no limited/floored cells.
  4. G4.1 mesh-refinement comparison (coarse vs medium) and the G4.3
     10-case robustness sweep -- the heavy evidence behind demo_report
     Sec P4 "supplementary analysis" (the coarse/medium surface-Cp
     sawtooth study + the shock-migration/lift-trend dashboard). These
     two solves are ~16 min and ~22 min, so they are OPT-IN behind
     PYFP3D_TRANSONIC_GATES=1 (same switch as the pytest gate); the
     committed figures/CSVs in results/ are the reference baseline and
     are what demo_report embeds, so a default run does not need them.

*** COST CAUTION (read before regenerating) ***
  Part 4 + part 5 together are ~40 min of Picard iterations. The figures
  and CSVs they produce are ALREADY COMMITTED in results/ and are what
  demo_report embeds. Do NOT re-run heavy mode casually -- only when the
  solver/mesh/reference actually changed in a way that would move these
  numbers, and you intend to commit the refreshed baseline. For a routine
  edit, the committed baseline is authoritative; verify the cheap coarse
  path (default run, ~4 min) instead. The pytest gate
  (tests/test_p4_transonic.py) is the correctness check; this heavy mode
  only refreshes the committed evidence PNGs.

Standalone + self-checking:  python cases/demo/p4_transonic/run_demo.py
  default (coarse only, ~4 min):  parts 1-3
  full evidence (~40 min):        PYFP3D_TRANSONIC_GATES=1 python ...
Outputs: cases/demo/p4_transonic/results/{*.png, *.csv, checks.csv}
Exit code 0 iff every acceptance check passes.
"""

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    CRITICAL, INK_2, MESH_DIR, MUTED, REFERENCE_DIR, S1_BLUE, S2_AQUA,
    S3_YELLOW, S4_ROSE, CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.kernels.jacobian import PicardOperator  # noqa: E402
from pyfp3d.kernels.upwind import UpwindOperator, upstream_elements  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field  # noqa: E402
from pyfp3d.post.section_cut import (  # noqa: E402
    cp_oscillation_metric,
    wall_cp_curve,
)
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients  # noqa: E402
from pyfp3d.solve.continuation import solve_transonic_lifting  # noqa: E402
from pyfp3d.solve.picard import solve_subsonic_lifting  # noqa: E402

OUT = Path(__file__).parent / "results"
M_INF, ALPHA = 0.80, 1.25


# ---------------------------------------------------------------------------
# Figure/CSV generators for the heavy G4.1 (coarse vs medium) + G4.3 evidence.
# These moved here from the P4 gate test (tests/test_p4_transonic.py) so the
# committed artifacts and the code that makes them live together in the demo
# (artifacts/ is gitignored, so demo_report cannot embed from there).
# ---------------------------------------------------------------------------
def _plot_cp_shock(curve, rep, level: str, out_name: str):
    """Matched-style Cp(x/c) + shock-marker figure for the coarse/medium
    refinement comparison (deliberately the same plain style at both levels
    so the sawtooth amplitude is visually comparable across meshes)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(curve["x_upper"], curve["cp_upper"], ".-", ms=3, lw=0.7,
            color=S4_ROSE, label="upper")
    ax.plot(curve["x_lower"], curve["cp_lower"], ".-", ms=3, lw=0.7,
            color=S1_BLUE, label="lower")
    ax.axhline(rep["cp_critical"], color=INK_2, lw=0.8, ls=":",
               label="Cp* (sonic)")
    for side, color in (("upper", S4_ROSE), ("lower", S1_BLUE)):
        if rep[side]["has_shock"]:
            ax.axvline(rep[side]["x_shock"], color=color, lw=0.8, ls="--")
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel("Cp")
    ax.set_title(f"V4.1 NACA0012 M={M_INF} alpha={ALPHA} ({level}): "
                 f"upper shock x/c={rep['upper']['x_shock']:.3f}")
    ax.legend()
    finish(fig, OUT, out_name)


def _g41_summary_rows(rep, r, forces, curve):
    rows = []
    for side in ("upper", "lower"):
        for k, v in rep[side].items():
            rows.append((f"{side}_{k}", v))
    # G6.1 baseline (roadmap P6): surface-Cp sawtooth metric on the supersonic
    # run. Upper is the dominant supersonic pocket and the primary gate target;
    # lower is reported (often a marginal pocket with a tiny Cp range).
    cp_star = rep["cp_critical"]
    osc_u = cp_oscillation_metric(curve["x_upper"], curve["cp_upper"], cp_star)
    osc_l = cp_oscillation_metric(curve["x_lower"], curve["cp_lower"], cp_star)
    rows += [
        ("cp_critical", rep["cp_critical"]),
        ("cl_pressure", forces["cl"]),
        ("cl_kj", 2.0 * float(r["gamma"][0])),
        ("gamma", float(r["gamma"][0])),
        ("mach_max", float(np.sqrt(r["mach2_max"]))),
        ("kutta_mismatch", r["kutta_mismatch"]),
        ("n_picard_total", r["n_picard_total"]),
        ("n_limited", r["n_limited"]),
        ("n_floored", r["n_floored"]),
        ("g61_cp_osc_upper", osc_u["metric"]),
        ("g61_n_super_upper", osc_u["n_super"]),
        ("g61_cp_osc_lower", osc_l["metric"]),
        ("g61_n_super_lower", osc_l["n_super"]),
    ]
    return rows


def _solve_case(mesh_path, m_inf, alpha):
    mesh = read_mesh(mesh_path)
    mc, wc = cut_wake(mesh)
    r = solve_transonic_lifting(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                max_gamma_evals=12, n_picard_eval=800,
                                verbose=True)
    dz = float(np.ptp(mc.nodes[:, 2]))
    curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf)
    rep = shock_report(curve, m_inf)
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=alpha, s_ref=dz, m_inf=m_inf,
    )
    return {"r": r, "curve": curve, "rep": rep, "forces": forces}


def part1_noop(cl: CheckList, mc, wc):
    print("\n[1/3] G4.2 subcritical no-op (M = 0.5)")
    r_p4 = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0, upwind_c=1.5)
    r_p3 = solve_subsonic_lifting(mc, wc, m_inf=0.5, alpha_deg=2.0, upwind_c=0.0)
    bits = bool(np.array_equal(r_p4["phi"], r_p3["phi"])
                and np.array_equal(r_p4["gamma"], r_p3["gamma"]))
    cl.add("G4.2", "nu_max at M=0.5", r_p4["nu_max"], "== 0.0",
           r_p4["nu_max"] == 0.0)
    cl.add("G4.2", "phi/Gamma vs P3 path", bits, "bitwise identical", bits)
    write_csv(OUT, "g42_noop.csv", "check,value",
              [("nu_max", r_p4["nu_max"]),
               ("bitwise_identical", bits),
               ("mach2_max", f"{r_p4['mach2_max']:.4f}")])


def part2_reach(cl: CheckList, mc, wc):
    print("\n[2/3] Upwind-reach evidence (dissipation starvation root cause)")
    r = solve_subsonic_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,
                               omega=0.9, tol_rho=1e-6, n_picard_max=300,
                               forcing=0.05)
    op = PicardOperator(mc.nodes, mc.elements)
    upw = UpwindOperator(mc.nodes, mc.elements)
    grad, q2 = op.velocities(r["phi"])

    # one-hop reach vs the shipped multi-hop walk, in units of the
    # element's own streamwise extent.
    u1 = np.empty(op.n_tets, dtype=np.int64)
    upstream_elements(upw.face_neighbors, upw.centroids, upw._nodes,
                      upw._elements, grad, u1, 1)
    um = np.empty(op.n_tets, dtype=np.int64)
    upstream_elements(upw.face_neighbors, upw.centroids, upw._nodes,
                      upw._elements, grad, um)
    vn = np.maximum(np.linalg.norm(grad, axis=1), 1e-30)
    vhat = grad / vn[:, None]
    pts = mc.nodes[mc.elements]
    proj = np.einsum("nkj,nj->nk", pts, vhat)
    extent = proj.max(axis=1) - proj.min(axis=1)
    d1 = -np.einsum("ij,ij->i", upw.centroids[u1] - upw.centroids, vhat) / extent
    dm = -np.einsum("ij,ij->i", upw.centroids[um] - upw.centroids, vhat) / extent
    act = mach_squared_field(q2, 0.70) > 1.0

    fig, ax = plt.subplots(figsize=(7, 4.6))
    bins = np.linspace(0, 2.0, 33)
    ax.hist(d1[act], bins=bins, alpha=0.65, color=S4_ROSE,
            label=f"single hop (median {np.median(d1[act]):.2f})")
    ax.hist(dm[act], bins=bins, alpha=0.65, color=S1_BLUE,
            label=f"multi-hop walk (median {np.median(dm[act]):.2f})")
    ax.axvline(1.0, color=CRITICAL, lw=1.2, ls="--",
               label="one streamwise cell (design intent)")
    ax.set_xlabel("upstream reach / element streamwise extent")
    ax.set_ylabel("supersonic elements")
    ax.set_title("why single-hop upwinding starved the scheme (M0.70 pocket)")
    ax.legend()
    finish(fig, OUT, "p4_upwind_reach.png")

    med1, medm = float(np.median(d1[act])), float(np.median(dm[act]))
    cl.add("P4-scheme", "single-hop median reach", f"{med1:.2f} extents",
           "documented deficiency (< 0.5)", med1 < 0.5)
    cl.add("P4-scheme", "multi-hop median reach", f"{medm:.2f} extents",
           ">= 0.8 (design intent ~1 cell)", medm >= 0.8)
    return r


def part3_transonic(cl: CheckList, mc, wc):
    print(f"\n[3/3] G4.1 evidence run: M = {M_INF}, alpha = {ALPHA} (coarse)")
    t0 = time.perf_counter()
    r = solve_transonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                max_gamma_evals=12, n_picard_eval=800,
                                verbose=True)
    t_solve = time.perf_counter() - t0
    dz = float(np.ptp(mc.nodes[:, 2]))
    curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=M_INF)
    rep = shock_report(curve, M_INF)
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF,
    )

    ref = {}
    with open(REFERENCE_DIR / "naca0012_m080" / "shock_reference.csv") as f:
        for row in csv.DictReader(f):
            ref[row["quantity"]] = (float(row["value"]), float(row["tolerance"]))

    up, lo = rep["upper"], rep["lower"]
    fig, ax = plt.subplots(figsize=(8.2, 6))
    ax.plot(curve["x_upper"], curve["cp_upper"], ".-", ms=3.5, lw=0.8,
            color=S4_ROSE, label="upper")
    ax.plot(curve["x_lower"], curve["cp_lower"], ".-", ms=3.5, lw=0.8,
            color=S1_BLUE, label="lower")
    ax.axhline(rep["cp_critical"], color=INK_2, lw=1.0, ls=":",
               label="Cp* (sonic)")
    x_ref, tol = ref["upper_shock_x_c"]
    ax.axvspan(x_ref - tol, x_ref + tol, color=S2_AQUA, alpha=0.15,
               label="reference shock band (see README)")
    for side, color in (("upper", S4_ROSE), ("lower", S1_BLUE)):
        if rep[side]["has_shock"]:
            ax.axvline(rep[side]["x_shock"], color=color, lw=1.0, ls="--")
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel("Cp")
    ax.set_title(f"V4.1 NACA0012 M={M_INF} alpha={ALPHA} (coarse): "
                 f"upper shock x/c = {up['x_shock']:.3f}, monotone, "
                 f"M_max = {np.sqrt(r['mach2_max']):.2f}")
    ax.legend()
    finish(fig, OUT, "g41_cp_shock.png")

    cl.add("G4.1", "converged (physical + Kutta)", r["converged"],
           "Gamma secant |F| < 2e-4, no limited cells",
           bool(r["converged"]), note=f"{t_solve:.0f} s coarse")
    cl.add("G4.1", "upper shock x/c", f"{up['x_shock']:.3f}",
           f"{x_ref} +/- {tol} (ref band)",
           abs(up["x_shock"] - x_ref) <= tol)
    cl.add("G4.1", "monotone shock, no expansion shock",
           f"monotone={up['monotone']}, exp={up['expansion_shock']}",
           "monotone and no expansion shock",
           up["monotone"] and not up["expansion_shock"])
    cl.add("G4.1", "shock sharpness", f"{up['n_cells']} station(s)",
           "<= 3 (2-3 cells)", up["n_cells"] <= 3)
    cl.add("G4.1", "lower weak shock x/c", f"{lo['x_shock']:.3f}",
           "~0.35 (reported, not gated)", True)
    cl.add("G4.1", "M_max", f"{np.sqrt(r['mach2_max']):.3f}",
           "physical (~1.3-1.5), below limiter cap", r["mach2_max"] < 9.0)
    write_csv(OUT, "g41_transonic.csv", "quantity,value",
              [("x_shock_upper", f"{up['x_shock']:.4f}"),
               ("x_shock_lower", f"{lo['x_shock']:.4f}"),
               ("cl_pressure", f"{forces['cl']:.5f}"),
               ("cl_kj", f"{2 * float(r['gamma'][0]):.5f}"),
               ("gamma", f"{float(r['gamma'][0]):.5f}"),
               ("mach_max", f"{np.sqrt(r['mach2_max']):.4f}"),
               ("kutta_mismatch", f"{r['kutta_mismatch']:.2e}"),
               ("n_picard_total", r["n_picard_total"]),
               ("t_solve_s", f"{t_solve:.0f}")])
    return {"r": r, "curve": curve, "rep": rep, "forces": forces}


def part4_refinement(cl: CheckList, coarse_case):
    """[heavy] G4.1 coarse-vs-medium refinement pair (demo_report Sec P4
    supplementary): re-emit the matched coarse figure from the part-3 solve
    (no re-solve), then run the full medium gate (~16 min) and emit its
    matched figure so the surface-Cp sawtooth is comparable across meshes."""
    print("\n[4/5] G4.1 mesh-refinement comparison (coarse vs medium) [heavy]")
    _plot_cp_shock(coarse_case["curve"], coarse_case["rep"], "coarse",
                   "g41_cp_shock_coarse.png")
    write_csv(OUT, "g41_summary_coarse.csv", "quantity,value",
              _g41_summary_rows(coarse_case["rep"], coarse_case["r"],
                                coarse_case["forces"], coarse_case["curve"]))

    t0 = time.perf_counter()
    med = _solve_case(MESH_DIR / "naca0012_2.5d" / "medium.msh", M_INF, ALPHA)
    t_solve = time.perf_counter() - t0
    r, curve, rep, forces = (med["r"], med["curve"], med["rep"], med["forces"])
    _plot_cp_shock(curve, rep, "medium", "g41_cp_shock_medium.png")
    write_csv(OUT, "g41_summary_medium.csv", "quantity,value",
              _g41_summary_rows(rep, r, forces, curve))

    ref = {}
    with open(REFERENCE_DIR / "naca0012_m080" / "shock_reference.csv") as f:
        for row in csv.DictReader(f):
            ref[row["quantity"]] = (float(row["value"]), float(row["tolerance"]))
    x_ref, tol = ref["upper_shock_x_c"]
    up = rep["upper"]
    cl.add("G4.1 medium", "converged (physical + Kutta)", r["converged"],
           "Kutta |F| < 2e-4, no limited cells", bool(r["converged"]),
           note=f"{t_solve:.0f} s medium")
    cl.add("G4.1 medium", "upper shock x/c", f"{up['x_shock']:.3f}",
           f"{x_ref} +/- {tol} (ref band)", abs(up["x_shock"] - x_ref) <= tol)
    cl.add("G4.1 medium", "M_max", f"{np.sqrt(r['mach2_max']):.3f}",
           "physical, below limiter cap", r["mach2_max"] < 9.0)
    cl.add("G4.1 medium", "no limited/floored cells",
           f"{r['n_limited']}/{r['n_floored']}", "0/0",
           r["n_limited"] == 0 and r["n_floored"] == 0)


def part5_sweep(cl: CheckList, mc, wc):
    """[heavy] G4.3 10-case robustness sweep + V4.3 dashboard (~22 min)."""
    print("\n[5/5] G4.3 robustness sweep (10 cases) [heavy]")
    dz = float(np.ptp(mc.nodes[:, 2]))
    ms = [0.74, 0.76, 0.78, 0.80, 0.82]
    header = ("alpha_deg,m_inf,converged,kutta_mismatch,mach_max,"
              "x_shock_upper,cl,n_limited")
    rows, results = [], {}
    for alpha in (0.0, 1.25):
        for m in ms:
            r = solve_transonic_lifting(mc, wc, m_inf=m, alpha_deg=alpha,
                                        max_gamma_evals=12, n_picard_eval=800)
            curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m)
            rep = shock_report(curve, m)
            forces = wall_force_coefficients(
                mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
                alpha_deg=alpha, s_ref=dz, m_inf=m,
            )
            results[(alpha, m)] = (r, rep, forces)
            rows.append((alpha, m, r["converged"],
                         f"{r['kutta_mismatch']:.2e}",
                         f"{np.sqrt(r['mach2_max']):.3f}",
                         f"{rep['upper']['x_shock']:.4f}",
                         f"{forces['cl']:.5f}", r["n_limited"]))
    write_csv(OUT, "g43_summary.csv", header, rows)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for alpha, mk, color in ((0.0, "o", S1_BLUE), (1.25, "s", S3_YELLOW)):
        xs = [results[(alpha, m)][1]["upper"]["x_shock"] for m in ms]
        cls = [results[(alpha, m)][2]["cl"] for m in ms]
        ax1.plot(ms, xs, mk + "-", color=color, label=f"alpha={alpha}")
        ax2.plot(ms, cls, mk + "-", color=color, label=f"alpha={alpha}")
    ax1.set_xlabel("M_inf"); ax1.set_ylabel("upper shock x/c"); ax1.legend()
    ax2.set_xlabel("M_inf"); ax2.set_ylabel("cl"); ax2.legend()
    ax1.set_title("V4.3 shock migration"); ax2.set_title("V4.3 lift trend")
    finish(fig, OUT, "g43_sweep_dashboard.png")

    n_conv = sum(r["converged"] for (r, _, _) in results.values())
    n_bad = sum(r["n_limited"] + r["n_floored"] for (r, _, _) in results.values())
    cl.add("G4.3", "all cases converged", f"{n_conv}/10", "10/10", n_conv == 10)
    cl.add("G4.3", "zero limited/floored across sweep", n_bad, "0", n_bad == 0)


def main():
    apply_style()
    cl = CheckList("P4 transonic artificial density (G4.1-G4.3)")
    mesh = read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)
    t0 = time.time()
    part1_noop(cl, mc, wc)
    part2_reach(cl, mc, wc)
    coarse_case = part3_transonic(cl, mc, wc)
    # Heavy mode = ~40 min. The results/ baseline is already committed; only
    # regenerate when a real solver/mesh/reference change moves the numbers
    # and you will commit the refresh (see the COST CAUTION in the docstring).
    if os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1":
        part4_refinement(cl, coarse_case)
        part5_sweep(cl, mc, wc)
    else:
        print("\n[4-5/5] medium G4.1 + G4.3 sweep skipped (~40 min); set "
              "PYFP3D_TRANSONIC_GATES=1 to regenerate. The committed\n"
              "        results/g41_cp_shock_{coarse,medium}.png, "
              "g43_sweep_dashboard.png and their *_summary CSVs\n"
              "        are the reference baseline embedded by demo_report.")
    print(f"\ntotal runtime {time.time() - t0:.0f} s")
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
