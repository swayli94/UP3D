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
     no limited/floored cells. (The medium-mesh gate run and the G4.3
     sweep run under PYFP3D_TRANSONIC_GATES=1 in pytest; their artifacts
     land in artifacts/G4.1, G4.3.)

Standalone + self-checking:  python cases/demo/p4_transonic/run_demo.py
Outputs: cases/demo/p4_transonic/results/{*.png, *.csv, checks.csv}
Exit code 0 iff every acceptance check passes. Runtime ~4 min.
"""

import csv
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
from pyfp3d.post.section_cut import wall_cp_curve  # noqa: E402
from pyfp3d.post.shock import shock_report  # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients  # noqa: E402
from pyfp3d.solve.continuation import solve_transonic_lifting  # noqa: E402
from pyfp3d.solve.picard import solve_subsonic_lifting  # noqa: E402

OUT = Path(__file__).parent / "results"
M_INF, ALPHA = 0.80, 1.25


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


def main():
    apply_style()
    cl = CheckList("P4 transonic artificial density (G4.1-G4.3)")
    mesh = read_mesh(MESH_DIR / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)
    t0 = time.time()
    part1_noop(cl, mc, wc)
    part2_reach(cl, mc, wc)
    part3_transonic(cl, mc, wc)
    print(f"\ntotal runtime {time.time() - t0:.0f} s")
    return cl.report(OUT)


if __name__ == "__main__":
    sys.exit(main())
