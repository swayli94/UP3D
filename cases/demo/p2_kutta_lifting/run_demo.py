"""
P2 demo -- wake cut, circulation and Kutta condition evidence
(gates G2.1-G2.5, all closed; NACA0012 alpha = 4 deg, incompressible).

What this shows, per docs/roadmap.md P2 and docs/demo_report.md:
  1. G2.1: the wake cut adds zero spurious residual -- phi = x with
     Gamma = 0 is preserved to machine zero on the cut mesh, including
     the folded wake-master rows.
  2. G2.2: a prescribed circulation is reproduced exactly -- the jump
     [phi] equals Gamma at every wake node pair to 1e-12.
  3. G2.3 = V3: the secant-accelerated Kutta loop converges in 2 updates
     (the plain relaxation map has slope b ~ 0.93, i.e. O(100) iterations),
     and the resulting cl matches the Hess-Smith panel reference within 2%.
  4. G2.4 = V6: lift computed two independent ways -- surface-pressure
     integration and Kutta-Joukowski from Gamma -- agrees to ~0.01%:
     the circulation the Kutta loop found IS the lift the pressure field
     carries; the wake machinery is physics, not tuning.
  5. Physical plausibility: the mid-span flow field leaves the trailing
     edge smoothly (Kutta enforced) and phi shows the branch cut across
     the wake with a constant jump.
  6. G2.5: the quasi-2D solution is spanwise-consistent -- p99 |w|/U_inf
     decays at clean 1st order under refinement, no wake stripe.

Standalone + self-checking:  python cases/demo/p2_kutta_lifting/run_demo.py
Outputs: cases/demo/p2_kutta_lifting/results/{*.png, summary.csv, checks.csv}
Exit code 0 iff every acceptance check passes. Runtime ~1-2 min.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cases.demo._common import (  # noqa: E402
    BASELINE, CRITICAL, DIV_BLUE_RED, GRID, INK, INK_2, MESH_DIR, MUTED,
    REFERENCE_DIR, REPO_ROOT, S1_BLUE, S2_AQUA, S3_YELLOW, S5_VIOLET,
    SEQ_BLUE, SURFACE, CheckList, apply_style, finish, write_csv,
)

import matplotlib.pyplot as plt  # noqa: E402

from pyfp3d.constraints.wake import WakeConstraint  # noqa: E402
from pyfp3d.kernels.residual import assemble_stiffness_matrix  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.post.section_cut import section_cut, wall_cp_curve  # noqa: E402
from pyfp3d.post.surface import (  # noqa: E402
    nodal_gradient_recovery, sectional_cl_from_gamma, wall_force_coefficients,
)
from pyfp3d.solve.picard import solve_laplace_lifting  # noqa: E402
from tests.mesh_utils import element_gradients_all  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
NACA_DIR = MESH_DIR / "naca0012_2.5d"
PANEL_DIR = REFERENCE_DIR / "naca0012_incompressible"
ALPHA_DEG = 4.0


def load_panel_cl(alpha_deg=ALPHA_DEG):
    with open(PANEL_DIR / "cl_reference.csv") as f:
        for row in csv.DictReader(f):
            if abs(float(row["alpha_deg"]) - alpha_deg) < 1e-9:
                return float(row["cl"])
    raise LookupError(f"alpha={alpha_deg} not in panel reference")


def load_panel_cp():
    x, cp, surf = [], [], []
    with open(PANEL_DIR / "cp_alpha4.csv") as f:
        for row in csv.DictReader(f):
            x.append(float(row["x_c"]))
            cp.append(float(row["cp"]))
            surf.append(row["surface"])
    x, cp = np.asarray(x), np.asarray(cp)
    return x, cp, np.asarray(surf) == "upper"


def lifting_case(level):
    mesh = read_mesh(NACA_DIR / f"{level}.msh")
    mesh_cut, wc = cut_wake(mesh)
    result = solve_laplace_lifting(mesh_cut, wc, alpha_deg=ALPHA_DEG)
    dz = float(np.ptp(mesh_cut.nodes[:, 2]))
    forces = wall_force_coefficients(
        mesh_cut.nodes, mesh_cut.elements, mesh_cut.boundary_faces["wall"],
        result["phi"], alpha_deg=ALPHA_DEG, s_ref=1.0 * dz)
    return {"level": level, "mesh_cut": mesh_cut, "wc": wc, "result": result,
            "forces": forces, "dz": dz,
            "n_tets": len(mesh_cut.elements)}


# ---------------------------------------------------------------------------
# 1+2. G2.1 freestream on the cut mesh / G2.2 prescribed jump
# ---------------------------------------------------------------------------
def demo_cut_exactness(cases, checks):
    print("[1/6] G2.1 freestream on cut mesh + G2.2 prescribed jump")
    rows, r_free, r_master = [], {}, {}
    for level in ("coarse", "medium"):
        mc, wc = cases[level]["mesh_cut"], cases[level]["wc"]
        A = assemble_stiffness_matrix(mc.nodes, mc.elements)
        con = WakeConstraint(A, wc)
        R = con.A_reduced @ mc.nodes[: con.n_reduced, 0]
        to_master = np.arange(len(mc.nodes))
        to_master[wc.slave_nodes] = wc.master_nodes
        check = np.ones(con.n_reduced, dtype=bool)
        for tag in ("wall", "farfield"):
            check[np.unique(to_master[np.unique(mc.boundary_faces[tag])])] = False
        masters_int = wc.master_nodes[check[wc.master_nodes]]
        r_free[level] = float(np.max(np.abs(R[check])))
        r_master[level] = float(np.max(np.abs(R[masters_int])))
        rows.append((level, f"{r_free[level]:.3e}", f"{r_master[level]:.3e}"))

    # G2.2: prescribed Gamma on the coarse mesh
    gamma = 0.3
    mc, wc = cases["coarse"]["mesh_cut"], cases["coarse"]["wc"]
    r = solve_laplace_lifting(mc, wc, alpha_deg=0.0, gamma_fixed=gamma)
    x_w = mc.nodes[wc.master_nodes, 0]
    order = np.argsort(x_w)
    phi_m = r["phi"][wc.master_nodes][order]
    phi_s = r["phi"][wc.slave_nodes][order]
    jump_err = np.abs(phi_s - phi_m - gamma)

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.0))
    ax = axes[0]
    labels, vals, colors = [], [], []
    for level in ("coarse", "medium"):
        labels += [f"{level}: free dofs", f"{level}: wake-master rows"]
        vals += [r_free[level], r_master[level]]
        colors += [S1_BLUE, S2_AQUA]
    y = np.arange(len(vals))[::-1]
    vals = np.array(vals)
    ax.hlines(y, vals.min() * 0.3, vals, color=GRID, linewidth=1.0, zorder=2)
    for yi, v, c in zip(y, vals, colors):
        ax.plot(v, yi, "o", color=c, markersize=9, markeredgecolor=SURFACE,
                markeredgewidth=2, zorder=3)
        ax.text(v * 2.2, yi, f"{v:.1e}", va="center", color=INK_2, fontsize=9)
    ax.axvline(1e-12, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax.set_ylim(-0.75, len(vals) - 0.4)
    ax.text(1.35e-12, -0.55, "gate: 1e-12", color=CRITICAL, fontsize=9, va="center")
    ax.set_xscale("log")
    ax.set_xlim(vals.min() * 0.3, 3e-11)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=INK_2)
    ax.set_xlabel("max |R(phi = x, Gamma = 0)|")
    ax.set_title("G2.1: the cut adds zero residual")
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    ax.plot(x_w[order], phi_s, "-", color=S1_BLUE, linewidth=1.8,
            label="phi+ (upper side)")
    ax.plot(x_w[order], phi_m, "-", color=S3_YELLOW, linewidth=1.8,
            label="phi- (lower side)")
    ax.set_xlabel("x along wake")
    ax.set_ylabel("phi")
    ax.set_title(f"G2.2: wake pair potentials, Gamma = {gamma}")
    ax.legend()

    ax = axes[2]
    ax.semilogy(x_w[order], np.maximum(jump_err, 1e-17), ".", color=S5_VIOLET,
                markersize=5)
    ax.axhline(1e-12, color=CRITICAL, linewidth=1.2, linestyle="--")
    ax.text(x_w.max(), 6.5e-13, "gate: 1e-12", color=CRITICAL, fontsize=9,
            ha="right", va="top")
    ax.set_xlabel("x along wake")
    ax.set_ylabel("| [phi] - Gamma |")
    ax.set_title("jump error at every wake node pair")
    fig.suptitle("Wake-cut machinery is exact: no spurious residual, "
                 "exact prescribed circulation", fontsize=12.5,
                 fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    finish(fig, OUT, "g21_g22_cut_exactness.png")
    write_csv(OUT, "g21_residuals.csv", "level,max_R_free,max_R_wake_master", rows)

    worst_free = max(r_free.values())
    worst_master = max(r_master.values())
    checks.add("G2.1", "max free-dof residual, phi=x on cut", f"{worst_free:.1e}",
               "< 1e-12", bool(worst_free < 1e-12))
    checks.add("G2.1", "max wake-master row residual", f"{worst_master:.1e}",
               "< 1e-13", bool(worst_master < 1e-13))
    checks.add("G2.2", "max |[phi] - Gamma|, prescribed jump",
               f"{jump_err.max():.1e}", "< 1e-12", bool(jump_err.max() < 1e-12))


# ---------------------------------------------------------------------------
# 3. G2.3 Kutta convergence + Cp vs panel
# ---------------------------------------------------------------------------
def demo_kutta_convergence(cases, checks):
    print("[2/6] Kutta loop convergence")
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    rows = []
    for level, color, marker, ms in (("coarse", S1_BLUE, "o", 11),
                                     ("medium", S2_AQUA, "s", 6.5)):
        res = cases[level]["result"]
        hist = np.asarray(res["gamma_history"], dtype=float).ravel()
        gamma_star = res["gamma"][0] if np.ndim(res["gamma"]) else res["gamma"]
        err = np.abs(hist - gamma_star)
        k = np.arange(len(hist))
        ax.semilogy(k, np.maximum(err, 1e-16), marker + "-", color=color,
                    markersize=ms, markeredgecolor=SURFACE, markeredgewidth=1.6,
                    label=f"{level}: {res['n_kutta_updates']} updates "
                          f"(Gamma* = {gamma_star:.4f})")
        rows += [(level, i, f"{g:.8f}", f"{e:.3e}") for i, (g, e) in
                 enumerate(zip(hist, err))]
    # what plain relaxation would look like at the measured map slope b=0.93
    k_ref = np.arange(0, 8)
    err0 = np.abs(np.asarray(cases["medium"]["result"]["gamma_history"],
                             dtype=float).ravel()[0]
                  - cases["medium"]["result"]["gamma"][0])
    ax.semilogy(k_ref, err0 * 0.93 ** k_ref, "--", color=MUTED, linewidth=1.4)
    ax.text(k_ref[-1], err0 * 0.93 ** k_ref[-1] * 1.3,
            "plain relaxation, b = 0.93\n(O(100) updates to converge)",
            color=MUTED, fontsize=9, ha="right")
    ax.set_xlabel("Kutta update")
    ax.set_ylabel("|Gamma_k - Gamma*|")
    ax.set_title("G2.3: secant-accelerated Kutta loop converges in 2 updates")
    ax.legend(loc="center right")
    finish(fig, OUT, "g23_kutta_convergence.png")
    write_csv(OUT, "g23_kutta_history.csv", "level,update,gamma,abs_err", rows)

    med = cases["medium"]["result"]
    checks.add("G2.3", "Kutta converged (medium)", bool(med["kutta_converged"]),
               "True", bool(med["kutta_converged"]))
    checks.add("G2.3", "Kutta updates (medium)", med["n_kutta_updates"],
               "< 20 (measured 2)", bool(med["n_kutta_updates"] < 20))


def demo_cp_vs_panel(cases, checks):
    print("[3/6] G2.3 Cp + cl vs panel reference")
    cl_ref = load_panel_cl()
    case = cases["medium"]
    curve = wall_cp_curve(case["mesh_cut"], case["result"]["phi"],
                          z=0.5 * case["dz"])
    x_ref, cp_ref, up = load_panel_cp()
    cl = case["forces"]["cl"]
    rel = 100 * (cl / cl_ref - 1)

    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.plot(x_ref[up], cp_ref[up], "-", color=INK_2, linewidth=1.4,
            label="Hess-Smith panel, upper (N=800)")
    ax.plot(x_ref[~up], cp_ref[~up], "--", color=INK_2, linewidth=1.4,
            label="Hess-Smith panel, lower")
    ax.plot(curve["x_upper"], curve["cp_upper"], ".", color=S1_BLUE,
            markersize=4.5, label="pyFP3D upper (mid-span cut)")
    ax.plot(curve["x_lower"], curve["cp_lower"], ".", color=S3_YELLOW,
            markersize=4.5, label="pyFP3D lower")
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel("Cp")
    ax.text(0.97, 0.05,
            f"cl = {cl:.5f} vs panel {cl_ref:.5f}  ({rel:+.2f}%, gate < 2%)",
            transform=ax.transAxes, ha="right", color=INK_2, fontsize=10)
    ax.set_title(f"G2.3 = V3: NACA0012 alpha = {ALPHA_DEG:.0f} deg, medium mesh "
                 f"({case['n_tets']/1e3:.0f}k tets)")
    ax.legend(loc="upper right")
    finish(fig, OUT, "g23_cp_vs_panel.png")
    write_csv(OUT, "g23_cp_midspan.csv", "x_c,cp,surface",
              [(f"{x:.6f}", f"{c:.6f}", "upper") for x, c in
               zip(curve["x_upper"], curve["cp_upper"])]
              + [(f"{x:.6f}", f"{c:.6f}", "lower") for x, c in
                 zip(curve["x_lower"], curve["cp_lower"])])

    checks.add("G2.3", "cl vs panel reference (medium)", f"{rel:+.2f}%",
               "|err| < 2%", bool(abs(rel) < 2.0))
    # physics sanity mirrors the gate test
    checks.add("G2.3", "suction peak on upper surface",
               f"{curve['cp_upper'].min():.2f}", "< -0.7",
               bool(curve["cp_upper"].min() < -0.7))
    return cl_ref


# ---------------------------------------------------------------------------
# 4. G2.4 lift cross-check
# ---------------------------------------------------------------------------
def demo_cl_crosscheck(cases, cl_ref, checks):
    print("[4/6] G2.4 lift cross-check")
    rows = []
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = np.arange(2)
    width = 0.3
    for i, (route, color) in enumerate((("pressure integration", S1_BLUE),
                                        ("Kutta-Joukowski 2*Gamma/(U c)", S2_AQUA))):
        vals = []
        for level in ("coarse", "medium"):
            case = cases[level]
            cl_p = case["forces"]["cl"]
            cl_g = float(sectional_cl_from_gamma(case["result"]["gamma"])[0])
            vals.append(cl_p if i == 0 else cl_g)
        bars = ax.bar(x + (i - 0.5) * width, vals, width * 0.92, color=color,
                      label=route)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.002, f"{v:.4f}",
                    ha="center", color=INK_2, fontsize=9.5)
    ax.axhline(cl_ref, color=INK_2, linewidth=1.4, linestyle="--")
    ax.text(-0.42, cl_ref + 0.003, f"panel reference {cl_ref:.4f}",
            color=INK_2, fontsize=9.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{lv} ({cases[lv]['n_tets']/1e3:.0f}k tets)"
                        for lv in ("coarse", "medium")], color=INK_2)
    ax.set_ylabel("sectional cl")
    ax.set_ylim(0.44, 0.50)
    ax.set_title("G2.4 = V6: two independent lift routes agree "
                 "-- circulation IS the pressure lift")
    ax.legend(loc="upper left")
    ax.grid(axis="x", visible=False)
    finish(fig, OUT, "g24_cl_crosscheck.png")

    for level in ("coarse", "medium"):
        case = cases[level]
        cl_p = case["forces"]["cl"]
        cl_g = float(sectional_cl_from_gamma(case["result"]["gamma"])[0])
        rows.append((level, f"{cl_p:.6f}", f"{cl_g:.6f}", f"{cl_ref:.6f}",
                     f"{100 * abs(cl_g - cl_p) / cl_p:.4f}"))
    write_csv(OUT, "g24_cl_crosscheck.csv",
              "level,cl_pressure,cl_gamma,cl_panel,gamma_vs_pressure_pct", rows)

    case = cases["medium"]
    cl_p = case["forces"]["cl"]
    cl_g = float(sectional_cl_from_gamma(case["result"]["gamma"])[0])
    rel = abs(cl_g - cl_p) / abs(cl_p)
    checks.add("G2.4", "cl(Gamma) vs cl(pressure), medium", f"{100*rel:.3f}%",
               "< 1%", bool(rel < 0.01))


# ---------------------------------------------------------------------------
# 5. lifting flow field at mid-span
# ---------------------------------------------------------------------------
def demo_flowfield(cases, checks):
    print("[5/6] lifting flow field (mid-span slice)")
    from pyfp3d.constraints.dirichlet import freestream_phi

    case = cases["medium"]
    mc = case["mesh_cut"]
    phi = case["result"]["phi"]
    grad = nodal_gradient_recovery(mc.nodes, mc.elements, phi)
    q = np.linalg.norm(grad, axis=1)
    phi_pert = phi - freestream_phi(mc.nodes, alpha_deg=ALPHA_DEG)
    sec = section_cut(mc, {"q": q, "phi_pert": phi_pert}, z=0.5 * case["dz"])

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0))
    ax = axes[0]
    t = ax.tricontourf(sec.points2d[:, 0], sec.points2d[:, 1], sec.triangles,
                       sec.fields["q"], levels=np.linspace(0.0, 1.6, 25),
                       cmap=SEQ_BLUE, extend="max")
    ax.set_xlim(-0.35, 1.6)
    ax.set_ylim(-0.7, 0.7)
    ax.set_aspect("equal")
    ax.set_title("speed: smooth flow off the TE (Kutta)")
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    ax.grid(False)
    cb = fig.colorbar(t, ax=ax, shrink=0.85, pad=0.02, label="|grad phi| / U_inf")
    cb.ax.tick_params(labelsize=9, colors=MUTED)
    cb.outline.set_edgecolor(BASELINE)
    ax.annotate("stagnation", xy=(0.01, -0.02), xytext=(-0.3, -0.4),
                color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-|>", color=INK_2, linewidth=1.2))

    ax = axes[1]
    gamma_val = float(case["result"]["gamma"][0])
    lim = 0.75 * gamma_val
    t = ax.tripcolor(sec.points2d[:, 0], sec.points2d[:, 1], sec.triangles,
                     sec.fields["phi_pert"], cmap=DIV_BLUE_RED,
                     vmin=-lim, vmax=lim, shading="gouraud")
    ax.set_xlim(-1.5, 5.0)
    ax.set_ylim(-2.2, 2.2)
    ax.set_aspect("equal")
    ax.set_title("phi - phi_freestream: branch cut carries "
                 f"[phi] = Gamma = {gamma_val:.3f}")
    ax.set_xlabel("x/c")
    ax.grid(False)
    cb = fig.colorbar(t, ax=ax, shrink=0.85, pad=0.02,
                      label="perturbation potential")
    cb.ax.tick_params(labelsize=9, colors=MUTED)
    cb.outline.set_edgecolor(BASELINE)
    fig.suptitle("Physical plausibility: lifting NACA0012 mid-span field "
                 "(medium mesh, alpha = 4 deg)",
                 fontsize=12.5, fontweight="semibold", color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    finish(fig, OUT, "lifting_flowfield.png")

    # Kutta physics check: the flow leaves the TE smoothly -- speed one cell
    # off the TE stays O(1), no point-vortex-like suction spike.
    te_zone = (np.abs(sec.points2d[:, 0] - 1.0) < 0.06) & \
              (np.abs(sec.points2d[:, 1]) < 0.03)
    q_te = sec.fields["q"][te_zone].max()
    checks.add("Kutta-sanity", "max speed in TE neighborhood", f"{q_te:.2f}",
               "< 1.6 (no TE singularity)", bool(q_te < 1.6))


# ---------------------------------------------------------------------------
# 6. G2.5 spanwise consistency
# ---------------------------------------------------------------------------
def demo_spanwise(cases, checks):
    print("[6/6] G2.5 spanwise consistency")
    stats, rows = {}, []
    for level in ("coarse", "medium"):
        mc = cases[level]["mesh_cut"]
        gz_free = element_gradients_all(mc.nodes, mc.elements,
                                        mc.nodes[:, 0].copy())[:, 2]
        w = np.abs(element_gradients_all(mc.nodes, mc.elements,
                                         cases[level]["result"]["phi"])[:, 2])
        stats[level] = {
            "free": float(np.max(np.abs(gz_free))),
            "p99": float(np.percentile(w, 99)),
            "rms": float(np.sqrt(np.mean(w**2))),
            "max": float(w.max()),
        }
        rows.append((level, f"{stats[level]['free']:.2e}",
                     f"{stats[level]['p99']:.3e}", f"{stats[level]['rms']:.3e}",
                     f"{stats[level]['max']:.3e}"))

    h = np.array([1.0, 0.5])  # medium halves the coarse spacing
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for key, color, label in (("p99", S1_BLUE, "p99 |w|/U_inf (gated)"),
                              ("rms", S2_AQUA, "RMS |w|/U_inf"),
                              ("max", MUTED, "max |w|/U_inf (LE peak, not gated)")):
        vals = np.array([stats["coarse"][key], stats["medium"][key]])
        ax.loglog(h, vals, "o-", color=color, markeredgecolor=SURFACE,
                  markeredgewidth=2, label=label)
    ref = stats["coarse"]["p99"] * h
    ax.loglog(h, ref, "--", color=MUTED, linewidth=1.2)
    ax.text(h[0] * 0.99, ref[0] * 0.78, "O(h) reference", color=MUTED,
            fontsize=9, ha="right")
    ratio = stats["coarse"]["p99"] / stats["medium"]["p99"]
    ax.text(0.5, 0.03, f"p99 decay ratio = {ratio:.2f} (gate >= 1.8, "
                       "clean 1st order)",
            transform=ax.transAxes, color=INK_2, fontsize=10, ha="center")
    ax.set_xlabel("relative mesh spacing h (coarse = 1)")
    ax.set_ylabel("spanwise velocity |w| / U_inf")
    ax.set_xlim(0.42, 1.18)
    ax.set_xticks([1.0, 0.5], labels=["coarse", "medium"])
    ax.minorticks_off()
    ax.set_title("G2.5: spanwise noise of the lifting solution decays "
                 "at 1st order")
    ax.legend(loc="center left", fontsize=9)
    ax.grid(True, which="major")
    finish(fig, OUT, "g25_spanwise_decay.png")
    write_csv(OUT, "g25_spanwise.csv",
              "level,max_w_freestream_interp,p99_w_lifting,rms_w_lifting,max_w_lifting",
              rows)

    worst_free = max(stats[lv]["free"] for lv in stats)
    checks.add("G2.5a", "spanwise gradient of interpolated phi=x",
               f"{worst_free:.1e}", "< 1e-12", bool(worst_free < 1e-12))
    checks.add("G2.5b", "p99 |w| decay ratio coarse/medium", f"{ratio:.2f}",
               ">= 1.8 (1st order at h-ratio 2)", bool(ratio >= 1.8))


def main():
    apply_style()
    t0 = time.time()
    checks = CheckList("P2 wake cut + Kutta (G2.1-G2.5, closed)")

    print("solving lifting cases (coarse + medium)...")
    cases = {level: lifting_case(level) for level in ("coarse", "medium")}

    demo_cut_exactness(cases, checks)
    demo_kutta_convergence(cases, checks)
    cl_ref = demo_cp_vs_panel(cases, checks)
    demo_cl_crosscheck(cases, cl_ref, checks)
    demo_flowfield(cases, checks)
    demo_spanwise(cases, checks)

    write_csv(OUT, "summary.csv", "metric,value",
              [(c["gate"] + " " + c["name"].replace(",", ";"),
                str(c["value"]) + " [" + c["status"] + "]")
               for c in checks.checks] + [("runtime_seconds", f"{time.time()-t0:.1f}")])
    code = checks.report(OUT)
    print(f"done in {time.time() - t0:.1f}s -> {OUT.relative_to(REPO_ROOT)}/")
    sys.exit(code)


if __name__ == "__main__":
    main()
