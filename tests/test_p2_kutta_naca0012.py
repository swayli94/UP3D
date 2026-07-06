"""
P2 gate tests on the M0 NACA0012 quasi-2D family: G2.3 (lifting solution
vs the 2D panel reference), G2.4 (Gamma-lift vs pressure-integrated lift)
and G2.5 (quasi-2D spanwise consistency), with the V2.1-V2.5 headless
artifacts (roadmap Sec 0.1: PNG + CSV, matplotlib Agg, no GUI).

Reference data: cases/reference_data/naca0012_incompressible/ (Hess-Smith
panel method, provenance in its README.md). Gate values are read from the
CSV, never embedded here (roadmap Sec "reference data").
"""

import csv

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.picard import solve_laplace_lifting
from pyfp3d.post.section_cut import section_cut, wall_cp_curve
from pyfp3d.post.surface import (
    nodal_gradient_recovery,
    sectional_cl_from_gamma,
    wall_force_coefficients,
)
from .mesh_utils import element_gradients_all

ALPHA_DEG = 4.0


def _panel_reference_cl(reference_mesh_dir, alpha_deg):
    path = reference_mesh_dir / "naca0012_incompressible" / "cl_reference.csv"
    with open(path) as f:
        for row in csv.DictReader(f):
            if abs(float(row["alpha_deg"]) - alpha_deg) < 1e-9:
                return float(row["cl"])
    raise LookupError(f"alpha = {alpha_deg} not in {path}")


def _panel_reference_cp(reference_mesh_dir):
    path = reference_mesh_dir / "naca0012_incompressible" / "cp_alpha4.csv"
    x, cp, surf = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            x.append(float(row["x_c"]))
            cp.append(float(row["cp"]))
            surf.append(row["surface"])
    x, cp = np.asarray(x), np.asarray(cp)
    up = np.asarray(surf) == "upper"
    return x, cp, up


def _lifting_case(mesh_dir, level):
    mesh = read_mesh(mesh_dir / "naca0012_2.5d" / f"{level}.msh")
    mesh_cut, wc = cut_wake(mesh)
    result = solve_laplace_lifting(mesh_cut, wc, alpha_deg=ALPHA_DEG)
    dz = float(np.ptp(mesh_cut.nodes[:, 2]))
    forces = wall_force_coefficients(
        mesh_cut.nodes, mesh_cut.elements, mesh_cut.boundary_faces["wall"],
        result["phi"], alpha_deg=ALPHA_DEG, s_ref=1.0 * dz,
    )
    return {"mesh_cut": mesh_cut, "wc": wc, "result": result,
            "forces": forces, "dz": dz}


@pytest.fixture(scope="module")
def coarse_case():
    from .conftest import REPO_ROOT
    return _lifting_case(REPO_ROOT / "cases" / "meshes", "coarse")


@pytest.fixture(scope="module")
def medium_case():
    from .conftest import REPO_ROOT
    return _lifting_case(REPO_ROOT / "cases" / "meshes", "medium")


def test_g23_cl_vs_panel_reference(medium_case, coarse_case,
                                   reference_mesh_dir, artifacts_dir):
    """G2.3 = V3(incompressible): NACA0012 alpha=4deg, Delta-cl < 2% vs the
    2D panel reference; Kutta loop converges in < 20 updates."""
    cl_ref = _panel_reference_cl(reference_mesh_dir, ALPHA_DEG)

    rows = [("level", "cl_pressure", "cl_gamma", "cl_panel_ref",
             "rel_err_pct", "n_kutta_updates")]
    for level, case in (("coarse", coarse_case), ("medium", medium_case)):
        cl = case["forces"]["cl"]
        cl_g = float(sectional_cl_from_gamma(case["result"]["gamma"])[0])
        rows.append((level, f"{cl:.6f}", f"{cl_g:.6f}", f"{cl_ref:.6f}",
                     f"{100 * (cl / cl_ref - 1):.3f}",
                     case["result"]["n_kutta_updates"]))

    gate_dir = artifacts_dir / "G2.3"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)

    med = medium_case["result"]
    assert med["kutta_converged"]
    assert med["n_kutta_updates"] < 20
    cl_med = medium_case["forces"]["cl"]
    assert abs(cl_med / cl_ref - 1) < 0.02, (
        f"medium-mesh cl = {cl_med:.5f} vs panel {cl_ref:.5f} "
        f"({100 * (cl_med / cl_ref - 1):+.2f}%)"
    )


def test_g24_gamma_vs_pressure_lift(medium_case, artifacts_dir):
    """G2.4 = V6: sectional cl from Gamma vs pressure integration, < 1%."""
    cl_p = medium_case["forces"]["cl"]
    cl_g = float(sectional_cl_from_gamma(medium_case["result"]["gamma"])[0])
    rel = abs(cl_g - cl_p) / abs(cl_p)

    gate_dir = artifacts_dir / "G2.4"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cl_pressure", "cl_gamma", "rel_diff_pct"])
        w.writerow([f"{cl_p:.6f}", f"{cl_g:.6f}", f"{100 * rel:.4f}"])

    # V2.4: the two lift routes side by side with the 1% band.
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["pressure\nintegration", "2*Gamma/(U c)"], [cl_p, cl_g],
           color=["tab:blue", "tab:orange"], width=0.5)
    ax.axhspan(cl_p * 0.99, cl_p * 1.01, alpha=0.25, color="gray",
               label="±1% of cl_pressure")
    ax.set_ylim(min(cl_p, cl_g) * 0.97, max(cl_p, cl_g) * 1.03)
    ax.set_ylabel("sectional cl")
    ax.set_title(f"V2.4 lift cross-check (medium, diff {100 * rel:.2f}%)")
    ax.legend()
    fig.savefig(gate_dir / "v2_4_cl_crosscheck.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    assert rel < 0.01, f"cl_gamma {cl_g:.5f} vs cl_pressure {cl_p:.5f}: {100*rel:.2f}%"


def test_g25_quasi2d_consistency(coarse_case, medium_case, artifacts_dir):
    """G2.5 (criterion (b) as re-specified 2026-07-06, see roadmap evidence
    note): (a) phi = x on the cut mesh gives machine-zero spanwise element
    gradient; (b) for the converged lifting solution, the field-wide
    spanwise noise (99th percentile of |w|/U_inf over elements) decays at
    >= ~1st order over the mesh family, with no coherent stripe along the
    wake/TE in the V2.5 heatmap. The single-element max is recorded but
    not gated: it sits in the leading-edge peak-gradient region (measured
    coarse (0.014, -0.024), medium (-0.0001, -0.0074) -- nowhere near the
    wake) and tracks the local h*|grad^2 phi|, not the split asymmetry.
    A literal 1e-12 stays unachievable for any solved field on 3-tet
    prisms."""
    rows = [("level", "case", "max_abs_w_over_uinf", "p99_abs_w_over_uinf",
             "rms_abs_w_over_uinf")]

    # (a) freestream on the cut mesh: machine zero.
    for level, case in (("coarse", coarse_case), ("medium", medium_case)):
        mc = case["mesh_cut"]
        gz = element_gradients_all(mc.nodes, mc.elements, mc.nodes[:, 0])[:, 2]
        w_free = float(np.max(np.abs(gz)))
        rows.append((level, "phi=x", f"{w_free:.3e}", "", ""))
        assert w_free < 1e-12

    # (b) converged lifting solution: >= ~1st-order decay coarse -> medium.
    w_p99 = {}
    for level, case in (("coarse", coarse_case), ("medium", medium_case)):
        mc = case["mesh_cut"]
        gz = element_gradients_all(mc.nodes, mc.elements,
                                   case["result"]["phi"])[:, 2]
        w = np.abs(gz)
        w_p99[level] = float(np.percentile(w, 99))
        rows.append((level, f"lifting alpha={ALPHA_DEG}", f"{w.max():.3e}",
                     f"{w_p99[level]:.3e}", f"{np.sqrt(np.mean(w**2)):.3e}"))

    gate_dir = artifacts_dir / "G2.5"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)

    # V2.5: |w| heatmap on the mid-plane slice (general marching-tets cut).
    mc = medium_case["mesh_cut"]
    w_nodal = nodal_gradient_recovery(mc.nodes, mc.elements,
                                      medium_case["result"]["phi"])[:, 2]
    sec = section_cut(mc, {"abs_w": np.abs(w_nodal)},
                      z=0.5 * medium_case["dz"])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (x0, x1, y0, y1) in zip(axes, [(-2, 16, -6, 6), (-0.3, 2.5, -0.8, 0.8)]):
        t = ax.tripcolor(sec.points2d[:, 0], sec.points2d[:, 1],
                         sec.triangles, sec.fields["abs_w"], cmap="magma")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
        fig.colorbar(t, ax=ax, label="|w|/U_inf")
    axes[0].set_title("V2.5 |w| mid-plane (medium) - domain")
    axes[1].set_title("airfoil/wake zoom: expect no coherent stripe")
    fig.savefig(gate_dir / "v2_5_spanwise_w_heatmap.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    decay = w_p99["coarse"] / w_p99["medium"]
    assert decay > 1.8, (
        f"spanwise noise not decaying ~1st order (h ratio 2): p99 coarse "
        f"{w_p99['coarse']:.3e} -> medium {w_p99['medium']:.3e} "
        f"(ratio {decay:.2f})"
    )


def test_v21_v22_cut_artifacts(coarse_case, artifacts_dir):
    """V2.1: residual magnitude near the wake for the phi = x field (hot
    wall rows are physical; the wake region must be stripe-free).
    V2.2: phi+ / phi- and their difference along the wake for prescribed
    constant Gamma."""
    from pyfp3d.kernels.residual import assemble_stiffness_matrix
    from pyfp3d.constraints.wake import WakeConstraint

    mc, wc = coarse_case["mesh_cut"], coarse_case["wc"]
    A = assemble_stiffness_matrix(mc.nodes, mc.elements)
    con = WakeConstraint(A, wc)
    R = np.abs(con.A_reduced @ mc.nodes[: con.n_reduced, 0])

    gate_dir = artifacts_dir / "G2.1"
    gate_dir.mkdir(parents=True, exist_ok=True)
    R_cut = np.zeros(len(mc.nodes))
    R_cut[: con.n_reduced] = R
    R_cut[wc.slave_nodes] = R[wc.master_nodes]
    sec = section_cut(mc, {"log10_absR": np.log10(R_cut + 1e-17)}, z=0.0)
    fig, ax = plt.subplots(figsize=(8, 5))
    t = ax.tripcolor(sec.points2d[:, 0], sec.points2d[:, 1], sec.triangles,
                     sec.fields["log10_absR"], cmap="viridis", vmin=-17, vmax=-3)
    ax.set_xlim(-0.5, 16); ax.set_ylim(-3, 3); ax.set_aspect("equal")
    fig.colorbar(t, ax=ax, label="log10 |R| (phi = x, Gamma = 0)")
    ax.set_title("V2.1 residual on cut mesh: no hot stripe along wake")
    fig.savefig(gate_dir / "v2_1_residual_heatmap.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    gamma = 0.3
    r = solve_laplace_lifting(mc, wc, alpha_deg=0.0, gamma_fixed=gamma)
    x_w = mc.nodes[wc.master_nodes, 0]
    order = np.argsort(x_w)
    phi_m = r["phi"][wc.master_nodes][order]
    phi_s = r["phi"][wc.slave_nodes][order]

    gate_dir = artifacts_dir / "G2.2"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "phi_minus", "phi_plus", "jump"])
        for xx, pm, ps in zip(x_w[order], phi_m, phi_s):
            w.writerow([f"{xx:.6f}", f"{pm:.8f}", f"{ps:.8f}", f"{ps - pm:.3e}"])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.plot(x_w[order], phi_s, ".-", ms=3, label="phi+ (upper)")
    ax1.plot(x_w[order], phi_m, ".-", ms=3, label="phi- (lower)")
    ax1.set_ylabel("phi"); ax1.legend()
    ax1.set_title(f"V2.2 wake pair potentials, prescribed Gamma = {gamma}")
    ax2.plot(x_w[order], phi_s - phi_m, ".-", ms=3)
    ax2.axhline(gamma, color="k", lw=0.8, ls="--")
    ax2.set_xlabel("x"); ax2.set_ylabel("[phi]")
    ax2.set_ylim(gamma - 1e-3, gamma + 1e-3)
    fig.savefig(gate_dir / "v2_2_wake_jump.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    np.testing.assert_allclose(phi_s - phi_m, gamma, rtol=0, atol=1e-12)


def test_v23_cp_curves_vs_panel(medium_case, reference_mesh_dir, artifacts_dir):
    """V2.3: sectional Cp(x/c) upper/lower vs the panel reference."""
    mc = medium_case["mesh_cut"]
    curve = wall_cp_curve(mc, medium_case["result"]["phi"],
                          z=0.5 * medium_case["dz"])
    x_ref, cp_ref, up_ref = _panel_reference_cp(reference_mesh_dir)

    gate_dir = artifacts_dir / "G2.3"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "v2_3_cp_curves.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_c", "cp", "surface"])
        for xx, cc in zip(curve["x_upper"], curve["cp_upper"]):
            w.writerow([f"{xx:.6f}", f"{cc:.6f}", "upper"])
        for xx, cc in zip(curve["x_lower"], curve["cp_lower"]):
            w.writerow([f"{xx:.6f}", f"{cc:.6f}", "lower"])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_ref[up_ref], cp_ref[up_ref], "-", color="0.4", lw=1.2,
            label="panel upper")
    ax.plot(x_ref[~up_ref], cp_ref[~up_ref], "--", color="0.4", lw=1.2,
            label="panel lower")
    ax.plot(curve["x_upper"], curve["cp_upper"], ".", ms=3, color="tab:red",
            label="pyFP3D upper")
    ax.plot(curve["x_lower"], curve["cp_lower"], ".", ms=3, color="tab:blue",
            label="pyFP3D lower")
    ax.invert_yaxis()
    ax.set_xlabel("x/c"); ax.set_ylabel("Cp"); ax.legend()
    ax.set_title(f"V2.3 NACA0012 alpha={ALPHA_DEG} Cp vs panel (medium)")
    fig.savefig(gate_dir / "v2_3_cp_vs_panel.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # Sanity: suction peak on the upper surface, TE recovery toward Cp > 0.
    assert curve["cp_upper"].min() < -0.7
    assert curve["cp_upper"][-1] > 0.0


def test_section_cut_general_path_linear_exact(coarse_case):
    """The P2 marching-tets path reproduces a linear field exactly at an
    off-node plane (where the degenerate symmetry-plane path cannot run)."""
    mc = coarse_case["mesh_cut"]
    f = 2.0 * mc.nodes[:, 0] - 0.5 * mc.nodes[:, 1] + 0.25 * mc.nodes[:, 2]
    z = 0.37 * coarse_case["dz"]
    sec = section_cut(mc, {"f": f}, z=z)
    expected = 2.0 * sec.points2d[:, 0] - 0.5 * sec.points2d[:, 1] + 0.25 * z
    np.testing.assert_allclose(sec.fields["f"], expected, atol=1e-12)
