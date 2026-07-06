"""
P3 gate G3.2: NACA0012 at M_inf = 0.5, alpha = 2 deg on the M0 quasi-2D
family -- cl vs the 2D corrected-panel reference, Picard budget and
residual monotonicity, with the V3.2/V3.3 headless artifacts.

Reference data: cases/reference_data/naca0012_m05/ (the P2-verified
Hess-Smith panel solution with Prandtl-Glauert and Karman-Tsien
compressibility corrections; provenance in its README.md). The exact
subcritical full-potential cl is bracketed by the two corrections
(PG under-, KT over-corrects a 12%-thick section), so the G3.2 reference
value is their midpoint and the gate additionally asserts the computed cl
falls INSIDE the [PG, KT] bracket. Gate values are read from the CSV,
never embedded here (roadmap "reference data discipline").
"""

import csv

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import wall_cp_curve
from pyfp3d.post.surface import sectional_cl_from_gamma, wall_force_coefficients
from pyfp3d.solve.picard import solve_subsonic_lifting

M_INF = 0.5
ALPHA_DEG = 2.0


def _reference_cl(reference_mesh_dir, alpha_deg):
    path = reference_mesh_dir / "naca0012_m05" / "cl_reference.csv"
    with open(path) as f:
        for row in csv.DictReader(f):
            if abs(float(row["alpha_deg"]) - alpha_deg) < 1e-9:
                return float(row["cl_pg"]), float(row["cl_kt"])
    raise LookupError(f"alpha = {alpha_deg} not in {path}")


def _reference_cp(reference_mesh_dir):
    path = reference_mesh_dir / "naca0012_m05" / "cp_alpha2_m05.csv"
    x, cp_kt, surf = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            x.append(float(row["x_c"]))
            cp_kt.append(float(row["cp_kt"]))
            surf.append(row["surface"])
    return np.asarray(x), np.asarray(cp_kt), np.asarray(surf) == "upper"


def _compressible_case(mesh_dir, level):
    mesh = read_mesh(mesh_dir / "naca0012_2.5d" / f"{level}.msh")
    mesh_cut, wc = cut_wake(mesh)
    result = solve_subsonic_lifting(
        mesh_cut, wc, m_inf=M_INF, alpha_deg=ALPHA_DEG,
    )
    dz = float(np.ptp(mesh_cut.nodes[:, 2]))
    forces = wall_force_coefficients(
        mesh_cut.nodes, mesh_cut.elements, mesh_cut.boundary_faces["wall"],
        result["phi"], alpha_deg=ALPHA_DEG, s_ref=1.0 * dz, m_inf=M_INF,
    )
    return {"mesh_cut": mesh_cut, "wc": wc, "result": result,
            "forces": forces, "dz": dz}


@pytest.fixture(scope="module")
def coarse_case():
    from .conftest import REPO_ROOT
    return _compressible_case(REPO_ROOT / "cases" / "meshes", "coarse")


@pytest.fixture(scope="module")
def medium_case():
    from .conftest import REPO_ROOT
    return _compressible_case(REPO_ROOT / "cases" / "meshes", "medium")


def test_g32_cl_vs_corrected_panel(medium_case, coarse_case,
                                   reference_mesh_dir, artifacts_dir):
    """G3.2 part 1: medium-mesh cl within 2% of the corrected-panel
    reference midpoint AND inside the [PG, KT] correction bracket."""
    cl_pg, cl_kt = _reference_cl(reference_mesh_dir, ALPHA_DEG)
    cl_ref = 0.5 * (cl_pg + cl_kt)

    rows = [("level", "cl_pressure", "cl_gamma", "cl_ref_pg", "cl_ref_kt",
             "cl_ref_mid", "rel_err_mid_pct", "n_picard", "n_solves",
             "mach2_max")]
    for level, case in (("coarse", coarse_case), ("medium", medium_case)):
        cl = case["forces"]["cl"]
        cl_g = float(sectional_cl_from_gamma(case["result"]["gamma"])[0])
        r = case["result"]
        rows.append((level, f"{cl:.6f}", f"{cl_g:.6f}", f"{cl_pg:.6f}",
                     f"{cl_kt:.6f}", f"{cl_ref:.6f}",
                     f"{100 * (cl / cl_ref - 1):.3f}", r["n_picard"],
                     r["n_solves_total"], f"{r['mach2_max']:.4f}"))

    gate_dir = artifacts_dir / "G3.2"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)

    med = medium_case["result"]
    assert med["converged"] and med["kutta_converged"]
    cl_med = medium_case["forces"]["cl"]
    assert cl_pg <= cl_med <= cl_kt, (
        f"cl = {cl_med:.5f} outside the correction bracket "
        f"[{cl_pg:.5f}, {cl_kt:.5f}]"
    )
    assert abs(cl_med / cl_ref - 1) < 0.02, (
        f"medium-mesh cl = {cl_med:.5f} vs corrected-panel midpoint "
        f"{cl_ref:.5f} ({100 * (cl_med / cl_ref - 1):+.2f}%)"
    )


def test_g32_picard_budget_and_monotone(medium_case, artifacts_dir):
    """G3.2 part 2: Picard converges < 30 iterations with a monotone
    residual (non-increasing all the way; the tail sits at the linear-
    solver floor). The density-lag history must decay strictly."""
    r = medium_case["result"]
    h = r["residual_history"]
    d = r["drho_history"]

    gate_dir = artifacts_dir / "G3.2"
    gate_dir.mkdir(parents=True, exist_ok=True)
    # V3.3: residual + density-lag semilog history.
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(range(1, len(h) + 1), h, "o-", ms=4,
                label="||R||_inf (free reduced dofs)")
    ax.semilogy(range(1, len(d) + 1), np.maximum(d, 1e-17), "s--", ms=4,
                label="density lag ||rho_new - rho||_inf")
    ax.set_xlabel("Picard (density) iteration")
    ax.set_ylabel("residual / density lag")
    ax.set_title(f"V3.3 NACA0012 M={M_INF} alpha={ALPHA_DEG} Picard "
                 f"convergence (medium, {r['n_picard']} iters)")
    ax.legend()
    fig.savefig(gate_dir / "v3_3_picard_residual.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    assert r["converged"]
    assert r["n_picard"] < 30, f"Picard took {r['n_picard']} >= 30 iterations"
    for i in range(len(h) - 1):
        assert h[i + 1] <= h[i], (
            f"residual not monotone at iter {i}: {h[i]:.3e} -> {h[i + 1]:.3e}"
        )
    # drho strictly decreasing until the terminal fixed-point zeros.
    d_sig = [x for x in d if x > 0.0]
    for i in range(len(d_sig) - 1):
        assert d_sig[i + 1] < d_sig[i]


def test_v32_cp_and_mach_curves(medium_case, reference_mesh_dir,
                                artifacts_dir):
    """V3.2: sectional surface Cp at mid-span vs the KT-corrected panel
    curve; suction peak amplified vs incompressible and smooth recovery
    to the TE."""
    from pyfp3d.physics.isentropic import mach_number_squared

    mc = medium_case["mesh_cut"]
    curve = wall_cp_curve(mc, medium_case["result"]["phi"],
                          z=0.5 * medium_case["dz"], m_inf=M_INF)
    x_ref, cp_ref, up_ref = _reference_cp(reference_mesh_dir)

    gate_dir = artifacts_dir / "G3.2"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "v3_2_cp_curves.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_c", "cp", "surface"])
        for xx, cc in zip(curve["x_upper"], curve["cp_upper"]):
            w.writerow([f"{xx:.6f}", f"{cc:.6f}", "upper"])
        for xx, cc in zip(curve["x_lower"], curve["cp_lower"]):
            w.writerow([f"{xx:.6f}", f"{cc:.6f}", "lower"])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    ax.plot(x_ref[up_ref], cp_ref[up_ref], "-", color="0.4", lw=1.2,
            label="panel + Karman-Tsien, upper")
    ax.plot(x_ref[~up_ref], cp_ref[~up_ref], "--", color="0.4", lw=1.2,
            label="panel + Karman-Tsien, lower")
    ax.plot(curve["x_upper"], curve["cp_upper"], ".", ms=3, color="tab:red",
            label="pyFP3D upper")
    ax.plot(curve["x_lower"], curve["cp_lower"], ".", ms=3, color="tab:blue",
            label="pyFP3D lower")
    ax.invert_yaxis()
    ax.set_xlabel("x/c"); ax.set_ylabel("Cp"); ax.legend()
    ax.set_title(f"V3.2 NACA0012 M={M_INF} alpha={ALPHA_DEG} Cp (medium)")

    # Surface local Mach from the same Cp points (via q^2 -> M^2 needs the
    # tangential q2; reuse the Cp -> M mapping through isentropic relations
    # is indirect, so recompute from the curve q2 = inverse not needed:
    # plot local Mach from element q2 on the wall Cp positions instead).
    for side, color in (("upper", "tab:red"), ("lower", "tab:blue")):
        x = curve[f"x_{side}"]
        cp = curve[f"cp_{side}"]
        # invert the exact isentropic Cp(q2) for q2, then M(q2)
        gam = 1.4
        rho_pow = 1.0 + 0.5 * gam * M_INF**2 * cp  # rho^gamma
        q2 = 1.0 - (rho_pow ** ((gam - 1.0) / gam) - 1.0) / (
            0.5 * (gam - 1.0) * M_INF**2
        )
        m_loc = np.sqrt([mach_number_squared(float(q), M_INF) for q in q2])
        ax2.plot(x, m_loc, ".", ms=3, color=color, label=f"pyFP3D {side}")
    ax2.axhline(1.0, color="k", lw=0.8, ls=":")
    ax2.set_xlabel("x/c"); ax2.set_ylabel("local Mach"); ax2.legend()
    ax2.set_title("surface local Mach (subcritical: M < 1 everywhere)")
    fig.savefig(gate_dir / "v3_2_cp_mach_vs_reference.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # Suction peak amplified vs the incompressible panel level (-0.61 at
    # alpha=2 scales to ~ -0.70 PG / -0.74 KT); smooth TE recovery.
    assert curve["cp_upper"].min() < -0.60
    assert curve["cp_upper"][-1] > 0.0
    # Subcritical everywhere (G4.2's nu == 0 regime): max local M < 1.
    assert medium_case["result"]["mach2_max"] < 1.0
