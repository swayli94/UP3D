"""
P4 gates G4.1 (transonic shock quality vs reference) and G4.3
(robustness sweep), with the V4.1/V4.3 headless artifacts.

Runtime policy: the always-on part of this module is a COARSE-mesh
transonic smoke check (~2.5 min, guards the whole P4 machinery in the
regression suite). The medium-mesh G4.1 gate run and the G4.3 sweep are
gated behind PYFP3D_TRANSONIC_GATES=1 (each is several minutes to tens
of minutes of Picard iterations -- run explicitly for gate closure;
their evidence lands in artifacts/G4.1 / artifacts/G4.3 and
docs/demo_report.md). P6 owns making these fast (Newton).

Reference: cases/reference_data/naca0012_m080/ (Euler anchor + documented
conservative-FP shift band; provenance in its README.md).

Convergence semantics at transonic (documented in
solve/picard.py::solve_subsonic_lifting and solve/continuation.py): the
pseudo-time-stabilized density iteration settles into an
engineering-converged regime -- physical M_max, no limited/floored
cells, Kutta mismatch below the evaluation-noise-matched tol_gamma --
rather than the subsonic 1e-10 residual; the sharp-shock residual tail
is a known, bounded limitation until P6's Newton.
"""

import csv
import os

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import wall_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.continuation import solve_transonic_lifting

M_INF = 0.80
ALPHA = 1.25

run_gates = pytest.mark.skipif(
    os.environ.get("PYFP3D_TRANSONIC_GATES", "0") != "1",
    reason="medium-mesh transonic gate runs take several minutes; "
           "set PYFP3D_TRANSONIC_GATES=1 for the gate-closure run",
)


def _reference(reference_mesh_dir):
    path = reference_mesh_dir / "naca0012_m080" / "shock_reference.csv"
    ref = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            ref[row["quantity"]] = (float(row["value"]),
                                    float(row["tolerance"]),
                                    row["gated"] == "yes")
    return ref


def _transonic_case(mesh_path, **kw):
    mesh = read_mesh(mesh_path)
    mc, wc = cut_wake(mesh)
    r = solve_transonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                max_gamma_evals=12, n_picard_eval=800, **kw)
    dz = float(np.ptp(mc.nodes[:, 2]))
    curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=M_INF)
    rep = shock_report(curve, M_INF)
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF,
    )
    return {"mesh_cut": mc, "wc": wc, "result": r, "curve": curve,
            "shock": rep, "forces": forces}


def _write_g41_artifacts(case, level, reference_mesh_dir, artifacts_dir):
    r, rep, curve = case["result"], case["shock"], case["curve"]
    gate_dir = artifacts_dir / "G4.1"
    gate_dir.mkdir(parents=True, exist_ok=True)

    with open(gate_dir / f"summary_{level}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value"])
        for side in ("upper", "lower"):
            for k, v in rep[side].items():
                w.writerow([f"{side}_{k}", v])
        w.writerow(["cp_critical", rep["cp_critical"]])
        w.writerow(["cl_pressure", case["forces"]["cl"]])
        w.writerow(["cl_kj", 2.0 * float(r["gamma"][0])])
        w.writerow(["gamma", float(r["gamma"][0])])
        w.writerow(["mach_max", float(np.sqrt(r["mach2_max"]))])
        w.writerow(["kutta_mismatch", r["kutta_mismatch"]])
        w.writerow(["n_picard_total", r["n_picard_total"]])
        w.writerow(["n_limited", r["n_limited"]])
        w.writerow(["n_floored", r["n_floored"]])

    # V4.1: Cp with Cp* line and shock markers.
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(curve["x_upper"], curve["cp_upper"], ".-", ms=3, lw=0.7,
            color="tab:red", label="upper")
    ax.plot(curve["x_lower"], curve["cp_lower"], ".-", ms=3, lw=0.7,
            color="tab:blue", label="lower")
    ax.axhline(rep["cp_critical"], color="k", lw=0.8, ls=":",
               label="Cp* (sonic)")
    for side, color in (("upper", "tab:red"), ("lower", "tab:blue")):
        if rep[side]["has_shock"]:
            ax.axvline(rep[side]["x_shock"], color=color, lw=0.8, ls="--")
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel("Cp")
    ax.set_title(f"V4.1 NACA0012 M={M_INF} alpha={ALPHA} ({level}): "
                 f"upper shock x/c={rep['upper']['x_shock']:.3f}")
    ax.legend()
    fig.savefig(gate_dir / f"v4_1_cp_shock_{level}.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)


def _assert_g41(case, reference_mesh_dir):
    rep = case["shock"]
    r = case["result"]
    ref = _reference(reference_mesh_dir)

    assert r["converged"], (
        f"transonic solve not converged: kutta mismatch "
        f"{r['kutta_mismatch']:.2e}, n_limited={r['n_limited']}"
    )
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    up = rep["upper"]
    assert up["has_shock"]
    assert up["monotone"], "non-monotone shock jump"
    assert not up["expansion_shock"], "expansion shock detected"
    assert up["n_cells"] <= 3, f"shock smeared over {up['n_cells']} stations"
    x_ref, tol, gated = ref["upper_shock_x_c"]
    assert gated
    assert abs(up["x_shock"] - x_ref) <= tol, (
        f"upper shock x/c = {up['x_shock']:.3f} vs reference "
        f"{x_ref} +/- {tol}"
    )


def test_g41_transonic_coarse_smoke(reference_mesh_dir, artifacts_dir):
    """Always-on coarse guard for the full P4 pipeline: continuation +
    Gamma root-find + shock quality checks (same asserts as the medium
    gate, coarse mesh)."""
    from .conftest import REPO_ROOT
    case = _transonic_case(
        REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
    _write_g41_artifacts(case, "coarse", reference_mesh_dir, artifacts_dir)
    _assert_g41(case, reference_mesh_dir)


@run_gates
def test_g41_transonic_medium_gate(reference_mesh_dir, artifacts_dir):
    """Gate G4.1 = V4 on the medium mesh."""
    from .conftest import REPO_ROOT
    case = _transonic_case(
        REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh")
    _write_g41_artifacts(case, "medium", reference_mesh_dir, artifacts_dir)
    _assert_g41(case, reference_mesh_dir)


@run_gates
def test_g43_robustness_sweep(artifacts_dir):
    """Gate G4.3: M in {0.74..0.82} x alpha in {0, 1.25} deg, ONE
    parameter set (TRANSONIC_DEFAULTS), all levels of both continuation
    passes converge with a physical field. V4.3 dashboard artifact."""
    from .conftest import REPO_ROOT
    from pyfp3d.constraints.wake import kutta_targets  # noqa: F401

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)
    dz = float(np.ptp(mc.nodes[:, 2]))

    rows = [("alpha_deg", "m_inf", "converged", "kutta_mismatch",
             "mach_max", "x_shock_upper", "cl", "n_limited")]
    results = {}
    for alpha in (0.0, 1.25):
        for m in (0.74, 0.76, 0.78, 0.80, 0.82):
            r = solve_transonic_lifting(mc, wc, m_inf=m, alpha_deg=alpha,
                                        max_gamma_evals=12,
                                        n_picard_eval=800)
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

    gate_dir = artifacts_dir / "G4.3"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)

    # V4.3 dashboard: shock x/c and cl trends vs Mach per alpha.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for alpha, mk in ((0.0, "o"), (1.25, "s")):
        ms = [0.74, 0.76, 0.78, 0.80, 0.82]
        xs = [results[(alpha, m)][1]["upper"]["x_shock"] for m in ms]
        cls = [results[(alpha, m)][2]["cl"] for m in ms]
        ax1.plot(ms, xs, mk + "-", label=f"alpha={alpha}")
        ax2.plot(ms, cls, mk + "-", label=f"alpha={alpha}")
    ax1.set_xlabel("M_inf"); ax1.set_ylabel("upper shock x/c"); ax1.legend()
    ax2.set_xlabel("M_inf"); ax2.set_ylabel("cl"); ax2.legend()
    ax1.set_title("V4.3 shock migration"); ax2.set_title("V4.3 lift trend")
    fig.savefig(gate_dir / "v4_3_sweep_dashboard.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    for (alpha, m), (r, rep, _f) in results.items():
        assert r["converged"], f"alpha={alpha} M={m} did not converge"
        assert r["n_limited"] == 0 and r["n_floored"] == 0
        assert not rep["upper"]["expansion_shock"]
        # smooth trend: shock exists for all supercritical lifting cases
        if alpha > 0 or m >= 0.78:
            assert rep["upper"]["has_shock"], f"no shock at alpha={alpha} M={m}"
