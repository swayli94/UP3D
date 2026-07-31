"""
P4 gates G4.1 (transonic shock quality vs reference) and G4.3
(robustness sweep). This module is the numeric GATE: it solves and
asserts (converged, shock position/monotonicity, no limited cells) and
logs a numeric summary CSV to the ephemeral artifacts/ dir.

The COMMITTED figures + CSVs that demo_report embeds (the coarse/medium
surface-Cp refinement pair and the V4.3 sweep dashboard) and the code
that generates them live in the demo instead --
cases/demo/p4_transonic/run_demo.py, heavy mode (PYFP3D_TRANSONIC_GATES=1)
-- because artifacts/ is gitignored. This test therefore no longer draws
figures; it verifies the numbers behind them.

Runtime policy: the always-on part of this module is a COARSE-mesh
transonic smoke check (~2.5 min, guards the whole P4 machinery in the
regression suite). The medium-mesh G4.1 gate run and the G4.3 sweep are
gated behind PYFP3D_TRANSONIC_GATES=1 (each is several minutes to tens
of minutes of Picard iterations -- run explicitly for gate closure).
P7 owns making these fast (Newton).

Reference: cases/reference_data/naca0012_m080/ (Euler anchor + documented
conservative-FP shift band; provenance in its README.md).

Convergence semantics at transonic (documented in
solve/picard.py::solve_subsonic_lifting and solve/continuation.py): the
pseudo-time-stabilized density iteration settles into an
engineering-converged regime -- physical M_max, no limited/floored
cells, Kutta mismatch below the evaluation-noise-matched tol_gamma --
rather than the subsonic 1e-10 residual; the sharp-shock residual tail
is a known, bounded limitation until P7's Newton.
"""

import csv
import os

import numpy as np
import pytest

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


def _write_g41_summary(case, level, artifacts_dir):
    """Numeric gate log to the ephemeral artifacts/ dir (the committed
    figure + this same summary are produced by the demo, see module
    docstring)."""
    r, rep = case["result"], case["shock"]
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
        # GS1b.8: the same writer now serves BOTH paths, and they report different
        # keys -- the Picard path has kutta_mismatch / n_picard_total, the coupled
        # Newton has n_newton and a field-residual history. Write whatever the
        # result actually carries rather than assuming one shape.
        for key in ("kutta_mismatch", "n_picard_total", "residual_final",
                    "engineering_converged", "n_newton", "converged"):
            if key in r:
                w.writerow([key, r[key]])
        if r.get("residual_history"):
            w.writerow(["residual_final_history", r["residual_history"][-1]])
        w.writerow(["n_limited", r["n_limited"]])
        w.writerow(["n_floored", r["n_floored"]])


def _transonic_case_newton(mesh_path, **kw):
    """G4.1's solve, re-spec'd onto the COUPLED NEWTON path (GS1b.8, user ruling
    2026-07-30).

    Why: GS1b.7 measured that the Picard-continuation terminal state at this very
    condition is NOT a solution of the discrete equations -- |R| = 2.20e-04 against
    the Newton's 5.46e-12, the Newton walks off it in six steps, and the Newton's
    residual evaluated at it is 2.198e-04. Its shock sits 0.054c from the converged
    one, 2.7x the +-0.02 a product criterion asks of a shock position. So the gate
    now runs on the path whose state IS a solution: converged to ~1e-12,
    reproducible from a heterogeneous seed, unclamped.
    """
    from pyfp3d.solve.newton import solve_newton_lifting
    mesh = read_mesh(mesh_path)
    mc, wc = cut_wake(mesh)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, upwind_c=1.5,
                             m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
                             precond="direct", direct_refactor_every=4,
                             n_newton_max=80, **kw)
    dz = float(np.ptp(mc.nodes[:, 2]))
    curve = wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=M_INF)
    rep = shock_report(curve, M_INF)
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF,
    )
    return {"mesh_cut": mc, "wc": wc, "result": r, "curve": curve,
            "shock": rep, "forces": forces}


def _assert_g41(case, reference_mesh_dir):
    rep = case["shock"]
    r = case["result"]
    ref = _reference(reference_mesh_dir)

    # `converged` is asserted by the CALLER (its meaning differs between the two
    # paths: the Newton's is a field-residual convergence, the Picard's is the
    # engineering regime -- GS1b.8 split the flags so they cannot impersonate each
    # other).
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


def test_p4_picard_path_historical_regression(reference_mesh_dir, artifacts_dir):
    """★ NOT A PHYSICS GATE -- a drift lock on the Picard-continuation path.

    GS1b.7 measured that this path's terminal state is not a solution of the
    discrete equations (|R| 2.20e-04 versus the Newton's 5.46e-12; the Newton walks
    off it; the Newton's residual at it is 2.198e-04), and its shock sits 0.054c
    from the converged one. So this test asserts only what the path actually
    achieves -- `engineering_converged`, the explicitly named old semantics -- plus
    the shock-quality checks, and it deliberately does NOT compare the position
    against the Euler-anchored reference. The physics gate moved onto the Newton
    path in test_g41_transonic_coarse_newton.
    """
    from .conftest import REPO_ROOT
    # ★ GS1b.11: entropy_correction is now the DEFAULT, and this test exists to lock the
    # HISTORICAL (isentropic Picard) behaviour, so it pins the flag OFF explicitly. With
    # the default it reads 0.5633 instead of 0.6041 -- a correct move, not a drift, and
    # exactly the use the switch was kept for (the user's 2026-07-30 ruling that the
    # correction stays switchable).
    case = _transonic_case(
        REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh",
        entropy_correction=False)
    _write_g41_summary(case, "coarse", artifacts_dir)
    r, up = case["result"], case["shock"]["upper"]
    assert r["engineering_converged"], r.get("not_converged_reason")
    assert not r["converged"], (
        "the Picard path now reports true convergence at this condition -- if that "
        "is real, GS1b.7's attribution needs revisiting and this test should be "
        "promoted back to a physics gate")
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert up["has_shock"] and up["monotone"] and not up["expansion_shock"]
    # drift lock on the (non-converged) historical value, NOT a physics claim
    assert abs(up["x_shock"] - 0.6041) < 0.01, (
        f"Picard-path shock drifted to {up['x_shock']:.4f} from the recorded "
        f"0.6041 (GS1b.7)")


#: ★ The GS1b.8 xfail is GONE: it said "this becomes a PASS when entropy_correction
#: becomes the default", and after GS1b.11 flipped it that prediction held --
#: x_shock 0.6073 at coarse, inside the Euler-anchored band 0.62 +- 0.03 (-0.0127 from
#: its centre), where the isentropic value 0.6581 fell outside. Recorded because a
#: prediction that is checked is worth more than one that is quietly dropped.
def test_g41_transonic_coarse_newton(reference_mesh_dir, artifacts_dir):
    """Gate G4.1 on the Newton path: the shock position of a state that IS a
    solution, against the Euler-anchored reference band."""
    from .conftest import REPO_ROOT
    case = _transonic_case_newton(
        REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
    r = case["result"]
    assert r["converged"] and not r["clamped"], (
        f"Newton path not converged: |R| = {r['residual_history'][-1]:.2e}")
    _assert_g41(case, reference_mesh_dir)


@run_gates
#: ★ Also un-xfailed by GS1b.11, but NOT smoothly -- worth recording. The GS1b.8 reason
#: predicted a pass, then GS1b.9's sigma polish (added afterwards) moved this leg to
#: x_shock 0.7031, OUT of the band, and the A3 prediction check caught it. Removing the
#: polish -- which had never converged, so its apparent benefit was a coincidence --
#: brings this leg to x_shock 0.6006, INSIDE the band (-0.0194 from centre), and the
#: isentropic leg here does not converge at all.
def test_g41_transonic_medium_gate(reference_mesh_dir, artifacts_dir):
    """Gate G4.1 = V4 on the medium mesh, on the Newton path (GS1b.8)."""
    from .conftest import REPO_ROOT
    case = _transonic_case_newton(
        REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh")
    r = case["result"]
    assert r["converged"] and not r["clamped"], (
        f"Newton path not converged: |R| = {r['residual_history'][-1]:.2e}")
    _write_g41_summary(case, "medium", artifacts_dir)
    _assert_g41(case, reference_mesh_dir)


@run_gates
def test_g43_robustness_sweep(artifacts_dir):
    """Gate G4.3: M in {0.74..0.82} x alpha in {0, 1.25} deg, ONE parameter set
    (TRANSONIC_DEFAULTS) -- a DRIVER-ROBUSTNESS sweep: does the Picard continuation
    survive the envelope with a physical field? V4.3 dashboard artifact.

    ★ GS1b.8 (2026-07-30): the assertion is on `engineering_converged`, the
    explicitly named old semantics (physical + Kutta), because
    solve_transonic_lifting's `converged` now also requires the FIELD RESIDUAL and
    the transonic points sit on the Picard shock plateau (|R| ~ 2e-04; GS1b.7
    measured that such a state is not a solution of the discrete equations). That
    is the faithful reading of what this gate always measured -- it asks whether
    the DRIVER survives the envelope, not whether the answer is right -- so the
    true residual is now RECORDED per point in the summary CSV instead of being
    implied by a flag that meant something weaker than its name.

    The answer-quality gate lives on the Newton path
    (test_g41_transonic_coarse_newton). A Newton-path version of this sweep is a
    registered follow-up, not done here (ten more coupled solves)."""
    from .conftest import REPO_ROOT
    from pyfp3d.constraints.wake import kutta_targets  # noqa: F401

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh")
    mc, wc = cut_wake(mesh)
    dz = float(np.ptp(mc.nodes[:, 2]))

    rows = [("alpha_deg", "m_inf", "engineering_converged", "converged",
             "residual_final", "kutta_mismatch",
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
            rows.append((alpha, m, r["engineering_converged"],
                         r["converged"], f"{r['residual_final']:.3e}",
                         f"{r['kutta_mismatch']:.2e}",
                         f"{np.sqrt(r['mach2_max']):.3f}",
                         f"{rep['upper']['x_shock']:.4f}",
                         f"{forces['cl']:.5f}", r["n_limited"]))

    gate_dir = artifacts_dir / "G4.3"
    gate_dir.mkdir(parents=True, exist_ok=True)
    with open(gate_dir / "summary.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)
    # The committed V4.3 dashboard is produced by the demo (see docstring).

    for (alpha, m), (r, rep, _f) in results.items():
        assert r["engineering_converged"], (
            f"alpha={alpha} M={m}: {r.get('not_converged_reason')}")
        assert r["n_limited"] == 0 and r["n_floored"] == 0
        assert not rep["upper"]["expansion_shock"]
        # smooth trend: shock exists for all supercritical lifting cases
        if alpha > 0 or m >= 0.78:
            assert rep["upper"]["has_shock"], f"no shock at alpha={alpha} M={m}"
