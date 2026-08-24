"""E01 —— G8.1 / G8.2 的已提交物理锚点（E 类：回归/不变性）。

★ 这三条比的是**代码自己提交过的值**：G8.1 coarse (cl 0.41006 / x_shock 0.6196 / M_max 1.3946)、
G8.1 medium (0.523 / 0.674 / 1.404)、G8.2 M6 medium (0.268691 / 1.99687 / 三站激波)。
★★ 有了 CFL3D Euler 参照之后**这三条升格为 D 类** —— E 是「等参照到位就搬空」的中转站。
★ `bench/run_g82_anchor_check.py` 与 `bench/run_m3_budget.py` 从这里取 `_m6_case`。
"""
import os
import time
from pathlib import Path
import numpy as np
import pytest
import scipy.sparse.linalg as spla
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.solve.picard import solve_subsonic_lifting
from bench.recipes import NEWTON_M6_RECIPE  # noqa: E402  (裁决 a, 2026-08-24)
from tests.conftest import REPO_ROOT


run_gates = pytest.mark.skipif(
    os.environ.get("PYFP3D_TRANSONIC_GATES", "0") != "1",
    reason="transonic Newton gate runs take minutes; "
           "set PYFP3D_TRANSONIC_GATES=1 for the gate-closure run",
)


NEWTON_TRANSONIC_RECIPE = dict(
    dm=0.025, dm_min=0.003, freeze_tol=1e-6,
    newton_kw=dict(freeze_refresh_max=8, precond="direct",
                   n_newton_max=60),
)


def _assert_terminal_quadratic(r):
    """G8.1 protocol: the final level must end below 1e-9 AND contain a
    consecutive PAIR of >= 1.5-digit residual collapses (only a quadratic
    tail does that at these levels; cf. Lopez Table 4.9). The pair is
    searched over the level history rather than read off the last two
    entries because the freeze-refresh honesty re-evaluations interleave
    live-residual jumps with the (quadratic) frozen phases -- a final
    frozen phase can legitimately converge in ONE step after a refresh."""
    h = r["residual_history"]
    assert h[-1] < 1e-9, f"final residual {h[-1]:.3e}"
    drops = [h[i + 1] / h[i] for i in range(len(h) - 1)]
    pairs = [(drops[i], drops[i + 1]) for i in range(len(drops) - 1)]
    assert any(a < 3e-2 and b < 3e-2 for a, b in pairs), (
        f"no consecutive quadratic-like drop pair in {drops}")


def _transonic_case(mesh_file, m_inf):
    from tests.conftest import REPO_ROOT
    from pyfp3d.post.section_cut import wall_cp_curve
    from pyfp3d.post.shock import shock_report
    from pyfp3d.solve.newton import solve_newton_transonic

    mesh = read_mesh(REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
                     / mesh_file)
    mc, wc = cut_wake(mesh)
    r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=1.25,
                               **NEWTON_TRANSONIC_RECIPE)
    dz = float(np.ptp(mc.nodes[:, 2]))
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz, m_inf=m_inf),
                       m_inf)
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=1.25, s_ref=dz, m_inf=m_inf)
    return r, rep, forces


@run_gates
def test_g81_terminal_quadratic_coarse_m080():
    """G8.1: coarse M0.80/alpha1.25 -- terminal quadratic convergence to the TRUE
    discrete solution.

    ★ RE-ANCHORED 2026-07-31 (GS1b.11) to entropy-ON values, the correction having
    become the default. Measured at the runner-default 16 threads with this file's
    NEWTON_TRANSONIC_RECIPE: |R| 4.1e-13, x_shock 0.6196, cl 0.41006, M_max 1.3946.

    ANCHORED TO WHAT: unusually for this project, this one has an EXTERNAL reference
    and passes it -- x_shock 0.6196 sits inside the Euler-anchored band 0.62 +- 0.03
    (cases/reference_data/naca0012_m080/), 0.0004 from its centre, where the
    isentropic value 0.658 fell OUTSIDE it. cl and M_max remain drift locks.

    SUPERSEDED isentropic values, kept per discipline #11: shock 0.658, cl 0.459,
    M_max 1.408. Still true of both: this is the Newton answer, dissipation-scan
    robust and continuation-path independent, and NOT the P4 Picard stall state
    0.604/0.334 whose Newton residual is 2.2e-4 (GS1b.7 measured that state's
    residual at 2.198e-04 and watched the Newton walk off it in six steps)."""
    r, rep, forces = _transonic_case("coarse.msh", 0.80)
    assert r["converged"]
    _assert_terminal_quadratic(r)
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert r["F_history"][-1] < 1e-12                   # coupled Kutta
    # per-level instrumentation (capability demo): ascending Mach ramp
    # with one full history set per level
    lvls = r["level_results"]
    assert [lr["m"] for lr in lvls] == sorted(lr["m"] for lr in lvls)
    assert all(len(lr["gamma_history"]) == len(lr["residual_history"])
               for lr in lvls)
    assert lvls[-1]["residual_history"][-1] == r["residual_history"][-1]
    assert rep["upper"]["has_shock"] and rep["upper"]["monotone"]
    assert abs(rep["upper"]["x_shock"] - 0.6196) < 0.012
    assert abs(forces["cl"] - 0.41006) < 0.01
    assert abs(float(np.sqrt(r["mach2_max"])) - 1.3946) < 0.02


@run_gates
@pytest.mark.xfail(strict=True, reason=(
    "GS1b.11: the entropy correction is now the default, and at MEDIUM its answer is "
    "not yet anchorable. Measured (GS1b.9): the sigma factor is FROZEN over a Newton "
    "step, and where it freezes depends on the recipe, so at medium M0.7875 two "
    "converged Newton recipes disagreed by 0.118 c of shock position with the "
    "correction on (0.6029 versus 0.4852) where the isentropic answer agrees to four "
    "decimals (0.6738 both). A sigma self-consistency POLISH cut that to 0.0029 c and "
    "was therefore adopted in GS1b.9 -- then REMOVED in GS1b.11 once measured to be a "
    "coincidence: it never converged, and at medium M0.80 it moved the shock 0.6146 -> "
    "0.7031, out of the Euler band. So the 0.118 c recipe spread stands as the reason "
    "this lock cannot be anchored: anchoring it would fix an arbitrary recipe. Root "
    "(GS1b.4): sigma inherits the churn of the upwind DONOR MAP, whose selection "
    "changes for 0.34 % of medium elements between adjacent Mach steps. Five "
    "discriminator families were tried and excluded (GS1b.10 sec 9), and the "
    "obstruction is now understood: a badly-shaped tet produces the same local Mach "
    "signature as a shock, so no local statistic can separate them. Re-anchor this "
    "lock when the SHOCK LOCALISATION problem is solved -- its shared fix is S2's "
    "gradient reconstruction, which is also G1.6's root."))
def test_g81_terminal_quadratic_medium_m07875():
    """G8.1 (re-specced): medium M0.7875/alpha1.25 -- the strongest
    medium-mesh condition with a reachable isolated solution (M0.80 sits
    at the fold, recorded finding). Regression lock around the measured
    Newton solution: shock 0.674, cl 0.523, M_max 1.404; the frozen
    finish reaches ~8e-11 with the live assignment-discontinuity floor
    reported (~1.3e-7 -- the C0 walk flux's intrinsic floor on this
    tie-degenerate prism mesh)."""
    r, rep, forces = _transonic_case("medium.msh", 0.7875)
    assert r["converged"]
    _assert_terminal_quadratic(r)
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert rep["upper"]["has_shock"] and rep["upper"]["monotone"]
    assert abs(rep["upper"]["x_shock"] - 0.674) < 0.012
    assert abs(forces["cl"] - 0.523) < 0.01
    assert abs(float(np.sqrt(r["mach2_max"])) - 1.404) < 0.02
    if r["residual_unfrozen"] is not None:
        assert r["residual_unfrozen"] < 1e-5


def _m6_case(level, m_inf=0.84, alpha=3.06):
    from pyfp3d.constraints.wake import tip_taper_factors
    from pyfp3d.meshgen.wing3d import B_SEMI
    from pyfp3d.post.section_cut import section_cp_curve
    from pyfp3d.post.shock import shock_report as _sr
    from pyfp3d.post.surface import cl_kj_3d, planform_area
    from pyfp3d.solve.newton import solve_newton_transonic

    mesh_dir = REPO_ROOT / 'cases' / 'meshes' / 'onera_m6'   #: ★ 深度无关
    p = mesh_dir / f"{level}.msh"
    if not p.exists():
        pytest.skip(f"onera_m6/{level}.msh not generated; run "
                    "cases/meshes/onera_m6/generate_onera_m6.py")
    t0 = time.perf_counter()
    mc, wc = cut_wake(read_mesh(p))
    #: ★★ RE-ANCHORED 2026-08-05/06, user ruling. Full record in
    #: phases/p2/docs/dev_phase_two/20260806-0600-g82-reanchor.md, written before this edit as
    #: roadmap sec 5 requires. Why the production taper is now part of the case:
    #:
    #: onera_m6/medium.msh was regenerated on 2026-08-04 when the base level names
    #: flipped flat -> round, and this recipe's anchors were measured on the FLAT cap.
    #: Verified: medium_flat.msh reproduces the committed anchors to 5-6 digits (cl_p
    #: 0.263882 vs 0.263888, M_max 2.10697 vs 2.10709) while the round mesh does not
    #: converge -- so the failure was 100 % the mesh, not the solver.
    #:
    #: Measured on the round mesh at target 0.84, same recipe otherwise:
    #:     no taper   -> NOT converged, |R| 3.13e-06, M_max 2.92867, highest converged
    #:                   LEVEL 0.82 (0.70/0.75/0.80/0.82 pass, 0.84 fails twice), 2949 s
    #:     with taper -> converged, |R| 6.155e-15, M_max 1.99687, 0/0 clamps, 341 s
    #: so the ruling's "re-anchor to the highest Mach the round tip reaches" needs NO
    #: reduction in Mach: the production configuration still reaches 0.84.
    #:
    #: ★ The shock position -- the only externally anchored quantity here (Euler band
    #: 0.60-0.63) -- barely moves: eta 0.44 goes 0.59582 -> 0.59632, five ten-thousandths
    #: of a chord. So this re-anchor moves MODEL quantities (cl, M_max), not the physical
    #: location.
    #:
    #: ★★ And it resolves the concern this file registered right below: with no taper the
    #: walk found a pre-shock Mach near 2 that was "very likely the P13 wing-tip free-edge
    #: singularity rather than a physical shock", i.e. the entropy correction possibly
    #: charging entropy to a NUMERICAL artefact. With the taper M_max drops to 1.99687
    #: against the flat anchor's 2.10709. (A mitigation of the registered item, NOT a
    #: claim that the tip singularity is solved -- P13 itself is untouched.)
    #:
    #: ⚠ cl_p rises +1.82 % vs the flat anchor while the taper ALONE lowers cl (-1.3 %
    #: measured on the wing-body), so two effects are mixed -- mesh flat->round and
    #: taper off->on -- and this round did NOT separate them, because the third leg it
    #: would need (round, no taper, CONVERGED) does not exist. Recorded as a composite,
    #: not attributed.
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", 0.05 * B_SEMI)
    kw = dict(NEWTON_M6_RECIPE)
    kw["newton_kw"] = dict(kw["newton_kw"], tip_taper=taper)
    r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha, **kw)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    forces = wall_force_coefficients(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
        alpha_deg=alpha, s_ref=s_ref, m_inf=m_inf)
    shocks = {}
    for eta in (0.44, 0.65, 0.90):
        c = section_cp_curve(mc, r["phi"], eta=eta, b_semi=B_SEMI,
                             m_inf=m_inf)
        shocks[eta] = _sr(c, m_inf)["upper"]["x_shock"]
    wall = time.perf_counter() - t0
    return r, forces, shocks, wall


@run_gates
def test_g82_m6_medium_newton_end_to_end():
    """G8.2: ONERA M6 medium (63k nodes / 351k tets), M0.84/alpha3.06 --
    Newton end to end (mesh read -> cut -> continuation -> forces+shocks)
    < 5 min single node (P5 Picard: 4539 s). Physics bands are regression
    locks around the measured Newton solution (2026-07-11; the P5 Picard
    numbers cl_p 0.2453 / M_max 1.995 carry the P4-erratum-in-kind caveat:
    Newton residual at that state is 7.6e-6 with Kutta |F| 5.8e-4, and the
    true discrete solution sits ~+7.9% in cl_p with the same shocks).

    Timing protocol (the CLAUDE.md 16-thread cap, quantified here): run
    with NUMBA_NUM_THREADS=16 AND OMP_NUM_THREADS=16 AND
    OPENBLAS_NUM_THREADS=16 on the 16C/32T box -- without the BLAS/OMP
    caps the same run measures ~333 s (oversubscription costs ~33%,
    measured A/B 2026-07-11); with them 252 s (G8.2 closure run), ~145 s
    since the G10.2 intermediate_tol promotion (same final-level
    solution, recipe comment above).

    ★ RE-ANCHORED 2026-07-31 to the ENTROPY-CORRECTED default, and the GS1b.11
    strict xfail REMOVED. That xfail recorded a real defect -- entropy ON gave
    converged=False, |R| stalled at 2.49e-06, 57 floored cells and sigma_min = 0.0
    EXACTLY -- which is now fixed at its root: refresh_sigma was handing the
    m_cap-LIMITED speed field to the entropy walk, so a limited cell's recovered
    "pre-shock Mach" was the cap itself (m1_max read exactly 2.9999999999999996) and
    sigma_RH(3.0) = 0.32834 entered the chain product. Pre-registered and measured in
    docs/dev_phase_two/20260731-2000 / -2200; the guard is locked by
    tests/A/test_A06_entropy_correction.py's mcap group.

    That xfail also left a named-but-unrun check -- cycle versus a very long chain --
    and it is now RUN: both routes are measured and locked. Separate capped shocks
    along an acyclic chain multiply one 0.32834 factor each (twelve give 1.570e-06,
    transport still converged); a capped shock inside a DONOR CYCLE gives exactly 0.0
    with converged=False, since pointer doubling squares the product every round. M6
    medium read exactly 0.0, so it was the CYCLE route -- the name was right.

    Measured post-fix (16 threads, 142 s, 14 Newton steps): converged, |R| 6.63e-15,
    Kutta |F| 2.08e-16, 0 limited / 0 floored, cl_p 0.263888, M_max 2.10709, shocks
    0.59582 / 0.53914 / 0.34225, sigma_min 0.73247, m1_max 1.97516 over 2901 shock
    cells. The bands below are re-centred on those; the SUPERSEDED isentropic anchors
    are cl_p 0.2646, M_max 2.134, shocks 0.596 / 0.541 / 0.362 (the eta = 0.90 one
    had only 0.0002 of margin left inside its own tolerance, which is why this is a
    re-anchor and not a "still passes").

    ⚠ REGISTERED, not resolved: m1_max 1.975 at M_inf 0.84 means the walk finds a
    pre-shock Mach near 2 somewhere, and this recipe carries no tip_taper, so that is
    very likely the P13 wing-tip free-edge singularity rather than a physical shock --
    i.e. the correction may be charging entropy to a NUMERICAL artefact. That is the
    registered shock-localisation problem resurfacing; it does not affect this fix
    (which strictly removes a non-physical input) but it must not be read as
    validation of the correction's magnitude here."""
    r, forces, shocks, wall = _m6_case("medium")
    assert r["converged"]
    _assert_terminal_quadratic(r)
    assert r["n_limited"] == 0 and r["n_floored"] == 0
    assert r["F_history"][-1] < 1e-12                   # coupled Kutta
    #: re-anchored (round mesh + production taper). SUPERSEDED flat-cap anchors, kept
    #: per discipline #11 rather than overwritten: cl 0.263888, M_max 2.10709,
    #: shock(0.44) 0.59582 -- reproducible today on medium_flat.msh to 5-6 digits.
    assert abs(forces["cl"] - 0.268691) < 0.005            # regression lock
    assert abs(float(np.sqrt(r["mach2_max"])) - 1.99687) < 0.05
    #: ★ the re-anchored shock positions, and their SPANWISE PATTERN is a consistency
    #: check on the re-spec rather than three independent numbers: the change grows
    #: monotonically outboard, which is exactly what a TIP model should do --
    #:     eta 0.44   0.59582 -> 0.59632   (+0.0005 c)
    #:     eta 0.65   0.53914 -> 0.54020   (+0.0011 c)
    #:     eta 0.90   0.34225 -> 0.37144   (+0.0292 c)   <- nearest the tip, moves most
    #: The inboard stations staying put is what makes "the physical location is
    #: preserved" a measurement rather than a hope; eta 0.90 exceeding the old 0.02
    #: tolerance is the taper doing its job on the tip loading, not drift.
    assert abs(shocks[0.44] - 0.59632) < 0.02
    assert abs(shocks[0.65] - 0.54020) < 0.02
    assert abs(shocks[0.90] - 0.37144) < 0.02
    #: ★★ The budget is asserted LAST, on purpose, and it was moved here on 2026-08-10.
    #: It used to sit above the physics anchors, and that ordering turned a machine-load
    #: spike into a report with the physics UNMEASURED -- pytest stops at the first
    #: failing assert, so a 588 s reading hid whether cl, M_max and the three shock
    #: positions were still in band. b7 and b9 sprang exactly that trap on 2026-08-09
    #: (both had later assertions hidden behind the first). Nothing about the threshold
    #: changed; only the order, so a slow machine can no longer suppress the physics.
    #:
    #: ★ And the reading that forced this: the SAME solve measured 588 s and then 109 s
    #: within half an hour on this box (load average 22 -> 15), with the answer
    #: bit-reproducing the committed anchors to six decimals -- cl 0.268691, M_max
    #: 1.996867, shocks 0.596316 / 0.540203 / 0.371440
    #: (bench/gate_results/g82_anchor_check.csv). A 5.4x spread on one machine means a
    #: fixed wall-clock bound is a CALIBRATION of the machine, not a guarantee about the
    #: solver -- the same lesson as the EW forcing, the taper r_c and the descent10
    #: threshold. It is kept as a gate because a real 2x regression must still fail, but
    #: a red here is a timing reading first and a capability claim only after the
    #: physics above has passed.
    #:
    #: budget raised 300 -> 450 s in the 2026-08-06 re-spec, and NOT because of the
    #: round tip: the FLAT leg measured 331 s in the same session and the round+taper leg
    #: 341 s, so both exceed the old bound and this machine is simply slower than the P8
    #: record's CI reference of 301.66 s. 450 s leaves ~30 % over the measured 341 s.
    assert wall < 450.0, (
        f"G8.2 budget: {wall:.0f} s >= 450 s. Physics anchors above all PASSED, so this "
        f"is a timing reading; check machine load before treating it as a regression.")
