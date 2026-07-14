"""
Track B / B15 gate GB15.1: frozen-selection (N5) on the LEVEL-SET per-side
upwind operator.

Why this exists (roadmap B15): the LS Newton re-selects the per-side upwind walk
live at every step, so on a shock the max(nu_e, nu_u) near-tie band churns, the
residual carries a moving discontinuity, and the 0-limited/0-floored convergence
gate (`newton_ls.py`) can never fire -- the recorded B6 medium M0-embedded
limit-cycle at 3e-6. The conforming path solved this at P8/N5 by FREEZING the
(upstream, branch) assignment; B15 ports that to the per-side masked-graph
operator.

The claims locked here:
  * GB15.1a  `freeze_side_state` captures a per-side (upstream, branch) pair.
  * GB15.1b  the frozen sweep REPRODUCES the live density BITWISE at the freeze
             point -- the property the whole method rests on (if this fails, the
             freeze silently changes the equations being solved).
  * GB15.1c  a clean freeze has n_floored == 0 BY DESIGN (no floor is re-applied
             on branches 0-2), which is exactly what unblocks the 0-clamped gate.
  * GB15.1d  ** THE FD GATE **: the Jacobian assembled under a frozen selection
             is the true derivative of the residual assembled under the SAME
             frozen selection. Both come from the shipped `LSNewtonSystem`, so
             this tests the solver's own code path, not a re-implementation.
  * GB15.1e  freeze_tol=None (the default) is bit-identical to the pre-B15 live
             solver.

FD protocol (inherited from P7/P8): within a frozen assignment the residual is
smooth EXCEPT at (i) the sonic threshold inside nu = C*max(0, 1 - Mc^2/M^2) and
(ii) the q^2 speed limiter, whose flat clamp gives a one-sided derivative. Probes
that straddle either locus read a branch AVERAGE, so those rows are excluded by
an epsilon-guard (measure-zero locus; the assignment is settled near a solution).
The field is noise-broken first -- separable fields on a prism-split mesh park
whole slabs exactly ON the tie (the P7 trap).
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.physics.isentropic import limit_q2_field, mach_squared_field
from pyfp3d.solve.newton_ls import LSNewtonSystem, solve_multivalued_newton
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
ALPHA = 2.0
M_INF = 0.75          # transonic: a real supersonic pocket, so Terms 2/3 bite
UPWIND_C, M_CRIT = 1.5, 0.95


def _mesh(level="coarse"):
    path = M0_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _mvop(mesh, alpha=ALPHA):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _transonic_state(mesh, mvop, seed=25):
    """A settled transonic state with a real supersonic pocket (nu > 0), then
    noise-broken so no slab of elements sits exactly on the max(nu_e,nu_u) tie."""
    r = solve_multivalued_lifting(
        mesh=mesh, mvop=mvop, m_inf=M_INF, alpha_deg=ALPHA,
        farfield="neumann", upwind_c=UPWIND_C, m_crit=M_CRIT,
        damping_theta=0.2, omega_rho=0.5, n_outer_max=seed, tol_residual=None)
    phi = r["phi_ext"].copy()
    rng = np.random.default_rng(20260715)
    phi += 1e-6 * rng.standard_normal(phi.shape)
    return phi


def _kink_rows(mvop, sysm, phi, frozen, band=5e-3):
    """Extended DOFs whose residual row is near a frozen-system kink: the sonic
    threshold (M^2 ~ Mc^2) or the q^2 limiter. Also pulls in elements whose
    UPSTREAM is a kink element (rho_tilde_e depends on q^2_{u(e)})."""
    bad = np.zeros(mvop.n_total, dtype=bool)
    up, lo = mvop.newton_side_data(
        phi, M_INF, UPWIND_C, M_CRIT, 1.4, 1.0, 3.0, 0.05, frozen=frozen)
    mc2 = M_CRIT ** 2
    for side in (up, lo):
        grad = side["grad"]
        q2 = np.einsum("ij,ij->i", grad, grad)
        q2l = limit_q2_field(q2, M_INF, 3.0, 1.4)
        m2 = mach_squared_field(q2l, M_INF, 1.4)
        near_sonic = np.abs(m2 - mc2) < band * mc2
        limited = ~side["lim_mask"]
        kink = near_sonic | limited
        # an element sees its upstream's m2 through s_u
        kink |= kink[side["upstream"]]
        dofs = side["dofvec"][kink & side["keep"]]
        if dofs.size:
            bad[np.unique(dofs)] = True
    return bad


# ---------------------------------------------------------------------------
# GB15.1a/b/c -- the freeze capture and its defining properties.
# ---------------------------------------------------------------------------

def test_freeze_side_state_shapes():
    mesh = _mesh()
    mvop = _mvop(mesh)
    phi = _transonic_state(mesh, mvop)
    frozen = mvop.freeze_side_state(phi, M_INF, UPWIND_C, M_CRIT)
    assert len(frozen) == 2                      # (upper, lower)
    n_tets = mvop.op.n_tets
    for upstream, branch in frozen:
        assert upstream.shape == (n_tets,) and upstream.dtype == np.int64
        assert branch.shape == (n_tets,) and branch.dtype == np.int8
        assert set(np.unique(branch)).issubset({0, 1, 2, 3})


def test_frozen_reproduces_live_density_bitwise_at_freeze_point():
    """THE load-bearing property: at the state where it was captured, the frozen
    sweep must reproduce the live sweep BIT FOR BIT on every consumed (own-side)
    entry. Otherwise freezing silently perturbs the equations."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    phi = _transonic_state(mesh, mvop)
    live_u, live_l = mvop.newton_side_data(phi, M_INF, UPWIND_C, M_CRIT)
    assert mvop.nu_max > 0.0, "no supersonic pocket: the test is vacuous"

    frozen = mvop.freeze_side_state(phi, M_INF, UPWIND_C, M_CRIT)
    fz_u, fz_l = mvop.newton_side_data(phi, M_INF, UPWIND_C, M_CRIT,
                                       frozen=frozen)
    for live, fz in ((live_u, fz_u), (live_l, fz_l)):
        keep = live["keep"]
        assert np.array_equal(live["rho_tilde"][keep], fz["rho_tilde"][keep])
        # the frozen selection IS the live walk's selection
        assert np.array_equal(live["upstream"][keep], fz["upstream"][keep])


def test_clean_freeze_has_zero_floored_by_design():
    """A freeze taken at a 0-floored state stays 0-floored under the frozen
    sweep (no floor is re-applied on branches 0-2) -- this is what unblocks the
    0-clamped convergence gate at a shock."""
    mesh = _mesh()
    mvop = _mvop(mesh)
    phi = _transonic_state(mesh, mvop)
    mvop.newton_side_data(phi, M_INF, UPWIND_C, M_CRIT)
    if mvop.n_floored != 0:
        pytest.skip("seed state is floored; the freeze trigger would not arm")

    frozen = mvop.freeze_side_state(phi, M_INF, UPWIND_C, M_CRIT)
    for _, branch in frozen:
        assert np.count_nonzero(branch == 3) == 0
    mvop.newton_side_data(phi, M_INF, UPWIND_C, M_CRIT, frozen=frozen)
    assert mvop.n_floored == 0


# ---------------------------------------------------------------------------
# GB15.1d -- THE FD GATE. J(frozen) is the derivative of R(frozen).
# ---------------------------------------------------------------------------

def test_frozen_jacobian_matches_fd():
    mesh = _mesh()
    mvop = _mvop(mesh)
    phi = _transonic_state(mesh, mvop)

    # b_base is phi-independent, so it cancels in the central difference and in
    # the Jacobian: the FD gate does not need the far-field RHS.
    sysm = LSNewtonSystem(mvop, M_INF, upwind_c=UPWIND_C, m_crit=M_CRIT)
    frozen = sysm.freeze(phi)

    A, R0, up, lo = sysm.residual(phi, frozen)
    assert mvop.nu_max > 0.0, "no supersonic pocket: the test is vacuous"
    J = sysm.jacobian(A, up, lo, phi)

    ff = np.unique(mesh.boundary_faces["farfield"])
    is_dir = np.zeros(mvop.n_total, dtype=bool)
    is_dir[ff] = True
    free = np.flatnonzero(~is_dir)

    rng = np.random.default_rng(7)
    d = rng.standard_normal(mvop.n_total)
    d[is_dir] = 0.0                       # stay on the admissible manifold
    d /= np.linalg.norm(d)

    # eps=1e-5 is the central-difference sweet spot here: below it the error is
    # round-off dominated and grows ~1/eps (measured 6.7e-9 / 5.8e-8 / 6.0e-7 at
    # 1e-5 / 1e-6 / 1e-7 -- that clean scaling is itself the evidence that this
    # is a true derivative and not a coincidental match).
    eps = 1e-5
    _, Rp, _, _ = sysm.residual(phi + eps * d, frozen)
    _, Rm, _, _ = sysm.residual(phi - eps * d, frozen)
    fd = (Rp - Rm) / (2.0 * eps)
    jv = J @ d

    bad = _kink_rows(mvop, sysm, phi, frozen)
    rows = free[~bad[free]]
    assert rows.size > 0.90 * free.size, (       # measured: 96.9% kept
        f"epsilon-guard excluded too much: {rows.size}/{free.size}")

    scale = max(np.max(np.abs(fd[rows])), np.max(np.abs(jv[rows])))
    rel = np.max(np.abs(fd[rows] - jv[rows])) / scale
    assert rel < 1e-7, f"frozen Jacobian vs FD: rel={rel:.3e}"   # measured 6.7e-9


# ---------------------------------------------------------------------------
# GB15.1e -- the default is the pre-B15 solver.
# ---------------------------------------------------------------------------

def test_live_newton_limit_cycles_without_freeze():
    """GB15.2 (the motivation, locked): with the LIVE per-side selection the LS
    Newton does NOT converge on the coarse transonic case -- it parks in a
    genuine limit cycle (measured period-6: 3.2e-7, 2.8e-7, 2.7e-7, 1.3e-6,
    8.6e-7, 4.3e-7, repeating) at |R| ~ 3e-7, three orders above tol, with 0
    limited/floored. This is the assignment churn the freeze exists to remove;
    if this test ever starts PASSING convergence, B15's premise has changed."""
    mesh = _mesh()
    live = solve_multivalued_newton(
        mvop=_mvop(mesh), mesh=mesh, m_inf=M_INF, alpha_deg=ALPHA,
        farfield="neumann", upwind_c=UPWIND_C, m_crit=M_CRIT, n_seed=25,
        n_newton_max=60, tol_residual=1e-10)
    assert not live["converged"]
    assert live["n_limited"] == 0 and live["n_floored"] == 0   # a CLEAN stall
    assert live["residual_history"][-1] > 1e-9                 # measured ~3e-7


def test_freeze_cures_the_limit_cycle_and_keeps_the_same_solution():
    """GB15.2 (the win): arming the freeze converts that limit cycle into a
    converged solve -- measured 22 steps to |R| 8.5e-13 (vs live: 60 steps
    stuck at 2.7e-7) -- WITHOUT moving the answer (gamma 0.218809 vs the live
    cycle's 0.218804).

    The freeze is a convergence AID, never a change of equations: if it moved
    the solution it would be worthless.
    """
    mesh = _mesh()
    common = dict(mesh=mesh, m_inf=M_INF, alpha_deg=ALPHA, farfield="neumann",
                  upwind_c=UPWIND_C, m_crit=M_CRIT, n_seed=25,
                  n_newton_max=60, tol_residual=1e-10)
    live = solve_multivalued_newton(mvop=_mvop(mesh), **common)
    frz = solve_multivalued_newton(mvop=_mvop(mesh), freeze_tol=1e-6, **common)

    assert frz["converged"] and frz["froze"]
    assert frz["n_newton"] < live["n_newton"]          # 22 vs 60
    assert frz["residual_history"][-1] < 1e-10         # measured 8.5e-13
    assert frz["n_limited"] == 0 and frz["n_floored"] == 0
    assert frz["n_freeze_reverts"] == 0
    # SAME solution as the (stuck) live iteration -- the freeze only removes the
    # churn, it does not select a different state.
    assert abs(frz["gamma"] - live["gamma"]) < 1e-4, (frz["gamma"],
                                                      live["gamma"])
    # HONEST acceptance: either the live residual is tight at the frozen
    # solution ("tol"), or we accepted the intrinsic assignment-discontinuity
    # floor ("assignment_cycle") -- and in the latter case the LIVE residual is
    # reported, because a frozen finish shows 0 floored BY DESIGN and so can
    # never be its own evidence.
    assert frz["accept_reason"] in ("tol", "assignment_cycle")
    assert frz["residual_unfrozen"] is not None


def test_freeze_trigger_is_tol_only_not_stall():
    """The freeze must NOT arm on a plateau-looking-but-still-DESCENDING trace.

    Regression for a measured B15 defect: the conforming driver also freezes on
    `live_stalled` (newton.py:573), and porting that verbatim made the LS solver
    freeze a still-MOVING assignment -> the frozen step diverges -> revert ->
    re-arm: 3 reverts and NO convergence on medium M0.75, a case the untouched
    live path converges (54 steps, 7.5e-12). With the stall trigger removed the
    same solve freezes late and converges (53 steps, 2.1e-12, 0 reverts, exactly
    the live gamma). So the trigger is `res < freeze_tol` ALONE; `stall_window`
    now only feeds the opt-in accept_on_stall route.
    """
    import inspect
    src = inspect.getsource(solve_multivalued_newton)
    trigger = src.split("--- freeze trigger")[1].split("frozen = sysm.freeze")[0]
    code = "\n".join(ln.split("#")[0] for ln in trigger.splitlines())
    assert "res < freeze_tol" in code
    assert "live_stalled" not in code, (
        "the stall trigger re-entered the freeze ARMING condition")
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["freeze_max_reverts"].default == 3      # fail-safe: disarm


def test_freeze_off_by_default_and_bit_identical():
    import inspect
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["freeze_tol"].default is None          # machinery disarmed
    assert p["tol_residual_loose"].default is None
    assert p["tol_residual_rel"].default is None
    assert p["accept_on_stall"].default is False

    mesh = _mesh()
    common = dict(mesh=mesh, m_inf=0.7, alpha_deg=ALPHA, farfield="neumann",
                  n_seed=20, n_newton_max=25)
    a = solve_multivalued_newton(mvop=_mvop(mesh), **common)
    b = solve_multivalued_newton(mvop=_mvop(mesh), freeze=False, **common)
    assert a["converged"] and b["converged"]
    assert np.array_equal(a["phi_ext"], b["phi_ext"])
    assert not a["froze"] and a["n_freeze_refresh"] == 0
    assert a["accept_reason"] == "tol"


# ---------------------------------------------------------------------------
# GB15.4 errata -- the four traps, locked so they cannot silently return.
# ---------------------------------------------------------------------------

def test_zero_te_nodes_raises_instead_of_returning_nan():
    """Erratum 1. A wake level-set that matches NO TE wall node carries no Kutta
    condition: Gamma is unpinned and the solve returns NaN. Measured on M6
    medium with a TE polyline off by 2e-4 (hand-rolled constants instead of
    meshgen.wing3d.x_te): 0 TE nodes -> 340k limited cells -> gamma = mean([]) =
    NaN, and the solver passed SILENTLY. It must raise instead."""
    from pyfp3d.solve.picard_ls import solve_multivalued_lifting

    mesh = _mesh()
    z = mesh.nodes[:, 2]
    # a TE polyline displaced well off the wall -> matches no TE node
    wls = WakeLevelSet(
        np.array([[5.0, 3.0, z.min()], [5.0, 3.0, z.max()]]),
        direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    if len(cm.te_nodes) != 0:
        pytest.skip("the displaced polyline still matched TE nodes")
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)

    for solver, kw in ((solve_multivalued_lifting, {}),
                       (solve_multivalued_newton, {})):
        with pytest.raises(ValueError, match="0 TE nodes"):
            solver(mvop=mvop, mesh=mesh, m_inf=0.5, alpha_deg=ALPHA, **kw)


def test_freeze_max_clamped_default_is_the_conforming_rule():
    """Erratum 4 / the P9-G9.1 wall. `freeze_max_clamped` defaults to 0 -- the
    conforming N5 rule (freeze only from a 0-limited/0-floored state), so the
    default path is bit-identical. It exists because at M6 medium M0.70 a SINGLE
    persistently-floored cell of 330k blocks the freeze at ANY freeze_tol, while
    the frozen sweep represents a clamped cell exactly (branch 3)."""
    import inspect
    p = inspect.signature(solve_multivalued_newton).parameters
    assert p["freeze_max_clamped"].default == 0


def test_m6_recipe_arms_freeze_above_the_churn_floor():
    """Errata 2 + 4, locked in the shipped recipe. The M6 recipe MUST use a
    freeze_tol above the measured churn floor (which RISES with Mach: <1e-6 at
    M0.60, 8.6e-6 at M0.65, 2.7e-4 at M0.70) and MUST tolerate the lone clamped
    cell -- with freeze_tol=1e-6 / freeze_max_clamped=0 the ramp DIES at M~0.66
    (measured)."""
    from pyfp3d.solve.newton_ls import B_NEWTON_M6_DEFAULTS as R

    assert R["freeze_tol"] >= 1e-3, (
        "freeze_tol must sit ABOVE the M0.70 churn floor (2.7e-4)")
    assert R["freeze_max_clamped"] >= 1, (
        "a single floored cell otherwise blocks the freeze at any freeze_tol")
