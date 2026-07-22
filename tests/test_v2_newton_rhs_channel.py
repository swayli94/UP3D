"""Track V V2 RHS channels in the two Newton drivers (binding:
docs/roadmap/track_v.md "V2 -- Transpiration channel through all three
drivers", gates GV2.1(b)/(c)).

Covers, on the NACA0012 coarse 2.5D case (M0.5, alpha 2; subcritical, so
the walk selection is inert and the raw residual is FD-safe -- the
test_p8 protocol):

  - GV2.1(b) conforming Newton: external_rhs=zeros is bit-identical to
    the channel-absent default (phi, gamma, residual history);
  - GV2.1(c) Jacobian stays EXACT under a lagged external_rhs:
    (i) the assembled coupled Jacobian at the same state is bit-identical
        with and without b_ext (b_ext never enters assemble_coupled);
    (ii) the FD oracle of eval_residual -- which INCLUDES the b_ext
        subtraction -- still matches the analytic J_ff (random directions)
        and B (Gamma columns) to the test_p8/test_b31 tolerances;
    (iii) the algebraic identity R'_free = R_free - (T^T b_ext)[free];
  - a nonzero transpiration load solves consistently (driver converges to
    T^T (R - b) = 0 with a shifted solution);
  - GV2.1(b) LS leg: wall_rhs=zeros into the existing b_base slot is
    bit-identical (phi_ext), and a nonzero load is live and converges.

Runs in both lanes: default JIT and PYFP3D_NOJIT=1.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.viscous.transpiration import assemble_transpiration_rhs
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
NACA_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"

M_INF = 0.5
ALPHA = 2.0
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05
CASE_ARGS = dict(upwind_c=UPWIND_C, m_crit=M_CRIT, m_cap=M_CAP,
                 rho_floor=RHO_FLOOR)


def _wall_load(nodes, wall_faces, m_dot_value):
    """Manufactured transpiration load: uniform blowing m_dot_value at the
    wall nodes (the viscous/transpiration.py channel input)."""
    mdot = np.zeros(len(nodes))
    mdot[np.unique(wall_faces)] = m_dot_value
    return assemble_transpiration_rhs(nodes, wall_faces, mdot)


# ---------------------------------------------------------------------------
# conforming Newton (solve/newton.py external_rhs)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def naca_coarse_cut():
    mesh = read_mesh(NACA_DIR / "coarse.msh")
    return cut_wake(mesh)


@pytest.fixture(scope="module")
def b_ext(naca_coarse_cut):
    mc, _ = naca_coarse_cut
    return _wall_load(mc.nodes, mc.boundary_faces["wall"], 0.01)


@pytest.fixture(scope="module")
def baseline(naca_coarse_cut):
    mc, wc = naca_coarse_cut
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             **CASE_ARGS)
    assert r["converged"]
    return r


def _state_of(ws, r):
    phi_free = np.asarray(r["phi"])[: ws.n_red][ws.free].copy()
    gamma = np.asarray(r["gamma"], dtype=np.float64).copy()
    return phi_free, gamma


def test_newton_zero_external_rhs_bit_identical(naca_coarse_cut, baseline):
    """GV2.1(b): external_rhs=zeros is bit-identical to the absent
    channel."""
    mc, wc = naca_coarse_cut
    r0 = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                              external_rhs=np.zeros(len(mc.nodes)),
                              **CASE_ARGS)
    assert np.array_equal(baseline["phi"], r0["phi"])
    assert np.array_equal(baseline["gamma"], r0["gamma"])
    assert baseline["residual_history"] == r0["residual_history"]


def test_newton_external_rhs_jacobian_invariant(naca_coarse_cut, b_ext,
                                                baseline):
    """GV2.1(c) structural clause: at the same state the coupled Jacobian
    assembled with and without the lagged b_ext is BIT-IDENTICAL (b_ext
    enters only the residual's RHS, never assemble_coupled)."""
    mc, wc = naca_coarse_cut
    ws0 = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws0.set_mach(M_INF)
    wsb = NewtonWorkspace(mc, wc, alpha_deg=ALPHA, external_rhs=b_ext)
    wsb.set_mach(M_INF)
    phi_free, gamma = _state_of(ws0, baseline)
    _, _, state0 = ws0.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                     M_CAP, RHO_FLOOR)
    _, _, stateb = wsb.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                     M_CAP, RHO_FLOOR)
    J0, B0 = ws0.assemble_coupled(state0, UPWIND_C, M_CRIT, RHO_FLOOR)
    Jb, Bb = wsb.assemble_coupled(stateb, UPWIND_C, M_CRIT, RHO_FLOOR)
    assert np.array_equal(J0.data, Jb.data)
    assert np.array_equal(B0.data, Bb.data)


def test_newton_external_rhs_residual_identity(naca_coarse_cut, b_ext,
                                               baseline):
    """GV2.1(c) algebra: R'_free = R_free - (T^T b_ext)[free] at the same
    state (the channel subtracts b_ext inside eval_residual, nothing
    else)."""
    mc, wc = naca_coarse_cut
    ws0 = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws0.set_mach(M_INF)
    wsb = NewtonWorkspace(mc, wc, alpha_deg=ALPHA, external_rhs=b_ext)
    wsb.set_mach(M_INF)
    phi_free, gamma = _state_of(ws0, baseline)
    R0, F0, _ = ws0.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                  M_CAP, RHO_FLOOR)
    Rb, Fb, _ = wsb.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                  M_CAP, RHO_FLOOR)
    expected = R0 - (ws0.con.T.T @ b_ext)[ws0.free]
    assert np.max(np.abs(Rb - expected)) < 1e-12
    assert np.array_equal(F0, Fb)


def test_newton_external_rhs_jacobian_fd(naca_coarse_cut, b_ext, baseline):
    """GV2.1(c) FD clause: with the nonzero lagged b_ext ACTIVE, the
    analytic coupled Jacobian still matches central FD of eval_residual
    -- J_ff on random phi directions and B on the Gamma columns (the
    test_p8 gamma-column protocol, rel err < 1e-6)."""
    mc, wc = naca_coarse_cut
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA, external_rhs=b_ext)
    ws.set_mach(M_INF)
    phi_free, gamma = _state_of(ws, baseline)
    _, _, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT,
                                   M_CAP, RHO_FLOOR)
    J_ff, B = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)

    rng = np.random.default_rng(11)
    eps = 1e-5
    for _ in range(3):
        delta = rng.standard_normal(ws.n_free)
        delta /= np.abs(delta).max()
        R_p, _, _ = ws.eval_residual(phi_free + eps * delta, gamma,
                                     UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        R_m, _, _ = ws.eval_residual(phi_free - eps * delta, gamma,
                                     UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
        fd = (R_p - R_m) / (2.0 * eps)
        exact = J_ff @ delta
        scale = max(np.abs(fd).max(), 1e-30)
        rel = np.abs(exact - fd).max() / scale
        assert rel < 1e-6, f"J_ff FD rel err {rel:.3e}"

    for j in range(ws.n_st):
        dg = np.zeros(ws.n_st)
        dg[j] = eps
        R_p, _, _ = ws.eval_residual(phi_free, gamma + dg, UPWIND_C,
                                     M_CRIT, M_CAP, RHO_FLOOR)
        R_m, _, _ = ws.eval_residual(phi_free, gamma - dg, UPWIND_C,
                                     M_CRIT, M_CAP, RHO_FLOOR)
        fd = (R_p - R_m) / (2.0 * eps)
        col = np.asarray(B[:, j].todense()).ravel()
        scale = max(np.abs(fd).max(), 1e-30)
        rel = np.abs(col - fd).max() / scale
        assert rel < 1e-6, f"Gamma column {j}: rel err {rel:.3e}"


def test_newton_external_rhs_solves_consistently(naca_coarse_cut, b_ext,
                                                 baseline):
    """The nonzero transpiration load is a consistent RHS: the driver
    converges to T^T (R - b_ext) = 0 with a shifted (phi, gamma)."""
    mc, wc = naca_coarse_cut
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             external_rhs=b_ext, **CASE_ARGS)
    assert r["converged"]
    assert r["residual_history"][-1] < 1e-10
    assert not np.array_equal(baseline["phi"], r["phi"])
    assert np.max(np.abs(baseline["gamma"] - r["gamma"])) > 1e-8


def test_newton_external_rhs_forwarded_with_workspace_guard(naca_coarse_cut,
                                                            b_ext):
    """Plumbing guard: an externally built workspace owns its b_ext;
    forwarding a DIFFERENT vector alongside it raises (the kutta_estimator
    guard precedent)."""
    mc, wc = naca_coarse_cut
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA, external_rhs=b_ext)
    ws.set_mach(M_INF)
    with pytest.raises(ValueError, match="external_rhs"):
        solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             workspace=ws,
                             external_rhs=np.zeros(len(mc.nodes)),
                             **CASE_ARGS)


# ---------------------------------------------------------------------------
# LS Newton (solve/newton_ls.py b_base slot)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ls_case():
    mesh = read_mesh(NACA_DIR / "coarse.msh")
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return mesh, mvop


LS_ARGS = dict(m_inf=0.3, alpha_deg=2.0, farfield="neumann",
               n_seed=10, n_newton_max=15)


def test_ls_wall_rhs_zero_bit_identical(ls_case):
    """GV2.1(b) LS leg: wall_rhs=zeros into the b_base slot is
    bit-identical to the channel-absent default."""
    mesh, mvop = ls_case
    a = solve_multivalued_newton(mvop=mvop, mesh=mesh, **LS_ARGS)
    b = solve_multivalued_newton(mvop=mvop, mesh=mesh,
                                 wall_rhs=np.zeros(len(mesh.nodes)),
                                 **LS_ARGS)
    assert np.array_equal(a["phi_ext"], b["phi_ext"])
    assert a["residual_history"] == b["residual_history"]


def test_ls_wall_rhs_nonzero_live(ls_case):
    """A nonzero wall_rhs rides b_base into the solve: shifted solution,
    still converged."""
    mesh, mvop = ls_case
    rhs = _wall_load(mesh.nodes, mesh.boundary_faces["wall"], 0.005)
    a = solve_multivalued_newton(mvop=mvop, mesh=mesh, **LS_ARGS)
    c = solve_multivalued_newton(mvop=mvop, mesh=mesh, wall_rhs=rhs,
                                 **LS_ARGS)
    assert c["converged"]
    assert c["residual_history"][-1] < 1e-8
    assert not np.array_equal(a["phi_ext"], c["phi_ext"])
