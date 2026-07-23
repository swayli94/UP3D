"""Track V V5 shared case/state builders (GV5.1).

The 2.5-D NACA0012 coarse wake-cut strip configuration and the loose-loop
k=1 state point, factored out of tests/test_v5_tight_jacobian.py so the
Stage-1 (J_phi,BL) and Stage-2 (J_BL,phi, tests/test_v5_tight_edge.py)
gates probe the SAME pre-registered state (the k=1 builder is verbatim
the Stage-1 module's k1_state fixture body; mirrors coupling.py:689-806).
"""

from pathlib import Path

import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import density_field, mach_squared_field
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from pyfp3d.viscous import closures as C
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    _lam_seed,
    _turb_seed,
    build_airfoil_case,
)
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.transpiration import edge_velocity_per_zone

REPO_ROOT = Path(__file__).parent.parent
NACA_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"

M_INF, ALPHA, RE = 0.5, 2.0, 3.0e6
# the GV3.1 Newton-leg settings (cases/analysis/v3_loose_coupling/run.py:67)
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05
CASE_ARGS = dict(
    upwind_c=UPWIND_C,
    m_crit=M_CRIT,
    m_cap=M_CAP,
    rho_floor=RHO_FLOOR,
    tol_residual=1e-10,
)


def build_naca_case():
    """The 2.5-D coarse strip + IBL case wiring (the test_v3_coupling.py
    fixture, lines 56-63)."""
    mc, wc = cut_wake(read_mesh(str(NACA_DIR / "coarse.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg
    )
    return mc, wc, cfg, case


def build_k1_state(naca_case):
    """The loose-loop k=1 state point (coupling.py:689-806 mirrored
    step-by-step, seeds included via coupling.py's own private helpers so
    the mirror cannot drift): inviscid-converged (phi, gamma) + ONE IBL
    Newton solve + the closure packet. The PRE_REGISTRATION risks section
    records the IBL floor on harsh k=1 states -- the solve's convergence
    is recorded, not asserted (the FD gate needs a smooth state point of
    the closure map, not a converged one)."""
    mc, wc, cfg, case = naca_case
    sm = case.sm
    mu = 1.0 / cfg.re_chord
    n_cut = len(mc.nodes)

    # inviscid baseline (coupling.py:689; the GV3.1 Newton recipe)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **CASE_ARGS)
    assert r["converged"]
    phi = np.asarray(r["phi"], dtype=np.float64)
    gamma = np.asarray(r["gamma"], dtype=np.float64)

    # u_e recovery + edge data (coupling.py:710-726)
    le_mask_vol = np.zeros(n_cut, dtype=bool)
    le_mask_vol[sm.volume_node_of[case.le_band_surf]] = True
    ue_vol = edge_velocity_per_zone(
        mc.nodes,
        case.wall_faces,
        phi,
        elements=case.elements,
        le_band_mask=le_mask_vol,
        n_smooth_passes=cfg.n_smooth_passes,
    )
    ue_surf = ue_vol[sm.volume_node_of]
    assert np.all(np.isfinite(ue_surf))
    q2 = np.sum(ue_surf ** 2, axis=1)
    q = np.sqrt(q2)
    rho_e = density_field(q2, cfg.m_inf, cfg.gamma_air)
    mach_e = np.sqrt(mach_squared_field(q2, cfg.m_inf, cfg.gamma_air))

    # inflow Dirichlet band, frozen at this state (coupling.py:738-763,
    # airfoil branch: stations exist, seed_kind None -> all laminar seeds)
    xc_n = case.stations.xc[case.stations.station_of]
    inflow_mask = xc_n <= cfg.inflow_band_x
    idx_in = np.where(inflow_mask)[0]
    inflow_state = np.stack([
        _lam_seed(case.seed_fetch[i], max(float(q[i]), 1.0e-8), 1.0, mu)
        for i in idx_in
    ])

    # U0 seed (coupling.py:784-792)
    q_floor = 0.02 * max(float(np.max(q)), 1.0e-12)
    U0 = np.zeros((sm.n_node, 6), dtype=np.float64)
    for i in range(sm.n_node):
        qq = max(q[i], q_floor)
        if case.turbulent_flags[i]:
            U0[i] = _turb_seed(case.seed_fetch[i], qq, 1.0, mu)
        else:
            U0[i] = _lam_seed(case.seed_fetch[i], qq, 1.0, mu)

    # ONE IBL solve (coupling.py:772-795)
    solver = IBL3Solver(
        sm,
        ue_surf,
        rho_e,
        mu,
        mach_e,
        case.turbulent_flags,
        inflow_mask,
        inflow_state,
        eps_diff=cfg.eps_diff,
        eps_diff_s=cfg.eps_diff_s,
    )
    U, ibl_info = solver.solve(U0, tol=cfg.ibl_tol, max_iter=cfg.ibl_max_iter)
    assert np.all(np.isfinite(U))

    # closure packet (coupling.py:800-806)
    outs = np.empty((sm.n_node, C.N_OUT), dtype=np.float64)
    douts = np.empty((sm.n_node, C.N_OUT, 6), dtype=np.float64)
    douts_e = np.empty((sm.n_node, C.N_OUT, 2), dtype=np.float64)
    C.closure_all(
        U, q, rho_e, np.full(sm.n_node, mu), mach_e,
        case.turbulent_flags, C.C_L_DEFAULT, outs, douts, douts_e,
    )

    # the workspace the FD oracle runs on; external_rhs is set
    # post-construction by the tests (the delta* = 0 identity test first)
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(M_INF)
    phi_free = phi[: ws.n_red][ws.free].copy()
    return {
        "mc": mc, "wc": wc, "cfg": cfg, "case": case, "sm": sm,
        "mu": mu, "n_cut": n_cut, "phi": phi, "gamma": gamma,
        "le_mask_vol": le_mask_vol, "ue_surf": ue_surf, "q": q,
        "rho_e": rho_e, "mach_e": mach_e, "U": U, "outs": outs,
        "douts": douts, "douts_e": douts_e, "ibl_info": ibl_info,
        "ws": ws, "phi_free": phi_free, "solver": solver,
        "inflow_mask": inflow_mask, "inflow_state": inflow_state,
    }
