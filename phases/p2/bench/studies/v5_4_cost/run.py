"""GV5.4 augmented-step cost on M6 medium with the block preconditioner.

Binding text: phases/p1/docs/roadmap/track_v.md GV5.4 ("augmented step wall-time
<= ~2x the inviscid Newton step on M6 medium with the block
preconditioner working; measured number recorded either way");
pre-registered: bench/studies/v5_4_cost/PRE_REGISTRATION.md (committed
a0c2a5b BEFORE the first code change).

  ONERA M6, TEST 2308 verbatim (M 0.8395 / alpha 3.06, Re_MAC 11.72e6,
  x_tr/c 0.05, tip band z > 0.95*b_semi). The GV5.1 augmented
  (phi, Gamma, U) system on the GV5.0 wing case (~124k DOFs medium),
  PROBE Kutta (D1: the only wired F_Gamma row; identical size/sparsity;
  the inviscid anchor measured on the same probe branch). Seed = the
  newton_tight semantics verbatim (D2): the A1 conf_newton chain
  VERBATIM (addendum 2026-07-25 #3: ONE solve_newton_transonic with
  NEWTON_M6_RECIPE, the ramp Picard-seeding level 0, the estimator
  default = probe -- the committed probe anchor belongs to this chain)
  + ONE standalone IBL solve (its ~1e-6 floor expected). N = 5 augmented
  Newton steps, the linear step swapped to GMRES with the
  block-preconditioner ladder on the row+col equilibrated system (D4):
  rung 1 block-Jacobi (AMG-phi / exact-Gamma / ILU-BL), rung 2 exact-BL
  Schur (the B14 pattern; escalation on GMRES non-convergence or the
  1800 s linsolve cap). The inviscid anchor = the seed ramp's
  final-level step_records in-session (D6).

  (a) cost ratio RECORDED vs <= ~2x (either way); (b) the block
  preconditioner WORKING adjudicated (D5: GMRES info=0 within budget on
  every measured step of the rung + every step accepted + W3 passed);
  (c) diagnostics RECORDED. Wiring guards W1 (probe seed cl_p vs the
  committed P14 probe lock 0.2646, 1.5% -- addendum #4) / W2 (DOF
  counts) / W3 (60-row sampled FD ladder) raise = recipe error, not
  verdicts. Exits 1 iff (b) FAIL.

Run:  python bench/studies/v5_4_cost/run.py [--levels coarse medium]
"""

import argparse
import os
import resource
import sys
import time

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le
from pyfp3d.physics.isentropic import density_field, mach_squared_field
from pyfp3d.post.surface import planform_area, wall_force_coefficients
from pyfp3d.solve.linear import (
    build_amg_preconditioner,
    build_ilu_preconditioner,
    solve_gmres,
)
from pyfp3d.solve.newton import (
    NewtonWorkspace,
    solve_newton_transonic,
)
from pyfp3d.viscous import tight_driver as td
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    _lam_seed,
    _turb_seed,
    build_wing_case,
)
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.transpiration import edge_velocity_per_zone

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..")))
from tests.test_p8_newton import NEWTON_M6_RECIPE  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
M6_DIR = os.path.join(REPO, "cases", "meshes", "onera_m6")

M_INF, ALPHA = 0.8395, 3.06            # TEST 2308 dataset label verbatim
MAC = 0.64607
RE_MAC = 11.72e6
RE_CHORD = RE_MAC / MAC                # per meter (meshes in NASA meters)
X_TR = 0.05                            # forced transition, both sides
TIP_FRAC = 0.05                        # tip mask = production tip_taper radius

# committed probe-branch anchor (addendum 2026-07-25 #4: the P14
# cross-model table's "conforming probe (G8.2 lock)" row at M0.84 medium,
# cross_model_medium_m084.csv -- the most recent committed probe reading;
# the A1-era 0.26918 is stale, superseded. The 1.5% guard absorbs the
# dM = 0.0005 label difference + the anchor's 4-digit precision)
PROBE_LOCK_CL_P = 0.2646

N_STEPS = 5                            # pre-registered measured steps
GMRES_RTOL, GMRES_RESTART, GMRES_MAXITER = 1.0e-8, 60, 5
LINSOLVE_CAP = 1800.0                  # s; a trip escalates the ladder

SUMMARY = []  # (gate, metric, band, measured, verdict)


def _record(gate, metric, band, measured, ok=None):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append((gate, metric, band, measured, verdict))
    print(f"  [{verdict:8s}] {gate} {metric}: measured={measured} "
          f"(band: {band})", flush=True)


def _write_csv(name, header, rows):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"  wrote {path}", flush=True)


# ---------------------------------------------------------------------------
# the probe-Kutta inviscid seed (D1/D2): the GV5.3 cold chain with the
# estimator default = probe at every level; NO early return on a
# non-converged seed (the GV5.3 addendum-#1 lesson)
# ---------------------------------------------------------------------------

def make_probe_seed(mc, wc):
    """The A1 conf_newton chain VERBATIM (addendum 2026-07-25 #3): ONE
    solve_newton_transonic with NEWTON_M6_RECIPE -- the ramp
    Picard-seeds level 0 itself (no separate M0.70 Newton solve, no
    phi_init handoff); the estimator default = probe at every level.
    The committed probe anchor (W1) belongs to THIS chain."""
    t0 = time.perf_counter()
    r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                               **NEWTON_M6_RECIPE)
    print(f"    A1 probe ramp: converged={r['converged']} levels="
          f"{[(lv[0], lv[1]) for lv in r['level_history']]} "
          f"({time.perf_counter() - t0:.0f}s)", flush=True)
    return r


# ---------------------------------------------------------------------------
# the tight-pack state at the k=0 inviscid state (D2): the wing branch of
# the GV5.1 state builder -- the closed-body inflow mask + per-node regime
# seeds (coupling.py:784-835 verbatim), then ONE standalone IBL solve
# ---------------------------------------------------------------------------

def build_wing_k0_state(mc, wc, cfg, case, phi, gamma):
    sm = case.sm
    mu = 1.0 / cfg.re_chord
    n_cut = len(mc.nodes)
    phi = np.asarray(phi, dtype=np.float64)
    gamma = np.asarray(gamma, dtype=np.float64)

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

    # the closed-body inflow branch (coupling.py:784-805 verbatim): the
    # LE-band candidates | the tip pin band, per-node regime seeds
    inflow_mask = case.inflow_candidates.copy()
    inflow_mask |= case.outflow_pin_surf
    idx_in = np.where(inflow_mask)[0]
    q_ref = float(np.percentile(q[idx_in], 75))
    q_floor = 0.05 * max(q_ref, 1.0e-12)
    inflow_state = np.stack([
        (_turb_seed(case.seed_fetch[i], max(float(q[i]), q_floor), 1.0, mu)
         if case.turbulent_flags[i] else
         _lam_seed(case.seed_fetch[i], max(float(q[i]), 1.0e-8), 1.0, mu))
        for i in idx_in
    ])

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

    # the per-node seed field (coupling.py:829-835 verbatim), then ONE
    # standalone IBL solve (D2; converged=False at the committed ~1e-6
    # floor EXPECTED)
    q_floor0 = 0.02 * max(float(np.max(q)), 1.0e-12)
    U0 = np.zeros((sm.n_node, 6), dtype=np.float64)
    for i in range(sm.n_node):
        qq = max(float(q[i]), q_floor0)
        if case.turbulent_flags[i]:
            U0[i] = _turb_seed(case.seed_fetch[i], qq, 1.0, mu)
        else:
            U0[i] = _lam_seed(case.seed_fetch[i], qq, 1.0, mu)
    t0 = time.perf_counter()
    U, ibl_info = solver.solve(U0)
    ibl_info["wall_s"] = time.perf_counter() - t0
    print(f"    standalone IBL solve: converged={ibl_info['converged']} "
          f"n_iter={ibl_info['n_iter']} final={ibl_info['final_residual']:.3e} "
          f"({ibl_info['wall_s']:.0f}s)", flush=True)

    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA,
                         farfield_spanwise_gamma=True,
                         kutta_estimator="probe")
    ws.set_mach(M_INF)
    phi_free = phi[: ws.n_red][ws.free].copy()
    return {
        "mc": mc, "wc": wc, "cfg": cfg, "case": case, "sm": sm,
        "mu": mu, "n_cut": n_cut, "phi": phi, "gamma": gamma,
        "le_mask_vol": le_mask_vol, "ue_surf": ue_surf, "q": q,
        "rho_e": rho_e, "mach_e": mach_e,
        "U": np.ascontiguousarray(U, dtype=np.float64),
        "ws": ws, "phi_free": phi_free, "solver": solver,
    }, ibl_info


# ---------------------------------------------------------------------------
# the block-preconditioner ladder (D4): one callback; rung 1 block-Jacobi,
# rung 2 exact-BL Schur; escalation on GMRES non-convergence, a setup
# failure, or the linsolve wall cap
# ---------------------------------------------------------------------------

class GmresLadder:
    def __init__(self, pack):
        self.pack = pack
        self.rung = 1
        self.step = 0
        self.rows = []        # per GMRES call
        self.escalations = []  # (step, reason)

    def _escalate(self, reason):
        if self.rung == 1:
            print(f"    [ladder] step {self.step}: rung 1 -> 2 ({reason})",
                  flush=True)
            self.escalations.append((self.step, reason))
            self.rung = 2

    def __call__(self, A_csc, b):
        self.step += 1
        n_f, n_g = self.pack.n_free, self.pack.n_st
        n_b = 6 * self.pack.n_s
        nh = n_f + n_g
        A = A_csc.tocsr()
        err = ""
        try:
            if self.rung == 1:
                x, n_inner, info, t_setup, t_lin, notes = \
                    self._rung1(A, b, n_f, n_g)
            else:
                x, n_inner, info, t_setup, t_lin, notes = \
                    self._rung2(A, b, nh)
        except Exception as exc:  # a SETUP failure: escalate and retry
            err = f"{type(exc).__name__}: {exc}"
            x, n_inner, info, t_setup, t_lin, notes = (None, -1, -1,
                                                       0.0, 0.0, err)
        if x is None and self.rung == 1:
            self._escalate(f"setup failed ({err})")
            x, n_inner, info, t_setup, t_lin, notes = \
                self._rung2(A, b, nh)
        rel_res = float(np.linalg.norm(b - A @ x)
                        / max(np.linalg.norm(b), 1.0e-300))
        cap = t_lin > LINSOLVE_CAP
        self.rows.append({
            "step": self.step, "rung": self.rung, "t_setup": t_setup,
            "t_linsolve": t_lin, "n_inner": n_inner, "info": info,
            "rel_res": rel_res, "cap_tripped": cap, "notes": notes,
        })
        print(f"    [gmres] step {self.step} rung {self.rung}: "
              f"setup {t_setup:.1f}s + linsolve {t_lin:.1f}s, "
              f"n_inner={n_inner} info={info} rel_res={rel_res:.2e}",
              flush=True)
        if info != 0:
            self._escalate(f"GMRES info={info}")
        if cap:
            self._escalate(f"linsolve {t_lin:.0f}s > {LINSOLVE_CAP:.0f}s")
        return x

    def _rung1(self, A, b, n_f, n_g):
        """Block-Jacobi: AMG on the (phi,phi) block, the exact inverse of
        the (Gamma,Gamma) block, ILU on the (BL,BL) block."""
        t0 = time.perf_counter()
        Jff = A[:n_f, :n_f]
        Jgg = A[n_f:n_f + n_g, n_f:n_f + n_g].toarray()
        Jbb = A[n_f + n_g:, n_f + n_g:]
        _, M_phi = build_amg_preconditioner(Jff)
        t_amg = time.perf_counter() - t0
        M_gam = np.linalg.inv(Jgg)
        t1 = time.perf_counter()
        M_bl = build_ilu_preconditioner(Jbb)
        t_ilu = time.perf_counter() - t1
        t_setup = time.perf_counter() - t0

        def pmv(x):
            y = np.empty_like(x)
            y[:n_f] = M_phi @ x[:n_f]
            y[n_f:n_f + n_g] = M_gam @ x[n_f:n_f + n_g]
            y[n_f + n_g:] = M_bl @ x[n_f + n_g:]
            return y

        M = spla.LinearOperator(A.shape, matvec=pmv)
        t2 = time.perf_counter()
        x, n_inner, info = solve_gmres(
            A, b, M=M, rtol=GMRES_RTOL, restart=GMRES_RESTART,
            maxiter=GMRES_MAXITER, on_fail="return")
        t_lin = time.perf_counter() - t2
        return x, n_inner, info, t_setup, t_lin, \
            f"t_amg={t_amg:.1f},t_ilu={t_ilu:.1f}"

    def _rung2(self, A, b, nh):
        """Exact-BL Schur (the B14 pattern): splu(J_BL,BL) once (the
        'measure before Schur' number), the matrix-free reduced operator
        on (phi,Gamma), AMG-phi/exact-Gamma on the reduced
        preconditioner, exact back-substitution."""
        n_f = self.pack.n_free
        t0 = time.perf_counter()
        Jhh = A[:nh, :nh]
        Jhb = A[:nh, nh:]
        Jbh = A[nh:, :nh]
        lu = spla.splu(A[nh:, nh:].tocsc())
        t_lu = time.perf_counter() - t0
        _, M_phi = build_amg_preconditioner(Jhh[:n_f, :n_f])
        M_gam = np.linalg.inv(Jhh[n_f:, n_f:].toarray())
        t_setup = time.perf_counter() - t0

        def mv(xh):
            return Jhh @ xh - Jhb @ lu.solve(Jbh @ xh)

        def pmv(xh):
            y = np.empty_like(xh)
            y[:n_f] = M_phi @ xh[:n_f]
            y[n_f:] = M_gam @ xh[n_f:]
            return y

        K = spla.LinearOperator((nh, nh), matvec=mv)
        M = spla.LinearOperator((nh, nh), matvec=pmv)
        rhs = b[:nh] - Jhb @ lu.solve(b[nh:])
        t2 = time.perf_counter()
        xh, n_inner, info = solve_gmres(
            K, rhs, M=M, rtol=GMRES_RTOL, restart=GMRES_RESTART,
            maxiter=GMRES_MAXITER, on_fail="return")
        t_lin = time.perf_counter() - t2
        x = np.concatenate([xh, lu.solve(b[nh:] - Jbh @ xh)])
        return x, n_inner, info, t_setup, t_lin, f"t_lu={t_lu:.1f}"


# ---------------------------------------------------------------------------
# the W3 sampled FD guard (the GV5.1 fd_gate ladder idiom on a 60-row
# sample: 24 phi (8 wake-adjacent + 16 spread), 6 Gamma (incl. the tip
# station), 30 BL (5 nodes x 6: tip-pinned / LE-band / max-q / min-q /
# mid); per-block median relative error < 1e-6)
# ---------------------------------------------------------------------------

def fd_guard(pack, mc, case, q, x):
    ws = pack.ws
    wc = ws.wc
    n_f, n_g, n_b = pack.n_free, pack.n_st, 6 * pack.n_s
    n = n_f + n_g + n_b

    # -- the row sample -------------------------------------------------------
    slave = np.asarray(wc.slave_nodes, dtype=np.int64)
    mark = np.isin(mc.elements, slave).any(axis=1)
    near_cut = np.union1d(np.unique(mc.elements[mark]), slave)
    touched = np.unique(ws.con.T[near_cut].indices)
    wake_free = touched[np.isin(touched, ws.free)]
    take = np.linspace(0, len(wake_free) - 1,
                       min(8, len(wake_free))).astype(int)
    phi_rows = list(wake_free[take])
    phi_rows += list(ws.free[np.linspace(0, n_f - 1, 16).astype(int)])
    gam_rows = sorted(set(min(g, n_g - 1)
                          for g in (0, 33, 66, 99, 132, n_g - 1)))
    pinned = np.where(case.outflow_pin_surf)[0]
    notpin = np.where(~case.outflow_pin_surf)[0]
    nodes = [int(pinned[0]), int(np.where(case.le_band_surf)[0][0]),
             int(np.argmax(q)),
             int(notpin[np.argmin(q[notpin])]), int(pack.n_s // 2)]
    bl_rows = sorted(6 * i + k for i in nodes for k in range(6))
    rows = (np.array(phi_rows),
            n_f + np.array(gam_rows),
            n_f + n_g + np.array(bl_rows))

    # the GV5.1 fallback row-mask idiom: drop the 6 rows of every
    # q <= 1e-12 fallback node and of every element-mate
    bad = np.where(q <= 1.0e-12)[0]
    keep_bl = np.ones(n_b, dtype=bool)
    if bad.size:
        drop = set(int(i) for i in bad)
        for e in range(len(case.sm.triangles)):
            if np.isin(case.sm.triangles[e], bad).any():
                drop.update(int(i) for i in case.sm.triangles[e])
        for i in drop:
            keep_bl[6 * i: 6 * i + 6] = False
    dropped = int((~keep_bl[rows[2] - n_f - n_g]).sum())
    rows = (rows[0], rows[1],
            rows[2][keep_bl[rows[2] - n_f - n_g]])
    print(f"    W3 sample: {len(rows[0])} phi + {len(rows[1])} gamma + "
          f"{len(rows[2])} BL rows ({dropped} fallback-dropped)",
          flush=True)

    J = td.augmented_jacobian(pack, x)
    rng = np.random.default_rng(97)
    ndir = 6
    out_rows = []
    errs = {"phi": [], "gamma": [], "BL": []}
    for d in range(ndir):
        v = rng.standard_normal(n)
        v /= np.max(np.abs(v))
        an = J @ v
        best = np.full(n, np.inf)
        ladder = []
        for e in (1.0e-5, 1.0e-6, 1.0e-7):
            ff = (td.augmented_residual(pack, x + e * v)
                  - td.augmented_residual(pack, x - e * v)) / (2.0 * e)
            dev = np.abs(ff - an)
            best = np.minimum(best, dev)
            ladder.append(float(np.max(dev)))
        for name, r in (("phi", rows[0]), ("gamma", rows[1]),
                        ("BL", rows[2])):
            # the GV5.1 scaled max-norm idiom, restricted to the sample:
            # the deviation scaled by the block's matvec scale
            scale = max(float(np.max(np.abs(an[r]))), 1.0e-300)
            rel = best[r] / scale
            errs[name].append(float(np.median(rel)))
        print(f"    W3 dir {d}: ladder max {ladder[0]:.2e}/"
              f"{ladder[1]:.2e}/{ladder[2]:.2e}; sampled medians "
              f"phi {errs['phi'][-1]:.2e} gamma {errs['gamma'][-1]:.2e} "
              f"BL {errs['BL'][-1]:.2e}", flush=True)
        out_rows.append((d, f"{ladder[0]:.3e}", f"{ladder[1]:.3e}",
                         f"{ladder[2]:.3e}", f"{errs['phi'][-1]:.3e}",
                         f"{errs['gamma'][-1]:.3e}",
                         f"{errs['BL'][-1]:.3e}", dropped))
    med = {k: float(np.median(v)) for k, v in errs.items()}
    for k, val in med.items():
        if val > 1.0e-6:
            raise RuntimeError(
                f"W3 wiring guard: the {k} block's median relative FD "
                f"error {val:.3e} > 1e-6 on the wing augmented J -- "
                "recipe error, not a verdict")
    print(f"    [W3 ok] per-block median rel FD error: phi "
          f"{med['phi']:.2e}, gamma {med['gamma']:.2e}, BL "
          f"{med['BL']:.2e}", flush=True)
    return out_rows, med, J


# ---------------------------------------------------------------------------
# one level (coarse shakedown / medium binding)
# ---------------------------------------------------------------------------

def run_level(level):
    print(f"--- GV5.4 {level}: ONERA M6 M={M_INF} alpha={ALPHA} "
          f"Re_MAC={RE_MAC:.3e}, probe Kutta (D1) ---", flush=True)
    mc, wc = cut_wake(read_mesh(os.path.join(M6_DIR, f"{level}.msh")))
    wall = mc.boundary_faces["wall"]
    cfg = CouplingConfig(re_chord=RE_CHORD, m_inf=M_INF, alpha_deg=ALPHA,
                         x_tr_upper=X_TR, x_tr_lower=X_TR)
    case = build_wing_case(mc.nodes, mc.elements, wall, cfg,
                           x_le=x_le, chord_at=chord_at,
                           tip_mask_frac=TIP_FRAC)
    s_ref = planform_area(mc.nodes, wall)
    print(f"    IBL surface: {case.sm.n_node} nodes / {case.sm.n_tri} "
          f"tris; tip-masked {int(case.outflow_pin_surf.sum())}",
          flush=True)

    # -- the inviscid probe seed (D1/D2, addendum #3: A1 verbatim) + W1 --
    t_seed0 = time.perf_counter()
    ramp = make_probe_seed(mc, wc)
    t_seed = time.perf_counter() - t_seed0
    f = wall_force_coefficients(mc.nodes, mc.elements, wall, ramp["phi"],
                                alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    cl_p = float(f["cl"])
    if not ramp["converged"]:
        raise RuntimeError(
            f"W1 wiring guard: the {level} probe ramp did not converge "
            "the final level -- recipe error, not a verdict")
    rel = abs(cl_p / PROBE_LOCK_CL_P - 1.0)
    # addendum 2026-07-25 #1: the cl_p tolerance binds MEDIUM only (the
    # A1 anchor is a medium number; the coarse deviation = the committed
    # pressure mesh effect -5.35%, measured here -4.97%); the coarse
    # cl_p is RECORDED with the cross-check quoted
    if level == "medium" and rel > 0.015:
        raise RuntimeError(
            f"W1 wiring guard: {level} probe seed cl_p {cl_p:.6f} vs "
            f"the committed P14 probe lock {PROBE_LOCK_CL_P:.5f} (rel "
            f"{rel:.3%} > 1.5%) -- recipe error, not a verdict")
    tag = (", medium-binding)"
           if level == "medium"
           else ", coarse RECORDED; cf. the committed pressure mesh "
                "effect -5.35%)")
    print(f"    [W1 ok] probe seed cl_p {cl_p:.5f} (anchor "
          f"{PROBE_LOCK_CL_P:.5f}, rel {rel:.3%}{tag}", flush=True)

    # -- the inviscid anchor (D6): the final level's step_records -------------
    sr = ramp["step_records"]
    walls = np.diff([0.0] + [float(r["wall_cum_s"]) for r in sr])
    a_rows = []
    for r, w in zip(sr, walls):
        a_rows.append((int(r["i"]), f"{w:.3f}",
                       *(f"{float(r.get(f't_{p}', 0.0)):.3f}"
                         for p in ("seed", "assembly", "precond",
                                   "linsolve", "residual", "kutta")),
                       int(r["n_lin_iters"]), int(r["n_refactor"]),
                       f"{float(r['residual']):.3e}"))
    _write_csv(f"inviscid_anchor_{level}.csv",
               "i,wall_step,t_seed,t_assembly,t_precond,t_linsolve,"
               "t_residual,t_kutta,n_lin_iters,n_refactor,residual",
               a_rows)
    inv_mean = float(np.mean(walls))
    inv_asm = float(np.mean([float(r.get("t_assembly", 0.0))
                             for r in sr]))
    inv_pre = float(np.mean([float(r.get("t_precond", 0.0)) for r in sr]))
    inv_lin = float(np.mean([float(r.get("t_linsolve", 0.0)) for r in sr]))
    n_refac = int(sum(int(r["n_refactor"]) for r in sr))
    print(f"    inviscid anchor: {len(sr)} final-level steps, mean "
          f"{inv_mean:.2f}s/step (assembly {inv_asm:.2f} + precond "
          f"{inv_pre:.2f} + linsolve {inv_lin:.2f}; {n_refac} "
          f"refactors)", flush=True)

    # -- the wing k=0 state + ONE standalone IBL solve (D2) -------------------
    st, ibl_info = build_wing_k0_state(mc, wc, cfg, case,
                                       ramp["phi"], ramp["gamma"])

    # -- the pack + W2 ---------------------------------------------------------
    pack = td.build_tight_pack(st)
    n_f, n_g, n_b = pack.n_free, pack.n_st, 6 * pack.n_s
    print(f"    augmented system: n_free={n_f} + n_st={n_g} + "
          f"6*n_s={n_b} = {n_f + n_g + n_b} DOFs", flush=True)
    if level == "medium":
        if n_g != 166 or n_b != 61230 or not (60000 <= n_f <= 64000):
            raise RuntimeError(
                f"W2 wiring guard: medium DOFs n_free={n_f} "
                f"(band 60k-64k), n_st={n_g} (=166?), 6*n_s={n_b} "
                "(=61230?) -- recipe error, not a verdict")
    print(f"    [W2 ok] DOF counts recorded", flush=True)

    # -- F at the seed + the assembly/eval timing samples (D7) ----------------
    x0 = pack.x_base()
    t0 = time.perf_counter()
    F0 = td.augmented_residual(pack, x0)
    t_feval = time.perf_counter() - t0
    bm0 = [float(np.max(np.abs(b))) for b in pack.split_F(F0)]
    print(f"    seed |F|: phi {bm0[0]:.3e} gamma {bm0[1]:.3e} BL "
          f"{bm0[2]:.3e} (F eval {t_feval:.1f}s)", flush=True)
    t0 = time.perf_counter()
    td.augmented_jacobian(pack, x0)
    t_asm1 = time.perf_counter() - t0
    print(f"    assembly sample 1: {t_asm1:.1f}s", flush=True)

    # -- W3 FD guard ------------------------------------------------------------
    fd_rows, fd_med, _ = fd_guard(pack, mc, case, st["q"], x0)
    _write_csv(f"fd_guard_{level}.csv",
               "dir,max_1e5,max_1e6,max_1e7,med_phi,med_gamma,med_BL,"
               "fallback_dropped", fd_rows)

    # -- the measurement: N augmented Newton steps with the ladder (D4/D7) -----
    ladder = GmresLadder(pack)
    t0 = time.perf_counter()
    res = td.newton_tight(
        pack, max_iter=N_STEPS, scaling="rowcol", lm_damping=False,
        floor_stop=False, step_solve=ladder, verbose=True)
    t_newton = time.perf_counter() - t0
    t0 = time.perf_counter()
    td.augmented_jacobian(pack, res["x"])
    t_asm2 = time.perf_counter() - t0
    t_asm = 0.5 * (t_asm1 + t_asm2)
    print(f"    newton: termination={res['termination']} "
          f"n_iter={res['n_iter']} ({t_newton:.0f}s); assembly sample 2 "
          f"{t_asm2:.1f}s -> assigned {t_asm:.1f}s/step", flush=True)

    # -- per-step accounting ----------------------------------------------------
    hist = res["history"][1:]
    ws_steps = np.diff([0.0] + [float(h["wall_s"]) for h in
                                res["history"]])
    s_rows = []
    for k, h in enumerate(hist):
        g = ladder.rows[k] if k < len(ladder.rows) else {
            "rung": "", "t_setup": 0.0, "t_linsolve": 0.0, "n_inner": -1,
            "info": -1, "rel_res": float("nan"), "cap_tripped": ""}
        bm = h["block_max"]
        resid_share = float(ws_steps[k + 1]) - t_asm \
            - g["t_setup"] - g["t_linsolve"]
        s_rows.append((int(h["iter"]), g["rung"],
                       f"{float(ws_steps[k + 1]):.2f}", f"{t_asm:.2f}",
                       f"{g['t_setup']:.2f}", f"{g['t_linsolve']:.2f}",
                       f"{resid_share:.2f}", g["n_inner"], g["info"],
                       f"{g['rel_res']:.3e}", f"{float(h['lam']):.4f}",
                       bool(h["accepted"]), f"{h['merit']:.6e}",
                       f"{bm[0]:.3e}", f"{bm[1]:.3e}", f"{bm[2]:.3e}",
                       f"{h['ds_change']:.3e}"))
    _write_csv(f"steps_{level}.csv",
               "iter,rung,wall_step,t_assembly,t_setup,t_linsolve,"
               "t_residual_share,n_inner,info,rel_res,lam,accepted,"
               "merit,f_phi_max,f_gamma_max,f_bl_max,ds_change", s_rows)
    _write_csv(f"gmres_{level}.csv",
               "step,rung,t_setup,t_linsolve,n_inner,info,rel_res,"
               "cap_tripped,notes",
               [(r["step"], r["rung"], f"{r['t_setup']:.2f}",
                 f"{r['t_linsolve']:.2f}", r["n_inner"], r["info"],
                 f"{r['rel_res']:.3e}", r["cap_tripped"], r["notes"])
                for r in ladder.rows])

    # -- the adjudication (D5) + the ratio (D8) ---------------------------------
    aug_mean = float(np.mean(ws_steps[1:len(hist) + 1])) if hist else \
        float("nan")
    accepted = [bool(h["accepted"]) for h in hist]
    working_rung = None
    for rung in (1, 2):
        idx = [k for k, r in enumerate(ladder.rows) if r["rung"] == rung]
        if not idx:
            continue
        gmres_ok = all(ladder.rows[k]["info"] == 0
                       and ladder.rows[k]["t_linsolve"] <= LINSOLVE_CAP
                       for k in idx)
        acc_ok = all(accepted[k] for k in idx if k < len(accepted))
        print(f"    rung {rung}: steps {len(idx)}, gmres_ok={gmres_ok}, "
              f"accepted_all={acc_ok}", flush=True)
        if gmres_ok and acc_ok:
            working_rung = rung
            break
    ratio = aug_mean / inv_mean

    _record("(a)", f"augmented step wall {level}", "<= ~2x the "
            "inviscid step (RECORDED either way)",
            f"{aug_mean:.2f}s / {inv_mean:.2f}s = {ratio:.2f}x")
    # (b) is adjudicated PASS/FAIL on the binding medium level only; the
    # coarse shakedown's reading is RECORDED (pre-registration section 2)
    _record("(b)", f"block preconditioner working {level}",
            "GMRES info=0 within budget on every step + accepted (D5)",
            f"rung {working_rung}" if working_rung else "NOT-WORKING",
            ok=(working_rung is not None) if level == "medium" else None)
    _record("(c)", f"DOFs {level}", "RECORDED",
            f"n_free {n_f} + n_st {n_g} + 6n_s {n_b}")
    _record("(c)", f"seed blocks {level}", "RECORDED",
            f"|F_phi| {bm0[0]:.2e} |F_gam| {bm0[1]:.2e} |F_BL| "
            f"{bm0[2]:.2e}; IBL floor {ibl_info['final_residual']:.2e} "
            f"(converged={ibl_info['converged']})")
    _record("(c)", f"inviscid anchor {level}", "RECORDED",
            f"{len(sr)} steps mean {inv_mean:.2f}s (asm {inv_asm:.2f} "
            f"pre {inv_pre:.2f} lin {inv_lin:.2f}, {n_refac} refactors)")
    _setup_mu = float(np.mean([r["t_setup"] for r in ladder.rows]))
    _lin_mu = float(np.mean([r["t_linsolve"] for r in ladder.rows]))
    _record("(c)", f"augmented phases {level}", "RECORDED",
            f"asm {t_asm:.1f}s + setup {_setup_mu:.1f}s + "
            f"linsolve {_lin_mu:.1f}s per step; termination "
            f"{res['termination']}")
    _record("(c)", f"FD guard medians {level}", "< 1e-6 (W3)",
            f"phi {fd_med['phi']:.2e} gamma {fd_med['gamma']:.2e} "
            f"BL {fd_med['BL']:.2e}")
    _maxrss_gb = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                  / 1.0e6)
    _record("(c)", f"memory {level}", "RECORDED",
            f"ru_maxrss {_maxrss_gb:.1f} GB")
    _record("(c)", f"wall {level} @8 threads", "RECORDED",
            f"seed {t_seed:.0f}s + IBL {ibl_info['wall_s']:.0f}s + "
            f"newton {t_newton:.0f}s (non-comparable)")

    # -- the figure ----------------------------------------------------------
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    it = np.arange(1, len(hist) + 1)
    asm = np.full(len(hist), t_asm)
    setup = np.array([r["t_setup"] for r in ladder.rows])
    lin = np.array([r["t_linsolve"] for r in ladder.rows])
    rest = np.maximum(ws_steps[1:len(hist) + 1] - asm - setup - lin, 0.0)
    ax0.bar(it, asm, label="assembly (sampled)")
    ax0.bar(it, setup, bottom=asm, label="precond setup")
    ax0.bar(it, lin, bottom=asm + setup, label="GMRES")
    ax0.bar(it, rest, bottom=asm + setup + lin, label="residual+backtrack")
    ax0.axhline(2.0 * inv_mean, color="k", ls="--",
                label=f"2x inviscid mean ({2 * inv_mean:.1f}s)")
    ax0.axhline(inv_mean, color="k", ls=":",
                label=f"inviscid mean ({inv_mean:.1f}s)")
    ax0.set_xlabel("augmented Newton step")
    ax0.set_ylabel("wall s @8 threads")
    ax0.set_title(f"GV5.4 {level}: step cost "
                  f"(ratio {ratio:.2f}x, rung {working_rung})")
    ax0.legend(fontsize=7)
    ax1.semilogy(it, [r["n_inner"] for r in ladder.rows], "o-",
                 label="GMRES n_inner")
    ax1.semilogy(it, [max(r["rel_res"], 1e-16) for r in ladder.rows],
                 "s--", label="rel_res")
    ax1.semilogy(it, [h["merit"] for h in hist], "^:", label="merit")
    ax1.set_xlabel("augmented Newton step")
    ax1.set_title("GMRES effort + step quality")
    ax1.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, f"cost_breakdown_{level}.png"),
                dpi=150)
    plt.close(fig)
    print(f"  wrote results/cost_breakdown_{level}.png", flush=True)

    return {"aug_mean": aug_mean, "inv_mean": inv_mean, "ratio": ratio,
            "working_rung": working_rung, "res": res, "ladder": ladder}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse", "medium"],
                    choices=["coarse", "medium"])
    args = ap.parse_args()
    t_all = time.perf_counter()
    out = {}
    for level in args.levels:
        out[level] = run_level(level)
    _write_csv("summary.csv", "gate,metric,band,measured,verdict",
               SUMMARY)
    fail = any(v == "FAIL" for _, _, _, _, v in SUMMARY)
    print(f"--- GV5.4 done in {time.perf_counter() - t_all:.0f}s; "
          f"verdicts: "
          f"{sum(1 for *_, v in SUMMARY if v == 'PASS')} PASS / "
          f"{sum(1 for *_, v in SUMMARY if v == 'FAIL')} FAIL / "
          f"{sum(1 for *_, v in SUMMARY if v == 'RECORDED')} RECORDED ---",
          flush=True)
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
