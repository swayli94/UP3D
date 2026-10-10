"""GV5.6 Schur-aware reduced-space preconditioner for the augmented step.

Binding text: phases/p1/docs/roadmap/track_v.md GV5.6 (the GV5.4 registered
follow-up, opening user-adjudicated 2026-07-25); pre-registered:
cases/analysis/v5_6_schur_prec/PRE_REGISTRATION.md (committed BEFORE the
first code change).

  ONERA M6, TEST 2308 verbatim (M 0.8395 / alpha 3.06, Re_MAC 11.72e6,
  x_tr/c 0.05, tip band z > 0.95*b_semi). The GV5.4 system/seed/protocol
  VERBATIM: the 124,216-DOF augmented (phi, Gamma, U) system on the GV5.0
  wing case, PROBE Kutta, the A1 conf_newton seed chain verbatim + ONE
  standalone IBL solve, scaling="rowcol", mu == 0, N = 5 measured
  Newton steps, GMRES budget rtol 1e-8 / restart 60 / maxiter 5 (<= 300
  it), on_fail="return", linsolve cap 1800 s, no stale reuse, the
  inviscid anchor = the seed ramp's final-level step_records in-session.
  The ladder (D2/D3):

  - rung 3 (binding, first): the GV5.4 rung-2 exact-BL Schur operator
    verbatim with the reduced-space preconditioner bdiag(AMG(S_hat_ff),
    M_Gamma), S_hat_ff = J_ff - C_hat_ff, C_hat = J_hB D_BB^-1 J_Bh
    assembled explicitly (D_BB = the per-node 6x6 block-diagonal of
    J_BL,BL via the (6,6)-BSR view; rcond guard 1e-12 -> zero-safe
    diagonal fallback, the count recorded) -- the quasi-simultaneous
    local BL response inside the AMG matrix.
  - rung 4 (escalation): full-system GMRES with the block
    upper-triangular preconditioner (y_B = lu.solve(r_B); y_h =
    P_hh(r_h - J_hB y_B), P_hh = bdiag(AMG(S_hat_ff), M_Gamma)).

  Escalation = the GV5.4 rule (GMRES info != 0 within budget, a setup
  failure, or the cap). (a) cost ratio RECORDED either way vs <= ~2x;
  (b) the preconditioner WORKING adjudicated (the GV5.4 D5 verbatim:
  GMRES info=0 within budget on every measured step of the rung + every
  step accepted + W3 passed); (c) diagnostics RECORDED. Wiring guards W1
  (probe seed cl_p vs the P14 probe lock 0.2646, 1.5% medium-binding) /
  W2 (DOF counts) / W3 (60-row sampled FD ladder) raise = recipe error,
  not verdicts; W4 = the pyfp3d tree untouched (no library change). Exits
  1 iff (b) FAIL.

Run:  python cases/analysis/v5_6_schur_prec/run.py [--levels coarse medium]
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

# committed probe-branch anchor (the GV5.4 addendum-#4 anchor: the P14
# cross-model table's "conforming probe (G8.2 lock)" row at M0.84 medium,
# cross_model_medium_m084.csv; the 1.5% guard absorbs the dM = 0.0005
# label difference + the anchor's 4-digit precision)
PROBE_LOCK_CL_P = 0.2646

N_STEPS = 5                            # pre-registered measured steps
GMRES_RTOL, GMRES_RESTART, GMRES_MAXITER = 1.0e-8, 60, 5
LINSOLVE_CAP = 1800.0                  # s; a trip escalates the ladder
RCOND_GUARD = 1.0e12                   # cond > this -> diagonal fallback

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
# the probe-Kutta inviscid seed (D1): the A1 conf_newton chain VERBATIM
# ---------------------------------------------------------------------------

def make_probe_seed(mc, wc):
    """The A1 conf_newton chain VERBATIM (GV5.4 addendum 2026-07-25 #3):
    ONE solve_newton_transonic with NEWTON_M6_RECIPE -- the ramp
    Picard-seeds level 0 itself (no separate M0.70 Newton solve, no
    phi_init handoff); the estimator default = probe at every level."""
    t0 = time.perf_counter()
    r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                               **NEWTON_M6_RECIPE)
    print(f"    A1 probe ramp: converged={r['converged']} levels="
          f"{[(lv[0], lv[1]) for lv in r['level_history']]} "
          f"({time.perf_counter() - t0:.0f}s)", flush=True)
    return r


# ---------------------------------------------------------------------------
# the tight-pack state at the k=0 inviscid state (D1 = GV5.4 D2 verbatim)
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
    # standalone IBL solve (converged=False at the committed ~1e-6 floor
    # EXPECTED)
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
# the Schur-aware ladder (D2/D3): rung 3 = exact-BL Schur operator with
# bdiag(AMG(S_hat_ff), M_Gamma), S_hat_ff = J_ff - (J_hB D_BB^-1 J_Bh)_ff;
# rung 4 = full-system GMRES with the block upper-triangular
# preconditioner. Escalation = the GV5.4 rule.
# ---------------------------------------------------------------------------

def _per_node_block_inv(Jbb, n_s):
    """D_BB^-1: the per-node 6x6 block-diagonal of J_BL,BL (the node-major
    6i+k layout) inverted node by node, via the (6,6)-BSR view. A block
    with cond > RCOND_GUARD (or not finite) falls back to its zero-safe
    diagonal inverse; the pinned identity rows invert exactly. Returns
    (D_inv csr, n_fallback)."""
    bsr = Jbb.tocsr().tobsr(blocksize=(6, 6))
    indptr, indices, data = bsr.indptr, bsr.indices, bsr.data
    blocks = np.zeros((n_s, 6, 6), dtype=np.float64)
    for i in range(n_s):
        s, e = indptr[i], indptr[i + 1]
        hits = np.nonzero(indices[s:e] == i)[0]  # NOT sorted-assumed
        if hits.size:
            blocks[i] = data[s + hits[0]]
        # else: a structurally empty diagonal block stays zero -> the
        # diagonal fallback below (cond of a zero block = inf)
    conds = np.linalg.cond(blocks)
    bad = ~np.isfinite(conds) | (conds > RCOND_GUARD)
    safe = np.where(bad[:, None, None],
                    np.broadcast_to(np.eye(6), (n_s, 6, 6)), blocks)
    inv = np.linalg.inv(safe)
    n_fallback = int(bad.sum())
    if n_fallback:
        d = np.diagonal(blocks[bad], axis1=1, axis2=2)
        dz = np.where(d != 0.0, d, 1.0)
        dinv = np.where(d != 0.0, 1.0 / dz, 0.0)
        eye = np.zeros((n_fallback, 6, 6), dtype=np.float64)
        k = np.arange(6)
        eye[:, k, k] = dinv
        inv[bad] = eye
    D_inv = sp.bsr_matrix(
        (inv, np.arange(n_s), np.arange(n_s + 1)), shape=Jbb.shape)
    return D_inv.tocsr(), n_fallback


class SchurLadder:
    def __init__(self, pack):
        self.pack = pack
        self.rung = 3
        self.step = 0
        self.rows = []        # per GMRES call
        self.escalations = []  # (step, reason)

    def _escalate(self, reason):
        if self.rung == 3:
            print(f"    [ladder] step {self.step}: rung 3 -> 4 ({reason})",
                  flush=True)
            self.escalations.append((self.step, reason))
            self.rung = 4

    def _setup(self, A, nh):
        """The shared per-step setup (fresh every step -- no stale reuse):
        splu(J_BB), D_BB^-1 + C_hat (t_corr), AMG(S_hat_ff), M_Gamma."""
        n_f = self.pack.n_free
        t0 = time.perf_counter()
        lu = spla.splu(A[nh:, nh:].tocsc())
        t_lu = time.perf_counter() - t0
        t0 = time.perf_counter()
        D_inv, n_fb = _per_node_block_inv(A[nh:, nh:], self.pack.n_s)
        C_hat = (A[:nh, nh:] @ D_inv) @ A[nh:, :nh]
        t_corr = time.perf_counter() - t0
        t0 = time.perf_counter()
        S_hat_ff = (A[:n_f, :n_f] - C_hat[:n_f, :n_f]).tocsr()
        _, M_phi = build_amg_preconditioner(S_hat_ff)
        M_gam = np.linalg.inv(A[n_f:nh, n_f:nh].toarray())
        t_amg = time.perf_counter() - t0
        t_setup = t_lu + t_corr + t_amg
        notes = (f"t_lu={t_lu:.1f},t_corr={t_corr:.1f},t_amg={t_amg:.1f},"
                 f"nnz_chat={C_hat.nnz},nnz_shat={S_hat_ff.nnz},"
                 f"n_fallback={n_fb}")
        return lu, C_hat, M_phi, M_gam, t_setup, notes

    def __call__(self, A_csc, b):
        self.step += 1
        n_f, n_g = self.pack.n_free, self.pack.n_st
        nh = n_f + n_g
        A = A_csc.tocsr()
        err = ""
        try:
            x, n_inner, info, t_setup, t_lin, notes = (
                self._rung3(A, b, nh) if self.rung == 3
                else self._rung4(A, b, nh))
        except Exception as exc:  # a SETUP failure: escalate and retry
            err = f"{type(exc).__name__}: {exc}"
            x, n_inner, info, t_setup, t_lin, notes = (None, -1, -1,
                                                       0.0, 0.0, err)
        if x is None and self.rung == 3:
            self._escalate(f"setup failed ({err})")
            x, n_inner, info, t_setup, t_lin, notes = \
                self._rung4(A, b, nh)
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

    def _rung3(self, A, b, nh):
        """Rung 3 (D2): the GV5.4 rung-2 exact-BL Schur operator verbatim
        with bdiag(AMG(S_hat_ff), M_Gamma) -- the AMG matrix now carries
        the sparsified Schur correction C_hat = J_hB D_BB^-1 J_Bh."""
        n_f = self.pack.n_free
        Jhh = A[:nh, :nh]
        Jhb = A[:nh, nh:]
        Jbh = A[nh:, :nh]
        lu, _, M_phi, M_gam, t_setup, notes = self._setup(A, nh)

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
        return x, n_inner, info, t_setup, t_lin, notes

    def _rung4(self, A, b, nh):
        """Rung 4 (D3): full-system GMRES with the block upper-triangular
        preconditioner y_B = lu.solve(r_B); y_h = P_hh(r_h - J_hB y_B),
        P_hh = bdiag(AMG(S_hat_ff), M_Gamma)."""
        n_f = self.pack.n_free
        Jhb = A[:nh, nh:]
        lu, _, M_phi, M_gam, t_setup, notes = self._setup(A, nh)

        def pmv(r):
            y_b = lu.solve(r[nh:])
            rh = r[:nh] - Jhb @ y_b
            y_h = np.empty(nh)
            y_h[:n_f] = M_phi @ rh[:n_f]
            y_h[n_f:] = M_gam @ rh[n_f:]
            return np.concatenate([y_h, y_b])

        M = spla.LinearOperator(A.shape, matvec=pmv)
        t2 = time.perf_counter()
        x, n_inner, info = solve_gmres(
            A, b, M=M, rtol=GMRES_RTOL, restart=GMRES_RESTART,
            maxiter=GMRES_MAXITER, on_fail="return")
        t_lin = time.perf_counter() - t2
        return x, n_inner, info, t_setup, t_lin, notes


# ---------------------------------------------------------------------------
# the W3 sampled FD guard (the GV5.4 guard verbatim: 24 phi (8 wake-adjacent
# + 16 spread), 6 Gamma (incl. the tip station), 30 BL (5 nodes x 6); the
# h = 1e-5/1e-6/1e-7 ladder + the fallback row mask; per-block median < 1e-6)
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
    print(f"--- GV5.6 {level}: ONERA M6 M={M_INF} alpha={ALPHA} "
          f"Re_MAC={RE_MAC:.3e}, probe Kutta (GV5.4 D1 verbatim) ---",
          flush=True)
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

    # -- the inviscid probe seed (D1: A1 verbatim) + W1 -----------------------
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
    # medium-binding (GV5.4 addendum #1); coarse recorded with the
    # committed -5.35% mesh-effect cross-check
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

    # -- the inviscid anchor (the GV5.4 D6 protocol) ----------------------------
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

    # -- the wing k=0 state + ONE standalone IBL solve --------------------------
    st, ibl_info = build_wing_k0_state(mc, wc, cfg, case,
                                       ramp["phi"], ramp["gamma"])

    # -- the pack + W2 ----------------------------------------------------------
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

    # -- F at the seed + the assembly/eval timing samples -----------------------
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

    # -- the measurement: N augmented Newton steps with the ladder --------------
    ladder = SchurLadder(pack)
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

    # -- the adjudication (D4 = the GV5.4 D5 verbatim) + the ratio (D5) ---------
    aug_mean = float(np.mean(ws_steps[1:len(hist) + 1])) if hist else \
        float("nan")
    accepted = [bool(h["accepted"]) for h in hist]
    working_rung = None
    for rung in (3, 4):
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
    _record("(c)", f"wall {level} @16 threads", "RECORDED",
            f"seed {t_seed:.0f}s + IBL {ibl_info['wall_s']:.0f}s + "
            f"newton {t_newton:.0f}s")

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
    ax0.set_ylabel("wall s @16 threads")
    ax0.set_title(f"GV5.6 {level}: step cost "
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
    print(f"--- GV5.6 done in {time.perf_counter() - t_all:.0f}s; "
          f"verdicts: "
          f"{sum(1 for *_, v in SUMMARY if v == 'PASS')} PASS / "
          f"{sum(1 for *_, v in SUMMARY if v == 'FAIL')} FAIL / "
          f"{sum(1 for *_, v in SUMMARY if v == 'RECORDED')} RECORDED ---",
          flush=True)
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
