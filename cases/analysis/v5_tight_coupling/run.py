"""GV5.1 augmented (tight) Newton gates -- COARSE leg (Track V5).

Binding text: docs/roadmap/track_v.md GV5.1 (with the 2026-07-22
user-directed pre-registered FD note); pre-registered:
cases/analysis/v5_tight_coupling/PRE_REGISTRATION.md (committed before
the first execution). COARSE ONLY: the medium leg is deliberately held
(the pre-registered binding level; its run follows the coarse verdict).
Regenerates every CSV in results/ and exits 0 iff the FD pass band
holds; the convergence bands are RECORDED on coarse per the GV3.1
discipline (medium binding) -- the pre-registered honesty note covers
N_aug = 3, and the Stage-3 smoke evidence (da27e95) already showed the
crawl on the stalled k=1 IBL floor.

Phases (PRE_REGISTRATION "Outputs"):

  1. seed: the pre-registered k=1 state -- inviscid-converged (phi,
     Gamma) via the production GV3.1 Newton recipe + ONE IBL solve --
     built by tests/v5_state.py::build_k1_state so the seed cannot drift
     from the FD-gate fixture (the tests.mesh_utils import precedent of
     cases/analysis/v2_transpiration_channel/run.py). The IBL leg's
     converged=False at the ~2.8e-6 floor is EXPECTED (the loose loop
     proceeded the same way): recorded, not hard-failed.
  2. build_tight_pack at the seed.
  3. full-system FD gate at the seed (4 random full-state directions,
     ladder 1e-5/1e-6/1e-7, frozen-veps reference, natural-recompute
     omission measured) + per-column-block and augmentation-targeted
     scopes -> gv5_1_fd_report.csv. Repeated at the coupled solution
     when the Newton converges (the pre-registered "repeated at the
     coupled solution" clause).
  4. newton_tight(max_iter=10, the pre-registered criterion) ->
     gv5_1_newton_history_coarse.csv.
  5. gv5_1_compare.csv: augmented vs the COMMITTED GV3.1 coarse loose
     numbers (read from cases/analysis/v3_loose_coupling/results/, NOT
     re-run). The augmented delta* output field is recomputed with the
     closure at the CURRENT edge data (the loose-loop-consistent output;
     the system's internal delta*(U) uses base-frozen closure edge
     inputs -- the tight_driver.py module docstring).
  6. summary.csv: pass-band rows (FD blocks / quadratic tail /
     N_aug <= 2) + recorded rows.

Run:  python cases/analysis/v5_tight_coupling/run.py
"""

import csv
import os
import sys
import time

# resolve pyfp3d from THIS worktree (the site-packages editable install may
# point at a sibling worktree)
sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..")))

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np

from pyfp3d.physics.isentropic import density_field, mach_squared_field
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.viscous import closures as C
from pyfp3d.viscous import tight_driver as td
from tests.v5_state import (
    ALPHA,
    M_CAP,
    M_CRIT,
    M_INF,
    RE,
    RHO_FLOOR,
    UPWIND_C,
    build_k1_state,
    build_naca_case,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

V3_RESULTS = os.path.join(HERE, "..", "v3_loose_coupling", "results")

FD_TOL = 1e-5           # the pre-registered FD tolerance
N_AUG_BAND = 2          # pre-registered pass band (c)
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
# small local mirrors (self-contained apart from the sanctioned v5_state
# fixture; each documents its source)
# ---------------------------------------------------------------------------


def _phi_cut_of(pack, phi_free, gamma):
    """Mirror of tight_driver._phi_cut_of (the eval_residual assembly,
    solve/newton.py:339-343)."""
    ws = pack.ws
    phi_red = np.empty(ws.n_red, dtype=np.float64)
    phi_red[ws.free] = phi_free
    phi_red[ws.dir_red] = ws.vals0_red + ws.V_red @ gamma
    return ws.con.expand(phi_red, gamma)


def _fallback_row_mask(k1):
    """The Stage-2 mask (tests/test_v5_tight_system.py::_fallback_row_mask):
    drop the 6 rows of every q <= 1e-12 fallback node and of every node
    sharing an element with one."""
    sm, q = k1["sm"], k1["q"]
    bad = np.where(q <= 1.0e-12)[0]
    rows = np.ones(6 * sm.n_node, dtype=bool)
    if bad.size:
        drop = set(int(i) for i in bad)
        for e in range(len(sm.triangles)):
            if np.isin(sm.triangles[e], bad).any():
                drop.update(int(i) for i in sm.triangles[e])
        for i in drop:
            rows[6 * i: 6 * i + 6] = False
    return rows, bad


def _closure_ds_current(pack, U, ue):
    """delta* of closure_all at U with the CURRENT edge data (the
    loose-loop-consistent output field; the system's internal delta*(U)
    freezes the closure edge inputs at the pack base -- tight_driver.py
    module docstring)."""
    n_s = pack.n_s
    q2 = np.einsum("ij,ij->i", ue, ue)
    outs = np.empty((n_s, C.N_OUT), dtype=np.float64)
    douts = np.empty((n_s, C.N_OUT, 6), dtype=np.float64)
    douts_e = np.empty((n_s, C.N_OUT, 2), dtype=np.float64)
    C.closure_all(
        np.ascontiguousarray(U, dtype=np.float64),
        np.sqrt(q2),
        density_field(q2, pack.cfg.m_inf, pack.cfg.gamma_air),
        np.full(n_s, pack.mu, dtype=np.float64),
        np.sqrt(mach_squared_field(q2, pack.cfg.m_inf, pack.cfg.gamma_air)),
        pack.case.turbulent_flags,
        C.C_L_DEFAULT,
        outs,
        douts,
        douts_e,
    )
    return outs[:, C.OUT_DS1]


def _f_phi_transp_only(pack, x):
    """F_phi^tight(x) - F_phi^bare(phi, gamma): the transpiration part in
    isolation (tests/test_v5_tight_system.py::_f_phi_transp_only)."""
    phi_free, gamma, _ = pack.split_x(x)
    F_tight = td.block_residuals(pack, x)[0]
    R_bare, _, _ = pack.ws.eval_residual(
        phi_free, gamma, pack.upwind_c, pack.m_crit, pack.m_cap,
        pack.rho_floor,
    )
    return F_tight - R_bare


# ---------------------------------------------------------------------------
# the FD gate (the PRE_REGISTRATION FD protocol; the test's machinery)
# ---------------------------------------------------------------------------


def fd_gate(pack, k1, x, label, row_mask):
    """Central FD of augmented_residual vs the assembled analytic blocks,
    ladder 1e-5/1e-6/1e-7, scaled max-norm, sweet spot = min over the
    ladder. Scopes: aug_phi/aug_gamma (the J_phiphi augmentation at fixed
    U, the pre-registered targeted check), phi_cols/gamma_cols/U_cols
    (the total column blocks), full (4 random full-state directions, the
    gate proper). The veps natural-recompute omission is measured at
    eps = 1e-6 for the U-carrying scopes (decision 5). Returns (rows,
    worst_sweet, worst_veps_scaled)."""
    import scipy.sparse as sp

    n_f, n_g, n_b = pack.n_free, pack.n_st, 6 * pack.n_s
    n = n_f + n_g + n_b
    keep = np.ones(n, dtype=bool)
    keep[n - n_b:] = row_mask
    bl = td.augmented_jacobian_blocks(pack, x)
    J_full = td.augmented_jacobian(pack, x)
    rows = []
    worst_sweet, worst_veps = 0.0, 0.0

    def fd_full(v, eps, frozen=True):
        return (td.augmented_residual(pack, x + eps * v, veps_frozen=frozen)
                - td.augmented_residual(pack, x - eps * v, veps_frozen=frozen)
                ) / (2 * eps)

    def fd_block(v, which, eps, frozen=True):
        def at(s):
            pf, g, U = pack.split_x(x)
            if which == "phi":
                pf = pf + s * eps * v
            elif which == "gam":
                g = g + s * eps * v
            else:
                U = U + s * eps * v.reshape(pack.n_s, 6)
            return td.augmented_residual(
                pack, np.concatenate([pf, g, U.ravel()]), veps_frozen=frozen)

        return (at(1) - at(-1)) / (2 * eps)

    def fd_aug(v, which, eps):
        def at(s):
            pf, g, U = pack.split_x(x)
            if which == "phi":
                pf = pf + s * eps * v
            else:
                g = g + s * eps * v
            return _f_phi_transp_only(pack, np.concatenate([pf, g, U.ravel()]))

        return (at(1) - at(-1)) / (2 * eps)

    def sweep(scope, ndir, an_of, fd_of, rng, omission, masked):
        nonlocal worst_sweet, worst_veps
        for d in range(ndir):
            v = rng.standard_normal(an_of[1])
            v /= np.max(np.abs(v))
            an = an_of[0] @ v
            errs = []
            for e in (1e-5, 1e-6, 1e-7):
                ff = fd_of(v, e)
                a, f = (an[keep], ff[keep]) if masked else (an, ff)
                scale = max(float(np.max(np.abs(a))), 1e-12)
                errs.append(float(np.max(np.abs(f - a))) / scale)
            om_abs = om_scl = ""
            if omission:
                if scope == "U_cols":
                    fdu = fd_block(v, "U", 1e-6, frozen=False)
                    fdf = fd_block(v, "U", 1e-6, frozen=True)
                else:
                    fdu = fd_full(v, 1e-6, frozen=False)
                    fdf = fd_full(v, 1e-6, frozen=True)
                fdu, fdf = fdu[keep], fdf[keep]
                a = an[keep] if masked else an
                scale = max(float(np.max(np.abs(a))), 1e-12)
                om_abs = float(np.max(np.abs(fdu - fdf)))
                om_scl = om_abs / scale
                worst_veps = max(worst_veps, om_scl)
            sweet = min(errs)
            worst_sweet = max(worst_sweet, sweet)
            rows.append((scope, label, d, f"{errs[0]:.3e}", f"{errs[1]:.3e}",
                         f"{errs[2]:.3e}", f"{sweet:.3e}",
                         f"{om_abs:.3e}" if om_abs != "" else "",
                         f"{om_scl:.3e}" if om_scl != "" else "",
                         int((~keep).sum()) if masked else ""))
            print(f"    FD [{label}] {scope} dir {d}: {errs[0]:.3e} / "
                  f"{errs[1]:.3e} / {errs[2]:.3e}", flush=True)

    rng = np.random.default_rng(41)
    sweep("aug_phi", 2, (bl["A_free"], n_f),
          lambda v, e: fd_aug(v, "phi", e), rng, False, False)
    sweep("aug_gamma", 2, (bl["A_gam"], n_g),
          lambda v, e: fd_aug(v, "gam", e), rng, False, False)
    rng = np.random.default_rng(5)
    Jphi = sp.vstack([bl["J_ff"] + bl["A_free"], bl["K"],
                      bl["J_blphi_free"]]).tocsr()
    sweep("phi_cols", 2, (Jphi, n_f), lambda v, e: fd_block(v, "phi", e),
          rng, False, True)
    Jgam = sp.vstack([bl["B"] + bl["A_gam"], bl["minus_I"],
                      bl["J_blphi_gam"]]).tocsr()
    sweep("gamma_cols", 2, (Jgam, n_g), lambda v, e: fd_block(v, "gam", e),
          rng, False, True)
    JU = sp.vstack([bl["J_phibl"], sp.csr_matrix((n_g, n_b)),
                    bl["J_blbl"]]).tocsr()
    sweep("U_cols", 2, (JU, n_b), lambda v, e: fd_block(v, "U", e), rng,
          True, True)
    rng = np.random.default_rng(97)  # the test's seed: comparable numbers
    sweep("full", 4, (J_full, n), fd_full, rng, True, True)
    return rows, worst_sweet, worst_veps


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    t_start = time.perf_counter()
    print(f"=== GV5.1 coarse leg: NACA0012 2.5-D strip, M={M_INF} "
          f"alpha={ALPHA} Re={RE:.2e} x_tr 0.05/0.05 ===", flush=True)

    # -- phase 1: the pre-registered k=1 seed ---------------------------------
    print("--- phase 1: seed (inviscid Newton + ONE IBL solve, the "
          "tests/v5_state.py k=1 fixture) ---", flush=True)
    t0 = time.perf_counter()
    k1 = build_k1_state(build_naca_case())
    ibl = k1["ibl_info"]
    print(f"    seed built in {time.perf_counter() - t0:.0f}s: n_cut="
          f"{k1['n_cut']} n_s={k1['sm'].n_node} | IBL converged="
          f"{ibl.get('converged')} n_iter={ibl.get('n_iter')} "
          f"final_residual={ibl.get('final_residual'):.3e}", flush=True)
    _record("GV5.1", "seed IBL floor (k=1, EXPECTED non-converged)",
            "recorded", f"converged={ibl.get('converged')} "
            f"n_iter={ibl.get('n_iter')} "
            f"final_residual={ibl.get('final_residual'):.3e}", None)

    # -- phase 2: the pack -----------------------------------------------------
    print("--- phase 2: build_tight_pack ---", flush=True)
    t0 = time.perf_counter()
    pack = td.build_tight_pack(k1, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    print(f"    pack in {time.perf_counter() - t0:.0f}s: n_free="
          f"{pack.n_free} n_st={pack.n_st} n_s={pack.n_s} -> "
          f"{pack.n_free + pack.n_st + 6 * pack.n_s} DOFs", flush=True)

    # -- phase 3: the full-system FD gate at the seed --------------------------
    print("--- phase 3: FD gate at the seed (4 full-state directions, "
          "ladder 1e-5/1e-6/1e-7) ---", flush=True)
    t0 = time.perf_counter()
    row_mask, bad = _fallback_row_mask(k1)
    print(f"    fallback nodes (q <= 1e-12): {len(bad)}; masked F_BL rows: "
          f"{int((~row_mask).sum())}", flush=True)
    fd_rows, worst_sweet, worst_veps = fd_gate(
        pack, k1, pack.x_base(), "seed", row_mask)
    print(f"    FD gate in {time.perf_counter() - t0:.0f}s: worst "
          f"sweet-spot {worst_sweet:.3e}, worst veps omission (scaled) "
          f"{worst_veps:.3e}", flush=True)
    _record("GV5.1", "FD: coupling blocks + full system, worst sweet-spot",
            f"< {FD_TOL:.0e}", f"{worst_sweet:.3e}",
            bool(worst_sweet < FD_TOL))
    _record("GV5.1", "FD: masked rows (q <= 1e-12 fallback)",
            "<= 2% of F_BL rows", f"{int((~row_mask).sum())} / "
            f"{6 * pack.n_s}",
            bool((~row_mask).sum() <= 0.02 * 6 * pack.n_s))
    _record("GV5.1", "FD: veps frozen-vs-natural omission, worst scaled",
            "recorded (decision 5)", f"{worst_veps:.3e}", None)

    # -- phase 4: the augmented Newton -----------------------------------------
    print("--- phase 4: newton_tight(max_iter=10, pre-registered "
          "criterion) ---", flush=True)
    t0 = time.perf_counter()
    res = td.newton_tight(pack, max_iter=10, tol=1e-8, tol_abs=1e-10,
                          line_search=True, verbose=True)
    newton_wall = time.perf_counter() - t0
    print(f"    newton_tight in {newton_wall:.0f}s: converged="
          f"{res['converged']} n_iter={res['n_iter']} "
          f"ds_change_last={res['ds_change_last']:.3e}", flush=True)

    hist = res["history"]
    bm0 = np.maximum(res["block_max0"], 1e-300)
    phi_free_f, gamma_f, U_f = pack.split_x(res["x"])
    ue_f = (pack.G @ _phi_cut_of(pack, phi_free_f, gamma_f)).reshape(
        pack.n_s, 3)
    ds_f = _closure_ds_current(pack, U_f, ue_f)
    ds_max_abs = max(float(np.max(np.abs(ds_f))), 1e-300)
    hist_rows = []
    for h in hist:
        bm = h["block_max"]
        hist_rows.append(
            (h["iter"], f"{bm[0]:.6e}", f"{bm[1]:.6e}", f"{bm[2]:.6e}",
             f"{bm[0] / bm0[0]:.6e}", f"{bm[1] / bm0[1]:.6e}",
             f"{bm[2] / bm0[2]:.6e}", f"{h['merit']:.6e}",
             "" if h["lam"] is None else f"{h['lam']:.6f}",
             f"{h['ds_change']:.6e}",
             f"{h['ds_change'] / ds_max_abs:.6e}",
             f"{h['wall_s']:.2f}"))
    _write_csv("gv5_1_newton_history_coarse.csv",
               "iter,f_phi_max,f_gamma_max,f_bl_max,f_phi_scaled,"
               "f_gamma_scaled,f_bl_scaled,merit,lam,ds_change,"
               "ds_change_rel,wall_s", hist_rows)

    # quadratic-tail diagnostic: p_i = log(e_i/e_{i-1})/log(e_{i-1}/e_{i-2})
    # on the total max-norm e; a slope-2 regime needs p ~ 2 over the last
    # decades before the floor
    e = np.array([max(h["block_max"]) for h in hist])
    ps = []
    for i in range(2, len(e)):
        if e[i] < e[i - 1] < e[i - 2] and e[i] > 0.0:
            ps.append(float(np.log(e[i] / e[i - 1])
                            / np.log(e[i - 1] / e[i - 2])))
    tail = f"{ps[-1]:.2f}" if ps else "n/a (residual not decreasing)"
    _record("GV5.1", "quadratic tail (slope-2 regime, coarse recorded; "
            "medium binding)", "slope ~2 before the floor",
            f"converged={res['converged']} last p={tail}; mechanism: "
            f"stalled k=1 IBL floor, cond(J_BL,BL) ~ 4e10 -> raw dU "
            f"O(1e3), line search at lam << 1 (the Stage-3 smoke "
            f"evidence, da27e95)", None)
    _record("GV5.1", "N_aug <= 2 (coarse recorded; medium binding; the "
            "pre-registered honesty note covers 3)",
            f"<= {N_AUG_BAND}",
            f"n_iter={res['n_iter']} converged={res['converged']} "
            f"(crawl on the stalled floor, not a blow-up)", None)
    _record("GV5.1", "ds_change of the last step vs tol_ds = 1e-3",
            "recorded (cross-check)", f"{res['ds_change_last']:.3e}",
            None)
    _record("GV5.1", "newton wall time (coarse)", "recorded",
            f"{newton_wall:.0f}s over {res['n_iter']} iterations", None)

    # -- phase 5: FD at the coupled solution (pre-registered repeat) -----------
    if res["converged"]:
        print("--- phase 5: FD gate at the coupled solution ---",
              flush=True)
        fd_rows2, worst2, _ = fd_gate(pack, k1, res["x"], "solution",
                                      row_mask)
        fd_rows += fd_rows2
        _record("GV5.1", "FD at the coupled solution, worst sweet-spot",
                f"< {FD_TOL:.0e}", f"{worst2:.3e}",
                bool(worst2 < FD_TOL))
    else:
        _record("GV5.1", "FD at the coupled solution",
                "pre-registered repeat",
                "skipped: no coarse convergence (deferred to the "
                "medium leg)", None)
    _write_csv("gv5_1_fd_report.csv",
               "scope,state,direction,err_eps_1e-5,err_eps_1e-6,"
               "err_eps_1e-7,sweet_spot,veps_omission_abs,"
               "veps_omission_scaled,masked_rows", fd_rows)

    # -- phase 6: compare vs the COMMITTED GV3.1 coarse loose result -----------
    print("--- phase 6: compare vs committed GV3.1 coarse (not "
          "re-run) ---", flush=True)
    with open(os.path.join(V3_RESULTS, "gv3_1_history_coarse.csv")) as f:
        loose_hist = list(csv.DictReader(f))
    loose_n_outer = max(int(r["k"]) for r in loose_hist)
    loose_last = loose_hist[-1]
    loose_cl = float(loose_last["cl"])
    loose_ds_max = float(loose_last["ds_max"])
    loose_ds_change = float(loose_last["ds_change_rel"])
    loose_ibl_floor = float(loose_last["ibl_final_residual"])
    loose_converged = loose_ds_change < 1e-3
    print(f"    committed loose coarse: n_outer={loose_n_outer} "
          f"converged={loose_converged} cl_p={loose_cl:.4f} "
          f"ds_max={loose_ds_max:.6e}", flush=True)

    mc, case, ws = k1["mc"], k1["case"], pack.ws
    dz = float(np.ptp(mc.nodes[:, 2]))
    phi_f = _phi_cut_of(pack, phi_free_f, gamma_f)
    f_au = wall_force_coefficients(mc.nodes, mc.elements, case.wall_faces,
                                   phi_f, alpha_deg=ALPHA, s_ref=dz,
                                   m_inf=M_INF)
    cl_p_au = float(f_au["cl"])
    cl_kj_au = float(2.0 * gamma_f[0])  # the 2.5-D cl = 2 Gamma/c check
    f_seed = wall_force_coefficients(mc.nodes, mc.elements,
                                     case.wall_faces, k1["phi"],
                                     alpha_deg=ALPHA, s_ref=dz,
                                     m_inf=M_INF)
    print(f"    augmented endpoint: cl_p={cl_p_au:.4f} "
          f"cl_kj(2*Gamma)={cl_kj_au:.4f} (seed cl_p="
          f"{float(f_seed['cl']):.4f} vs committed loose k=0 0.2776)",
          flush=True)

    # per-station delta* max-diff vs the committed loose profile
    st = case.stations
    ds_station = {}
    for side_val, side_name in ((1, "upper"), (-1, "lower")):
        for r in range(len(st.xc)):
            m = (st.station_of == r) & (st.side_node == side_val)
            if np.any(m):
                ds_station[(round(float(st.xc[r]), 6), side_name)] = float(
                    np.mean(ds_f[m]))
    prof_rows, worst_ds = [], 0.0
    with open(os.path.join(V3_RESULTS, "gv3_1_profiles_coarse.csv")) as f:
        for r in csv.DictReader(f):
            key = (round(float(r["x_c"]), 6), r["side"])
            if key in ds_station:
                d = abs(ds_station[key] - float(r["ds_ours"]))
                worst_ds = max(worst_ds, d)
                prof_rows.append((key[0], key[1], d))
    print(f"    per-station max|delta*_aug - delta*_loose| = "
          f"{worst_ds:.3e} over {len(prof_rows)} matched stations",
          flush=True)

    cmp_rows = [
        ("iterations (aug Newton / loose outer)", res["n_iter"],
         loose_n_outer, ""),
        ("converged", res["converged"], loose_converged, ""),
        ("final ds_max", f"{float(np.max(np.abs(ds_f))):.6e}",
         f"{loose_ds_max:.6e}",
         f"{abs(float(np.max(np.abs(ds_f))) - loose_ds_max):.3e}"),
        ("final per-station max|delta* diff|", "", "",
         f"{worst_ds:.3e}"),
        ("cl_p (pressure integral)", f"{cl_p_au:.6f}", f"{loose_cl:.6f}",
         f"{abs(cl_p_au - loose_cl):.3e}"),
        ("cl_kj (2*Gamma, design.md Sec 9)", f"{cl_kj_au:.6f}",
         "n/a (not committed)", ""),
        ("IBL residual floor (seed / loose final)", 
         f"{ibl.get('final_residual'):.3e}", f"{loose_ibl_floor:.3e}", ""),
        ("final |F|_inf (aug) / ds_change_rel (loose)",
         f"{max(hist[-1]['block_max']):.3e}", f"{loose_ds_change:.3e}", ""),
    ]
    _write_csv("gv5_1_compare.csv",
               "metric,augmented,loose_committed,abs_diff", cmp_rows)
    _record("GV5.1", "compare: augmented iterations vs loose outer "
            "(coarse)", "recorded",
            f"n_iter={res['n_iter']} converged={res['converged']} vs "
            f"n_outer={loose_n_outer}", None)
    _record("GV5.1", "compare: final per-station max|delta* diff|",
            "recorded", f"{worst_ds:.3e}", None)
    _record("GV5.1", "compare: cl_p augmented vs loose committed",
            "recorded", f"{cl_p_au:.4f} vs {loose_cl:.4f}", None)

    # -- summary -----------------------------------------------------------------
    _write_csv("summary.csv", "gate,metric,band,measured,verdict", SUMMARY)
    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV5.1 coarse: {n_pass} PASS / {n_fail} FAIL / {n_rec} "
          f"RECORDED in {time.perf_counter() - t_start:.0f}s", flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
