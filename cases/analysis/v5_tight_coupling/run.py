"""GV5.1 augmented (tight) Newton gates -- AMENDED protocol (Track V5).

Binding text: docs/roadmap/track_v.md GV5.1; pre-registered:
cases/analysis/v5_tight_coupling/PRE_REGISTRATION.md Addendum 2
(2026-07-23, user-adjudicated seed amendment, committed at 606b149
BEFORE this execution). Trigger: the k=1-seed coarse leg (71df59a,
protocol=k1seed rows kept in every output CSV) passed band (a) FD
exactness but not (b)/(c) -- the k=1 state sits in a near-null
direction of the steady BL block (cond(J_BL,BL) ~ 4e10), NOT a
coupling defect. Amended protocol per Addendum 2:

  1. SEED: the loose loop's CONVERGED state, regenerated with the
     committed GV3.1 recipe (make_picard_lifting_driver +
     run_loose_coupling, CouplingConfig(re_chord=3e6, m_inf=0.5,
     alpha_deg=2) defaults: omega=1.0, n_outer_max=10, tol_ds=1e-3).
     The regen is the Addendum-2 sanctioned exception to "never
     recompute the loose numbers" (the committed artifacts hold CSVs,
     not arrays); a cross-check against the committed GV3.1 history
     guards the recipe. Reframed per the 2026-07-23 user decision
     after the medium-seed diagnosis
     (results/gv5_1_medium_seed_diagnosis.md): the guard is the
     WIRING check (converged + |dcl_k0| <= 1e-8 on the inviscid
     baseline -- a mismatch stops the leg as a recipe error); the
     pointwise final-cl/ds_max offsets are the loose trajectory's
     scatter across the IBL-floor manifold (three code/env
     trajectories -> three fixed points; the committed GV3.1 medium
     fixed point is not reproducible in this environment; the coarse
     one reproduces to 1.5e-4) -- a recorded FINDING, not a recipe
     error. The HEAD regen itself is bit-identical run-to-run.
  2. build_tight_pack at the loose-converged state (phi = res.phi; the
     edge data recomputed from res.phi; U = res.U; the inflow band
     frozen at the converged res.U values -- pre-registration
     decision 4).
  3. FD gate at the seed (the Stage-3 machinery, unchanged).
  4. newton_tight(max_iter=10) as a POLISH -> gv5_1_newton_history_
     {level}.csv; band (b) read on the polish history.
  5. FD gate again at the polish endpoint (the converged solution, or
     the last iterate -- flagged).
  6. gv5_1_compare.csv + summary.csv: N_polish and the honest total
     N_loose + N_polish vs the committed loose baseline (4 coarse /
     5 medium, READ from cases/analysis/v3_loose_coupling/results/,
     never recomputed); band (c) standalone N_aug <= 2 stays recorded
     NOT met (the k1seed rows persist; no further standalone retries).

Merge discipline: the committed k1seed CSVs are read, tagged
protocol=k1seed (level=coarse), kept, and the amended rows appended;
summary.csv/fd_report/compare carry a protocol (+level) column. The
committed gv5_1_newton_history_coarse.csv was renamed to
gv5_1_newton_history_coarse_k1seed.csv (both protocols' records kept).

Run:  python cases/analysis/v5_tight_coupling/run.py --levels coarse
      python cases/analysis/v5_tight_coupling/run.py --levels medium
"""

import argparse
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

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.physics.isentropic import density_field, mach_squared_field
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import NewtonWorkspace
from pyfp3d.viscous import closures as C
from pyfp3d.viscous import tight_driver as td
from pyfp3d.viscous.coupling import (
    CouplingConfig,
    build_airfoil_case,
    make_picard_lifting_driver,
    run_loose_coupling,
)
from pyfp3d.viscous.ibl3 import IBL3Solver
from pyfp3d.viscous.transpiration import edge_velocity_per_zone
from tests.v5_state import (
    ALPHA,
    M_CAP,
    M_CRIT,
    M_INF,
    NACA_DIR,
    RE,
    RHO_FLOOR,
    UPWIND_C,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

V3_RESULTS = os.path.join(HERE, "..", "v3_loose_coupling", "results")

FD_TOL = 1e-5           # the pre-registered FD tolerance
# the recipe cross-check, reframed per the 2026-07-23 user decision
# (gv5_1_medium_seed_diagnosis.md): the wiring guard = converged + the
# inviscid k=0 cl match; the pointwise final-cl/ds_max offsets are the
# IBL-floor trajectory scatter (a recorded finding, not a recipe error)
RECIPE_TOL_K0 = 1e-8    # |dcl| at the inviscid k=0 (wiring discriminator)
SUMMARY = []            # (protocol, gate, metric, band, measured, verdict)
SUMMARY_WRITTEN = 0     # how many SUMMARY rows are already on disk


def _record(gate, metric, band, measured, ok=None):
    verdict = "RECORDED" if ok is None else ("PASS" if ok else "FAIL")
    SUMMARY.append(("amended", gate, metric, band, measured, verdict))
    print(f"  [{verdict:8s}] {gate} {metric}: measured={measured} "
          f"(band: {band})", flush=True)


def _write_lines(name, header, lines):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as f:
        f.write(header + "\n")
        for ln in lines:
            f.write(ln + "\n")
    print(f"  wrote {path}", flush=True)


def _read_body(path):
    """Existing CSV body lines (header dropped), [] when absent."""
    if not os.path.exists(path):
        return []
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    return lines[1:] if lines else []


def _merge_tagged(name, new_header, tag_legacy, drop_levels, fresh_lines):
    """Keep the committed rows (tagged protocol=k1seed when the file
    predates the protocol column), drop the amended rows of the levels
    being re-run, append the fresh rows."""
    kept = []
    for ln in _read_body(os.path.join(RESULTS, name)):
        if ln.startswith("amended,"):
            flds = ln.split(",")
            if len(flds) > 1 and flds[1] in drop_levels:
                continue
            kept.append(ln)
        elif ln.startswith("k1seed,"):
            kept.append(ln)
        else:
            kept.append(tag_legacy + ln)
    _write_lines(name, new_header, kept + fresh_lines)


def _write_summary(level):
    """summary.csv: protocol,gate,metric,band,measured,verdict. Legacy
    5-column k1seed rows are kept verbatim with a k1seed tag (their
    metric/band/measured fields contain free text with commas -- the
    committed file's own convention); amended rows carry the level as a
    '<level> ' metric prefix and comma-free gate/metric/band fields."""
    global SUMMARY_WRITTEN
    kept = []
    for ln in _read_body(os.path.join(RESULTS, "summary.csv")):
        if ln.startswith("amended,"):
            flds = ln.split(",")
            if len(flds) >= 3 and flds[2].startswith(level + " "):
                continue
            kept.append(ln)
        elif ln.startswith("k1seed,"):
            kept.append(ln)
        else:
            kept.append("k1seed," + ln)
    fresh = [",".join(str(x) for x in t) for t in SUMMARY[SUMMARY_WRITTEN:]]
    SUMMARY_WRITTEN = len(SUMMARY)
    _write_lines("summary.csv", "protocol,gate,metric,band,measured,verdict",
                 kept + fresh)


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


def _fallback_row_mask(st):
    """The Stage-2 mask (tests/test_v5_tight_system.py::_fallback_row_mask):
    drop the 6 rows of every q <= 1e-12 fallback node and of every node
    sharing an element with one."""
    sm, q = st["sm"], st["q"]
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
# the FD gate (the PRE_REGISTRATION FD protocol; the Stage-3 machinery,
# verbatim from the k1seed leg -- only the unused state parameter renamed)
# ---------------------------------------------------------------------------


def fd_gate(pack, st, x, label, row_mask):
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
# the amended seed: the loose loop's converged state (Addendum 2)
# ---------------------------------------------------------------------------


def build_loose_converged_state(mc, wc, cfg, case, res):
    """The tight-pack state dict at the loose loop's CONVERGED result
    (the tests/v5_state.py::build_k1_state keys, same construction
    rules): phi = res.phi (the state AFTER the final FP re-solve,
    coupling.py:836-864); the edge data recomputed from res.phi (the
    loose loop's res.ue_surf belongs to the PREVIOUS phi); U = res.U;
    the inflow band's Dirichlet state frozen at the converged res.U
    values (pre-registration decision 4: frozen at the seed -- the
    loop's k=1 laminar seeds would kick F_BL at the band)."""
    sm = case.sm
    mu = 1.0 / cfg.re_chord
    n_cut = len(mc.nodes)
    phi = np.asarray(res.phi, dtype=np.float64)
    gamma = np.asarray(res.gamma, dtype=np.float64)

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

    # the airfoil inflow band (coupling.py:738-740), state frozen at the
    # converged U (decision 4)
    inflow_mask = case.stations.xc[case.stations.station_of] <= cfg.inflow_band_x
    inflow_state = np.array(res.U[inflow_mask], dtype=np.float64).copy()

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
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(M_INF)
    phi_free = phi[: ws.n_red][ws.free].copy()
    return {
        "mc": mc, "wc": wc, "cfg": cfg, "case": case, "sm": sm,
        "mu": mu, "n_cut": n_cut, "phi": phi, "gamma": gamma,
        "le_mask_vol": le_mask_vol, "ue_surf": ue_surf, "q": q,
        "rho_e": rho_e, "mach_e": mach_e,
        "U": np.ascontiguousarray(res.U, dtype=np.float64),
        "ws": ws, "phi_free": phi_free, "solver": solver,
        "inflow_mask": inflow_mask, "inflow_state": inflow_state,
    }


# ---------------------------------------------------------------------------
# one amended leg (one mesh level)
# ---------------------------------------------------------------------------


def run_amended_leg(level):
    t_leg = time.perf_counter()
    walls = {}
    print(f"=== GV5.1 amended leg [{level}]: NACA0012 2.5-D strip, "
          f"M={M_INF} alpha={ALPHA} Re={RE:.2e} x_tr 0.05/0.05; seed = "
          f"the loose loop's CONVERGED state (Addendum 2) ===",
          flush=True)

    # -- phase 1: the loose-converged seed (the committed GV3.1 recipe) ------
    print("--- phase 1: loose-loop regeneration (committed GV3.1 recipe: "
          "picard driver, omega=1.0, n_outer_max=10, tol_ds=1e-3) ---",
          flush=True)
    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, f"{level}.msh")))
    cfg = CouplingConfig(re_chord=RE, m_inf=M_INF, alpha_deg=ALPHA)
    case = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz  # chord x span thickness (the v3 convention)

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        print(f"    loose k={k}: cl_p={f['cl']:.8f} "
              f"cl_kj(2*Gamma)={2.0 * gamma[0]:.8f} "
              f"cd_p={f['cd_pressure']:.8f}", flush=True)
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    t0 = time.perf_counter()
    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)
    res = run_loose_coupling(driver, case, cfg, probe=probe)
    walls["loose"] = time.perf_counter() - t0
    print(f"    loose regen in {walls['loose']:.0f}s: converged="
          f"{res.converged} n_outer={res.n_outer}", flush=True)

    # the recipe cross-check against the COMMITTED GV3.1 history (read,
    # never recomputed). Reframed per the 2026-07-23 user decision after
    # the medium-seed diagnosis (gv5_1_medium_seed_diagnosis.md): the
    # recipe-guard function is the WIRING check -- converged + the
    # inviscid k=0 cl (identical mesh/FP path; |dcl_k0| <= 1e-8 stops
    # the leg as a recipe error). The pointwise final-cl/ds_max offsets
    # are the loose trajectory's scatter across the IBL-floor manifold
    # (the diagnosis: three code/env trajectories -> three fixed
    # points), RECORDED with the caveat, not verdicted.
    import csv
    with open(os.path.join(V3_RESULTS, f"gv3_1_history_{level}.csv")) as f:
        loose_hist = list(csv.DictReader(f))
    loose_n_outer = max(int(r["k"]) for r in loose_hist)
    loose_last = loose_hist[-1]
    loose_cl = float(loose_last["cl"])
    loose_ds_max = float(loose_last["ds_max"])
    loose_ds_change = float(loose_last["ds_change_rel"])
    loose_ibl_floor = float(loose_last["ibl_final_residual"])
    loose_converged = loose_ds_change < 1e-3
    loose_k0_cl = float(loose_hist[0]["cl"])
    regen_cl = float(res.history[-1]["cl"])
    regen_ds_max = float(res.history[-1]["ds_max"])
    d_cl = abs(regen_cl - loose_cl)
    d_ds = abs(regen_ds_max - loose_ds_max)
    d_k0 = abs(float(res.history[0]["cl"]) - loose_k0_cl)
    recipe_ok = bool(res.converged) and d_k0 <= RECIPE_TOL_K0
    print(f"    cross-check vs committed GV3.1 {level}: n_outer "
          f"{res.n_outer} vs {loose_n_outer} | wiring guard: converged="
          f"{res.converged}, dcl_k0={d_k0:.3e} (<= {RECIPE_TOL_K0:.0e}) "
          f"| recorded offsets: dcl_final={d_cl:.3e} "
          f"dds_max_final={d_ds:.3e} (rel {d_ds / loose_ds_max:.3e})",
          flush=True)
    # the reproducibility evidence: committed vs regen, per outer k
    repro_lines = []
    for h_comm, h_reg in zip(loose_hist, res.history):
        repro_lines.append(",".join([
            str(h_comm["k"]),
            h_comm.get("cl", ""), f"{h_reg.get('cl', ''):.16g}",
            h_comm.get("ds_max", ""),
            (f"{h_reg['ds_max']:.16g}"
             if h_reg.get("ds_max") is not None else ""),
            h_comm.get("ibl_n_iter", ""), str(h_reg.get("ibl_n_iter", ""))]))
    _write_lines(
        f"gv5_1_seed_reproducibility_{level}.csv",
        "k,cl_committed,cl_regen,ds_max_committed,ds_max_regen,"
        "ibl_n_iter_committed,ibl_n_iter_regen", repro_lines)
    _record("GV5.1", f"{level} loose-regen recipe cross-check vs "
            "committed GV3.1 (Addendum 2; wiring guard per the "
            "2026-07-23 user decision)",
            f"converged; |dcl_k0| <= {RECIPE_TOL_K0:.0e}",
            f"converged={res.converged}; dcl_k0={d_k0:.3e}; "
            f"n_outer {res.n_outer} vs committed {loose_n_outer}",
            bool(recipe_ok))
    _record("GV5.1", f"{level} seed pointwise offsets vs committed "
            "GV3.1 (loose-trajectory scatter across the IBL-floor "
            "manifold; a recorded FINDING, not a recipe error -- "
            "gv5_1_medium_seed_diagnosis.md)", "recorded (scatter caveat)",
            f"dcl_final={d_cl:.3e}; dds_max={d_ds:.3e} rel "
            f"{d_ds / loose_ds_max:.3e}; n_outer {res.n_outer} vs "
            f"{loose_n_outer}", None)
    if level == "medium":
        _record("GV5.1", "medium loose trajectory scatter (FINDING, "
                "2026-07-23 diagnosis): three code/env trajectories -> "
                "three fixed points", "recorded",
                "committed cl=0.2719/ds=6.841e-3 (n_outer=5); HEAD "
                "regen cl=0.2814/ds=3.454e-3 (n_outer=3); c2dc325 "
                "regen cl=0.2217/ds=9.728e-3 (n_outer=6). Mechanism: "
                "1e-12-level douts perturbations amplified on the "
                "cond(J_BL,BL)~4e10 near-null manifold at the "
                "100-iter-capped IBL solves; HEAD bit-identical "
                "run-to-run; coarse well-conditioned (0.14% k=1 ds)",
                None)
    if level == "coarse":
        _record("GV5.1", "coarse seed reproducibility diagnostic "
                "(2026-07-23, why the cross-check band is not bit-exact)",
                "recorded",
                "HEAD regen bit-identical run-to-run (|dphi|=|dU|=0, "
                "same IBL iters); pre-Stage-2 code (c2dc325, same env) "
                "gives k=0 cl drift 1.6e-12 and n_outer=3 vs committed "
                "4 -> the committed artifacts predate an environment "
                "shift; the IBL floor stopping amplifies 1e-12 to "
                "~1e-4 in cl. Band set ~5x above the measured drift, "
                "~10x below any wiring-error signature", None)
    _record("GV5.1", f"{level} seed state (loose-converged)",
            "recorded",
            f"n_outer={res.n_outer} converged={res.converged} "
            f"cl={regen_cl:.8f} ds_max={regen_ds_max:.6e} "
            f"ds_change_rel={float(res.history[-1]['ds_change_rel']):.3e} "
            f"wall={walls['loose']:.0f}s", None)
    ibl = res.ibl_info
    _record("GV5.1", f"{level} seed IBL floor (loose-final, EXPECTED "
            "non-converged)", "recorded",
            f"converged={ibl.get('converged')} n_iter={ibl.get('n_iter')} "
            f"final_residual={ibl.get('final_residual'):.3e}", None)
    if not recipe_ok:
        print("RECIPE ERROR: the wiring guard failed (the regenerated "
              "loose loop did not converge, or the inviscid k=0 cl "
              "mismatches the committed GV3.1 baseline) -- STOP (only "
              "the summary rows and the reproducibility CSV written)",
              flush=True)
        _write_summary(level)
        sys.exit(2)

    # -- phase 2: the pack at the loose-converged state -----------------------
    print("--- phase 2: build_tight_pack at the loose-converged state "
          "---", flush=True)
    t0 = time.perf_counter()
    st = build_loose_converged_state(mc, wc, cfg, case, res)
    pack = td.build_tight_pack(st, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    walls["pack"] = time.perf_counter() - t0
    print(f"    pack in {walls['pack']:.0f}s: n_free={pack.n_free} "
          f"n_st={pack.n_st} n_s={pack.n_s} -> "
          f"{pack.n_free + pack.n_st + 6 * pack.n_s} DOFs", flush=True)

    # -- phase 3: the full-system FD gate at the seed -------------------------
    print("--- phase 3: FD gate at the seed (4 full-state directions, "
          "ladder 1e-5/1e-6/1e-7) ---", flush=True)
    t0 = time.perf_counter()
    row_mask, bad = _fallback_row_mask(st)
    print(f"    fallback nodes (q <= 1e-12): {len(bad)}; masked F_BL "
          f"rows: {int((~row_mask).sum())}", flush=True)
    fd_rows, worst_sweet, worst_veps = fd_gate(
        pack, st, pack.x_base(), "seed", row_mask)
    walls["fd_seed"] = time.perf_counter() - t0
    print(f"    FD gate in {walls['fd_seed']:.0f}s: worst sweet-spot "
          f"{worst_sweet:.3e}, worst veps omission (scaled) "
          f"{worst_veps:.3e}", flush=True)
    _record("GV5.1", f"{level} FD at the seed: coupling blocks + full "
            "system, worst sweet-spot", f"< {FD_TOL:.0e}",
            f"{worst_sweet:.3e}", bool(worst_sweet < FD_TOL))
    _record("GV5.1", f"{level} FD at the seed: masked rows "
            "(q <= 1e-12 fallback)", "<= 2% of F_BL rows",
            f"{int((~row_mask).sum())} / {6 * pack.n_s}",
            bool((~row_mask).sum() <= 0.02 * 6 * pack.n_s))
    _record("GV5.1", f"{level} FD at the seed: veps frozen-vs-natural "
            "omission, worst scaled", "recorded (decision 5)",
            f"{worst_veps:.3e}", None)

    # -- phase 4: the Newton polish -------------------------------------------
    print("--- phase 4: newton_tight(max_iter=10) as a POLISH from the "
          "loose-converged seed ---", flush=True)
    t0 = time.perf_counter()
    res_n = td.newton_tight(pack, max_iter=10, tol=1e-8, tol_abs=1e-10,
                            line_search=True, verbose=True)
    walls["newton"] = time.perf_counter() - t0
    print(f"    newton_tight in {walls['newton']:.0f}s: converged="
          f"{res_n['converged']} n_iter={res_n['n_iter']} "
          f"ds_change_last={res_n['ds_change_last']:.3e}", flush=True)

    hist = res_n["history"]
    bm0 = np.maximum(res_n["block_max0"], 1e-300)
    phi_free_f, gamma_f, U_f = pack.split_x(res_n["x"])
    ue_f = (pack.G @ _phi_cut_of(pack, phi_free_f, gamma_f)).reshape(
        pack.n_s, 3)
    ds_f = _closure_ds_current(pack, U_f, ue_f)
    ds_max_abs = max(float(np.max(np.abs(ds_f))), 1e-300)
    hist_lines = []
    for h in hist:
        bm = h["block_max"]
        hist_lines.append(",".join([
            str(h["iter"]), f"{bm[0]:.6e}", f"{bm[1]:.6e}", f"{bm[2]:.6e}",
            f"{bm[0] / bm0[0]:.6e}", f"{bm[1] / bm0[1]:.6e}",
            f"{bm[2] / bm0[2]:.6e}", f"{h['merit']:.6e}",
            "" if h["lam"] is None else f"{h['lam']:.6f}",
            f"{h['ds_change']:.6e}",
            f"{h['ds_change'] / ds_max_abs:.6e}",
            f"{h['wall_s']:.2f}"]))
    _write_lines(
        f"gv5_1_newton_history_{level}.csv",
        "iter,f_phi_max,f_gamma_max,f_bl_max,f_phi_scaled,"
        "f_gamma_scaled,f_bl_scaled,merit,lam,ds_change,"
        "ds_change_rel,wall_s", hist_lines)

    # band (b): the quadratic-tail read on the polish history --
    # p_i = log(e_i/e_{i-1})/log(e_{i-1}/e_{i-2}) on the total max-norm
    # e; a slope-2 regime shows p ~ 2 over the last decades before the
    # floor; if the BL floor pins, record where the tail flattens
    e = np.array([max(h["block_max"]) for h in hist])
    ps = []
    for i in range(2, len(e)):
        if e[i] < e[i - 1] < e[i - 2] and e[i] > 0.0:
            ps.append(float(np.log(e[i] / e[i - 1])
                            / np.log(e[i - 1] / e[i - 2])))
    tail_ps = ("; ".join(f"p={p:.2f}" for p in ps[-3:])
               if ps else "n/a (residual not decreasing)")
    flatten = (f"e_last={e[-1]:.3e} f_bl_last="
               f"{hist[-1]['block_max'][2]:.3e}")
    if level == "coarse":
        tail_ok = None  # coarse recorded per the GV3.1 discipline
    else:
        # medium binding: a converged polish whose last contraction
        # ratios sit in the slope-2 band (vacuously true when the
        # polish converged in < 3 measured contractions)
        tail_ok = bool(res_n["converged"]) and (
            not ps or 1.5 <= float(np.median(ps[-3:])) <= 2.5)
    _record("GV5.1", f"{level} quadratic tail read on the polish "
            "history (slope-2 regime over the last decades before the "
            "floor)", "slope ~2 before the floor" +
            (" (coarse recorded; medium binding)" if level == "coarse"
             else " (binding; PASS = converged and median(last 3 p) "
             "in [1.5, 2.5], vacuous when < 3 contractions)"),
            f"converged={res_n['converged']} n_iter={res_n['n_iter']} "
            f"{tail_ps} | {flatten}", tail_ok)

    # band (c): the honest counts (the standalone N_aug <= 2 stays
    # recorded NOT met on the k=1 seed -- the k1seed rows persist)
    n_polish = res_n["n_iter"]
    _record("GV5.1", f"{level} N_polish and N_total = N_loose + "
            "N_polish vs the committed loose baseline (Addendum 2 (c); "
            "standalone N_aug <= 2 recorded NOT met on the k=1 seed, "
            "no further standalone retries)", "recorded (honest counts)",
            f"N_polish={n_polish} converged={res_n['converged']} "
            f"N_total={res.n_outer + n_polish} vs committed loose "
            f"n_outer={loose_n_outer}", None)
    _record("GV5.1", f"{level} ds_change of the last polish step vs "
            "tol_ds = 1e-3", "recorded (cross-check)",
            f"{res_n['ds_change_last']:.3e}", None)
    _record("GV5.1", f"{level} newton wall time", "recorded",
            f"{walls['newton']:.0f}s over {n_polish} iterations", None)

    # -- phase 5: FD gate at the polish endpoint ------------------------------
    endpoint_label = "solution" if res_n["converged"] else "last-iterate"
    print(f"--- phase 5: FD gate at the polish endpoint "
          f"({endpoint_label}) ---", flush=True)
    t0 = time.perf_counter()
    fd_rows2, worst2, _ = fd_gate(pack, st, res_n["x"], endpoint_label,
                                  row_mask)
    walls["fd_solution"] = time.perf_counter() - t0
    fd_rows += fd_rows2
    flag = "" if res_n["converged"] else " (LAST ITERATE, not converged)"
    _record("GV5.1", f"{level} FD at the coupled solution, worst "
            f"sweet-spot{flag}", f"< {FD_TOL:.0e}", f"{worst2:.3e}",
            bool(worst2 < FD_TOL))
    _merge_tagged(
        "gv5_1_fd_report.csv",
        "protocol,level,scope,state,direction,err_eps_1e-5,err_eps_1e-6,"
        "err_eps_1e-7,sweet_spot,veps_omission_abs,veps_omission_scaled,"
        "masked_rows",
        "k1seed,coarse,", [level],
        [",".join(str(x) for x in ("amended", level) + r)
         for r in fd_rows])

    # -- phase 6: compare vs the COMMITTED GV3.1 loose result -----------------
    print(f"--- phase 6: compare vs committed GV3.1 {level} (the loose "
          "numbers read, never recomputed) ---", flush=True)
    phi_f = _phi_cut_of(pack, phi_free_f, gamma_f)
    f_au = wall_force_coefficients(mc.nodes, mc.elements, case.wall_faces,
                                   phi_f, alpha_deg=ALPHA, s_ref=s_ref,
                                   m_inf=M_INF)
    cl_p_au = float(f_au["cl"])
    cl_kj_au = float(2.0 * gamma_f[0])  # the 2.5-D cl = 2 Gamma/c check
    f_seed = wall_force_coefficients(mc.nodes, mc.elements,
                                     case.wall_faces, res.phi,
                                     alpha_deg=ALPHA, s_ref=s_ref,
                                     m_inf=M_INF)
    cl_p_seed = float(f_seed["cl"])
    print(f"    polish endpoint: cl_p={cl_p_au:.6f} "
          f"cl_kj(2*Gamma)={cl_kj_au:.6f} | seed cl_p={cl_p_seed:.6f} "
          f"(committed loose final {loose_cl:.6f}, k=0 "
          f"{loose_k0_cl:.6f})", flush=True)

    # per-station delta* max-diff vs the committed loose profile
    sta = case.stations
    ds_station = {}
    for side_val, side_name in ((1, "upper"), (-1, "lower")):
        for r in range(len(sta.xc)):
            m = (sta.station_of == r) & (sta.side_node == side_val)
            if np.any(m):
                ds_station[(round(float(sta.xc[r]), 6), side_name)] = \
                    float(np.mean(ds_f[m]))
    prof_rows, worst_ds = [], 0.0
    with open(os.path.join(V3_RESULTS,
                           f"gv3_1_profiles_{level}.csv")) as f:
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
        ("N_polish (aug polish iterations)", n_polish, "", ""),
        ("N_total = N_loose + N_polish vs loose outer",
         res.n_outer + n_polish, loose_n_outer, ""),
        ("converged", res_n["converged"], loose_converged, ""),
        ("final |F_phi|_max", f"{hist[-1]['block_max'][0]:.6e}", "", ""),
        ("final |F_gamma|_max", f"{hist[-1]['block_max'][1]:.6e}", "", ""),
        ("final |F_BL|_max", f"{hist[-1]['block_max'][2]:.6e}", "", ""),
        ("ds_change_last (polish) / ds_change_rel (loose final)",
         f"{res_n['ds_change_last']:.6e}", f"{loose_ds_change:.6e}", ""),
        ("final ds_max", f"{float(np.max(np.abs(ds_f))):.6e}",
         f"{loose_ds_max:.6e}",
         f"{abs(float(np.max(np.abs(ds_f))) - loose_ds_max):.3e}"),
        ("final per-station max|delta* diff|", "", "", f"{worst_ds:.3e}"),
        ("cl_p (pressure integral)", f"{cl_p_au:.6f}", f"{loose_cl:.6f}",
         f"{abs(cl_p_au - loose_cl):.3e}"),
        ("cl_kj (2*Gamma, design.md Sec 9)", f"{cl_kj_au:.6f}",
         "n/a (not committed)", ""),
        ("seed cl_p (loose-converged) vs committed loose final",
         f"{cl_p_seed:.6f}", f"{loose_cl:.6f}",
         f"{abs(cl_p_seed - loose_cl):.3e}"),
        ("IBL residual floor (seed / loose final)",
         f"{ibl.get('final_residual'):.3e}", f"{loose_ibl_floor:.3e}", ""),
        ("final |F|_inf (aug) / ds_change_rel (loose)",
         f"{max(hist[-1]['block_max']):.3e}", f"{loose_ds_change:.3e}", ""),
    ]
    caveat = (" (seed scatter caveat: the committed medium fixed point "
              "is one IBL-floor trajectory sample -- "
              "gv5_1_medium_seed_diagnosis.md)" if level == "medium"
              else "")
    if caveat:
        cmp_rows.append(
            ("seed-vs-committed offset = IBL-floor trajectory scatter "
             "(finding, not a defect)", "see summary.csv",
             "see gv5_1_medium_seed_diagnosis.md", ""))
    _merge_tagged(
        "gv5_1_compare.csv",
        "protocol,level,metric,augmented,loose_committed,abs_diff",
        "k1seed,coarse,", [level],
        [",".join(str(x) for x in ("amended", level) + r)
         for r in cmp_rows])
    _record("GV5.1", f"{level} compare: N_polish and N_total vs "
            "committed loose outer", "recorded",
            f"N_polish={n_polish} converged={res_n['converged']} "
            f"N_total={res.n_outer + n_polish} vs "
            f"n_outer={loose_n_outer}", None)
    _record("GV5.1", f"{level} compare: final per-station max|delta* "
            f"diff|{caveat}", "recorded", f"{worst_ds:.3e}", None)
    _record("GV5.1", f"{level} compare: cl_p augmented vs loose "
            f"committed{caveat}", "recorded",
            f"{cl_p_au:.6f} vs {loose_cl:.6f}", None)

    walls["leg"] = time.perf_counter() - t_leg
    _record("GV5.1", f"{level} wall times per phase", "recorded",
            f"loose={walls['loose']:.0f}s pack={walls['pack']:.0f}s "
            f"fd_seed={walls['fd_seed']:.0f}s newton="
            f"{walls['newton']:.0f}s fd_solution="
            f"{walls['fd_solution']:.0f}s leg_total={walls['leg']:.0f}s",
            None)
    _write_summary(level)
    print(f"GV5.1 amended [{level}] done in {walls['leg']:.0f}s",
          flush=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["coarse"],
                    choices=["coarse", "medium"],
                    help="mesh level(s) to run (default: coarse)")
    args = ap.parse_args()

    # the committed k1seed history is kept under its own name (both
    # protocols' records); idempotent
    src = os.path.join(RESULTS, "gv5_1_newton_history_coarse.csv")
    dst = os.path.join(RESULTS, "gv5_1_newton_history_coarse_k1seed.csv")
    if os.path.exists(src) and not os.path.exists(dst):
        os.rename(src, dst)
        print(f"renamed {src} -> {dst} (k1seed record kept)", flush=True)

    for level in args.levels:
        run_amended_leg(level)

    n_fail = sum(1 for *_, v in SUMMARY if v == "FAIL")
    n_pass = sum(1 for *_, v in SUMMARY if v == "PASS")
    n_rec = sum(1 for *_, v in SUMMARY if v == "RECORDED")
    print(f"\nGV5.1 amended: {n_pass} PASS / {n_fail} FAIL / {n_rec} "
          f"RECORDED", flush=True)
    if n_fail:
        print("HONEST FAIL -- see summary.csv", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
