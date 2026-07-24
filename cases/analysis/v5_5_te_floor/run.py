#!/usr/bin/env python
"""GV5.5 — TE-band (B,delta) formula-level floor breaking.

Pre-registered route (a), variant V1: TE outflow-row replacement — for each
TE node i with upstream partner up (pyfp3d.viscous.coupling.te_outflow_pairs),
the delta-carrier row 6i+0 and H-carrier row 6i+2 of the IBL system are
replaced by first-order extrapolation

    R[6i+0] = U[6i+0] - U[6*up+0]   (delta)
    R[6i+2] = U[6i+1] - U[6*up+1]   (H)

with matching exact Jacobian rows (out-of-pattern entries raise).  Flag
default OFF and bit-identical when off (tests/test_v5_te_outflow.py).

Protocol per level (coarse=crash-stop, medium=binding):
  1. rebuild the loose state and amended GV5.1 seeds verbatim
     (gv5_1b runner build_loose_state);
  2. V0 control: flag OFF, diagnostic Q7 protocol
     (solver.solve(U_seed, tol=ibl_tol, max_iter=ibl_max_iter)) — must
     reproduce the committed seed floor within 1 %, else the (b) read
     switches to the seed's own flag-OFF floor (medium 8-thread scatter
     clause) and the VERDICT marks it;
  3. band (a) live FD: jacobian-vector action of the VARIANT system at the
     seed vs central FD, 4 random directions, rel err < 1e-5;
  4. V1: flag ON, same protocol from the same seed.
       m1 = max|R_variant(U_V1)|   (variant system's own floor — recorded)
       m2 = max|R_orig(U_V1)|      (original system on V1 terminal — BINDING)
     PASS = m2 <= 0.5 * committed floor on BOTH levels + guards green;
     0.5-0.9x / >=0.9x / worse => RECORDED bands;
  5. secondary read: tight polish (GV5.1b recipe, floor_stop) on the variant
     solver, terminal F_BL vs the committed GV5.1b final;
  6. guards (c): GV1.1 plate H band with flag ON (lam 2.55-2.75,
     turb 1.2-2.0), loose smoke flag ON (<=10 outer, cl within committed
     cl +/-2.5 %).  The tight fleet + full suite flag-OFF run is executed
     separately (pytest).

Threading default 16 (agent-rules); this session ran with 8 (temporary).

Run:  python cases/analysis/v5_5_te_floor/run.py
"""
import importlib.util
import os
import sys
from pathlib import Path

for _var in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_var, "16")

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pyfp3d.viscous.coupling import te_outflow_pairs  # noqa: E402
from pyfp3d.viscous.ibl3 import IBL3Solver  # noqa: E402
import pyfp3d.viscous.closures as C  # noqa: E402

from tests.v5_state import UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR  # noqa: E402


def _load(mod_name, rel):
    spec = importlib.util.spec_from_file_location(mod_name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gb = _load("gv5_1b_run", "cases/analysis/v5_1b_scaled_newton/run.py")
v1r = _load("v1_ibl3_run", "cases/analysis/v1_ibl3_standalone/run.py")
from pyfp3d.viscous import tight_driver as td  # noqa: E402


def _gv51b_final(level):
    """The GV5.1b committed final row (f_bl_max, merit) -- read, never
    recomputed (same source as the gv5_1c runner helper)."""
    import csv

    path = (ROOT / "cases/analysis/v5_1b_scaled_newton/results"
            / f"newton_history_{level}.csv")
    with open(path) as f:
        last = list(csv.DictReader(f))[-1]
    return float(last["f_bl_max"]), float(last["merit"])

RESULTS = Path(__file__).resolve().parent / "results"

LEVELS = ("coarse", "medium")

# GV1.1 plate protocol constants (runner v1 defaults)
PLATE_UE = lambda x: np.array([v1r.Q, 0.0, 0.0])  # noqa: E731
PLATE_BANDS = {"lam": (2.55, 2.75), "turb": (1.2, 2.0)}


def _write_lines(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_summary(rows, name="summary.csv"):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    lines = [",".join(keys)]
    for r in rows:
        lines.append(",".join(str(r.get(k, "")) for k in keys))
    _write_lines(RESULTS / name, lines)


def _record(**kw):
    return kw


def _fd_check(solver, U, n_dirs=4, eps=1e-6, seed=1234):
    """band (a) live FD: jacobian action vs central differences (the
    tests/test_v5_te_outflow.py idiom, (n, 6) state shape)."""
    R0, J = solver.residual_jacobian(U)
    rng = np.random.default_rng(seed)
    errs = []
    for _ in range(n_dirs):
        v = rng.standard_normal(U.shape)
        v /= np.max(np.abs(v))
        fd = (solver.residual(U + eps * v)
              - solver.residual(U - eps * v)).ravel() / (2.0 * eps)
        exact = J @ v.ravel()
        errs.append(float(np.max(np.abs(exact - fd))
                          / max(np.max(np.abs(exact)), 1.0e-12)))
    return errs


def _build_variant_solver(st):
    cfg = st["cfg"]
    case = st["case"]
    pairs = te_outflow_pairs(case)
    solver1 = IBL3Solver(
        st["sm"], st["ue_surf"], st["rho_e"], st["mu"], st["mach_e"],
        case.turbulent_flags, st["inflow_mask"], st["inflow_state"],
        eps_diff=cfg.eps_diff, eps_diff_s=cfg.eps_diff_s,
        te_pairs=pairs, te_extrapolate=True,
    )
    return solver1, pairs


def _plate_guard(te):
    """GV1.1 plate H bands with the flag forced on; returns dict of rows."""
    out = {}
    lam = v1r.run_fe(160, 32, False, PLATE_UE,
                     lambda x: C.blasius_seed(x, q=v1r.Q, rho=v1r.RHO, mu=v1r.MU),
                     te=te)
    H_lam = v1r._centerline(lam, C.OUT_H1)[1]
    turb = v1r.run_fe(100, 16, True, PLATE_UE, v1r._turb_seed, te=te)
    H_turb = v1r._centerline(turb, C.OUT_H1)[1]
    out["lam_H_min"] = float(H_lam.min())
    out["lam_H_max"] = float(H_lam.max())
    out["turb_H_min"] = float(H_turb.min())
    out["turb_H_max"] = float(H_turb.max())
    lo, hi = PLATE_BANDS["lam"]
    out["lam_ok"] = bool(out["lam_H_min"] >= lo and out["lam_H_max"] <= hi)
    lo, hi = PLATE_BANDS["turb"]
    out["turb_ok"] = bool(out["turb_H_min"] >= lo and out["turb_H_max"] <= hi)
    return out


def _loose_smoke(st, level):
    """Flag-ON loose-loop smoke (the build_loose_state recipe verbatim but
    te_extrapolate=True): <=10 outer, converged, cl within the committed
    flag-OFF cl +/-2.5 %."""
    import dataclasses

    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.surface import wall_force_coefficients
    from pyfp3d.viscous.coupling import (
        build_airfoil_case,
        make_picard_lifting_driver,
        run_loose_coupling,
    )
    from tests.v5_state import ALPHA, M_INF, NACA_DIR

    mc, wc = cut_wake(read_mesh(os.path.join(NACA_DIR, f"{level}.msh")))
    cfg2 = dataclasses.replace(st["cfg"], te_extrapolate=True)
    case2 = build_airfoil_case(
        mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg2)
    dz = float(np.ptp(mc.nodes[:, 2]))
    s_ref = 1.0 * dz

    def probe(phi, gamma, k):
        f = wall_force_coefficients(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
            alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        return {"cl": f["cl"], "cd_p": f["cd_pressure"]}

    driver = make_picard_lifting_driver(mc, wc, M_INF, ALPHA)
    res2 = run_loose_coupling(driver, case2, cfg2, probe=probe)
    cl_ref = float(st["loose_res"].history[-1]["cl"])
    cl_new = float(res2.history[-1]["cl"])
    n_outer = int(res2.n_outer)
    rel = abs(cl_new - cl_ref) / max(abs(cl_ref), 1e-300)
    return {
        "loose_outer": n_outer,
        "loose_converged": bool(res2.converged),
        "cl_ref": cl_ref,
        "cl_te": cl_new,
        "cl_rel": rel,
        "loose_ok": bool(n_outer <= 10 and res2.converged and rel <= 0.025),
    }


def _tight_polish(st, solver):
    """GV5.1b tight polish on a given solver; returns terminal F_BL dict."""
    st1 = dict(st)
    st1["solver"] = solver
    pack = td.build_tight_pack(st1, UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR)
    res = td.newton_tight(pack, x0=pack.x_base(), max_iter=10, tol=1e-8,
                          tol_abs=1e-10, line_search=True, scaling="rowcol",
                          lm_damping=True, floor_stop=True, verbose=False)
    F = pack.F_BL(res["x"])
    return {
        "tight_n_iter": res["n_iter"],
        "tight_F_max": float(np.max(np.abs(F))),
        "tight_merit": float(0.5 * np.dot(F, F)),
        "tight_converged": bool(res["converged"]),
    }


def _anatomy(solver0, U, path):
    """Per-node original-system residual anatomy at U (diag-style CSV)."""
    R = solver0.residual(U)
    n = U.size // 6
    delta = np.abs(R[0::6][:n])
    H_row = np.abs(R[2::6][:n])
    ke = np.abs(R[1::6][:n])
    df = pd.DataFrame({
        "node": np.arange(n),
        "abs_delta_row": delta,
        "abs_H_row": H_row,
        "abs_ke_row": ke,
        "abs_max_row": np.maximum(np.maximum(delta, H_row), ke),
    })
    df.to_csv(path, index=False)


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    guard_rows = []
    for level in LEVELS:
        print(f"\n==================== GV5.5 level={level} ====================")
        st = gb.build_loose_state(level)
        cfg = st["cfg"]
        solver0 = st["solver"]
        seed = st["U"].copy()
        committed = gb._gv51_seed_ibl_floor(level)

        # --- V0 control: flag OFF reproduces committed floor ----------------
        U_v0, info0 = solver0.solve(seed, tol=cfg.ibl_tol,
                                    max_iter=cfg.ibl_max_iter, verbose=True)
        floor0 = float(np.max(np.abs(solver0.residual(U_v0))))
        control_rel = abs(floor0 - committed) / committed
        floor_ref = committed if control_rel <= 0.01 else floor0
        print(f"  V0 floor={floor0:.3e} committed={committed:.3e} "
              f"control_rel={control_rel:.2e} "
              f"{'OK' if control_rel <= 0.01 else 'SCATTER-CLAUSE'}")

        # --- V1 variant ------------------------------------------------------
        solver1, pairs = _build_variant_solver(st)
        fd_errs = _fd_check(solver1, seed)
        fd_ok = max(fd_errs) < 1e-5
        print(f"  band(a) FD errs={['%.2e' % e for e in fd_errs]} ok={fd_ok}")

        U_v1, info1 = solver1.solve(seed, tol=cfg.ibl_tol,
                                    max_iter=cfg.ibl_max_iter, verbose=True)
        m1 = float(np.max(np.abs(solver1.residual(U_v1))))
        m2 = float(np.max(np.abs(solver0.residual(U_v1))))
        ratio = m2 / floor_ref
        print(f"  V1 m1={m1:.3e} m2={m2:.3e} ratio_vs_floor={ratio:.3f}")

        # --- secondary tight polish -----------------------------------------
        tight_v1 = _tight_polish(st, solver1)
        committed_final, committed_merit = _gv51b_final(level)
        print(f"  tight V1: F_max={tight_v1['tight_F_max']:.3e} "
              f"(committed GV5.1b final {committed_final:.3e})")

        # --- guards -----------------------------------------------------------
        plate = _plate_guard(te=True)
        loose = _loose_smoke(st, level)
        guards_ok = bool(plate["lam_ok"] and plate["turb_ok"]
                         and loose["loose_ok"] and fd_ok)
        print(f"  guards: plate lam [{plate['lam_H_min']:.3f},{plate['lam_H_max']:.3f}] "
              f"turb [{plate['turb_H_min']:.3f},{plate['turb_H_max']:.3f}] "
              f"loose outer={loose['loose_outer']} cl_rel={loose['cl_rel']:.2e} "
              f"ok={guards_ok}")

        # --- artifacts --------------------------------------------------------
        n_hist = max(len(info0["residual_history"]), len(info1["residual_history"]))
        h0 = list(info0["residual_history"]) + [np.nan] * (n_hist - len(info0["residual_history"]))
        h1 = list(info1["residual_history"]) + [np.nan] * (n_hist - len(info1["residual_history"]))
        pd.DataFrame({"iter": range(n_hist), "V0_flag_off": h0,
                      "V1_te_extrapolate": h1}).to_csv(
            RESULTS / f"floor_probe_{level}.csv", index=False)
        _anatomy(solver0, U_v1, RESULTS / f"residual_anatomy_{level}.csv")

        rows.append(_record(
            level=level, mode="V0_control", floor=floor0,
            committed_floor=committed, control_rel=control_rel,
            converged=info0["converged"], n_iter=info0["n_iter"]))
        rows.append(_record(
            level=level, mode="V1_te_extrapolate", m1_variant_floor=m1,
            m2_orig_residual=m2, floor_ref=floor_ref, ratio_vs_floor=ratio,
            converged=info1["converged"], n_iter=info1["n_iter"],
            fd_max_err=max(fd_errs),
            tight_F_max=tight_v1["tight_F_max"],
            tight_merit=tight_v1["tight_merit"],
            tight_n_iter=tight_v1["tight_n_iter"],
            committed_gv51b_final=committed_final))
        guard_rows.append(_record(level=level, fd_ok=fd_ok, guards_ok=guards_ok,
                                  **plate, **loose))

        del st, solver0, solver1
        import gc
        gc.collect()

    _write_summary(rows)
    _write_summary(guard_rows, "guards.csv")

    print("\n--- summary ---")
    for r in rows:
        print(r)
    print("--- guards ---")
    for r in guard_rows:
        print(r)


if __name__ == "__main__":
    main()
