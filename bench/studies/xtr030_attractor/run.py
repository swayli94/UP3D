r"""Phase 7 S5: controlled xtr030 loose-coupling attractor attribution.

Pre-registration: docs/dev_phase_seven/20261103-0245-xtr030-attractor-prereg.md.
Every coupling leg calls the production run_loose_coupling. Seed/frozen-BC
injection and IBL recording are wrappers local to the worker processes.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WORK = HERE / "_work"
sys.path.insert(0, str(REPO))

MESH = REPO / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"
MAX_COUPLING_LEGS = 10
TIME_BUDGET_S = 75.0 * 60.0
GATE_NODES = (
    "tests/D/test_D08_rans_naca0012.py::TestCommittedEvidenceIsLoadBearing",
    "tests/D/test_D13_ibl_vs_xfoil.py::TestCommittedEvidenceIsLoadBearing",
)
SUMMARY_FIELDS = (
    "leg", "kind", "status", "level", "x_tr", "blas_threads", "repeat",
    "seed_source", "bc_mode", "tol_ds", "n_outer_max", "converged", "n_outer",
    "cl", "cd_p", "ds_change_rel", "ds_max", "mdot_max", "ibl_n_iter",
    "ibl_converged", "ibl_final_residual", "fp_converged", "fp_residual",
    "k0_cl", "k0_phi_sha256", "k0_gamma_sha256", "phi_sha256", "U_sha256",
    "delta_sha256", "wall_s", "n_threads", "mesh_sha256", "state_file",
    "first_file", "error",
)


def write_json(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def read_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def array_sha(a):
    import numpy as np
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def f12(x):
    return "" if x is None else f"{float(x):.12e}"


def read_csv(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n",
                           extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(rows)


def append_summary(row):
    path = HERE / "summary.csv"
    rows = [r for r in read_csv(path) if r["leg"] != row["leg"]]
    rows.append(row)
    rows.sort(key=lambda r: r["leg"])
    write_csv(path, SUMMARY_FIELDS, rows)


def require_blas(expected, reader):
    actual = reader()
    if actual != expected:
        raise RuntimeError(f"BLAS provenance failure: expected {expected}, got {actual}")
    return actual


def fixed_point_label(dcl, dds, ecl, eds, converged=True):
    if not converged:
        return "UNDEFINED_NONCONVERGED"
    if dcl <= 1e-6 and dds <= 1e-8 and ecl <= 1e-6 and eds <= 1e-8:
        return "SAME"
    distinct = (dcl > 1e-5 and dcl > 20 * ecl) or (dds > 1e-7 and dds > 20 * eds)
    return "DISTINCT" if distinct else "UNDEFINED_GREY"


def separated(dcl, dds, ecl, eds):
    return dcl > max(1e-4, 20 * ecl) or dds > max(1e-7, 20 * eds)


def make_recorder(original, captured, seed_U=None, frozen_bc=None):
    import numpy as np

    class Recorder(original):
        def __init__(self, smesh, u_e, rho, mu, mach, turbulent_flags,
                     inflow_mask, inflow_state, *args, **kwargs):
            if frozen_bc is not None:
                inflow_mask = np.asarray(frozen_bc["inflow_mask"], dtype=bool)
                inflow_state = np.asarray(frozen_bc["inflow_state"], dtype=float)
            if not captured.get("constructor"):
                captured["constructor"] = dict(
                    ue=np.array(u_e, copy=True), rho=np.array(rho, copy=True),
                    mu=np.array(mu, copy=True), mach=np.array(mach, copy=True),
                    flags=np.array(turbulent_flags, copy=True),
                    inflow_mask=np.array(inflow_mask, copy=True),
                    inflow_state=np.array(inflow_state, copy=True),
                )
            super().__init__(smesh, u_e, rho, mu, mach, turbulent_flags,
                             inflow_mask, inflow_state, *args, **kwargs)

        def solve(self, U0, *args, **kwargs):
            first = not captured.get("solved")
            if first and seed_U is not None:
                U0 = np.asarray(seed_U, dtype=float)
            if first:
                captured["U0"] = np.array(U0, copy=True)
            U, info = super().solve(U0, *args, **kwargs)
            if first:
                captured["U_first"] = np.array(U, copy=True)
                captured["info_first"] = dict(info)
                captured["solved"] = True
            return U, info

    return Recorder


def configure_worker(blas):
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ["NUMBA_NUM_THREADS"] = "16"
    os.environ["OMP_NUM_THREADS"] = "16"
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["OPENBLAS_NUM_THREADS"] = str(blas)
    if blas == 1:
        os.environ.pop("PYFP3D_ALLOW_BLAS_THREADS", None)
    else:
        os.environ["PYFP3D_ALLOW_BLAS_THREADS"] = "1"


def write_leg_tables(leg, res, case):
    hist_fields = (
        "leg", "k", "cl", "cd_p", "ds_max", "ds_change_rel", "ds_neg_floored",
        "mdot_max", "ibl_n_iter", "ibl_converged", "ibl_final_residual",
        "inflow_n_pinned", "inflow_q_ref",
    )
    rows = []
    for h in res.history:
        row = {"leg": leg}
        for key in hist_fields[1:]:
            value = h.get(key, "")
            row[key] = f12(value) if isinstance(value, (float, int)) and key not in (
                "k", "ibl_n_iter", "inflow_n_pinned") else value
        rows.append(row)
    write_csv(HERE / f"history_{leg}.csv", hist_fields, rows)

    st = case.stations
    xc = st.xc[st.station_of]
    with open(HERE / f"delta_{leg}.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["leg", "node", "x_c", "side", "delta_star"])
        for i in range(len(res.delta_star)):
            w.writerow([leg, i, f12(xc[i]), int(st.side_node[i]),
                        f12(res.delta_star[i])])


def coupling_worker(spec):
    configure_worker(spec["blas"])
    started = time.monotonic()
    result_path = WORK / f"result_{spec['leg']}.json"
    try:
        import numpy as np
        from pyfp3d import blas_thread_count
        require_blas(spec["blas"], blas_thread_count)
        from pyfp3d.mesh.reader import read_mesh
        from pyfp3d.mesh.wake_cut import cut_wake
        from pyfp3d.post.surface import wall_force_coefficients
        import pyfp3d.viscous.coupling as coupling

        leg = spec["leg"]
        mc, wc = cut_wake(read_mesh(str(MESH)))
        cfg = coupling.CouplingConfig(
            re_chord=3.0e6, m_inf=0.5, alpha_deg=2.0,
            x_tr_upper=spec["x_tr"], x_tr_lower=spec["x_tr"],
            n_outer_max=spec["n_outer_max"], tol_ds=spec["tol_ds"],
        )
        case = coupling.build_airfoil_case(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
        dz = float(np.ptp(mc.nodes[:, 2]))
        base_driver = coupling.make_picard_lifting_driver(mc, wc, 0.5, 2.0)
        seed = np.load(spec["seed_state"]) if spec.get("seed_state") else None
        frozen = np.load(spec["bc_first"]) if spec.get("bc_first") else None
        k0 = {}

        def driver(rhs, fp_seed):
            phi, gamma, info = base_driver(rhs, fp_seed)
            if rhs is None and fp_seed is None:
                if seed is not None:
                    phi = np.array(seed["phi"], copy=True)
                    gamma = np.array(seed["gamma"], copy=True)
                    info = dict(info)
                    info["instrument_seed_override"] = True
                k0["phi"] = np.array(phi, copy=True)
                k0["gamma"] = np.array(gamma, copy=True)
            return phi, gamma, info

        def probe(phi, gamma, k):
            f = wall_force_coefficients(
                mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                alpha_deg=2.0, u_inf=1.0, s_ref=dz, m_inf=0.5)
            return {"cl": float(f["cl"]), "cd_p": float(f["cd_pressure"])}

        captured = {}
        recorder = make_recorder(
            coupling.IBL3Solver, captured,
            seed_U=None if seed is None else np.asarray(seed["U"]),
            frozen_bc=frozen,
        )
        original = coupling.IBL3Solver
        coupling.IBL3Solver = recorder
        try:
            res = coupling.run_loose_coupling(driver, case, cfg, probe=probe)
        finally:
            coupling.IBL3Solver = original
        require_blas(spec["blas"], blas_thread_count)

        state_rel = f"_work/state_{leg}.npz"
        first_rel = f"_work/first_{leg}.npz"
        np.savez_compressed(
            HERE / state_rel, phi=res.phi, gamma=res.gamma, U=res.U,
            delta_star=res.delta_star, m_dot=res.m_dot, ue_surf=res.ue_surf,
            outs=res.outs, k0_phi=k0["phi"], k0_gamma=k0["gamma"],
        )
        first = captured["constructor"]
        np.savez_compressed(
            HERE / first_rel, **first, U0=captured["U0"],
            U_first=captured["U_first"],
        )
        write_leg_tables(leg, res, case)
        last = res.history[-1]
        info_first = captured["info_first"]
        result = dict(
            leg=leg, kind=spec["kind"], status="ok", level="medium",
            x_tr=spec["x_tr"], blas_threads=blas_thread_count(),
            repeat=spec.get("repeat", ""), seed_source=spec.get("seed_source", ""),
            bc_mode=spec.get("bc_mode", "native"), tol_ds=f12(cfg.tol_ds),
            n_outer_max=cfg.n_outer_max, converged=int(res.converged),
            n_outer=res.n_outer, cl=f12(last["cl"]), cd_p=f12(last["cd_p"]),
            ds_change_rel=f12(last["ds_change_rel"]), ds_max=f12(last["ds_max"]),
            mdot_max=f12(last["mdot_max"]),
            ibl_n_iter=info_first.get("n_iter", ""),
            ibl_converged=int(info_first["converged"]),
            ibl_final_residual=f12(info_first["final_residual"]),
            fp_converged=int(bool(res.driver_info.get("converged", True))),
            fp_residual=f12(res.driver_info.get("residual_norm")),
            k0_cl=f12(res.history[0]["cl"]),
            k0_phi_sha256=array_sha(k0["phi"]),
            k0_gamma_sha256=array_sha(k0["gamma"]),
            phi_sha256=array_sha(res.phi), U_sha256=array_sha(res.U),
            delta_sha256=array_sha(res.delta_star),
            wall_s=f12(time.monotonic() - started), n_threads=16,
            mesh_sha256=sha256_file(MESH), state_file=state_rel,
            first_file=first_rel, error="",
        )
        write_json(result_path, result)
        print("RESULT " + json.dumps(result), flush=True)
    except Exception as exc:
        failure = dict(leg=spec["leg"], kind=spec["kind"], status="failed",
                       error=f"{type(exc).__name__}: {exc}",
                       traceback=traceback.format_exc(),
                       wall_s=time.monotonic() - started)
        write_json(result_path, failure)
        raise


def replay_worker(spec):
    configure_worker(spec["blas"])
    started = time.monotonic()
    try:
        import numpy as np
        from pyfp3d import blas_thread_count
        require_blas(spec["blas"], blas_thread_count)
        from pyfp3d.mesh.reader import read_mesh
        from pyfp3d.mesh.wake_cut import cut_wake
        from pyfp3d.viscous.coupling import CouplingConfig, build_airfoil_case
        from pyfp3d.viscous.ibl3 import IBL3Solver

        mc, _wc = cut_wake(read_mesh(str(MESH)))
        cfg = CouplingConfig(re_chord=3.0e6, m_inf=0.5, alpha_deg=2.0,
                             x_tr_upper=0.30, x_tr_lower=0.30)
        case = build_airfoil_case(mc.nodes, mc.elements,
                                  mc.boundary_faces["wall"], cfg)
        f = np.load(WORK / "first_native_xtr030_b1_r1.npz")
        solver = IBL3Solver(
            case.sm, f["ue"], f["rho"], f["mu"], f["mach"], f["flags"],
            f["inflow_mask"], f["inflow_state"], eps_diff=cfg.eps_diff,
            eps_diff_s=cfg.eps_diff_s,
        )
        U, info = solver.solve(f["U0"], tol=cfg.ibl_tol,
                               max_iter=cfg.ibl_max_iter)
        np.savez_compressed(WORK / f"replay_b{spec['blas']}.npz", U=U)
        native = f["U_first"]
        result = dict(
            leg=f"replay_b{spec['blas']}", blas_threads=blas_thread_count(),
            status="ok", converged=int(info["converged"]),
            n_iter=info["n_iter"], final_residual=f12(info["final_residual"]),
            native_array_equal=bool(np.array_equal(U, native)),
            native_max_abs=f12(np.max(np.abs(U - native))),
            U_sha256=array_sha(U), wall_s=f12(time.monotonic() - started),
        )
        write_json(WORK / f"result_replay_b{spec['blas']}.json", result)
        print("REPLAY " + json.dumps(result), flush=True)
    except Exception as exc:
        write_json(WORK / f"result_replay_b{spec['blas']}.json",
                   dict(status="failed", error=f"{type(exc).__name__}: {exc}",
                        traceback=traceback.format_exc()))
        raise


def worker_env(blas):
    env = os.environ.copy()
    env.update(PYTHONNOUSERSITE="1", NUMBA_NUM_THREADS="16",
               OMP_NUM_THREADS="16", OPENBLAS_NUM_THREADS=str(blas),
               MPLBACKEND="Agg")
    if blas == 1:
        env.pop("PYFP3D_ALLOW_BLAS_THREADS", None)
    else:
        env["PYFP3D_ALLOW_BLAS_THREADS"] = "1"
    return env


def run_worker(spec, replay=False):
    WORK.mkdir(parents=True, exist_ok=True)
    leg = spec["leg"]
    spec_path = WORK / f"spec_{leg}.json"
    write_json(spec_path, spec)
    out_path = WORK / f"worker_{leg}.txt"
    cmd = [sys.executable, str(Path(__file__).resolve()),
           "--worker-replay" if replay else "--worker", str(spec_path)]
    started = time.monotonic()
    print(f"[run] {leg} (BLAS={spec['blas']})", flush=True)
    try:
        proc = subprocess.run(cmd, cwd=REPO, env=worker_env(spec["blas"]),
                              text=True, capture_output=True, timeout=TIME_BUDGET_S)
        out_path.write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
        print(proc.stdout, flush=True)
        if proc.returncode:
            print(proc.stderr, file=sys.stderr, flush=True)
    except subprocess.TimeoutExpired as exc:
        out_path.write_text(str(exc), encoding="utf-8")
        write_json(WORK / f"result_{leg}.json",
                   dict(leg=leg, kind=spec["kind"], status="failed",
                        error="worker timeout"))
    return time.monotonic() - started


def run_state():
    path = HERE / "run_state.json"
    if path.exists():
        return read_json(path)
    return dict(fresh_coupling_legs=0, cumulative_coupling_s=0.0,
                replay_s=0.0, max_coupling_legs=MAX_COUPLING_LEGS,
                time_budget_s=TIME_BUDGET_S, completed=[])


def save_state(state):
    write_json(HERE / "run_state.json", state)


def run_coupling_leg(spec, state):
    rows = {r["leg"]: r for r in read_csv(HERE / "summary.csv")}
    if spec["leg"] in rows and rows[spec["leg"]]["status"] == "ok":
        print(f"[skip] {spec['leg']} already complete", flush=True)
        return True
    if state["fresh_coupling_legs"] >= MAX_COUPLING_LEGS:
        raise RuntimeError("fresh coupling leg budget exhausted")
    if state["cumulative_coupling_s"] >= TIME_BUDGET_S:
        raise RuntimeError("coupling wall-time budget exhausted")
    elapsed = run_worker(spec)
    state["fresh_coupling_legs"] += 1
    state["cumulative_coupling_s"] += elapsed
    state["completed"].append(spec["leg"])
    save_state(state)
    result_path = WORK / f"result_{spec['leg']}.json"
    result = read_json(result_path)
    if result["status"] == "ok":
        append_summary({k: result.get(k, "") for k in SUMMARY_FIELDS})
    else:
        row = {k: "" for k in SUMMARY_FIELDS}
        row.update(leg=spec["leg"], kind=spec["kind"], status="failed",
                   level="medium", x_tr=spec["x_tr"], blas_threads=spec["blas"],
                   repeat=spec.get("repeat", ""),
                   seed_source=spec.get("seed_source", ""),
                   bc_mode=spec.get("bc_mode", "native"),
                   tol_ds=f12(spec["tol_ds"]), n_outer_max=spec["n_outer_max"],
                   wall_s=f12(elapsed), n_threads=16,
                   mesh_sha256=sha256_file(MESH), error=result.get("error", ""))
        append_summary(row)
    return result["status"] == "ok"


def run_replays(state):
    rows = {r["leg"]: r for r in read_csv(HERE / "replay.csv")}
    for blas in (1, 16):
        leg = f"replay_b{blas}"
        if rows.get(leg, {}).get("status") == "ok":
            continue
        spec = dict(leg=leg, kind="replay", blas=blas)
        elapsed = run_worker(spec, replay=True)
        state["replay_s"] += elapsed
        save_state(state)
        result = read_json(WORK / f"result_{leg}.json")
        rows[leg] = result
        write_csv(HERE / "replay.csv",
                  ("leg", "blas_threads", "status", "converged", "n_iter",
                   "final_residual", "native_array_equal", "native_max_abs",
                   "U_sha256", "wall_s"), list(rows.values()))
        if result["status"] != "ok":
            raise RuntimeError(f"{leg} failed: {result.get('error')}")
        if blas == 1 and not result["native_array_equal"]:
            raise RuntimeError("B1 replay does not reproduce the native first IBL solve; "
                               "the replay instrument is invalid")


def native_specs():
    return [dict(leg=f"native_xtr030_b{blas}_r{repeat}", kind="native",
                 blas=blas, x_tr=0.30, repeat=repeat, tol_ds=1e-3,
                 n_outer_max=10, bc_mode="native")
            for blas in (1, 16) for repeat in (1, 2)]


def control_specs():
    return [dict(leg=f"control_xtr005_b{blas}", kind="control", blas=blas,
                 x_tr=0.05, repeat=1, tol_ds=1e-3, n_outer_max=10,
                 bc_mode="native") for blas in (1, 16)]


def tight_specs():
    specs = []
    for source in ("b1", "b16"):
        for bc_mode in ("native", "commonB1"):
            bc_source = "b1" if bc_mode == "commonB1" else source
            specs.append(dict(
                leg=f"tight_{source}seed_{bc_mode}", kind="tight", blas=1,
                x_tr=0.30, repeat="", seed_source=source, bc_mode=bc_mode,
                tol_ds=1e-8, n_outer_max=30,
                seed_state=str(WORK / f"state_native_xtr030_{source}_r1.npz"),
                bc_first=str(WORK / f"first_native_xtr030_{bc_source}_r1.npz"),
            ))
    return specs


def state_distance(row_a, row_b):
    import numpy as np
    a = np.load(HERE / row_a["state_file"])
    b = np.load(HERE / row_b["state_file"])
    def rms(x):
        return float(np.sqrt(np.mean(np.asarray(x, dtype=float) ** 2)))
    return dict(
        d_cl=abs(float(row_a["cl"]) - float(row_b["cl"])),
        d_delta_max=float(np.max(np.abs(a["delta_star"] - b["delta_star"]))),
        d_phi_max=float(np.max(np.abs(a["phi"] - b["phi"]))),
        d_phi_rms=rms(a["phi"] - b["phi"]),
        d_gamma_max=float(np.max(np.abs(a["gamma"] - b["gamma"]))),
        k0_phi_equal=bool(np.array_equal(a["k0_phi"], b["k0_phi"])),
        k0_gamma_equal=bool(np.array_equal(a["k0_gamma"], b["k0_gamma"])),
        k0_d_phi_max=float(np.max(np.abs(a["k0_phi"] - b["k0_phi"]))),
        k0_d_gamma_max=float(np.max(np.abs(a["k0_gamma"] - b["k0_gamma"]))),
    )


def comparison_row(name, rows, converged=True):
    a, b = name.split("|")[1:3]
    out = dict(comparison=name, leg_a=a, leg_b=b)
    out.update(state_distance(rows[a], rows[b]))
    out["converged_pair"] = int(converged and rows[a]["converged"] == "1"
                                and rows[b]["converged"] == "1")
    return out


def replay_comparison():
    import numpy as np
    a_path = WORK / "replay_b1.npz"
    b_path = WORK / "replay_b16.npz"
    if not a_path.exists() or not b_path.exists():
        return {}
    a, b = np.load(a_path)["U"], np.load(b_path)["U"]
    return dict(replay_array_equal=bool(np.array_equal(a, b)),
                replay_d_U_max=float(np.max(np.abs(a - b))))


def analyze():
    rows = {r["leg"]: r for r in read_csv(HERE / "summary.csv")
            if r["status"] == "ok"}
    comparisons = []

    def add(name):
        row = comparison_row(name, rows)
        comparisons.append(row)
        return row

    repeat_rows = {}
    if all(f"native_xtr030_b1_r{i}" in rows for i in (1, 2)):
        repeat_rows["b1"] = add("repeat_b1|native_xtr030_b1_r1|native_xtr030_b1_r2")
    if all(f"native_xtr030_b16_r{i}" in rows for i in (1, 2)):
        repeat_rows["b16"] = add("repeat_b16|native_xtr030_b16_r1|native_xtr030_b16_r2")
    cross = None
    if all(f"native_xtr030_b{b}_r1" in rows for b in (1, 16)):
        cross = add("cross_native|native_xtr030_b1_r1|native_xtr030_b16_r1")
    tight_repeat = None
    if all(k in rows for k in ("tight_b1seed_native", "tight_b1seed_commonB1")):
        tight_repeat = add("repeat_tight|tight_b1seed_native|tight_b1seed_commonB1")
    tight_native = None
    tight_common = None
    if all(k in rows for k in ("tight_b1seed_native", "tight_b16seed_native")):
        tight_native = add("tight_native|tight_b1seed_native|tight_b16seed_native")
    if all(k in rows for k in ("tight_b1seed_commonB1", "tight_b16seed_commonB1")):
        tight_common = add("tight_common|tight_b1seed_commonB1|tight_b16seed_commonB1")

    comp_fields = (
        "comparison", "leg_a", "leg_b", "d_cl", "d_delta_max", "d_phi_max",
        "d_phi_rms", "d_gamma_max", "k0_phi_equal", "k0_gamma_equal",
        "k0_d_phi_max", "k0_d_gamma_max", "converged_pair", "repeat_e_cl",
        "repeat_e_ds", "separated", "fixed_point_label",
    )
    metrics = {}
    ecl = max((r["d_cl"] for r in repeat_rows.values()), default=0.0)
    eds = max((r["d_delta_max"] for r in repeat_rows.values()), default=0.0)
    if tight_repeat is not None:
        ecl = max(ecl, tight_repeat["d_cl"])
        eds = max(eds, tight_repeat["d_delta_max"])
    for row in comparisons:
        row["repeat_e_cl"] = f12(ecl)
        row["repeat_e_ds"] = f12(eds)
        if row["comparison"].startswith("cross_native"):
            row["separated"] = int(separated(row["d_cl"], row["d_delta_max"], ecl, eds))
        else:
            row["separated"] = ""
        if row["comparison"].startswith("tight_"):
            row["fixed_point_label"] = fixed_point_label(
                row["d_cl"], row["d_delta_max"], ecl, eds,
                converged=bool(row["converged_pair"]))
        else:
            row["fixed_point_label"] = ""
    write_csv(HERE / "comparisons.csv", comp_fields, comparisons)

    metrics["repeat_e_cl"] = ecl
    metrics["repeat_e_ds"] = eds
    if cross is not None:
        metrics.update(native_separated=separated(cross["d_cl"], cross["d_delta_max"], ecl, eds),
                       native_d_cl=cross["d_cl"], native_d_delta_max=cross["d_delta_max"],
                       k0_phi_equal=cross["k0_phi_equal"],
                       k0_gamma_equal=cross["k0_gamma_equal"])
    metrics.update(replay_comparison())
    if tight_native is not None:
        metrics["tight_native_label"] = fixed_point_label(
            tight_native["d_cl"], tight_native["d_delta_max"], ecl, eds,
            converged=bool(tight_native["converged_pair"]))
    if tight_common is not None:
        metrics["tight_common_label"] = fixed_point_label(
            tight_common["d_cl"], tight_common["d_delta_max"], ecl, eds,
            converged=bool(tight_common["converged_pair"]))
    verdict = "PENDING"
    if cross is not None:
        if not metrics["native_separated"]:
            verdict = "RECORDED_NO_B16_SEPARATION"
        elif tight_native is None or tight_common is None:
            verdict = "PENDING_TIGHT"
        elif metrics["tight_native_label"] == "SAME":
            verdict = "PASS_EARLY_STOP"
        elif (metrics["tight_native_label"] == "DISTINCT"
              and metrics["tight_common_label"] == "SAME"):
            verdict = "PASS_FROZEN_INFLOW_PATH"
        elif metrics["tight_common_label"] == "DISTINCT":
            verdict = "PASS_MULTIPLE_FIXED_POINTS"
        else:
            verdict = "UNDEFINED_TIGHT"
    metrics["verdict"] = verdict
    write_csv(HERE / "verdict_metrics.csv", ("metric", "value"),
              [dict(metric=k, value=v) for k, v in metrics.items()])
    return metrics


def self_test():
    import numpy as np
    rows = []
    def record(name, ok, detail=""):
        rows.append(dict(case=name, passed=bool(ok), detail=detail))

    record("fixed_point_same", fixed_point_label(1e-8, 1e-10, 0.0, 0.0) == "SAME")
    record("fixed_point_distinct", fixed_point_label(1e-3, 1e-5, 0.0, 0.0) == "DISTINCT")
    record("fixed_point_grey", fixed_point_label(3e-6, 1e-9, 0.0, 0.0) == "UNDEFINED_GREY")
    record("fixed_point_nonconverged",
           fixed_point_label(0.0, 0.0, 0.0, 0.0, False) == "UNDEFINED_NONCONVERGED")
    record("separation_teeth", separated(2e-4, 0.0, 0.0, 0.0)
           and not separated(1e-6, 0.0, 0.0, 0.0))
    try:
        require_blas(1, lambda: None)
        record("blas_missing_is_loud", False)
    except RuntimeError:
        record("blas_missing_is_loud", True)

    class Dummy:
        def __init__(self, *args):
            self.args = args
        def solve(self, U0, *args, **kwargs):
            return np.asarray(U0) + 1.0, dict(converged=True, n_iter=1,
                                              final_residual=0.0)

    captured = {}
    recorder = make_recorder(Dummy, captured, seed_U=np.full((2, 6), 7.0),
                             frozen_bc=dict(inflow_mask=np.array([True, False]),
                                            inflow_state=np.full((1, 6), 3.0)))
    args = (None, np.zeros((2, 3)), np.ones(2), 1e-6, np.ones(2),
            np.zeros(2), np.array([False, True]), np.zeros((1, 6)))
    solver = recorder(*args)
    U, _ = solver.solve(np.zeros((2, 6)))
    record("seed_override_has_teeth", np.array_equal(U, np.full((2, 6), 8.0)))
    record("frozen_bc_override_has_teeth",
           np.array_equal(captured["constructor"]["inflow_mask"], [True, False])
           and np.array_equal(captured["constructor"]["inflow_state"], np.full((1, 6), 3.0)))
    U_second, _ = solver.solve(np.zeros((2, 6)))
    record("seed_override_first_outer_only",
           np.array_equal(U_second, np.ones((2, 6))))
    write_csv(HERE / "instrument_self_test.csv", ("case", "passed", "detail"), rows)
    for row in rows:
        print(f"[self-test] {row['case']}: {row['passed']}")
    failed = [r["case"] for r in rows if not r["passed"]]
    if failed:
        raise AssertionError(f"instrument self-test failed: {failed}")


def write_environment():
    packages = {}
    for name in ("numpy", "scipy", "numba", "pyamg", "meshio", "pytest"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "missing"
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO,
                                     text=True).strip()
    write_json(HERE / "environment.json", dict(
        python=sys.version, platform=platform.platform(), processor=platform.processor(),
        packages=packages, commit=commit, mesh=str(MESH.relative_to(REPO)),
        mesh_sha256=sha256_file(MESH), numba_threads=16, omp_threads=16,
        blas_arms=(1, 16), cache_policy="no cache deletion; two gate warm/check runs first",
        max_coupling_legs=MAX_COUPLING_LEGS, coupling_time_budget_s=TIME_BUDGET_S,
        gate_time_budget_s=1800.0,
    ))


def require_self_test():
    rows = read_csv(HERE / "instrument_self_test.csv")
    if not rows or any(r["passed"] != "True" for r in rows):
        raise RuntimeError("run --self-test successfully before measurements")


def gate_checks():
    require_self_test()
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for pass_no in (1, 2):
        cmd = [sys.executable, "-m", "pytest", *GATE_NODES, "-q"]
        started = time.monotonic()
        print(f"[gate-check] pass {pass_no}", flush=True)
        try:
            proc = subprocess.run(cmd, cwd=REPO, env=worker_env(1), text=True,
                                  capture_output=True, timeout=900.0)
            output = proc.stdout + "\n" + proc.stderr
            code = proc.returncode
        except subprocess.TimeoutExpired as exc:
            output = str(exc)
            code = 124
        (WORK / f"gate_check_{pass_no}.txt").write_text(output, encoding="utf-8")
        print(output, flush=True)
        failed = re.findall(r"^FAILED\s+(\S+)", output, flags=re.MULTILINE)
        rows.append(dict(pass_no=pass_no, returncode=code,
                         wall_s=f12(time.monotonic() - started),
                         failed_nodes=";".join(failed), blas_expected=1,
                         command=" ".join(cmd)))
        write_csv(HERE / "gate_checks.csv",
                  ("pass_no", "returncode", "wall_s", "failed_nodes",
                   "blas_expected", "command"), rows)


def require_gate_checks():
    rows = read_csv(HERE / "gate_checks.csv")
    if len(rows) < 2:
        raise RuntimeError("run --gate-check (two passes) before measurements")


def plot_results():
    import csv as _csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    rows = {r["leg"]: r for r in read_csv(HERE / "summary.csv") if r["status"] == "ok"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for leg, row in sorted(rows.items()):
        hist = read_csv(HERE / f"history_{leg}.csv")
        if not hist:
            continue
        k = [int(h["k"]) for h in hist]
        cl = [float(h["cl"]) for h in hist]
        ds = [max(float(h["ds_change_rel"]), 1e-16) for h in hist[1:]]
        style = "-" if row["kind"] == "native" else ("--" if row["kind"] == "tight" else ":")
        label = f"{leg} (B{row['blas_threads']})"
        axes[0].plot(k, cl, style, marker="o", ms=3, label=label)
        axes[1].semilogy(k[1:], ds, style, marker="o", ms=3, label=label)
    axes[0].set(xlabel="outer iteration", ylabel="$c_l$", title="Loose-coupling trajectories")
    axes[1].set(xlabel="outer iteration", ylabel="$\Delta\delta^*$ relative change",
                title="Stopping history")
    axes[1].axhline(1e-3, color="0.4", ls=":", lw=1)
    axes[1].axhline(1e-8, color="0.4", ls="--", lw=1)
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "trajectories.png", dpi=150)
    plt.close(fig)

    selected = [leg for leg in rows if leg.startswith("native_") or leg.startswith("tight_")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for leg in selected:
        data = read_csv(HERE / f"delta_{leg}.csv")
        for side, col, ax in ((1, "upper", axes[0]), (-1, "lower", axes[1])):
            groups = {}
            for r in data:
                if int(r["side"]) == side:
                    groups.setdefault(float(r["x_c"]), []).append(float(r["delta_star"]))
            x = np.array(sorted(groups))
            y = np.array([np.mean(groups[float(v)]) for v in x])
            ax.plot(x, y, lw=1.0, label=leg)
    for ax, title in zip(axes, ("upper", "lower")):
        ax.set(xlabel="$x/c$", ylabel="$\delta^*/c$", title=f"Final delta*: {title}")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "delta_profiles.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--gate-check", action="store_true")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--worker-replay", type=Path)
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    if args.worker:
        coupling_worker(read_json(args.worker))
        return
    if args.worker_replay:
        replay_worker(read_json(args.worker_replay))
        return
    if args.self_test:
        self_test()
        write_environment()
        return
    if args.gate_check:
        gate_checks()
        return
    if args.analyze_only:
        analyze()
        plot_results()
        return

    require_self_test()
    require_gate_checks()
    state = run_state()
    for spec in native_specs() + control_specs():
        if not run_coupling_leg(spec, state):
            analyze()
            raise RuntimeError(f"baseline leg failed: {spec['leg']}")
    if args.baseline_only:
        analyze()
        plot_results()
        return
    run_replays(state)
    metrics = analyze()
    if metrics.get("native_separated"):
        for spec in tight_specs():
            run_coupling_leg(spec, state)
        analyze()
    plot_results()


if __name__ == "__main__":
    main()
