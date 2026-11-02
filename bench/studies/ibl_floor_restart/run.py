r"""Phase 7: fixed-state IBL residual-floor attribution and restart diagnosis.

Pre-registration: docs/dev_phase_seven/20261103-0420-ibl-floor-prereg.md.
All nonlinear solves call the production IBL3Solver; no coupling leg is run here.
"""
import argparse
import contextlib
import csv
import hashlib
import importlib.metadata
import io
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
S5 = REPO / "bench" / "studies" / "xtr030_attractor"
MESH = REPO / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"
sys.path.insert(0, str(REPO))

MAX_RESTARTS = 4
TIME_BUDGET_S = 20.0 * 60.0
VAR_NAMES = ("delta", "A", "B", "Psi", "Ct1", "Ct2")
LEG_SPECS = (
    dict(leg="native_xtr030_b1_r1", role="primary", x_tr=0.30),
    dict(leg="native_xtr030_b16_r1", role="primary", x_tr=0.30),
    dict(leg="control_xtr005_b1", role="control", x_tr=0.05),
    dict(leg="control_xtr005_b16", role="control", x_tr=0.05),
)
SUMMARY_FIELDS = (
    "leg", "role", "x_tr", "status", "error", "blas_threads", "n_threads",
    "state_file", "state_sha256", "first_file", "first_sha256", "mesh_sha256",
    "U_sha256", "phi_sha256", "delta_sha256", "provenance_guard",
    "recorded_residual", "native_residual", "native_rel_error", "native_guard",
    "global_max_node", "global_max_x_c", "global_max_side", "te_energy_share",
    "near_te_energy_share", "bd_joint_share", "top_two_equations",
    "delta_min_active", "re_d_min_active", "min_delta", "min_re_d",
    "residual_signature", "restart_converged", "restart_class", "restart_n_iter",
    "restart_initial_residual", "restart_initial_rel_error", "restart_final_residual",
    "restart_cfl_final", "restart_n_rejected", "restart_tail_rel_range",
    "restart_tol_threshold", "wall_s",
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


def read_csv(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n",
                                extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)


def append_summary(row):
    rows = [r for r in read_csv(HERE / "summary.csv") if r["leg"] != row["leg"]]
    rows.append(row)
    order = {s["leg"]: i for i, s in enumerate(LEG_SPECS)}
    rows.sort(key=lambda r: order.get(r["leg"], 99))
    write_csv(HERE / "summary.csv", SUMMARY_FIELDS, rows)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha(array):
    import numpy as np
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def f12(value):
    return "" if value is None else f"{float(value):.12e}"


def configure_worker():
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ["NUMBA_NUM_THREADS"] = "16"
    os.environ["OMP_NUM_THREADS"] = "16"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    os.environ.pop("PYFP3D_ALLOW_BLAS_THREADS", None)


def require_blas():
    from pyfp3d import blas_thread_count
    actual = blas_thread_count()
    if actual != 1:
        raise RuntimeError(f"BLAS provenance failure: expected 1, got {actual}")
    return actual


def native_residual_ok(actual, recorded):
    return abs(float(actual) - float(recorded)) <= 5.0e-12 * abs(float(recorded))


def require_digest(actual, expected, label):
    if actual != expected:
        raise RuntimeError(f"{label} digest mismatch: {actual} != {expected}")


def residual_signature(F, x_c, delta, re_d):
    import numpy as np
    from pyfp3d.viscous import closures as C

    F = np.asarray(F, dtype=float)
    x_c = np.asarray(x_c, dtype=float)
    delta = np.asarray(delta, dtype=float)
    re_d = np.asarray(re_d, dtype=float)
    fmax = np.max(np.abs(F), axis=1)
    top_node = int(np.argmax(fmax))
    norm = float(np.linalg.norm(F))
    eq_norm = np.linalg.norm(F, axis=0) / max(norm, 1.0e-300)
    top_two = set(np.argsort(eq_norm)[-2:].tolist())
    bd_joint = float(np.linalg.norm(F[:, [0, 2]]) / max(norm, 1.0e-300))
    delta_active = delta < C.DELTA_MIN
    re_d_active = re_d < C.RE_D_MIN
    total_energy = float(np.sum(F ** 2))
    signature = bool(
        x_c[top_node] >= 0.90
        and top_two == {0, 2}
        and bd_joint >= 0.90
        and not np.any(delta_active)
        and not np.any(re_d_active)
    )
    return dict(
        global_max_node=top_node,
        global_max_x_c=f12(x_c[top_node]),
        te_energy_share=f12(np.sum(F[x_c >= 0.90] ** 2) / max(total_energy, 1.0e-300)),
        near_te_energy_share=f12(np.sum(F[x_c >= 0.96] ** 2) / max(total_energy, 1.0e-300)),
        bd_joint_share=f12(bd_joint),
        top_two_equations="+".join(sorted(VAR_NAMES[i] for i in top_two)),
        delta_min_active=int(np.sum(delta_active)),
        re_d_min_active=int(np.sum(re_d_active)),
        min_delta=f12(np.min(delta)),
        min_re_d=f12(np.min(re_d)),
        residual_signature=int(signature),
        _delta_active=delta_active,
        _re_d_active=re_d_active,
        _eq_norm=eq_norm,
    )


def classify_restart(converged, residual_history, final_residual, cfl_final,
                     rejected, tol=1.0e-9):
    import numpy as np
    history = np.asarray(residual_history, dtype=float)
    if history.size == 0:
        raise ValueError("empty IBL residual history")
    threshold = float(tol * max(max(history[0], 1.0e-30), 1.0))
    final_ok = float(final_residual) < threshold
    tail = history[-min(10, history.size):]
    tail_range = float((np.max(tail) - np.min(tail)) /
                       max(float(np.max(np.abs(tail))), 1.0e-300))
    if converged:
        if not final_ok:
            raise RuntimeError("production converged=True but final residual fails its criterion")
        label = "RESTART_CONVERGED"
    elif (abs(float(cfl_final) - 1.0e-3) <= 1.0e-12 and rejected
          and tail_range <= 1.0e-3):
        label = "CONTROLLER_FLOOR"
    else:
        label = "OTHER_NONCONVERGED"
    return dict(restart_class=label, restart_tail_rel_range=f12(tail_range),
                restart_tol_threshold=f12(threshold))


ITER_RE = re.compile(r"^\[ibl3\] it=(\d+) cfl=([0-9.eE+-]+) \|R\|inf=([0-9.eE+-]+)$")
REJECT_RE = re.compile(r"^\[ibl3\]\s+step rejected, cfl <- ([0-9.eE+-]+)$")


def parse_verbose(text):
    records = []
    for line in text.splitlines():
        match = ITER_RE.match(line)
        if match:
            records.append(dict(record="iter", iter=int(match.group(1)),
                                cfl=float(match.group(2)), printed_r_inf=float(match.group(3))))
            continue
        match = REJECT_RE.match(line)
        if match:
            records.append(dict(record="reject", iter=len([r for r in records if r["record"] == "iter"]) - 1,
                                cfl=float(match.group(1)), printed_r_inf=None))
    return records


def restart_rows(verbose, info):
    parsed = parse_verbose(verbose)
    history = info["residual_history"]
    iter_records = [r for r in parsed if r["record"] == "iter"]
    if len(iter_records) != info["n_iter"]:
        raise RuntimeError(f"verbose parse found {len(iter_records)} iterations, "
                           f"production reports {info['n_iter']}")
    if [r["iter"] for r in iter_records] != list(range(info["n_iter"])):
        raise RuntimeError("verbose iteration indices are not contiguous")
    rows = []
    for record in parsed:
        row = dict(record=record["record"], iter=record["iter"],
                   cfl=f12(record["cfl"]), r_inf="")
        if record["record"] == "iter":
            row["r_inf"] = f12(history[record["iter"]])
        rows.append(row)
    rows.append(dict(record="final", iter=info["n_iter"], cfl=f12(info["cfl_final"]),
                     r_inf=f12(info["final_residual"])))
    return rows


def decide(rows):
    if len(rows) != 4 or any(r.get("status") != "ok" for r in rows):
        return "FAIL"
    if any(r.get("provenance_guard") not in (1, "1", True) for r in rows):
        return "FAIL"
    primary = [r for r in rows if r["role"] == "primary"]
    control = [r for r in rows if r["role"] == "control"]
    if any(r["restart_class"] == "OTHER_NONCONVERGED" for r in rows):
        return "UNDEFINED"
    if all(r["restart_class"] == "RESTART_CONVERGED" for r in primary):
        return "PASS_RESTART_CONVERGENCE"
    floor = lambda r: (r["restart_class"] == "CONTROLLER_FLOOR"
                       and r["residual_signature"] in (1, "1", True))
    if all(floor(r) for r in primary) and any(floor(r) for r in control):
        return "PASS_KNOWN_FORMULATION_FLOOR"
    return "UNDEFINED"


def self_test():
    import tempfile
    import numpy as np

    rows = []

    def check(name, expected, function):
        try:
            actual = function()
            passed = actual == expected
        except Exception as exc:
            actual = f"{type(exc).__name__}: {exc}"
            passed = False
        rows.append(dict(case=name, expected=expected, actual=actual, passed=int(passed)))

    a = np.arange(12, dtype=float)
    check("array_digest_same", True, lambda: array_sha(a) == array_sha(a.copy()))
    check("array_digest_mismatch", True, lambda: array_sha(a) != array_sha(a + 1.0))
    check("native_residual_roundoff", True, lambda: native_residual_ok(1.0 + 2.0e-12, 1.0))
    check("native_residual_reject", False, lambda: native_residual_ok(1.0 + 2.0e-11, 1.0))

    def file_guard():
        with tempfile.TemporaryDirectory(dir=WORK) as tmp:
            path = Path(tmp) / "input.bin"
            path.write_bytes(b"ibc")
            require_digest(sha256_file(path), sha256_file(path), "self-test")
            try:
                require_digest(sha256_file(path), "0" * 64, "self-test")
            except RuntimeError:
                return True
        return False

    check("file_digest_mismatch_has_teeth", True, file_guard)

    x = np.array([0.01, 0.5, 0.98, 0.99])
    F = np.zeros((4, 6))
    F[:, 0] = [0.0, 0.1, 0.4, 0.3]
    F[:, 2] = [0.0, 0.1, 0.8, 0.4]
    F[:, 1] = 0.05
    delta = np.full(4, 1.0e-5)
    re_d = np.full(4, 1.0)
    check("residual_signature_positive", 1,
          lambda: residual_signature(F, x, delta, re_d)["residual_signature"])
    F_mid = F.copy()
    F_mid[1, 2] = 2.0
    check("residual_signature_rejects_mid_max", 0,
          lambda: residual_signature(F_mid, x, delta, re_d)["residual_signature"])
    check("residual_signature_rejects_active_floor", 0,
          lambda: residual_signature(F, x, np.zeros(4), re_d)["residual_signature"])

    sample = ("[ibl3] it=0 cfl=1.000000e+00 |R|inf=2.000000e-06\n"
              "[ibl3]   step rejected, cfl <- 1.250000e-01\n"
              "[ibl3] it=1 cfl=1.250000e-01 |R|inf=2.000000e-06\n")
    check("verbose_parser", 3, lambda: len(parse_verbose(sample)))
    check("verbose_parser_reject_cfl", 0.125,
          lambda: parse_verbose(sample)[1]["cfl"])
    check("restart_converged", "RESTART_CONVERGED",
          lambda: classify_restart(True, [1.0e-2, 1.0e-12], 1.0e-12, 4.0, False)["restart_class"])
    check("restart_controller_floor", "CONTROLLER_FLOOR",
          lambda: classify_restart(False, [2.0e-6] * 21, 2.0e-6, 1.0e-3, True)["restart_class"])
    check("restart_budget_is_other", "OTHER_NONCONVERGED",
          lambda: classify_restart(False, np.geomspace(1.0e-2, 1.0e-6, 100),
                                   1.0e-6, 1.0e8, False)["restart_class"])

    def bad_convergence():
        try:
            classify_restart(True, [1.0e-2], 1.0e-6, 1.0, False)
        except RuntimeError:
            return True
        return False

    check("restart_convergence_guard_has_teeth", True, bad_convergence)

    base = [dict(role=s["role"], status="ok", provenance_guard=1,
                 residual_signature=1, restart_class="CONTROLLER_FLOOR")
            for s in LEG_SPECS]
    check("decision_known_floor", "PASS_KNOWN_FORMULATION_FLOOR", lambda: decide(base))
    restart = [dict(r, restart_class="RESTART_CONVERGED") for r in base]
    check("decision_restart", "PASS_RESTART_CONVERGENCE", lambda: decide(restart))
    mixed = [dict(r) for r in base]
    mixed[0]["restart_class"] = "RESTART_CONVERGED"
    check("decision_mixed_primary", "UNDEFINED", lambda: decide(mixed))
    other = [dict(r) for r in base]
    other[1]["restart_class"] = "OTHER_NONCONVERGED"
    check("decision_other", "UNDEFINED", lambda: decide(other))
    no_control = [dict(r) for r in base]
    no_control[2]["residual_signature"] = 0
    no_control[2]["restart_class"] = "RESTART_CONVERGED"
    no_control[3]["residual_signature"] = 0
    no_control[3]["restart_class"] = "RESTART_CONVERGED"
    check("decision_requires_control", "UNDEFINED", lambda: decide(no_control))
    broken = [dict(r, status="failed") for r in base]
    check("decision_instrument_fail", "FAIL", lambda: decide(broken))

    write_csv(HERE / "instrument_self_test.csv",
              ("case", "expected", "actual", "passed"), rows)
    for row in rows:
        print(f"[self-test] {row['case']}: {row['actual']}", flush=True)
    failed = [r["case"] for r in rows if not r["passed"]]
    if failed:
        raise AssertionError(f"self-test failed: {failed}")
    print(f"[self-test] {len(rows)}/{len(rows)} passed", flush=True)


def input_paths(spec):
    leg = spec["leg"]
    return (S5 / "_work" / f"state_{leg}.npz",
            S5 / "_work" / f"first_{leg}.npz")


def coupling_worker(spec):
    configure_worker()
    started = time.monotonic()
    result_path = WORK / f"result_{spec['leg']}.json"
    try:
        import numpy as np
        from pyfp3d.mesh.reader import read_mesh
        from pyfp3d.mesh.wake_cut import cut_wake
        from pyfp3d.viscous import closures as C
        import pyfp3d.viscous.coupling as coupling

        blas = require_blas()
        leg = spec["leg"]
        state_path, first_path = input_paths(spec)
        for path in (state_path, first_path, MESH):
            if not path.exists():
                raise FileNotFoundError(path)
        s5_rows = {r["leg"]: r for r in read_csv(S5 / "summary.csv")}
        history = read_csv(S5 / f"history_{leg}.csv")
        if leg not in s5_rows or not history:
            raise RuntimeError(f"missing committed S5 provenance for {leg}")
        source = s5_rows[leg]
        recorded = float(history[-1]["ibl_final_residual"])
        with np.load(state_path) as data:
            state = {k: data[k] for k in data.files}
        with np.load(first_path) as data:
            first = {k: data[k] for k in data.files}

        array_hashes = {
            "U_sha256": array_sha(state["U"]),
            "phi_sha256": array_sha(state["phi"]),
            "delta_sha256": array_sha(state["delta_star"]),
        }
        for key, actual in array_hashes.items():
            require_digest(actual, source[key], f"{leg}.{key}")
        mesh_sha = sha256_file(MESH)
        require_digest(mesh_sha, source["mesh_sha256"], f"{leg}.mesh")

        mc, _ = cut_wake(read_mesh(str(MESH)))
        cfg = coupling.CouplingConfig(
            re_chord=3.0e6, m_inf=0.5, alpha_deg=2.0,
            x_tr_upper=spec["x_tr"], x_tr_lower=spec["x_tr"],
        )
        case = coupling.build_airfoil_case(
            mc.nodes, mc.elements, mc.boundary_faces["wall"], cfg)
        if not np.array_equal(first["flags"], case.turbulent_flags):
            raise RuntimeError(f"{leg}: turbulent flags differ from production case")
        if float(first["mu"]) != 1.0 / cfg.re_chord:
            raise RuntimeError(f"{leg}: viscosity provenance mismatch")

        ue = np.asarray(state["ue_surf"], dtype=float)
        U = np.asarray(state["U"], dtype=float)
        q2 = np.sum(ue ** 2, axis=1)
        q = np.sqrt(q2)
        rho = coupling.density_field(q2, cfg.m_inf, cfg.gamma_air)
        mach = np.sqrt(coupling.mach_squared_field(q2, cfg.m_inf, cfg.gamma_air))
        solver = coupling.IBL3Solver(
            case.sm, ue, rho, float(first["mu"]), mach,
            case.turbulent_flags, first["inflow_mask"], first["inflow_state"],
            eps_diff=cfg.eps_diff, eps_diff_s=cfg.eps_diff_s,
            te_pairs=None, te_extrapolate=False,
        )
        F = solver.residual(U)
        native = float(np.max(np.abs(F)))
        native_ok = native_residual_ok(native, recorded)

        stations = case.stations
        x_c = stations.xc[stations.station_of]
        side = np.asarray(stations.side_node, dtype=int)
        inflow = np.asarray(first["inflow_mask"], dtype=bool)
        turb = np.asarray(case.turbulent_flags, dtype=bool)
        regime = np.where(inflow, "inflow_pin", np.where(turb, "turbulent", "laminar"))
        region = np.where(x_c <= cfg.inflow_band_x, "le_inflow",
                          np.where(x_c <= cfg.le_band_x, "le_band",
                                   np.where(x_c >= 0.90, "te_band", "mid")))
        delta = U[:, 0]
        re_d = rho * q * np.maximum(delta, C.DELTA_MIN) / float(first["mu"])
        signature = residual_signature(F, x_c, delta, re_d)
        signature["global_max_side"] = int(side[signature["global_max_node"]])
        anatomy_fields = ("node", "x_c", "side", "regime", "region") + tuple(
            f"F_{name}" for name in VAR_NAMES) + (
            "F_max", "delta", "re_d", "delta_min_active", "re_d_min_active")
        anatomy = []
        for i in range(case.sm.n_node):
            row = dict(node=i, x_c=f12(x_c[i]), side=int(side[i]),
                       regime=regime[i], region=region[i],
                       F_max=f12(np.max(np.abs(F[i]))), delta=f12(delta[i]), re_d=f12(re_d[i]),
                       delta_min_active=int(signature["_delta_active"][i]),
                       re_d_min_active=int(signature["_re_d_active"][i]))
            row.update({f"F_{name}": f12(F[i, j]) for j, name in enumerate(VAR_NAMES)})
            anatomy.append(row)
        write_csv(HERE / f"residual_anatomy_{leg}.csv", anatomy_fields, anatomy)

        base = dict(
            leg=leg, role=spec["role"], x_tr=f12(spec["x_tr"]), status="ok",
            blas_threads=blas, n_threads=16,
            state_file=str(state_path.relative_to(REPO)), state_sha256=sha256_file(state_path),
            first_file=str(first_path.relative_to(REPO)), first_sha256=sha256_file(first_path),
            mesh_sha256=mesh_sha, **array_hashes, provenance_guard=1,
            recorded_residual=f12(recorded), native_residual=f12(native),
            native_rel_error=f12(abs(native - recorded) / max(recorded, 1.0e-300)),
            native_guard=int(native_ok),
        )
        base.update({k: v for k, v in signature.items() if not k.startswith("_")})
        if not native_ok:
            raise RuntimeError(f"{leg}: rebuilt residual {native:.16e} does not match "
                               f"committed native residual {recorded:.16e}")

        capture = io.StringIO()
        with contextlib.redirect_stdout(capture):
            U_restart, info = solver.solve(
                U, tol=cfg.ibl_tol, max_iter=cfg.ibl_max_iter, verbose=True)
        initial = float(info["residual_history"][0])
        if not native_residual_ok(initial, recorded):
            raise RuntimeError(f"{leg}: restart initial residual does not match native")
        rows = restart_rows(capture.getvalue(), info)
        write_csv(HERE / f"restart_history_{leg}.csv",
                  ("record", "iter", "cfl", "r_inf"), rows)
        classification = classify_restart(
            info["converged"], info["residual_history"], info["final_residual"],
            info["cfl_final"], any(r["record"] == "reject" for r in rows), tol=cfg.ibl_tol)
        np.savez_compressed(WORK / f"restart_{leg}.npz", U=U_restart,
                            residual=solver.residual(U_restart))
        result = dict(
            **base, restart_converged=int(info["converged"]),
            restart_n_iter=info["n_iter"],
            restart_initial_residual=f12(initial),
            restart_initial_rel_error=f12(abs(initial - recorded) / max(recorded, 1.0e-300)),
            restart_final_residual=f12(info["final_residual"]),
            restart_cfl_final=f12(info["cfl_final"]),
            restart_n_rejected=sum(r["record"] == "reject" for r in rows),
            wall_s=f12(time.monotonic() - started), **classification,
        )
        write_json(result_path, result)
        print("RESULT " + json.dumps(result), flush=True)
    except Exception as exc:
        write_json(result_path, dict(
            leg=spec["leg"], role=spec["role"], status="failed",
            error=f"{type(exc).__name__}: {exc}", traceback=traceback.format_exc(),
            wall_s=f12(time.monotonic() - started)))
        raise


def worker_env():
    env = os.environ.copy()
    env.update(PYTHONNOUSERSITE="1", NUMBA_NUM_THREADS="16", OMP_NUM_THREADS="16",
               OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg")
    env.pop("PYFP3D_ALLOW_BLAS_THREADS", None)
    return env


def run_worker(spec):
    WORK.mkdir(parents=True, exist_ok=True)
    spec_path = WORK / f"spec_{spec['leg']}.json"
    write_json(spec_path, spec)
    output_path = WORK / f"worker_{spec['leg']}.txt"
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", str(spec_path)]
    print(f"[run] {spec['leg']}", flush=True)
    started = time.monotonic()
    try:
        process = subprocess.run(command, cwd=REPO, env=worker_env(), text=True,
                                 capture_output=True, timeout=TIME_BUDGET_S)
        output_path.write_text(process.stdout + "\n" + process.stderr, encoding="utf-8")
        print(process.stdout, flush=True)
        if process.returncode:
            print(process.stderr, file=sys.stderr, flush=True)
    except subprocess.TimeoutExpired as exc:
        output_path.write_text(str(exc), encoding="utf-8")
        write_json(WORK / f"result_{spec['leg']}.json",
                   dict(leg=spec["leg"], role=spec["role"], status="failed",
                        error="worker timeout"))
    return time.monotonic() - started


def run_state():
    path = HERE / "run_state.json"
    if path.exists():
        return read_json(path)
    return dict(fresh_restarts=0, cumulative_restart_s=0.0,
                max_restarts=MAX_RESTARTS, time_budget_s=TIME_BUDGET_S, completed=[])


def analyze():
    rows = read_csv(HERE / "summary.csv")
    verdict = decide(rows)
    by_leg = {r["leg"]: r for r in rows}
    metrics = {"verdict": verdict, "n_complete": len(rows)}
    for spec in LEG_SPECS:
        row = by_leg.get(spec["leg"], {})
        for key in ("native_residual", "residual_signature", "restart_class",
                    "restart_final_residual", "restart_cfl_final", "wall_s"):
            metrics[f"{spec['leg']}__{key}"] = row.get(key, "")
    write_csv(HERE / "verdict_metrics.csv", ("metric", "value"),
              [dict(metric=k, value=v) for k, v in metrics.items()])

    self_tests = read_csv(HERE / "instrument_self_test.csv")
    instruments_ok = bool(self_tests) and all(r["passed"] == "1" for r in self_tests)
    provenance_ok = len(rows) == 4 and all(r.get("provenance_guard") == "1" for r in rows)
    native_ok = len(rows) == 4 and all(r.get("native_guard") == "1" for r in rows)
    primary = [r for r in rows if r.get("role") == "primary" and r.get("status") == "ok"]
    controls = [r for r in rows if r.get("role") == "control" and r.get("status") == "ok"]
    criteria = [
        dict(criterion="I0", verdict="PASS" if instruments_ok and provenance_ok and native_ok else "FAIL",
             evidence="instrument_self_test.csv; provenance.csv; summary.csv",
             reading=f"self-test={sum(r['passed'] == '1' for r in self_tests)}/{len(self_tests)}, "
                     f"provenance={provenance_ok}, native_residual={native_ok}"),
        dict(criterion="C1", verdict="RECORDED" if len(rows) == 4 else "UNDEFINED",
             evidence="residual_anatomy_*.csv; residual_map.png; summary.csv",
             reading=f"old-floor residual signature on {sum(r.get('residual_signature') == '1' for r in rows)}/4 states"),
        dict(criterion="C2", verdict="RECORDED" if len(rows) == 4 else "UNDEFINED",
             evidence="restart_history_*.csv; restart_trajectories.png; summary.csv",
             reading="; ".join(f"{r['leg']}={r.get('restart_class', '')}" for r in rows)),
        dict(criterion="C3", verdict=verdict.split("_")[0] if verdict.startswith(("PASS", "FAIL")) else "UNDEFINED",
             evidence="criteria.csv; verdict_metrics.csv",
             reading=f"primary={len(primary)}, controls={len(controls)}, decision={verdict}"),
        dict(criterion="OVERALL", verdict=verdict.split("_")[0] if verdict.startswith(("PASS", "FAIL")) else "UNDEFINED",
             evidence="verdict_metrics.csv", reading=verdict),
    ]
    write_csv(HERE / "criteria.csv", ("criterion", "verdict", "evidence", "reading"), criteria)

    provenance_fields = ("leg", "state_file", "state_sha256", "first_file", "first_sha256",
                         "mesh_sha256", "U_sha256", "phi_sha256", "delta_sha256",
                         "blas_threads", "provenance_guard", "recorded_residual",
                         "native_residual", "native_rel_error", "native_guard")
    write_csv(HERE / "provenance.csv", provenance_fields, rows)
    write_csv(HERE / "residual_summary.csv",
              ("leg", "role", "native_residual", "global_max_node", "global_max_x_c",
               "global_max_side", "te_energy_share", "near_te_energy_share", "bd_joint_share",
               "top_two_equations", "delta_min_active", "re_d_min_active", "min_delta",
               "min_re_d", "residual_signature"), rows)
    write_csv(HERE / "restart_summary.csv",
              ("leg", "role", "restart_converged", "restart_class", "restart_n_iter",
               "restart_initial_residual", "restart_final_residual", "restart_cfl_final",
               "restart_n_rejected", "restart_tail_rel_range", "restart_tol_threshold", "wall_s"), rows)
    return metrics


def plot_results():
    import matplotlib.pyplot as plt
    import numpy as np

    rows = [r for r in read_csv(HERE / "summary.csv") if r.get("status") == "ok"]
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for row in rows:
        leg = row["leg"]
        anatomy = read_csv(HERE / f"residual_anatomy_{leg}.csv")
        x = np.array([float(r["x_c"]) for r in anatomy])
        side = np.array([int(r["side"]) for r in anatomy])
        fmax = np.array([float(r["F_max"]) for r in anatomy])
        for axis, value, title in zip(axes, (1, -1), ("upper", "lower")):
            mask = side == value
            axis.semilogy(x[mask], np.maximum(fmax[mask], 1.0e-16), ".", ms=3, label=leg)
            axis.axvline(0.90, color="gray", ls=":", lw=0.8)
            axis.set(xlabel="x/c", title=f"IBL residual anatomy: {title}")
            axis.grid(alpha=0.25)
    axes[0].set_ylabel(r"$\max_j |F_j|$")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "residual_map.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for row in rows:
        leg = row["leg"]
        history = read_csv(HERE / f"restart_history_{leg}.csv")
        iterations = [r for r in history if r["record"] == "iter"]
        k = np.array([int(r["iter"]) for r in iterations])
        residual = np.array([float(r["r_inf"]) for r in iterations])
        cfl = np.array([float(r["cfl"]) for r in iterations])
        axes[0].semilogy(k, residual, marker="o", ms=2, label=leg)
        axes[1].semilogy(k, cfl, marker="o", ms=2, label=leg)
    axes[0].set(xlabel="restart iteration", ylabel=r"$\|R\|_\infty$", title="Pure steady IBL residual")
    axes[1].set(xlabel="restart iteration", ylabel="CFL", title="Pseudo-time controller")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "restart_trajectories.png", dpi=150)
    plt.close(fig)


def write_environment():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    packages = {name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "numba", "meshio", "pyamg", "pytest")}
    write_json(HERE / "environment.json", dict(
        commit=commit, platform=platform.platform(), processor=platform.processor(),
        python=sys.version, packages=packages, blas_threads=1, numba_threads=16,
        omp_threads=16, cache_policy="no cache deletion; warm cache; fresh worker per state",
        mesh=str(MESH.relative_to(REPO)), mesh_sha256=sha256_file(MESH),
        max_restarts=MAX_RESTARTS, fresh_coupling_legs=0, dense_svd=0,
        time_budget_s=TIME_BUDGET_S,
    ))


def require_self_test():
    rows = read_csv(HERE / "instrument_self_test.csv")
    if not rows or any(r["passed"] != "1" for r in rows):
        raise RuntimeError("run --self-test successfully before formal measurements")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    if args.worker:
        coupling_worker(read_json(args.worker))
        return
    if args.self_test:
        configure_worker()
        self_test()
        write_environment()
        return
    if args.analyze_only:
        analyze()
        plot_results()
        return

    require_self_test()
    state = run_state()
    for spec in LEG_SPECS:
        rows = {r["leg"]: r for r in read_csv(HERE / "summary.csv")}
        if rows.get(spec["leg"], {}).get("status") == "ok":
            print(f"[skip] {spec['leg']} already complete", flush=True)
            continue
        if state["fresh_restarts"] >= MAX_RESTARTS:
            raise RuntimeError("restart budget exhausted")
        if state["cumulative_restart_s"] >= TIME_BUDGET_S:
            raise RuntimeError("diagnostic wall-time budget exhausted")
        elapsed = run_worker(spec)
        state["fresh_restarts"] += 1
        state["cumulative_restart_s"] += elapsed
        state["completed"].append(spec["leg"])
        write_json(HERE / "run_state.json", state)
        result = read_json(WORK / f"result_{spec['leg']}.json")
        if result["status"] == "ok":
            append_summary({k: result.get(k, "") for k in SUMMARY_FIELDS})
        else:
            row = {k: "" for k in SUMMARY_FIELDS}
            row.update(leg=spec["leg"], role=spec["role"], x_tr=f12(spec["x_tr"]),
                       status="failed", error=result.get("error", ""), wall_s=f12(elapsed))
            append_summary(row)
        analyze()
        if result["status"] != "ok":
            plot_results()
            raise RuntimeError(f"instrument failed for {spec['leg']}: {result.get('error')}")
    analyze()
    plot_results()


if __name__ == "__main__":
    main()
