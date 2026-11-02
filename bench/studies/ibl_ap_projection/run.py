r"""Phase 7: A/Psi residual Jacobian projection and accepted-step anatomy.

Pre-registration: docs/dev_phase_seven/20261103-0454-ibl-ap-projection-prereg.md.
Four fixed states only; no coupling solve and no production-library changes.
"""
import argparse
import csv
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import types

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WORK = HERE / "_work"
PREVIOUS = REPO / "bench" / "studies" / "ibl_floor_restart"
S5 = REPO / "bench" / "studies" / "xtr030_attractor"
MESH = REPO / "cases" / "meshes" / "naca0012_2.5d" / "medium.msh"
sys.path.insert(0, str(REPO))


def _load_previous():
    spec = importlib.util.spec_from_file_location(
        "ibl_floor_restart_run", PREVIOUS / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


floor = _load_previous()
read_csv = floor.read_csv
write_csv = floor.write_csv
read_json = floor.read_json
write_json = floor.write_json
sha256_file = floor.sha256_file
array_sha = floor.array_sha
f12 = floor.f12
configure_worker = floor.configure_worker
require_blas = floor.require_blas
native_residual_ok = floor.native_residual_ok
require_digest = floor.require_digest

LEG_SPECS = floor.LEG_SPECS
VAR_NAMES = floor.VAR_NAMES
AP_VARS = (1, 3)
BD_VARS = (2, 0)
N_NULL = 20
SOFT_REL = 1.0e-3
DOMINANT = 0.50
MAX_LEGS = 4
WORKER_BUDGET_S = 10.0 * 60.0
TIME_BUDGET_S = 30.0 * 60.0

SUMMARY_FIELDS = (
    "leg", "role", "x_tr", "status", "error", "blas_threads", "n_threads",
    "state_file", "state_sha256", "first_file", "first_sha256", "mesh_sha256",
    "U_sha256", "phi_sha256", "delta_sha256", "provenance_guard",
    "recorded_residual", "native_residual", "native_rel_error", "native_guard",
    "n_dof", "sigma_max", "sigma_min", "cond_equilibrated", "sigma20_rel",
    "n_lt_1e-6", "n_lt_1e-8", "n_lt_1e-10", "row_norm_range",
    "col_norm_range", "right_AP_share20", "right_BD_share20",
    "n_right_AP_dominant20", "capture20_full", "capture20_AP", "capture20_BD",
    "capture20_AP_dominant_modes", "dimension_baseline20", "AP_capture_over_baseline",
    "svd_reconstruction_rel_error", "projection_partition_error",
    "right_partition_error", "projection_guard", "projection_class",
    "restart_converged", "restart_n_iter", "restart_n_rejected_iterations",
    "restart_n_accepted", "restart_final_residual", "restart_cfl_final",
    "replay_u_equal", "replay_residual_equal", "replay_history_max_rel",
    "replay_cfl_max_rel", "replay_guard", "ap_steps_reduced", "ap_steps_amplified",
    "ap_steps_flat", "ap_accept_energy_before", "ap_accept_energy_after",
    "ap_accept_energy_rel_change", "accepted_step_class", "wall_svd_s",
    "wall_restart_s", "wall_s",
)


def append_summary(row):
    rows = [r for r in read_csv(HERE / "summary.csv") if r["leg"] != row["leg"]]
    rows.append(row)
    order = {s["leg"]: i for i, s in enumerate(LEG_SPECS)}
    rows.sort(key=lambda r: order.get(r["leg"], 99))
    write_csv(HERE / "summary.csv", SUMMARY_FIELDS, rows)


def block_mask(n_node, variables):
    import numpy as np
    mask = np.zeros((n_node, 6), dtype=float)
    mask[:, variables] = 1.0
    return mask.ravel()


def capture_energy(coefficients, bottom, target):
    import numpy as np
    norm2 = float(target @ target)
    if norm2 == 0.0:
        raise ValueError("projection target has zero energy")
    return float(np.sum(coefficients[bottom] ** 2) / norm2)


def dominant_ap_projection(right_ap_share20, capture20_ap, sigma20_rel):
    return int(right_ap_share20 >= DOMINANT and capture20_ap >= DOMINANT
               and sigma20_rel <= SOFT_REL)


def projection_decomposition(J, F, k=N_NULL):
    import numpy as np
    J = np.asarray(J, dtype=float)
    F = np.asarray(F, dtype=float).ravel()
    if J.ndim != 2 or J.shape[0] != J.shape[1] or J.shape[0] != F.size:
        raise ValueError("J must be square and match flattened F")
    n = J.shape[0]
    if n % 6 or not 0 < k < n:
        raise ValueError("invalid six-variable system or projection rank")
    row_norm = np.linalg.norm(J, axis=1)
    row_scale = 1.0 / np.maximum(row_norm, 1.0e-300)
    J1 = J * row_scale[:, None]
    col_norm = np.linalg.norm(J1, axis=0)
    col_scale = 1.0 / np.maximum(col_norm, 1.0e-300)
    Js = J1 * col_scale[None, :]
    left, sigma, vh = np.linalg.svd(Js, full_matrices=False)
    if not (np.all(np.isfinite(left)) and np.all(np.isfinite(sigma))
            and np.all(np.isfinite(vh))):
        raise RuntimeError("non-finite SVD output")
    f = F * row_scale
    n_node = n // 6
    targets = {
        "full": f,
        "AP": f * block_mask(n_node, AP_VARS),
        "BD": f * block_mask(n_node, BD_VARS),
    }
    coefficients = {name: left.T @ target for name, target in targets.items()}
    bottom = np.arange(n - k, n)
    captures = {name: capture_energy(coeff, bottom, targets[name])
                for name, coeff in coefficients.items()}
    complement = {
        name: capture_energy(coeff, np.arange(0, n - k), targets[name])
        for name, coeff in coefficients.items()
    }
    partition_error = max(abs(captures[name] + complement[name] - 1.0)
                          for name in targets)
    right_by_var = np.sum(vh.reshape(n, n_node, 6) ** 2, axis=1)
    right_partition_error = float(np.max(np.abs(np.sum(right_by_var, axis=1) - 1.0)))
    right_ap = np.sum(right_by_var[:, AP_VARS], axis=1)
    right_bd = np.sum(right_by_var[:, BD_VARS], axis=1)
    ap_dominant = bottom[right_ap[bottom] >= DOMINANT]
    reconstruction = (left * sigma[None, :]) @ vh
    reconstruction_error = float(
        np.linalg.norm(reconstruction - Js) / max(np.linalg.norm(Js), 1.0e-300))
    guard = bool(reconstruction_error <= 1.0e-10 and partition_error <= 1.0e-10
                 and right_partition_error <= 1.0e-10)
    metrics = dict(
        n_dof=n,
        sigma_max=f12(sigma[0]), sigma_min=f12(sigma[-1]),
        cond_equilibrated=f12(sigma[0] / max(sigma[-1], 1.0e-300)),
        sigma20_rel=f12(sigma[n - k] / sigma[0]),
        **{"n_lt_1e-6": int(np.sum(sigma < 1.0e-6 * sigma[0])),
           "n_lt_1e-8": int(np.sum(sigma < 1.0e-8 * sigma[0])),
           "n_lt_1e-10": int(np.sum(sigma < 1.0e-10 * sigma[0]))},
        row_norm_range=f12(np.max(row_norm) / max(np.min(row_norm), 1.0e-300)),
        col_norm_range=f12(np.max(col_norm) / max(np.min(col_norm), 1.0e-300)),
        right_AP_share20=f12(np.sum(right_ap[bottom]) / k),
        right_BD_share20=f12(np.sum(right_bd[bottom]) / k),
        n_right_AP_dominant20=int(ap_dominant.size),
        capture20_full=f12(captures["full"]), capture20_AP=f12(captures["AP"]),
        capture20_BD=f12(captures["BD"]),
        capture20_AP_dominant_modes=f12(
            np.sum(coefficients["AP"][ap_dominant] ** 2)
            / float(targets["AP"] @ targets["AP"])),
        dimension_baseline20=f12(k / n),
        AP_capture_over_baseline=f12(captures["AP"] / (k / n)),
        svd_reconstruction_rel_error=f12(reconstruction_error),
        projection_partition_error=f12(partition_error),
        right_partition_error=f12(right_partition_error),
        projection_guard=int(guard),
        projection_class=dominant_ap_projection(
            float(np.sum(right_ap[bottom]) / k), captures["AP"], sigma[n - k] / sigma[0]),
    )
    return dict(metrics=metrics, left=left, sigma=sigma, vh=vh,
                right_by_var=right_by_var, coefficients=coefficients,
                targets=targets, bottom=bottom, row_scale=row_scale,
                col_scale=col_scale, Js=Js, reconstruction_error=reconstruction_error,
                partition_error=partition_error,
                right_partition_error=right_partition_error, guard=guard)


def write_projection_modes(leg, decomposition):
    sigma = decomposition["sigma"]
    right = decomposition["right_by_var"]
    coefficients = decomposition["coefficients"]
    fields = (("mode", "sigma", "sigma_rel")
              + tuple(f"right_{name}" for name in VAR_NAMES)
              + ("right_AP", "right_BD", "proj_full", "proj_AP", "proj_BD"))
    rows = []
    for i in range(sigma.size):
        row = dict(mode=i, sigma=f12(sigma[i]), sigma_rel=f12(sigma[i] / sigma[0]),
                   right_AP=f12(right[i, AP_VARS].sum()),
                   right_BD=f12(right[i, BD_VARS].sum()),
                   proj_full=f12(coefficients["full"][i]),
                   proj_AP=f12(coefficients["AP"][i]),
                   proj_BD=f12(coefficients["BD"][i]))
        row.update({f"right_{name}": f12(right[i, j])
                    for j, name in enumerate(VAR_NAMES)})
        rows.append(row)
    write_csv(HERE / f"projection_modes_{leg}.csv", fields, rows)


def build_native(spec):
    import numpy as np
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    import pyfp3d.viscous.coupling as coupling

    blas = require_blas()
    leg = spec["leg"]
    state_path, first_path = floor.input_paths(spec)
    restart_path = PREVIOUS / "_work" / f"restart_{leg}.npz"
    for path in (state_path, first_path, restart_path, MESH):
        if not path.exists():
            raise FileNotFoundError(path)
    previous_rows = {r["leg"]: r for r in read_csv(PREVIOUS / "summary.csv")}
    previous = previous_rows.get(leg)
    if previous is None or previous.get("status") != "ok":
        raise RuntimeError(f"missing valid round-7 provenance for {leg}")
    if previous.get("provenance_guard") != "1" or previous.get("native_guard") != "1":
        raise RuntimeError(f"round-7 guards are not green for {leg}")
    s5_rows = {r["leg"]: r for r in read_csv(S5 / "summary.csv")}
    history = read_csv(S5 / f"history_{leg}.csv")
    if leg not in s5_rows or not history:
        raise RuntimeError(f"missing committed S5 provenance for {leg}")
    recorded = float(history[-1]["ibl_final_residual"])
    with np.load(state_path) as data:
        state = {k: data[k] for k in data.files}
    with np.load(first_path) as data:
        first = {k: data[k] for k in data.files}
    hashes = {
        "state_sha256": sha256_file(state_path),
        "first_sha256": sha256_file(first_path),
        "mesh_sha256": sha256_file(MESH),
        "U_sha256": array_sha(state["U"]),
        "phi_sha256": array_sha(state["phi"]),
        "delta_sha256": array_sha(state["delta_star"]),
    }
    for key, actual in hashes.items():
        require_digest(actual, previous[key], f"{leg}.{key}")
        if key in ("U_sha256", "phi_sha256", "delta_sha256", "mesh_sha256"):
            require_digest(actual, s5_rows[leg][key], f"{leg}.S5.{key}")

    mc, _ = cut_wake(read_mesh(str(MESH)))
    cfg = coupling.CouplingConfig(
        re_chord=3.0e6, m_inf=0.5, alpha_deg=2.0,
        x_tr_upper=spec["x_tr"], x_tr_lower=spec["x_tr"])
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
        te_pairs=None, te_extrapolate=False)
    F = solver.residual(U)
    native = float(np.max(np.abs(F)))
    native_ok = native_residual_ok(native, recorded)
    base = dict(
        leg=leg, role=spec["role"], x_tr=f12(spec["x_tr"]), status="ok",
        blas_threads=blas, n_threads=16,
        state_file=str(state_path.relative_to(REPO)),
        first_file=str(first_path.relative_to(REPO)), **hashes,
        provenance_guard=1, recorded_residual=f12(recorded),
        native_residual=f12(native),
        native_rel_error=f12(abs(native - recorded) / max(recorded, 1.0e-300)),
        native_guard=int(native_ok))
    if not native_ok:
        raise RuntimeError(f"{leg}: native residual mismatch {native} != {recorded}")
    return dict(solver=solver, U=U, F=F, cfg=cfg, base=base,
                restart_path=restart_path, previous_history=read_csv(
                    PREVIOUS / f"restart_history_{leg}.csv"))


def _cfl_from_dt_inv(solver, dt_inv):
    import numpy as np
    if dt_inv is None:
        return float("nan")
    q = np.asarray(solver._q)
    h = np.sqrt(np.maximum(solver.smesh.node_area, 1.0e-30))
    valid = q > 0.0
    ratio = np.asarray(dt_inv)[valid] * h[valid] / q[valid]
    if not np.all(np.isfinite(ratio)) or np.any(ratio <= 0.0):
        raise RuntimeError("cannot recover CFL from pseudo-time weight")
    cfl = float(1.0 / np.median(ratio))
    if float(np.max(np.abs(ratio * cfl - 1.0))) > 1.0e-10:
        raise RuntimeError("pseudo-time weight is not scalar CFL times q/h")
    return cfl


def solve_with_trace(solver, U0, **solve_kwargs):
    import numpy as np
    from pyfp3d.viscous.ibl3 import IBL3Solver

    original_rj = solver.residual_jacobian
    original_residual = solver.residual
    contexts = []
    events = []
    current = None
    sentinel = object()
    saved = {name: solver.__dict__.get(name, sentinel)
             for name in ("residual", "residual_jacobian")}

    def wrapped_rj(self, U, dt_inv=None, U_old=None):
        nonlocal current
        R, J = original_rj(U, dt_inv=dt_inv, U_old=U_old)
        current = dict(iter=len(contexts), U=np.array(U, copy=True),
                       U_old=None if U_old is None else np.array(U_old, copy=True),
                       dt_inv=None if dt_inv is None else np.array(dt_inv, copy=True),
                       cfl=_cfl_from_dt_inv(self, dt_inv), R=np.array(R, copy=True),
                       fpnorm=float(np.max(np.abs(R))), F_before=None, trials=[])
        contexts.append(current)
        return R, J

    def wrapped_residual(self, U, dt_inv=None, U_old=None):
        result = original_residual(U, dt_inv=dt_inv, U_old=U_old)
        if dt_inv is None:
            if (current is not None and current["F_before"] is None
                    and np.array_equal(U, current["U"])):
                current["F_before"] = np.array(result, copy=True)
            return result
        if current is None or current["F_before"] is None:
            raise RuntimeError("line-search trace has no base state")
        trial_index = len(current["trials"])
        if trial_index >= 12:
            raise RuntimeError("production line search exceeded 12 recorded trials")
        omega = 0.5 ** trial_index
        fnew = float(np.max(np.abs(result)))
        accepted = bool(fnew <= (1.0 - 1.0e-4 * omega) * current["fpnorm"])
        pure = np.array(original_residual(U), copy=True)
        trial = dict(iter=current["iter"], cfl=current["cfl"], omega=omega,
                     U=np.array(U, copy=True), R=np.array(result, copy=True),
                     F=pure, fnew=fnew, accepted=accepted)
        current["trials"].append(trial)
        if accepted:
            events.append(dict(
                iter=current["iter"], cfl=current["cfl"], omega=omega,
                U_before=current["U"].copy(), U_after=trial["U"].copy(),
                R_before=current["R"].copy(), R_after=trial["R"].copy(),
                F_before=current["F_before"].copy(), F_after=trial["F"].copy(),
                fpnorm_before=current["fpnorm"], fpnorm_after=fnew))
        return result

    solver.residual_jacobian = types.MethodType(wrapped_rj, solver)
    solver.residual = types.MethodType(wrapped_residual, solver)
    try:
        U, info = IBL3Solver.solve(solver, U0, **solve_kwargs)
    finally:
        for name, old in saved.items():
            if old is sentinel:
                delattr(solver, name)
            else:
                setattr(solver, name, old)
    return U, info, contexts, events


def classify_ap_step(before_l2, after_l2):
    before = float(before_l2)
    after = float(after_l2)
    if after < before * (1.0 - 1.0e-12):
        return "AP_REDUCED"
    if after > before * (1.0 + 1.0e-12):
        return "AP_AMPLIFIED"
    return "AP_FLAT"


def accepted_class(events):
    import numpy as np
    if not events:
        return dict(accepted_step_class="NO_ACCEPTED_STEPS",
                    ap_accept_energy_before=f12(0.0), ap_accept_energy_after=f12(0.0),
                    ap_accept_energy_rel_change="")
    before = float(sum(np.linalg.norm(e["F_before"][:, AP_VARS]) ** 2 for e in events))
    after = float(sum(np.linalg.norm(e["F_after"][:, AP_VARS]) ** 2 for e in events))
    label = classify_ap_step(np.sqrt(before), np.sqrt(after))
    return dict(accepted_step_class=label,
                ap_accept_energy_before=f12(before), ap_accept_energy_after=f12(after),
                ap_accept_energy_rel_change=f12((after - before) / max(before, 1.0e-300)))


def event_metrics(event):
    import numpy as np
    F0, F1 = event["F_before"], event["F_after"]
    R0, R1 = event["R_before"], event["R_after"]
    ap0 = float(np.linalg.norm(F0[:, AP_VARS]))
    ap1 = float(np.linalg.norm(F1[:, AP_VARS]))
    return dict(
        iter=event["iter"], cfl=f12(event["cfl"]), omega=f12(event["omega"]),
        pure_inf_before=f12(np.max(np.abs(F0))), pure_inf_after=f12(np.max(np.abs(F1))),
        pure_l2_before=f12(np.linalg.norm(F0)), pure_l2_after=f12(np.linalg.norm(F1)),
        pt_inf_before=f12(np.max(np.abs(R0))), pt_inf_after=f12(np.max(np.abs(R1))),
        ap_l2_before=f12(ap0), ap_l2_after=f12(ap1),
        ap_l2_ratio=f12(ap1 / max(ap0, 1.0e-300)),
        bd_l2_before=f12(np.linalg.norm(F0[:, BD_VARS])),
        bd_l2_after=f12(np.linalg.norm(F1[:, BD_VARS])),
        ap_step_class=classify_ap_step(ap0, ap1),
        dF_AP_l2=f12(np.linalg.norm((F1 - F0)[:, AP_VARS])),
        dU_AP_l2=f12(np.linalg.norm((event["U_after"] - event["U_before"])[:, AP_VARS])),
    )


def write_trace(leg, role, contexts, events):
    import numpy as np
    event_by_iter = {e["iter"]: e for e in events}
    iteration_rows = []
    for ctx in contexts:
        accepted = ctx["iter"] in event_by_iter
        event = event_by_iter.get(ctx["iter"])
        iteration_rows.append(dict(
            leg=leg, role=role, iter=ctx["iter"], cfl=f12(ctx["cfl"]),
            pure_inf="" if ctx["F_before"] is None else f12(np.max(np.abs(ctx["F_before"]))),
            pt_inf=f12(ctx["fpnorm"]), n_trials=len(ctx["trials"]),
            accepted=int(accepted), omega_accept="" if event is None else f12(event["omega"]),
            rejected_iteration=int(bool(ctx["trials"]) and not accepted)))
    write_csv(HERE / f"iteration_trace_{leg}.csv",
              ("leg", "role", "iter", "cfl", "pure_inf", "pt_inf", "n_trials",
               "accepted", "omega_accept", "rejected_iteration"), iteration_rows)

    step_fields = ("leg", "role", "iter", "cfl", "omega", "node", "equation",
                   "F_before", "F_after", "dF", "Rpt_before", "Rpt_after", "dU")
    step_rows = []
    summary_rows = []
    for event in events:
        for i in range(event["F_before"].shape[0]):
            for j, name in enumerate(VAR_NAMES):
                step_rows.append(dict(
                    leg=leg, role=role, iter=event["iter"], cfl=f12(event["cfl"]),
                    omega=f12(event["omega"]), node=i, equation=name,
                    F_before=f12(event["F_before"][i, j]),
                    F_after=f12(event["F_after"][i, j]),
                    dF=f12(event["F_after"][i, j] - event["F_before"][i, j]),
                    Rpt_before=f12(event["R_before"][i, j]),
                    Rpt_after=f12(event["R_after"][i, j]),
                    dU=f12(event["U_after"][i, j] - event["U_before"][i, j])))
        summary_rows.append(dict(leg=leg, role=role, **event_metrics(event)))
    write_csv(HERE / f"accepted_steps_{leg}.csv", step_fields, step_rows)
    fields = ("leg", "role", "iter", "cfl", "omega", "pure_inf_before",
              "pure_inf_after", "pure_l2_before", "pure_l2_after", "pt_inf_before",
              "pt_inf_after", "ap_l2_before", "ap_l2_after", "ap_l2_ratio",
              "bd_l2_before", "bd_l2_after", "ap_step_class", "dF_AP_l2", "dU_AP_l2")
    write_csv(HERE / f"accepted_event_metrics_{leg}.csv", fields, summary_rows)
    return iteration_rows, step_rows, summary_rows


def replay_metrics(U, final_F, info, contexts, restart_path, previous_history):
    import numpy as np
    with np.load(restart_path) as data:
        previous_U = data["U"]
        previous_F = data["residual"]
    iterations = [r for r in previous_history if r["record"] == "iter"]
    if len(iterations) != info["n_iter"] or len(contexts) != info["n_iter"]:
        raise RuntimeError("restart iteration count differs from round-7 evidence")
    recorded_r = np.array([float(r["r_inf"]) for r in iterations])
    actual_r = np.asarray(info["residual_history"], dtype=float)
    history_rel = np.abs(actual_r - recorded_r) / np.maximum(np.abs(recorded_r), 1.0)
    recorded_cfl = np.array([float(r["cfl"]) for r in iterations])
    actual_cfl = np.array([ctx["cfl"] for ctx in contexts])
    cfl_rel = np.abs(actual_cfl - recorded_cfl) / np.maximum(np.abs(recorded_cfl), 1.0e-300)
    u_equal = bool(np.array_equal(U, previous_U))
    residual_equal = bool(np.array_equal(final_F, previous_F))
    history_max = float(np.max(history_rel)) if history_rel.size else 0.0
    cfl_max = float(np.max(cfl_rel)) if cfl_rel.size else 0.0
    guard = bool(u_equal and residual_equal and history_max <= 5.0e-12
                 and cfl_max <= 5.0e-4)
    return dict(replay_u_equal=int(u_equal), replay_residual_equal=int(residual_equal),
                replay_history_max_rel=f12(history_max), replay_cfl_max_rel=f12(cfl_max),
                replay_guard=int(guard))


def coupling_worker(spec):
    configure_worker()
    started = time.monotonic()
    result_path = WORK / f"result_{spec['leg']}.json"
    try:
        native = build_native(spec)
        solver, U, F, cfg = native["solver"], native["U"], native["F"], native["cfg"]
        svd_started = time.monotonic()
        _, J = solver.residual_jacobian(U)
        decomposition = projection_decomposition(J.toarray(), F)
        write_projection_modes(spec["leg"], decomposition)
        if not decomposition["guard"]:
            raise RuntimeError(f"{spec['leg']}: SVD/projection numerical guard failed")
        projection_metrics = decomposition["metrics"]
        wall_svd = time.monotonic() - svd_started
        write_json(WORK / f"projection_{spec['leg']}.json", projection_metrics)
        del decomposition, J

        restart_started = time.monotonic()
        U_final, info, contexts, events = solve_with_trace(
            solver, U, tol=cfg.ibl_tol, max_iter=cfg.ibl_max_iter)
        write_trace(spec["leg"], spec["role"], contexts, events)
        final_F = solver.residual(U_final)
        replay = replay_metrics(U_final, final_F, info, contexts,
                                native["restart_path"], native["previous_history"])
        if not replay["replay_guard"]:
            raise RuntimeError(f"{spec['leg']}: production replay differs from round 7")
        wall_restart = time.monotonic() - restart_started
        labels = [event_metrics(e)["ap_step_class"] for e in events]
        result = dict(
            **native["base"], **projection_metrics,
            restart_converged=int(info["converged"]), restart_n_iter=info["n_iter"],
            restart_n_rejected_iterations=sum(
                bool(c["trials"]) and c["iter"] not in {e["iter"] for e in events}
                for c in contexts),
            restart_n_accepted=len(events),
            restart_final_residual=f12(info["final_residual"]),
            restart_cfl_final=f12(info["cfl_final"]), **replay,
            ap_steps_reduced=labels.count("AP_REDUCED"),
            ap_steps_amplified=labels.count("AP_AMPLIFIED"),
            ap_steps_flat=labels.count("AP_FLAT"), **accepted_class(events),
            wall_svd_s=f12(wall_svd), wall_restart_s=f12(wall_restart),
            wall_s=f12(time.monotonic() - started))
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
                                 capture_output=True, timeout=WORKER_BUDGET_S)
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
    return dict(fresh_svd=0, fresh_restarts=0, cumulative_wall_s=0.0,
                max_svd=MAX_LEGS, max_restarts=MAX_LEGS,
                worker_budget_s=WORKER_BUDGET_S, time_budget_s=TIME_BUDGET_S,
                completed=[])


def decide(rows):
    if len(rows) != MAX_LEGS:
        return "UNDEFINED"
    guards = ("provenance_guard", "native_guard", "projection_guard", "replay_guard")
    if any(r.get("status") != "ok" or any(r.get(g) not in (1, "1", True) for g in guards)
           for r in rows):
        return "FAIL_INSTRUMENT"
    primary = [r for r in rows if r["role"] == "primary"]
    if len(primary) != 2:
        return "FAIL_INSTRUMENT"
    classes = [int(r["projection_class"]) for r in primary]
    if classes == [1, 1]:
        return "PASS"
    if classes == [0, 0]:
        return "FAIL"
    return "UNDEFINED"


def _self_test_fake_solver(bad=False):
    import numpy as np
    import scipy.sparse as sp
    from pyfp3d.viscous.ibl3 import IBL3Solver

    class FakeSolver:
        solve = IBL3Solver.solve

        def __init__(self):
            self.n = 1
            self._q = np.ones(1)
            self.smesh = types.SimpleNamespace(node_area=np.ones(1))
            self.eps_diff = 0.005
            self.eps_diff_s = 0.02
            self._veps = 0.005
            self._veps_s = 0.02

        def residual(self, U, dt_inv=None, U_old=None):
            U = np.asarray(U, dtype=float).reshape(1, 6)
            R = U.copy()
            if dt_inv is not None:
                R = R + (self.smesh.node_area * dt_inv)[:, None] * (U - U_old)
            return R

        def residual_jacobian(self, U, dt_inv=None, U_old=None):
            U = np.asarray(U, dtype=float).reshape(1, 6)
            R = U.copy()
            if dt_inv is not None:
                R = R + (self.smesh.node_area * dt_inv)[:, None] * (U - U_old)
            if bad:
                diagonal = np.full(6, 1.0e-6)
            else:
                diagonal = np.ones(6)
                if dt_inv is not None:
                    diagonal = diagonal + self.smesh.node_area[0] * dt_inv
            return R, sp.csr_matrix(np.diag(diagonal))

    return FakeSolver()


def self_test():
    import numpy as np
    from pyfp3d.viscous.ibl3 import IBL3Solver

    rows = []

    def check(name, expected, function):
        try:
            actual = function()
            passed = actual == expected
        except Exception as exc:
            actual = f"{type(exc).__name__}: {exc}"
            passed = False
        rows.append(dict(case=name, expected=expected, actual=actual, passed=int(passed)))

    rng = np.random.default_rng(20261103)
    qu, _ = np.linalg.qr(rng.normal(size=(24, 24)))
    qv, _ = np.linalg.qr(rng.normal(size=(24, 24)))
    singular = np.geomspace(1.0, 1.0e-8, 24)
    base = qu @ np.diag(singular) @ qv.T
    J = (np.geomspace(1.0e-3, 1.0e3, 24)[:, None] * base
         * np.geomspace(1.0e2, 1.0e-2, 24)[None, :])
    F = rng.normal(size=24)
    decomposition = projection_decomposition(J, F, k=2)
    check("svd_reconstruction", True, lambda: decomposition["reconstruction_error"] <= 1.0e-10)
    check("projection_partition", True, lambda: decomposition["partition_error"] <= 1.0e-10)
    check("right_partition", True, lambda: decomposition["right_partition_error"] <= 1.0e-10)
    check("projection_guard", 1, lambda: decomposition["metrics"]["projection_guard"])

    row_scale = decomposition["row_scale"]
    bottom_target = decomposition["left"][:, -2:] @ np.array([0.3, -0.7])
    F_bottom = bottom_target / row_scale
    bottom_decomp = projection_decomposition(J, F_bottom, k=2)
    check("known_bottom_projection", True,
          lambda: abs(float(bottom_decomp["metrics"]["capture20_full"]) - 1.0) <= 1.0e-10)
    expected_mask = np.zeros(24)
    expected_mask[[1, 3, 7, 9, 13, 15, 19, 21]] = 1.0
    check("ap_mask_layout", True, lambda: np.array_equal(block_mask(4, AP_VARS), expected_mask))

    def zero_target():
        try:
            capture_energy(np.ones(4), np.arange(2, 4), np.zeros(4))
        except ValueError:
            return True
        return False

    check("zero_projection_target_has_teeth", True, zero_target)
    check("projection_class_positive", 1, lambda: dominant_ap_projection(0.8, 0.7, 1.0e-5))
    check("projection_class_rejects_capture", 0, lambda: dominant_ap_projection(0.8, 0.49, 1.0e-5))
    check("projection_class_rejects_right", 0, lambda: dominant_ap_projection(0.49, 0.7, 1.0e-5))
    check("projection_class_rejects_hard_spectrum", 0, lambda: dominant_ap_projection(0.8, 0.7, 1.0e-2))

    normal = _self_test_fake_solver()
    normal_baseline_U, normal_baseline_info = IBL3Solver.solve(
        normal, np.ones((1, 6)), tol=0.6, max_iter=4)
    normal_trace_U, normal_trace_info, normal_contexts, normal_events = solve_with_trace(
        normal, np.ones((1, 6)), tol=0.6, max_iter=4)
    check("trace_normal_replay", True,
          lambda: np.array_equal(normal_trace_U, normal_baseline_U)
          and normal_trace_info["residual_history"] == normal_baseline_info["residual_history"])
    check("trace_normal_one_accept", 1, lambda: len(normal_events))
    check("trace_normal_zero_residual", "AP_REDUCED",
          lambda: event_metrics(normal_events[0])["ap_step_class"])

    def export_normal_trace():
        paths = [HERE / f"{prefix}_selftest.csv" for prefix in (
            "iteration_trace", "accepted_steps", "accepted_event_metrics")]
        try:
            iterations, steps, summaries = write_trace(
                "selftest", "control", normal_contexts, normal_events)
            return (len(iterations) == 2 and len(steps) == 6 and len(summaries) == 1
                    and all(path.exists() for path in paths))
        finally:
            for path in paths:
                path.unlink(missing_ok=True)

    check("trace_export_normal", True, export_normal_trace)

    bad = _self_test_fake_solver(bad=True)
    bad_baseline_U, bad_baseline_info = IBL3Solver.solve(
        bad, np.ones((1, 6)), max_iter=20)
    bad_trace_U, bad_trace_info, bad_contexts, bad_events = solve_with_trace(
        bad, np.ones((1, 6)), max_iter=20)
    check("trace_bad_replay", True,
          lambda: np.array_equal(bad_trace_U, bad_baseline_U)
          and bad_trace_info["residual_history"] == bad_baseline_info["residual_history"])
    check("trace_bad_all_rejected", True,
          lambda: not bad_events and len(bad_contexts) == 11
          and all(len(c["trials"]) == 12 for c in bad_contexts))
    check("no_accepted_steps", "NO_ACCEPTED_STEPS", lambda: accepted_class([])["accepted_step_class"])
    check("ap_reduced", "AP_REDUCED", lambda: classify_ap_step(1.0, 0.5))
    check("ap_amplified", "AP_AMPLIFIED", lambda: classify_ap_step(1.0, 1.5))
    check("ap_flat", "AP_FLAT", lambda: classify_ap_step(1.0, 1.0))

    base_rows = [dict(leg=s["leg"], role=s["role"], status="ok", provenance_guard=1,
                      native_guard=1, projection_guard=1, replay_guard=1,
                      projection_class=1 if s["role"] == "primary" else 0)
                 for s in LEG_SPECS]
    check("decision_pass", "PASS", lambda: decide(base_rows))
    failed_projection = [dict(r, projection_class=0) for r in base_rows]
    check("decision_physical_fail", "FAIL", lambda: decide(failed_projection))
    mixed = [dict(r) for r in base_rows]
    mixed[1]["projection_class"] = 0
    check("decision_mixed_undefined", "UNDEFINED", lambda: decide(mixed))
    broken = [dict(r) for r in base_rows]
    broken[0]["replay_guard"] = 0
    check("decision_replay_has_teeth", "FAIL_INSTRUMENT", lambda: decide(broken))
    broken_svd = [dict(r) for r in base_rows]
    broken_svd[0]["projection_guard"] = 0
    check("decision_svd_has_teeth", "FAIL_INSTRUMENT", lambda: decide(broken_svd))

    write_csv(HERE / "instrument_self_test.csv",
              ("case", "expected", "actual", "passed"), rows)
    for row in rows:
        print(f"[self-test] {row['case']}: {row['actual']}", flush=True)
    failed = [r["case"] for r in rows if not r["passed"]]
    if failed:
        raise AssertionError(f"self-test failed: {failed}")
    print(f"[self-test] {len(rows)}/{len(rows)} passed", flush=True)


def analyze():
    import numpy as np
    rows = read_csv(HERE / "summary.csv")
    verdict = decide(rows)
    by_leg = {r["leg"]: r for r in rows}
    metrics = {"verdict": verdict, "n_complete": len(rows)}
    for spec in LEG_SPECS:
        row = by_leg.get(spec["leg"], {})
        for key in ("projection_class", "right_AP_share20", "capture20_AP",
                    "capture20_BD", "sigma20_rel", "cond_equilibrated",
                    "restart_n_accepted", "accepted_step_class",
                    "ap_accept_energy_rel_change", "wall_svd_s", "wall_restart_s"):
            metrics[f"{spec['leg']}__{key}"] = row.get(key, "")
    write_csv(HERE / "verdict_metrics.csv", ("metric", "value"),
              [dict(metric=k, value=v) for k, v in metrics.items()])

    self_tests = read_csv(HERE / "instrument_self_test.csv")
    instruments_ok = bool(self_tests) and all(r["passed"] == "1" for r in self_tests)
    guards_ok = len(rows) == MAX_LEGS and all(
        r.get("status") == "ok" and all(r.get(g) == "1" for g in (
            "provenance_guard", "native_guard", "projection_guard", "replay_guard"))
        for r in rows)
    primary = [r for r in rows if r.get("role") == "primary" and r.get("status") == "ok"]
    accepted_primary = sum(int(r.get("restart_n_accepted") or 0) for r in primary)
    c2_verdict = "RECORDED" if guards_ok and accepted_primary else "UNDEFINED"
    criteria = [
        dict(criterion="I0", verdict="PASS" if instruments_ok and guards_ok else "FAIL",
             evidence="instrument_self_test.csv; provenance.csv; summary.csv",
             reading=f"self-test={sum(r['passed'] == '1' for r in self_tests)}/{len(self_tests)}, guards={guards_ok}"),
        dict(criterion="C1", verdict=verdict if verdict != "FAIL_INSTRUMENT" else "FAIL",
             evidence="projection_modes_*.csv; projection_summary.csv; projection.png",
             reading="; ".join(f"{r['leg']}={r.get('projection_class', '')}" for r in primary)),
        dict(criterion="C2", verdict=c2_verdict,
             evidence="accepted_steps_*.csv; accepted_step_summary.csv; accepted_steps.png",
             reading=f"primary accepted events={accepted_primary}"),
        dict(criterion="C3", verdict=verdict,
             evidence="criteria.csv; verdict_metrics.csv", reading=verdict),
        dict(criterion="OVERALL", verdict=verdict,
             evidence="verdict_metrics.csv", reading=verdict),
    ]
    write_csv(HERE / "criteria.csv", ("criterion", "verdict", "evidence", "reading"), criteria)
    write_csv(HERE / "provenance.csv",
              ("leg", "state_file", "state_sha256", "first_file", "first_sha256",
               "mesh_sha256", "U_sha256", "phi_sha256", "delta_sha256", "blas_threads",
               "provenance_guard", "recorded_residual", "native_residual",
               "native_rel_error", "native_guard", "replay_u_equal",
               "replay_residual_equal", "replay_history_max_rel", "replay_cfl_max_rel",
               "replay_guard"), rows)
    write_csv(HERE / "projection_summary.csv",
              ("leg", "role", "n_dof", "sigma_max", "sigma_min", "cond_equilibrated",
               "sigma20_rel", "n_lt_1e-6", "n_lt_1e-8", "n_lt_1e-10",
               "row_norm_range", "col_norm_range", "right_AP_share20",
               "right_BD_share20", "n_right_AP_dominant20", "capture20_full",
               "capture20_AP", "capture20_BD", "capture20_AP_dominant_modes",
               "dimension_baseline20", "AP_capture_over_baseline",
               "svd_reconstruction_rel_error", "projection_partition_error",
               "right_partition_error", "projection_guard", "projection_class",
               "wall_svd_s"), rows)

    accepted_rows = []
    iteration_rows = []
    for spec in LEG_SPECS:
        accepted_rows.extend(read_csv(HERE / f"accepted_event_metrics_{spec['leg']}.csv"))
        iteration_rows.extend(read_csv(HERE / f"iteration_trace_{spec['leg']}.csv"))
    write_csv(HERE / "accepted_event_summary.csv",
              ("leg", "role", "iter", "cfl", "omega", "pure_inf_before",
               "pure_inf_after", "pure_l2_before", "pure_l2_after", "pt_inf_before",
               "pt_inf_after", "ap_l2_before", "ap_l2_after", "ap_l2_ratio",
               "bd_l2_before", "bd_l2_after", "ap_step_class", "dF_AP_l2", "dU_AP_l2"),
              accepted_rows)
    write_csv(HERE / "accepted_step_summary.csv",
              ("leg", "role", "restart_n_accepted", "restart_n_rejected_iterations",
               "ap_steps_reduced", "ap_steps_amplified", "ap_steps_flat",
               "ap_accept_energy_before", "ap_accept_energy_after",
               "ap_accept_energy_rel_change", "accepted_step_class", "wall_restart_s"), rows)
    write_csv(HERE / "iteration_summary.csv",
              ("leg", "role", "iter", "cfl", "pure_inf", "pt_inf", "n_trials",
               "accepted", "omega_accept", "rejected_iteration"), iteration_rows)
    return metrics


def plot_results():
    import matplotlib.pyplot as plt
    import numpy as np

    rows = [r for r in read_csv(HERE / "summary.csv") if r.get("status") == "ok"]
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for row in rows:
        modes = read_csv(HERE / f"projection_modes_{row['leg']}.csv")
        sigma_rel = np.array([float(r["sigma_rel"]) for r in modes])
        axes[0].semilogy(np.arange(sigma_rel.size), np.maximum(sigma_rel, 1.0e-300),
                         label=row["leg"], lw=1.0)
    axes[0].axhline(SOFT_REL, color="gray", ls=":", lw=0.8)
    axes[0].axvline(int(rows[0]["n_dof"]) - N_NULL, color="gray", ls="--", lw=0.8)
    axes[0].set(xlabel="descending singular-mode index", ylabel=r"$\sigma_i/\sigma_0$",
                title="Equilibrated Jacobian spectrum")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    labels = [r["leg"] for r in rows]
    x = np.arange(len(rows))
    width = 0.25
    for offset, key, name in ((-width, "right_AP_share20", "right A/Psi share"),
                              (0.0, "capture20_AP", "A/Psi residual capture"),
                              (width, "capture20_BD", "B/delta residual capture")):
        axes[1].bar(x + offset, [float(r[key]) for r in rows], width, label=name)
    axes[1].axhline(DOMINANT, color="gray", ls=":", lw=0.8)
    axes[1].set_xticks(x, labels, rotation=20, ha="right")
    axes[1].set(ylim=(0.0, 1.05), title="Bottom-20 projection signatures")
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "projection.png", dpi=150)
    plt.close(fig)

    events = read_csv(HERE / "accepted_event_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for row in rows:
        leg_events = [r for r in events if r["leg"] == row["leg"]]
        if not leg_events:
            continue
        iteration = np.array([int(r["iter"]) for r in leg_events])
        for axis, before, after, title in (
                (axes[0], "ap_l2_before", "ap_l2_after", "Accepted-step A/Psi residual"),
                (axes[1], "bd_l2_before", "bd_l2_after", "Accepted-step B/delta residual")):
            ratio = np.array([float(r[after]) / max(float(r[before]), 1.0e-300)
                              for r in leg_events])
            axis.semilogy(iteration, np.maximum(ratio, 1.0e-300), "o-", ms=3,
                          label=row["leg"])
    for axis in axes:
        axis.axhline(1.0, color="gray", ls=":", lw=0.8)
        axis.set(xlabel="production restart iteration", ylabel="after / before L2")
        axis.grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "accepted_steps.png", dpi=150)
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
        max_dense_svd=MAX_LEGS, dense_svd_per_leg=1, max_restarts=MAX_LEGS,
        fresh_coupling_legs=0, fresh_external_sa_sst=0, iterative_solvers=0,
        projection_rank=N_NULL, worker_budget_s=WORKER_BUDGET_S,
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
    write_environment()
    state = run_state()
    for spec in LEG_SPECS:
        rows = {r["leg"]: r for r in read_csv(HERE / "summary.csv")}
        if rows.get(spec["leg"], {}).get("status") == "ok":
            print(f"[skip] {spec['leg']} already complete", flush=True)
            continue
        if state["fresh_svd"] >= MAX_LEGS or state["fresh_restarts"] >= MAX_LEGS:
            raise RuntimeError("SVD/restart budget exhausted")
        if state["cumulative_wall_s"] >= TIME_BUDGET_S:
            raise RuntimeError("diagnostic wall-time budget exhausted")
        elapsed = run_worker(spec)
        state["fresh_svd"] += 1
        state["fresh_restarts"] += 1
        state["cumulative_wall_s"] += elapsed
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
