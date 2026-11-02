r"""Phase 7: accepted-step nonlinear remainder in the fixed left complement.

Pre-registration: docs/dev_phase_seven/20261103-0535-ibl-step-nonlinearity-prereg.md.
Only the 18 committed xtr030 B1 accepted events are replayed; no production solve.
"""
import argparse
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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WORK = HERE / "_work"
PREVIOUS = REPO / "bench" / "studies" / "ibl_ap_projection"
sys.path.insert(0, str(REPO))


def _load_previous():
    spec = importlib.util.spec_from_file_location(
        "ibl_ap_projection_run", PREVIOUS / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


previous = _load_previous()
floor = previous.floor
read_csv = floor.read_csv
write_csv = floor.write_csv
read_json = floor.read_json
write_json = floor.write_json
sha256_file = floor.sha256_file
f12 = floor.f12
configure_worker = floor.configure_worker

SPEC = previous.LEG_SPECS[0]
LEG = SPEC["leg"]
VAR_NAMES = previous.VAR_NAMES
REGIONS = ("FULL", "BOTTOM20", "COMPLEMENT")
EVENT_ITERS = (0, 1, 2, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 21, 23, 25, 27)
MAX_EVENTS = 18
MAX_SVD = 1
LINEAR_MAX = 0.10
DIRECTION_DEADBAND = 1.0e-12
WORKER_BUDGET_S = 10.0 * 60.0
TIME_BUDGET_S = 15.0 * 60.0

SUMMARY_FIELDS = (
    "leg", "role", "x_tr", "status", "error", "blas_threads", "n_threads",
    "state_file", "state_sha256", "first_file", "first_sha256", "mesh_sha256",
    "U_sha256", "phi_sha256", "delta_sha256", "provenance_guard",
    "recorded_residual", "native_residual", "native_rel_error", "native_guard",
    "n_events", "n_direction_replays", "n_svd", "replay_roundtrip_mismatches",
    "replay_u_equal", "replay_residual_equal", "replay_guard",
    "projection_sigma_mismatches", "projection_capture_mismatches",
    "svd_reconstruction_rel_error", "svd_projection_partition_error",
    "projection_guard", "response_partition_max_error", "source_partition_max_error",
    "energy_guard", "wall_s",
)

RAW_FIELDS = (
    "leg", "iter", "cfl", "omega", "node", "equation", "U_before", "U_after",
    "dU", "F_before", "J_dU", "F_pred", "F_after", "D_actual", "N_nonlinear",
    "Rpt_before", "Rpt_after",
)
METRIC_FIELDS = (
    "leg", "iter", "cfl", "omega", "equation", "region",
    "P_l2", "D_l2", "N_l2", "F0_l2", "Fpred_l2", "F1_l2",
    "P_energy", "D_energy", "N_energy", "F0_energy", "Fpred_energy", "F1_energy",
    "dot_PD", "cos_PD", "N_over_predactual_sum", "response_scale_S",
    "N_over_S", "N_over_P", "N_over_D", "response_class",
    "predicted_direction", "actual_direction", "direction_class",
)
CROSS_FIELDS = tuple(f"{name}_{suffix}" for name in ("P", "D", "N")
                     for suffix in ("source_energy_sum", "cross_energy"))
EVENT_FIELDS = tuple(
    field for field in METRIC_FIELDS if field != "equation") + CROSS_FIELDS


def _formatted(rows):
    return [{k: f12(v) if isinstance(v, float) else v for k, v in row.items()}
            for row in rows]


def _write_formatted(path, fields, rows):
    write_csv(path, fields, _formatted(rows))


def response_quantities(F0, F1, J, dU):
    import numpy as np
    F0 = np.asarray(F0, dtype=float)
    F1 = np.asarray(F1, dtype=float)
    dU = np.asarray(dU, dtype=float)
    if F0.shape != F1.shape or F0.shape != dU.shape or F0.ndim != 2 or F0.shape[1] != 6:
        raise ValueError("F0/F1/dU must have identical (n,6) shapes")
    if J.shape != (F0.size, F0.size):
        raise ValueError("J shape does not match the six-equation state")
    P = np.asarray(J @ dU.ravel()).reshape(F0.shape)
    D = F1 - F0
    Fpred = F0 + P
    N = F1 - Fpred
    return dict(P=P, D=D, N=N, F0=F0, Fpred=Fpred, F1=F1)


def classify_response(P_l2, D_l2, N_l2):
    p, d, n = float(P_l2), float(D_l2), float(N_l2)
    if min(p, d, n) < 0.0:
        raise ValueError("norms cannot be negative")
    if p + d == 0.0:
        if n != 0.0:
            raise ValueError("zero P/D with nonzero nonlinear remainder")
        return "EXACT_ZERO"
    if n <= LINEAR_MAX * p:
        return "LINEAR_DOMINATED"
    if n >= p and n >= d:
        return "NONLINEAR_DOMINATED"
    return "MIXED"


def norm_direction(before, after):
    before, after = float(before), float(after)
    if before < 0.0 or after < 0.0:
        raise ValueError("norms cannot be negative")
    if after < before * (1.0 - DIRECTION_DEADBAND):
        return "REDUCED"
    if after > before * (1.0 + DIRECTION_DEADBAND):
        return "AMPLIFIED"
    return "FLAT"


def direction_pair(before, predicted_after, actual_after):
    predicted = norm_direction(before, predicted_after)
    actual = norm_direction(before, actual_after)
    if predicted == actual:
        label = {"REDUCED": "BOTH_REDUCED", "AMPLIFIED": "BOTH_AMPLIFIED",
                 "FLAT": "BOTH_FLAT"}[predicted]
    elif predicted == "FLAT":
        label = "PRED_FLAT_ONLY"
    elif actual == "FLAT":
        label = "ACTUAL_FLAT_ONLY"
    else:
        label = "DIRECTION_MISMATCH"
    return predicted, actual, label


def project_regions(X, equation, row_scale, left20):
    import numpy as np
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != 6:
        raise ValueError("projection input must be (n,6)")
    if equation is not None:
        mask = np.zeros_like(X)
        mask[:, equation] = 1.0
        X = X * mask
    full = X.ravel() * np.asarray(row_scale, dtype=float)
    bottom = np.asarray(left20, dtype=float) @ (np.asarray(left20, dtype=float).T @ full)
    return full, bottom, full - bottom


def metric_values(vectors):
    import numpy as np
    norms = {name: float(np.linalg.norm(vectors[name]))
             for name in ("P", "D", "N", "F0", "Fpred", "F1")}
    p, d, n = norms["P"], norms["D"], norms["N"]
    dot = float(vectors["P"] @ vectors["D"])
    denominator = norms["Fpred"] + norms["F1"]
    scale = p + d
    predicted, actual, direction = direction_pair(norms["F0"], norms["Fpred"], norms["F1"])
    return dict(
        P_l2=p, D_l2=d, N_l2=n, F0_l2=norms["F0"],
        Fpred_l2=norms["Fpred"], F1_l2=norms["F1"],
        P_energy=p * p, D_energy=d * d, N_energy=n * n,
        F0_energy=norms["F0"] ** 2, Fpred_energy=norms["Fpred"] ** 2,
        F1_energy=norms["F1"] ** 2, dot_PD=dot,
        cos_PD=dot / (p * d) if p > 0.0 and d > 0.0 else "",
        N_over_predactual_sum=n / denominator if denominator > 0.0 else "",
        response_scale_S=scale, N_over_S=n / scale if scale > 0.0 else "",
        N_over_P=n / p if p > 0.0 else "", N_over_D=n / d if d > 0.0 else "",
        response_class=classify_response(p, d, n),
        predicted_direction=predicted, actual_direction=actual,
        direction_class=direction,
    )


def replay_path(solver, U0, event_specs):
    import numpy as np
    import scipy.sparse.linalg as spla
    U = np.array(U0, dtype=float, copy=True)
    U_old = U.copy()
    h = np.sqrt(np.maximum(solver.smesh.node_area, 1.0e-30))
    events = []
    for spec in event_specs:
        cfl, omega = float(spec["cfl"]), float(spec["omega"])
        dt_inv = solver._q / (cfl * h)
        U_before = U.copy()
        Rpt0, Jpt = solver.residual_jacobian(U_before, dt_inv=dt_inv, U_old=U_old)
        delta = spla.spsolve(Jpt, -Rpt0.ravel())
        if not np.all(np.isfinite(delta)):
            raise RuntimeError("non-finite offline Newton direction")
        dU = (omega * delta).reshape(U.shape)
        U_after = U_before + dU
        Rpt1 = solver.residual(U_after, dt_inv=dt_inv, U_old=U_old)
        F0 = solver.residual(U_before)
        F1 = solver.residual(U_after)
        events.append(dict(iter=int(spec["iter"]), cfl=cfl, omega=omega,
                           U_before=U_before, U_after=U_after, dU=dU,
                           F0=F0, F1=F1, Rpt0=Rpt0, Rpt1=Rpt1))
        U_old = U_before
        U = U_after
    return events, U


def committed_events():
    metrics = [r for r in read_csv(
        PREVIOUS / f"accepted_event_metrics_{LEG}.csv") if r["leg"] == LEG]
    if len(metrics) != MAX_EVENTS or tuple(int(r["iter"]) for r in metrics) != EVENT_ITERS:
        raise RuntimeError("committed accepted-event set differs from pre-registration")
    specs = []
    for row in metrics:
        if float(row["omega"]) > 2.0 ** -8:
            raise RuntimeError("committed event is outside the strong-backtracking domain")
        specs.append(dict(iter=int(row["iter"]), cfl=float(row["cfl"]),
                          omega=float(row["omega"])))
    rows = read_csv(PREVIOUS / f"accepted_steps_{LEG}.csv")
    grouped = {it: [] for it in EVENT_ITERS}
    for row in rows:
        grouped[int(row["iter"])].append(row)
    for it in EVENT_ITERS:
        event_rows = grouped[it]
        if len(event_rows) != 410 * 6:
            raise RuntimeError(f"event {it}: expected 2460 committed rows")
        ordered = sorted(event_rows, key=lambda r: (int(r["node"]), VAR_NAMES.index(r["equation"])))
        if [(int(r["node"]), VAR_NAMES.index(r["equation"])) for r in ordered] != [
                (node, eq) for node in range(410) for eq in range(6)]:
            raise RuntimeError(f"event {it}: committed node/equation layout is incomplete")
        grouped[it] = ordered
    return specs, grouped


def compare_replay(events, grouped, restart_path):
    import numpy as np
    mapping = (("F0", "F_before"), ("F1", "F_after"), ("Rpt0", "Rpt_before"),
               ("Rpt1", "Rpt_after"), ("dU", "dU"))
    rows = []
    total_mismatches = 0
    for event in events:
        committed = grouped[event["iter"]]
        row = dict(leg=LEG, iter=event["iter"], cfl=event["cfl"], omega=event["omega"])
        mismatches = 0
        for actual_name, committed_name in mapping + (("D", "dF"),):
            actual = (event["F1"] - event["F0"] if actual_name == "D"
                      else event[actual_name])
            max_abs = 0.0
            max_rel = 0.0
            count = 0
            for node in range(actual.shape[0]):
                for eq in range(6):
                    record = committed[node * 6 + eq]
                    expected = record[committed_name]
                    value = float(actual[node, eq])
                    difference = abs(value - float(expected))
                    max_abs = max(max_abs, difference)
                    max_rel = max(max_rel, difference / max(abs(float(expected)), 1.0e-300))
                    count += int(f12(value) != expected)
            mismatches += count
            row[f"{committed_name}_max_abs"] = max_abs
            row[f"{committed_name}_max_rel"] = max_rel
            row[f"{committed_name}_mismatches"] = count
        row["roundtrip_mismatches"] = mismatches
        total_mismatches += mismatches
        rows.append(row)
    with np.load(restart_path) as data:
        previous_U = data["U"]
        previous_F = data["residual"]
    final_U = events[-1]["U_after"]
    final_F = events[-1]["F1"]
    u_equal = bool(np.array_equal(final_U, previous_U))
    residual_equal = bool(np.array_equal(final_F, previous_F))
    guard = total_mismatches == 0 and u_equal and residual_equal
    return rows, dict(replay_roundtrip_mismatches=total_mismatches,
                      replay_u_equal=int(u_equal), replay_residual_equal=int(residual_equal),
                      replay_guard=int(guard))


def save_replay(events):
    import numpy as np
    keys = ("U_before", "U_after", "dU", "F0", "F1", "Rpt0", "Rpt1")
    arrays = {key: np.asarray([event[key] for event in events]) for key in keys}
    arrays.update(iter=np.array([event["iter"] for event in events]),
                  cfl=np.array([event["cfl"] for event in events]),
                  omega=np.array([event["omega"] for event in events]))
    np.savez(WORK / "accepted_path_replay.npz", **arrays)


def projection_guard(native):
    solver, U, F = native["solver"], native["U"], native["F"]
    _, J = solver.residual_jacobian(U)
    decomposition = previous.projection_decomposition(J.toarray(), F)
    modes = read_csv(PREVIOUS / f"projection_modes_{LEG}.csv")
    if len(modes) != decomposition["sigma"].size:
        raise RuntimeError("round-8 spectrum length mismatch")
    sigma_mismatches = sum(f12(value) != row["sigma"]
                           for value, row in zip(decomposition["sigma"], modes))
    summary = {r["leg"]: r for r in read_csv(PREVIOUS / "projection_summary.csv")}[LEG]
    capture_mismatches = sum(
        f12(decomposition["metrics"][key]) != summary[key]
        for key in ("capture20_full", "capture20_AP", "capture20_BD"))
    guard = bool(decomposition["guard"] and sigma_mismatches == 0
                 and capture_mismatches == 0)
    row = dict(leg=LEG, n_svd=1, projection_sigma_mismatches=sigma_mismatches,
               projection_capture_mismatches=capture_mismatches,
               svd_reconstruction_rel_error=decomposition["reconstruction_error"],
               svd_projection_partition_error=decomposition["partition_error"],
               right_partition_error=decomposition["right_partition_error"],
               projection_guard=int(guard))
    _write_formatted(HERE / "projection_guard.csv", tuple(row), [row])
    import numpy as np
    np.savez(WORK / "fixed_left_projector.npz", row_scale=decomposition["row_scale"],
             left20=decomposition["left"][:, -previous.N_NULL:],
             sigma=decomposition["sigma"])
    if not guard:
        raise RuntimeError("fixed left projector does not reproduce round-8 evidence")
    decomposition["guard_metrics"] = row
    return decomposition


def _partition_error(full, bottom, complement):
    import numpy as np
    full_energy = float(full @ full)
    other = float(bottom @ bottom) + float(complement @ complement)
    if full_energy == 0.0:
        return 0.0 if other == 0.0 else float("inf")
    return abs(full_energy - other) / full_energy


def measure_events(solver, events, decomposition, output_dir=HERE):
    import numpy as np
    row_scale = decomposition["row_scale"]
    left20 = decomposition["left"][:, -previous.N_NULL:]
    raw_rows = []
    metric_rows = []
    event_rows = []
    guard_rows = []
    for event in events:
        _, J = solver.residual_jacobian(event["U_before"])
        quantities = response_quantities(event["F0"], event["F1"], J, event["dU"])
        for node in range(event["dU"].shape[0]):
            for eq, name in enumerate(VAR_NAMES):
                raw_rows.append(dict(
                    leg=LEG, iter=event["iter"], cfl=event["cfl"], omega=event["omega"],
                    node=node, equation=name, U_before=event["U_before"][node, eq],
                    U_after=event["U_after"][node, eq], dU=event["dU"][node, eq],
                    F_before=event["F0"][node, eq], J_dU=quantities["P"][node, eq],
                    F_pred=quantities["Fpred"][node, eq], F_after=event["F1"][node, eq],
                    D_actual=quantities["D"][node, eq], N_nonlinear=quantities["N"][node, eq],
                    Rpt_before=event["Rpt0"][node, eq], Rpt_after=event["Rpt1"][node, eq]))
        response_error = 0.0
        source_error = 0.0
        source_energies = {name: {region: np.zeros(6) for region in REGIONS}
                           for name in quantities}
        source_vectors = {name: {region: np.zeros(row_scale.size) for region in REGIONS}
                          for name in quantities}
        for eq, name in enumerate(VAR_NAMES):
            projected = {}
            for quantity, values in quantities.items():
                regions = project_regions(values, eq, row_scale, left20)
                projected[quantity] = regions
                response_error = max(response_error, _partition_error(*regions))
                for region, vector in zip(REGIONS, regions):
                    source_energies[quantity][region][eq] = float(vector @ vector)
                    source_vectors[quantity][region] += vector
            for region_index, region in enumerate(REGIONS):
                vectors = {quantity: projected[quantity][region_index]
                           for quantity in quantities}
                metric_rows.append(dict(
                    leg=LEG, iter=event["iter"], cfl=event["cfl"], omega=event["omega"],
                    equation=name, region=region, **metric_values(vectors)))
        direct = {}
        for quantity, values in quantities.items():
            regions = project_regions(values, None, row_scale, left20)
            direct[quantity] = regions
            response_error = max(response_error, _partition_error(*regions))
            for region_index, region in enumerate(REGIONS):
                direct_norm = float(np.linalg.norm(regions[region_index]))
                difference = float(np.linalg.norm(
                    regions[region_index] - source_vectors[quantity][region]))
                if direct_norm == 0.0:
                    error = 0.0 if difference == 0.0 else float("inf")
                else:
                    error = difference / direct_norm
                source_error = max(source_error, error)
        for region_index, region in enumerate(REGIONS):
            vectors = {quantity: direct[quantity][region_index] for quantity in quantities}
            cross = {}
            for quantity in ("P", "D", "N"):
                source_energy = float(np.sum(source_energies[quantity][region]))
                cross[f"{quantity}_source_energy_sum"] = source_energy
                cross[f"{quantity}_cross_energy"] = (
                    float(vectors[quantity] @ vectors[quantity]) - source_energy)
            event_rows.append(dict(
                leg=LEG, iter=event["iter"], cfl=event["cfl"], omega=event["omega"],
                region=region, **metric_values(vectors), **cross))
        guard_rows.append(dict(leg=LEG, iter=event["iter"],
                               response_partition_error=response_error,
                               source_partition_error=source_error,
                               energy_guard=int(response_error <= 1.0e-10
                                                and source_error <= 1.0e-10)))
    _write_formatted(output_dir / f"accepted_nonlinearity_{LEG}.csv", RAW_FIELDS, raw_rows)
    _write_formatted(output_dir / "event_equation_metrics.csv", METRIC_FIELDS, metric_rows)
    _write_formatted(output_dir / "event_summary.csv", EVENT_FIELDS, event_rows)
    _write_formatted(output_dir / "energy_guard.csv",
                     ("leg", "iter", "response_partition_error", "source_partition_error",
                      "energy_guard"), guard_rows)
    response_max = max(r["response_partition_error"] for r in guard_rows)
    source_max = max(r["source_partition_error"] for r in guard_rows)
    guard = response_max <= 1.0e-10 and source_max <= 1.0e-10
    return dict(response_partition_max_error=response_max,
                source_partition_max_error=source_max, energy_guard=int(guard))


def write_provenance(native):
    files = {
        "accepted_event_file": PREVIOUS / f"accepted_event_metrics_{LEG}.csv",
        "accepted_step_file": PREVIOUS / f"accepted_steps_{LEG}.csv",
        "projection_mode_file": PREVIOUS / f"projection_modes_{LEG}.csv",
        "projection_summary_file": PREVIOUS / "projection_summary.csv",
        "restart_file": native["restart_path"],
    }
    row = dict(native["base"])
    for name, path in files.items():
        row[name] = str(path.relative_to(REPO))
        row[f"{name}_sha256"] = sha256_file(path)
    fields = tuple(native["base"]) + tuple(
        item for name in files for item in (name, f"{name}_sha256"))
    write_csv(HERE / "provenance.csv", fields, [row])


def decide(summary, metric_rows):
    guards = ("provenance_guard", "native_guard", "replay_guard", "projection_guard",
              "energy_guard")
    if summary.get("status") != "ok" or any(summary.get(g) not in (1, "1", True) for g in guards):
        return "FAIL_INSTRUMENT"
    if int(summary.get("n_events", 0)) == 0 or not metric_rows:
        return "UNDEFINED"
    if all(r["response_class"] == "EXACT_ZERO" for r in metric_rows):
        return "UNDEFINED"
    return "RECORDED"


def _stats(values):
    import numpy as np
    values = [float(v) for v in values if v != ""]
    if not values:
        return "", "", ""
    return float(np.median(values)), float(np.min(values)), float(np.max(values))


def analyze():
    import numpy as np
    summary_rows = read_csv(HERE / "summary.csv")
    summary = summary_rows[0] if summary_rows else {"status": "missing"}
    metrics = read_csv(HERE / "event_equation_metrics.csv")
    equation_rows = []
    for equation in VAR_NAMES:
        for region in REGIONS:
            rows = [r for r in metrics if r["equation"] == equation and r["region"] == region]
            if not rows:
                continue
            energies = {name: sum(float(r[f"{name}_energy"]) for r in rows)
                        for name in ("P", "D", "N")}
            pooled = classify_response(*(float(np.sqrt(energies[name])) for name in ("P", "D", "N")))
            row = dict(leg=LEG, equation=equation, region=region, n_events=len(rows),
                       P_energy=energies["P"], D_energy=energies["D"], N_energy=energies["N"],
                       P_l2=float(np.sqrt(energies["P"])), D_l2=float(np.sqrt(energies["D"])),
                       N_l2=float(np.sqrt(energies["N"])), pooled_response_class=pooled)
            for label in ("LINEAR_DOMINATED", "NONLINEAR_DOMINATED", "MIXED", "EXACT_ZERO"):
                row[f"n_{label}"] = sum(r["response_class"] == label for r in rows)
            for direction in ("BOTH_REDUCED", "BOTH_AMPLIFIED", "BOTH_FLAT",
                              "DIRECTION_MISMATCH", "PRED_FLAT_ONLY", "ACTUAL_FLAT_ONLY"):
                row[f"n_{direction}"] = sum(r["direction_class"] == direction for r in rows)
            for key in ("N_over_S", "N_over_P", "cos_PD", "N_over_predactual_sum"):
                median, minimum, maximum = _stats([r[key] for r in rows])
                row[f"{key}_median"] = median
                row[f"{key}_min"] = minimum
                row[f"{key}_max"] = maximum
            for name in ("P", "D", "N"):
                row[f"max_single_{name}_energy_share"] = max(
                    float(r[f"{name}_energy"]) for r in rows) / max(energies[name], 1.0e-300)
            equation_rows.append(row)
    fields = ("leg", "equation", "region", "n_events", "P_l2", "D_l2", "N_l2",
              "P_energy", "D_energy", "N_energy", "pooled_response_class")
    fields += tuple(f"n_{label}" for label in (
        "LINEAR_DOMINATED", "NONLINEAR_DOMINATED", "MIXED", "EXACT_ZERO"))
    fields += tuple(f"n_{label}" for label in (
        "BOTH_REDUCED", "BOTH_AMPLIFIED", "BOTH_FLAT", "DIRECTION_MISMATCH",
        "PRED_FLAT_ONLY", "ACTUAL_FLAT_ONLY"))
    fields += tuple(f"{key}_{stat}" for key in (
        "N_over_S", "N_over_P", "cos_PD", "N_over_predactual_sum")
        for stat in ("median", "min", "max"))
    fields += tuple(f"max_single_{name}_energy_share" for name in ("P", "D", "N"))
    _write_formatted(HERE / "equation_summary.csv", fields, equation_rows)

    self_tests = read_csv(HERE / "instrument_self_test.csv")
    instruments_ok = bool(self_tests) and all(r["passed"] == "1" for r in self_tests)
    verdict = decide(summary, metrics) if instruments_ok else "FAIL_INSTRUMENT"
    complete = len(metrics) == MAX_EVENTS * 6 * 3 and len(equation_rows) == 18
    i0 = verdict != "FAIL_INSTRUMENT" and instruments_ok
    criteria = [
        dict(criterion="I0", verdict="PASS" if i0 else "FAIL",
             evidence="instrument_self_test.csv; provenance.csv; replay_guard.csv; projection_guard.csv; energy_guard.csv",
             reading=f"self-test={sum(r['passed'] == '1' for r in self_tests)}/{len(self_tests)}; guards verdict={verdict}"),
        dict(criterion="C1", verdict="RECORDED" if i0 and complete else verdict,
             evidence="accepted_nonlinearity_*.csv; event_equation_metrics.csv; equation_summary.csv; nonlinearity.png",
             reading=f"event-equation-region rows={len(metrics)}/324"),
        dict(criterion="C2", verdict="RECORDED" if i0 and complete else verdict,
             evidence="event_equation_metrics.csv; event_summary.csv; equation_summary.csv",
             reading="predicted/actual direction classes recorded"),
        dict(criterion="C3", verdict=verdict, evidence="criteria.csv; verdict_metrics.csv",
             reading=verdict),
        dict(criterion="OVERALL", verdict=verdict, evidence="verdict_metrics.csv", reading=verdict),
    ]
    write_csv(HERE / "criteria.csv", ("criterion", "verdict", "evidence", "reading"), criteria)
    verdict_metrics = {"verdict": verdict, "n_complete_metrics": len(metrics),
                       "n_events": summary.get("n_events", 0)}
    for key in SUMMARY_FIELDS:
        if key.endswith("guard") or key in ("wall_s", "replay_roundtrip_mismatches"):
            verdict_metrics[key] = summary.get(key, "")
    for row in equation_rows:
        prefix = f"{row['equation']}__{row['region']}"
        verdict_metrics[f"{prefix}__pooled_response_class"] = row["pooled_response_class"]
        verdict_metrics[f"{prefix}__N_over_P"] = (
            float(row["N_l2"]) / float(row["P_l2"]) if float(row["P_l2"]) > 0.0 else "")
        verdict_metrics[f"{prefix}__n_DIRECTION_MISMATCH"] = row["n_DIRECTION_MISMATCH"]
    _write_formatted(HERE / "verdict_metrics.csv", ("metric", "value"),
                     [dict(metric=key, value=value) for key, value in verdict_metrics.items()])
    return verdict


def plot_results():
    import matplotlib.pyplot as plt
    import numpy as np
    rows = read_csv(HERE / "equation_summary.csv")
    if len(rows) != 18:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.0))
    complement = [r for r in rows if r["region"] == "COMPLEMENT"]
    x = np.arange(6)
    width = 0.25
    for offset, key, label in ((-width, "P_energy", "linear J dU"),
                               (0.0, "D_energy", "actual response"),
                               (width, "N_energy", "nonlinear remainder")):
        axes[0].bar(x + offset, [max(float(r[key]), 1.0e-300) for r in complement],
                    width, label=label)
    axes[0].set_xticks(x, VAR_NAMES)
    axes[0].set(yscale="log", ylabel="pooled response energy",
                title="Fixed left complement: source-equation decomposition")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=8)
    metrics = read_csv(HERE / "event_equation_metrics.csv")
    for eq, name in enumerate(VAR_NAMES):
        selected = [r for r in metrics if r["equation"] == name and r["region"] == "COMPLEMENT"]
        selected.sort(key=lambda r: EVENT_ITERS.index(int(r["iter"])))
        ratios = [float(r["N_over_P"]) if r["N_over_P"] else np.nan for r in selected]
        axes[1].semilogy([int(r["iter"]) for r in selected], ratios, "o-", ms=3, label=name)
    axes[1].axhline(LINEAR_MAX, color="gray", ls=":", lw=0.8, label="0.1 linear bound")
    axes[1].axhline(1.0, color="black", ls="--", lw=0.8, label="N = P")
    axes[1].set(xlabel="production restart iteration", ylabel=r"$||N||_2/||JdU||_2$",
                title="Strong-backtracking events in the left complement")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "nonlinearity.png", dpi=150)
    plt.close(fig)


def self_test():
    import numpy as np
    import scipy.sparse as sp
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
    linear = rng.normal(size=(12, 12))
    quadratic = rng.normal(size=12) * 0.05
    U = rng.normal(size=(2, 6))
    dU = rng.normal(size=(2, 6)) * 0.01

    def residual(value):
        flat = value.ravel()
        return (linear @ flat + 0.5 * quadratic * flat ** 2).reshape(2, 6)

    F0, F1 = residual(U), residual(U + dU)
    J = linear + np.diag(quadratic * U.ravel())
    quantities = response_quantities(F0, F1, sp.csr_matrix(J), dU)
    exact_N = 0.5 * quadratic * dU.ravel() ** 2
    check("quadratic_remainder_oracle", True,
          lambda: np.allclose(quantities["N"].ravel(), exact_N, rtol=1.0e-11, atol=1.0e-13))
    linear_q = response_quantities(F0, F0 + (J @ dU.ravel()).reshape(2, 6),
                                   sp.csr_matrix(J), dU)
    check("linear_response_zero_N", 0.0, lambda: float(np.linalg.norm(linear_q["N"])))
    check("class_linear", "LINEAR_DOMINATED", lambda: classify_response(1.0, 1.0, 0.0))
    check("class_nonlinear", "NONLINEAR_DOMINATED", lambda: classify_response(0.0, 1.0, 1.0))
    check("class_mixed", "MIXED", lambda: classify_response(1.0, 1.5, 0.5))
    check("class_exact_zero", "EXACT_ZERO", lambda: classify_response(0.0, 0.0, 0.0))

    def bad_zero():
        try:
            classify_response(0.0, 0.0, 1.0)
        except ValueError:
            return True
        return False

    check("zero_guard_has_teeth", True, bad_zero)
    q, _ = np.linalg.qr(rng.normal(size=(24, 24)))
    left20 = q[:, -2:]
    scales = np.geomspace(1.0e-2, 1.0e2, 24)
    X = rng.normal(size=(4, 6))

    def projection_oracle():
        direct = project_regions(X, None, scales, left20)
        errors = [_partition_error(*direct)]
        for region_index in range(3):
            source_vector = np.zeros(24)
            for eq in range(6):
                regions = project_regions(X, eq, scales, left20)
                errors.append(_partition_error(*regions))
                source_vector += regions[region_index]
            errors.append(float(np.linalg.norm(direct[region_index] - source_vector))
                          / float(np.linalg.norm(direct[region_index])))
        return max(errors) <= 1.0e-10

    def projection_cross_term_has_teeth():
        direct = project_regions(X, None, scales, left20)[2]
        direct_energy = float(direct @ direct)
        source_energy = sum(float(project_regions(X, eq, scales, left20)[2] @
                                  project_regions(X, eq, scales, left20)[2])
                            for eq in range(6))
        return abs(direct_energy - source_energy) > 1.0e-6 * direct_energy

    check("known_projection_and_source_partition", True, projection_oracle)
    check("projected_energy_sum_is_false", True, projection_cross_term_has_teeth)
    check("wrong_projection_has_teeth", True, lambda: _partition_error(
        X.ravel() * scales, np.zeros(24), np.zeros(24)) > 1.0e-10)
    check("mask_layout", True, lambda: np.array_equal(
        project_regions(X, 1, np.ones(24), np.zeros((24, 1)))[0].reshape(4, 6)[:, [0, 2, 3, 4, 5]],
        np.zeros((4, 5))))
    check("direction_both_reduced", "BOTH_REDUCED", lambda: direction_pair(1.0, 0.5, 0.6)[2])
    check("direction_both_amplified", "BOTH_AMPLIFIED", lambda: direction_pair(1.0, 1.5, 1.2)[2])
    check("direction_both_flat", "BOTH_FLAT", lambda: direction_pair(1.0, 1.0, 1.0)[2])
    check("direction_mismatch", "DIRECTION_MISMATCH", lambda: direction_pair(1.0, 0.5, 1.2)[2])
    check("direction_zero_before", "BOTH_AMPLIFIED", lambda: direction_pair(0.0, 1.0, 2.0)[2])
    check("direction_zero_pair", "BOTH_FLAT", lambda: direction_pair(0.0, 0.0, 0.0)[2])

    fake = previous._self_test_fake_solver()
    baseline_U, baseline_info = IBL3Solver.solve(fake, np.ones((1, 6)), tol=0.6, max_iter=4)
    traced_U, traced_info, contexts, events = previous.solve_with_trace(
        fake, np.ones((1, 6)), tol=0.6, max_iter=4)
    specs = [dict(iter=e["iter"], cfl=e["cfl"], omega=e["omega"]) for e in events]
    replayed, final_U = replay_path(fake, np.ones((1, 6)), specs)
    check("offline_replay_normal", True, lambda: len(replayed) == 1 and np.array_equal(
        replayed[0]["U_after"], events[0]["U_after"]) and np.array_equal(final_U, baseline_U))
    changed = [dict(specs[0], omega=specs[0]["omega"] * 0.5)]
    changed_replay, changed_final = replay_path(fake, np.ones((1, 6)), changed)
    check("offline_replay_omega_has_teeth", True,
          lambda: not np.array_equal(changed_final, baseline_U))

    def export_normal_measurements():
        paths = [WORK / f"accepted_nonlinearity_{LEG}.csv",
                 WORK / "event_equation_metrics.csv", WORK / "event_summary.csv",
                 WORK / "energy_guard.csv"]
        try:
            decomposition = dict(row_scale=np.ones(6), left=np.eye(6))
            guard = measure_events(fake, replayed, decomposition, output_dir=WORK)
            return (guard["energy_guard"] == 1
                    and len(read_csv(paths[0])) == 6 and len(read_csv(paths[1])) == 18
                    and len(read_csv(paths[2])) == 3 and len(read_csv(paths[3])) == 1)
        finally:
            for path in paths:
                path.unlink(missing_ok=True)

    check("measurement_export_normal", True, export_normal_measurements)

    committed = {}
    for node in range(1):
        for eq, name in enumerate(VAR_NAMES):
            committed[node * 6 + eq] = dict(
                F_before=f12(events[0]["F_before"][node, eq]),
                F_after=f12(events[0]["F_after"][node, eq]),
                Rpt_before=f12(events[0]["R_before"][node, eq]),
                Rpt_after=f12(events[0]["R_after"][node, eq]),
                dU=f12(events[0]["U_after"][node, eq] - events[0]["U_before"][node, eq]),
                dF=f12(events[0]["F_after"][node, eq] - events[0]["F_before"][node, eq]))

    restart_path = WORK / "selftest_restart.npz"
    np.savez(restart_path, U=events[0]["U_after"], residual=events[0]["F_after"])
    replay_rows, replay_summary = compare_replay(replayed, {0: committed}, restart_path)
    check("roundtrip_guard_normal", 1, lambda: replay_summary["replay_guard"])
    committed_bad = [dict(row) for row in committed.values()]
    committed_bad[0]["dU"] = f12(float(committed_bad[0]["dU"]) + 1.0)
    _, bad_summary = compare_replay(replayed, {0: committed_bad}, restart_path)
    check("roundtrip_cell_has_teeth", 0, lambda: bad_summary["replay_guard"])
    restart_path.unlink(missing_ok=True)

    good_summary = dict(status="ok", n_events=18, **{key: 1 for key in (
        "provenance_guard", "native_guard", "replay_guard", "projection_guard", "energy_guard")})
    metric_rows = [dict(response_class="MIXED")]
    check("decision_recorded", "RECORDED", lambda: decide(good_summary, metric_rows))
    check("decision_zero_undefined", "UNDEFINED", lambda: decide(
        good_summary, [dict(response_class="EXACT_ZERO")]))
    check("decision_no_events_undefined", "UNDEFINED", lambda: decide(
        dict(good_summary, n_events=0), []))
    check("decision_guard_has_teeth", "FAIL_INSTRUMENT", lambda: decide(
        dict(good_summary, replay_guard=0), metric_rows))

    write_csv(HERE / "instrument_self_test.csv",
              ("case", "expected", "actual", "passed"), rows)
    for row in rows:
        print(f"[self-test] {row['case']}: {row['actual']}", flush=True)
    failed = [r["case"] for r in rows if not r["passed"]]
    if failed:
        raise AssertionError(f"self-test failed: {failed}")
    print(f"[self-test] {len(rows)}/{len(rows)} passed", flush=True)


def write_environment():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    packages = {name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "numba", "meshio", "pyamg", "pytest")}
    write_json(HERE / "environment.json", dict(
        commit=commit, platform=platform.platform(), processor=platform.processor(),
        python=sys.version, packages=packages, blas_threads=1, numba_threads=16,
        omp_threads=16, cache_policy="no cache deletion; warm cache; fresh worker",
        mesh=str(previous.MESH.relative_to(REPO)), mesh_sha256=sha256_file(previous.MESH),
        max_dense_svd=MAX_SVD, dense_svd_scope="native xtr030 B1 fixed projector only",
        offline_direction_replays=MAX_EVENTS, production_solve_restarts=0,
        fresh_coupling_legs=0, fresh_external_sa_sst=0, other_leg_svd=0,
        projection_rank=previous.N_NULL, worker_budget_s=WORKER_BUDGET_S,
        time_budget_s=TIME_BUDGET_S,
    ))


def worker():
    configure_worker()
    started = time.monotonic()
    try:
        native = previous.build_native(SPEC)
        write_provenance(native)
        specs, grouped = committed_events()
        events, final_U = replay_path(native["solver"], native["U"], specs)
        save_replay(events)
        replay_rows, replay = compare_replay(events, grouped, native["restart_path"])
        replay_fields = ["leg", "iter", "cfl", "omega"]
        for name in ("F_before", "F_after", "Rpt_before", "Rpt_after", "dU", "dF"):
            replay_fields += [f"{name}_max_abs", f"{name}_max_rel", f"{name}_mismatches"]
        replay_fields.append("roundtrip_mismatches")
        replay_fields += ["replay_u_equal", "replay_residual_equal", "replay_guard"]
        for row in replay_rows:
            row.update(replay)
        _write_formatted(HERE / "replay_guard.csv", tuple(replay_fields), replay_rows)
        if not replay["replay_guard"]:
            raise RuntimeError("offline accepted-path replay differs from committed evidence")
        decomposition = projection_guard(native)
        energy = measure_events(native["solver"], events, decomposition)
        if not energy["energy_guard"]:
            raise RuntimeError("response/source projection energy guard failed")
        result = dict(native["base"], n_events=len(events), n_direction_replays=len(events),
                      **replay, **decomposition["guard_metrics"], **energy,
                      wall_s=time.monotonic() - started)
        write_csv(HERE / "summary.csv", SUMMARY_FIELDS,
                  [{key: result.get(key, "") for key in SUMMARY_FIELDS}])
        write_json(WORK / "result.json", dict(result, wall_s=f12(result["wall_s"])))
    except Exception as exc:
        row = {key: "" for key in SUMMARY_FIELDS}
        row.update(leg=LEG, role=SPEC["role"], x_tr=f12(SPEC["x_tr"]), status="failed",
                   error=f"{type(exc).__name__}: {exc}", wall_s=time.monotonic() - started)
        write_csv(HERE / "summary.csv", SUMMARY_FIELDS, [row])
        write_json(WORK / "result.json", dict(row, traceback=traceback.format_exc()))
        raise


def require_self_test():
    rows = read_csv(HERE / "instrument_self_test.csv")
    if not rows or any(r["passed"] != "1" for r in rows):
        raise RuntimeError("run --self-test successfully before formal measurements")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    if args.worker:
        worker()
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
    state_path = HERE / "run_state.json"
    if state_path.exists() and read_json(state_path).get("completed"):
        raise RuntimeError("formal measurement already completed; use --analyze-only")
    output_path = WORK / "worker.txt"
    started = time.monotonic()
    try:
        process = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker"],
                                 cwd=REPO, env=previous.worker_env(), text=True,
                                 capture_output=True, timeout=WORKER_BUDGET_S)
        output_path.write_text(process.stdout + "\n" + process.stderr, encoding="utf-8")
        print(process.stdout, flush=True)
        if process.returncode:
            print(process.stderr, file=sys.stderr, flush=True)
    except subprocess.TimeoutExpired as exc:
        output_path.write_text(str(exc), encoding="utf-8")
        raise RuntimeError("formal worker timeout") from exc
    elapsed = time.monotonic() - started
    write_json(state_path, dict(completed=process.returncode == 0, cumulative_wall_s=elapsed,
                                fresh_svd=MAX_SVD, offline_direction_replays=MAX_EVENTS,
                                production_solve_restarts=0, max_svd=MAX_SVD,
                                worker_budget_s=WORKER_BUDGET_S, time_budget_s=TIME_BUDGET_S))
    analyze()
    plot_results()
    if process.returncode:
        raise RuntimeError(f"formal worker failed; see {output_path}")


if __name__ == "__main__":
    main()
