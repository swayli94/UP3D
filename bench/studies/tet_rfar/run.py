r"""Phase 7 S4: NACA0012 tet r_far 15/20/30 non-convergence attribution.

Pre-registration: docs/dev_phase_seven/20261102-2030-tet-rfar-prereg.md.
The production library, meshes, gates and committed gate evidence are read-only.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WORK = HERE / "_work"
PRODUCTION_MESH = REPO / "cases" / "meshes" / "naca0012_2.5d"
sys.path.insert(0, str(REPO))

LEVELS = ("coarse", "medium")
RADII = (15, 20, 30)
M_INF, ALPHA = 0.80, 1.25
RECIPE = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
              freeze_refresh_max=8, precond="direct", direct_refactor_every=4,
              n_newton_max=400, n_picard_seed=5)
TIME_BUDGET_S = 60.0 * 60.0
SUMMARY_FIELDS = (
    "leg", "mode", "level", "r_far", "converged", "classification",
    "period", "period_score", "accept_reason", "residual", "n_newton",
    "n_points", "n_limited", "n_floored", "n_gmres_total",
    "n_gmres_stalled", "n_refactor", "froze", "n_freeze_refresh",
    "n_freeze_reverts", "sigma_frozen", "sigma_churn", "last_sigma_delta",
    "tail_n_shock", "cl", "cd", "x_shock", "mach_max", "m_final",
    "target_reached", "wall_s", "cumulative_solve_s", "n_threads",
    "blas_threads", "mesh_sha256", "history_csv",
)


def f12(value):
    return "" if value is None else f"{float(value):.12e}"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def classify(residual, converged, n_limited=0, n_floored=0):
    r = np.asarray(residual, dtype=float)
    if n_limited + n_floored > 0:
        return dict(classification="clamped", period="", period_score="")
    if converged:
        return dict(classification="converged", period="", period_score="")
    if len(r) < 201 or not np.all(np.isfinite(r)) or r[-1] <= 1e-10:
        return dict(classification="nonconverged_other", period="",
                    period_score="")
    tail = r[-200:]
    if np.any(tail <= 0.0):
        return dict(classification="nonconverged_other", period="",
                    period_score="")
    ratio = float(np.median(tail[-20:]) / np.median(tail[:20]))
    spread = float(tail.max() / tail.min())
    bounded = 0.5 <= ratio <= 2.0 and spread < 10.0
    if not bounded:
        return dict(classification="nonconverged_other", period="",
                    period_score=f12(spread))
    scale = max(float(np.median(tail)), np.finfo(float).tiny)
    if np.ptp(tail) <= 1e-14 * scale:
        return dict(classification="fixed_point_stall", period=1,
                    period_score=f12(spread))
    best_p, best_score = None, np.inf
    for p in range(2, 9):
        folded = np.concatenate([tail[j::p] - np.mean(tail[j::p])
                                 for j in range(p)])
        score = float(np.sqrt(np.mean(folded ** 2)) / np.ptp(tail))
        if score < best_score - 1e-12:
            best_p, best_score = p, score
    label = ("periodic_limit_cycle" if best_score <= 0.25
             else "bounded_stall")
    return dict(classification=label, period=best_p,
                period_score=f12(best_score))


def self_test():
    i = np.arange(400, dtype=float)
    cases = (
        ("converged", np.geomspace(1.0, 1e-12, 400), True, 0, 0,
         "converged"),
        ("clamped", 2e-7 * (1 + 0.05 * np.sin(i)), False, 1, 0,
         "clamped"),
        ("period4", 2e-7 * (1 + 0.05 * np.sin(2 * np.pi * i / 4)),
         False, 0, 0, "periodic_limit_cycle"),
        ("bounded_nonperiodic",
         2e-7 * (1 + 0.04 * np.sin(0.123 * i) + 0.02 * np.sin(0.317 * i)),
         False, 0, 0, "bounded_stall"),
        ("fixed_point", np.full(400, 2e-7), False, 0, 0,
         "fixed_point_stall"),
    )
    rows = []
    for name, hist, converged, limited, floored, expected in cases:
        got = classify(hist, converged, limited, floored)
        actual = got["classification"]
        rows.append(dict(case=name, expected=expected, actual=actual,
                         period=got["period"], period_score=got["period_score"],
                         passed=actual == expected))
    with open(HERE / "classifier_self_test.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, lineterminator="\n", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"[self-test] {row['case']}: {row['actual']} "
              f"(period={row['period']}, score={row['period_score']})")
    failed = [r["case"] for r in rows if not r["passed"]]
    if failed:
        raise AssertionError(f"classifier self-test failed: {failed}")


def production_generator():
    path = PRODUCTION_MESH / "generate_naca0012.py"
    spec = importlib.util.spec_from_file_location("_naca0012_generator", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def generate_one(level, radius, repeat):
    from pyfp3d.mesh.reader import write_mesh
    from pyfp3d.meshgen.extrude import extrude_single_layer
    from pyfp3d.meshgen.planar import naca0012_wake_2d

    p = production_generator().LEVELS[level]
    points, triangles, edges, interior = naca0012_wake_2d(
        r_far=float(radius), h_wall=p["h_wall"], h_far=p["h_far"],
        h_wake=p["h_wake"], dist_min=p["dist_min"], dist_max=p["dist_max"],
        wake_dist_max=p["wake_dist_max"], n_half=p["n_half"])
    mesh = extrude_single_layer(
        points, triangles, edges, interior_edge_groups=interior, dz=p["dz"],
        z0=0.0, name=f"naca0012_2.5d_{level}_r{radius}")
    path = WORK / f"{level}_r{radius}_{repeat}.msh"
    write_mesh(mesh, path)
    counts = {name: len(faces) for name, faces in mesh.boundary_faces.items()}
    return path, dict(n_nodes=len(mesh.nodes), n_tets=len(mesh.elements),
                      n_tris_farfield=counts.get("farfield", 0),
                      n_tris_wall=counts.get("wall", 0),
                      n_tris_wake=counts.get("wake", 0),
                      n_tris_symmetry=counts.get("symmetry", 0))


def generate_matrix():
    import gmsh

    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for level in LEVELS:
        for radius in RADII:
            print(f"[mesh] {level} r={radius}", flush=True)
            a, counts = generate_one(level, radius, "a")
            b, counts_b = generate_one(level, radius, "b")
            sha_a, sha_b = sha256(a), sha256(b)
            production = PRODUCTION_MESH / f"{level}.msh"
            production_sha = sha256(production) if radius == 15 else ""
            row = dict(level=level, r_far=radius, **counts,
                       bytes=a.stat().st_size, sha256_a=sha_a, sha256_b=sha_b,
                       repeat_match=sha_a == sha_b and counts == counts_b,
                       production_sha256=production_sha,
                       production_match=(radius == 15 and sha_a == production_sha),
                       gmsh=gmsh.__version__, mesh_a=str(a.relative_to(HERE)))
            rows.append(row)
    fields = ("level", "r_far", "n_nodes", "n_tets", "n_tris_farfield",
              "n_tris_wall", "n_tris_wake", "n_tris_symmetry", "bytes",
              "sha256_a", "sha256_b", "repeat_match", "production_sha256",
              "production_match", "gmsh", "mesh_a")
    with open(HERE / "mesh_manifest.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    bad_repeat = [f"{r['level']}:r{r['r_far']}" for r in rows
                  if not r["repeat_match"]]
    bad_reference = [f"{r['level']}:r{r['r_far']}" for r in rows
                     if r["r_far"] == 15 and not r["production_match"]]
    if bad_repeat or bad_reference:
        raise RuntimeError(f"mesh provenance mismatch: repeat={bad_repeat}, "
                           f"r15={bad_reference}; stop before flow solves")
    return rows


def read_rows(name):
    path = HERE / name
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def append_row(path, fields, row):
    new = not path.exists()
    with open(path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, lineterminator="\n", fieldnames=fields)
        if new:
            writer.writeheader()
        writer.writerow(row)


def write_history(leg, result):
    name = f"history_{leg}.csv"
    steps = result.get("step_records") or []
    residual = result["residual_history"]
    with open(HERE / name, "w", newline="") as fh:
        fields = ("step", "residual", "F", "gamma_mean", "gamma_root",
                  "n_limited", "n_floored", "lam", "eta", "frozen",
                  "n_lin_solves", "n_lin_iters", "n_refactor", "sigma_min",
                  "sigma_delta")
        writer = csv.DictWriter(fh, lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        for i, r in enumerate(residual):
            rec = steps[i] if i < len(steps) else {}
            writer.writerow(dict(
                step=i, residual=f12(r),
                F=f12(result["F_history"][i]), gamma_mean=f12(rec.get("gamma_mean")),
                gamma_root=f12(rec.get("gamma_root")),
                n_limited=rec.get("n_limited", ""), n_floored=rec.get("n_floored", ""),
                lam=f12(rec.get("lam")), eta=f12(rec.get("eta")),
                frozen=rec.get("frozen", ""), n_lin_solves=rec.get("n_lin_solves", ""),
                n_lin_iters=rec.get("n_lin_iters", ""),
                n_refactor=rec.get("n_refactor", ""),
                sigma_min=f12(rec.get("sigma_min")),
                sigma_delta=f12(rec.get("sigma_delta"))))
    return name


def write_cp(leg, curve):
    with open(HERE / "cp_curves.csv", "a", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        if fh.tell() == 0:
            writer.writerow(["leg", "side", "x", "cp"])
        for side in ("upper", "lower"):
            x = np.asarray(curve[f"x_{side}"])
            cp = np.asarray(curve[f"cp_{side}"])
            order = np.argsort(x)
            x, cp = x[order], cp[order]
            keep = np.concatenate(([True], np.diff(x) > 0))
            for xx, yy in zip(x[keep], cp[keep]):
                writer.writerow([leg, side, f12(xx), f12(yy)])


def solve_leg(mode, level, radius, solve_seconds):
    from pyfp3d import blas_thread_count
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.section_cut import section_cp_curve
    from pyfp3d.post.shock import shock_report
    from pyfp3d.post.surface import wall_force_coefficients
    from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic

    assert blas_thread_count() == 1, f"BLAS must be 1, got {blas_thread_count()}"
    leg = f"{mode}_{level}_r{radius}"
    mesh_path = WORK / f"{level}_r{radius}_a.msh"
    mesh, wake = cut_wake(read_mesh(mesh_path))
    kwargs = dict(RECIPE)
    if mode == "entropy_off":
        kwargs["entropy_correction"] = False
    print(f"[solve] {leg}; cumulative solver time {solve_seconds:.1f}s", flush=True)
    start = time.perf_counter()
    if mode == "continuation":
        result = solve_newton_transonic(
            mesh, wake, m_inf=M_INF, alpha_deg=ALPHA, m_start=0.70, dm=0.05,
            dm_min=0.01, newton_kw=kwargs)
    else:
        result = solve_newton_lifting(mesh, wake, m_inf=M_INF, alpha_deg=ALPHA,
                                      **kwargs)
    wall = time.perf_counter() - start
    phi = np.asarray(result["phi"])
    m_state = result.get("m_final", M_INF)
    if m_state is None:
        raise RuntimeError(f"{leg}: no final Mach state for post-processing")
    m_state = float(m_state)
    dz = float(np.ptp(mesh.nodes[:, 2]))
    force = wall_force_coefficients(mesh.nodes, mesh.elements,
                                    mesh.boundary_faces["wall"], phi,
                                    alpha_deg=ALPHA, u_inf=1.0, s_ref=dz,
                                    m_inf=m_state)
    curve = section_cp_curve(mesh, phi, z=float(np.mean(mesh.nodes[:, 2])),
                             smooth_passes=1, m_inf=m_state)
    upper = shock_report(curve, m_state)["upper"]
    freeze = result.get("sigma_freeze_report") or {}
    history = write_history(leg, result)
    write_cp(leg, curve)
    state = classify(result["residual_history"], result["converged"],
                     result["n_limited"], result["n_floored"])
    row = dict(
        leg=leg, mode=mode, level=level, r_far=radius,
        converged=int(bool(result["converged"])), **state,
        accept_reason=result.get("accept_reason") or "",
        residual=f12(result["residual_history"][-1]), n_newton=result["n_newton"],
        n_points=len(result["residual_history"]), n_limited=result["n_limited"],
        n_floored=result["n_floored"], n_gmres_total=result["n_gmres_total"],
        n_gmres_stalled=result["n_gmres_stalled"], n_refactor=result["n_refactor"],
        froze=int(bool(result["froze"])), n_freeze_refresh=result["n_freeze_refresh"],
        n_freeze_reverts=result["n_freeze_reverts"],
        sigma_frozen=int(bool(freeze.get("frozen_in_transient", False))),
        sigma_churn=int(bool(freeze.get("selection_churn", False))),
        last_sigma_delta=f12(freeze.get("last_sigma_delta")),
        tail_n_shock=freeze.get("tail_n_shock") or "", cl=f12(force["cl"]),
        cd=f12(force["cd_pressure"]),
        x_shock=(f12(upper["x_shock"]) if upper.get("has_shock") else ""),
        mach_max=f12(np.sqrt(result["mach2_max"])),
        m_final=f12(result.get("m_final", M_INF)),
        target_reached=int(bool(result.get("target_reached", result["converged"]))),
        wall_s=f12(wall), cumulative_solve_s=f12(solve_seconds + wall),
        n_threads=os.environ.get("NUMBA_NUM_THREADS", "unset"),
        blas_threads=blas_thread_count(), mesh_sha256=sha256(mesh_path),
        history_csv=history)
    append_row(HERE / "summary.csv", SUMMARY_FIELDS, row)
    if mode == "continuation":
        for i, level_result in enumerate(result["level_results"]):
            append_row(HERE / "probe_levels.csv",
                       ("leg", "level_index", "mach", "converged", "accept_reason",
                        "n_newton", "residual", "n_lin_iters", "n_refactor", "wall_s"),
                       dict(leg=leg, level_index=i, mach=f12(level_result["m"]),
                            converged=int(level_result["converged"]),
                            accept_reason=level_result["accept_reason"] or "",
                            n_newton=level_result["n_newton"],
                            residual=f12(level_result["residual_history"][-1]),
                            n_lin_iters=level_result["n_lin_iters"],
                            n_refactor=level_result["n_refactor"],
                            wall_s=f12(level_result["wall_s"])))
    print(f"[done] {leg}: {row['classification']}, converged={row['converged']}, "
          f"|R|={row['residual']}, cl={row['cl']}, wall={wall:.1f}s", flush=True)
    return wall


def write_environment():
    from pyfp3d import blas_thread_count

    versions = {name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "numba", "llvmlite", "gmsh")}
    data = dict(python=sys.version, platform=platform.platform(),
                processor=platform.processor(), packages=versions,
                numba_threads=os.environ.get("NUMBA_NUM_THREADS"),
                omp_threads=os.environ.get("OMP_NUM_THREADS"),
                openblas_threads_env=os.environ.get("OPENBLAS_NUM_THREADS"),
                blas_threads_actual=blas_thread_count(),
                python_no_user_site=os.environ.get("PYTHONNOUSERSITE"),
                cache_policy="existing warm cache; one process for all legs",
                recipe={k: str(v) for k, v in RECIPE.items()}, mach=M_INF,
                alpha_deg=ALPHA, time_budget_s=TIME_BUDGET_S)
    with open(HERE / "environment.json", "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
        fh.write("\n")


def plot_results():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = read_rows("summary.csv")
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    for level, ax in zip(LEVELS, axes):
        for row in rows:
            if row["level"] != level:
                continue
            hist = read_rows(row["history_csv"])
            residual = [float(r["residual"]) for r in hist]
            style = "-" if row["mode"] == "baseline" else "--"
            ax.semilogy(residual, style, lw=1.0,
                        label=f"{row['mode']} r{row['r_far']}: {row['classification']}")
        ax.set_title(f"{level} residual history")
        ax.set_xlabel("Newton state index")
        ax.set_ylabel("|R|_inf")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "residual_histories.png", dpi=150)
    plt.close(fig)

    curves = {}
    for row in read_rows("cp_curves.csv"):
        curves.setdefault(row["leg"], {}).setdefault(row["side"], []).append(
            (float(row["x"]), float(row["cp"])))
    medium = [r for r in rows if r["level"] == "medium"]
    fig, axes = plt.subplots(1, max(1, len(RADII)), figsize=(16, 4.8), sharey=True)
    for radius, ax in zip(RADII, np.atleast_1d(axes)):
        for row in medium:
            if int(row["r_far"]) != radius:
                continue
            for side, line in (("upper", "-"), ("lower", "--")):
                points = np.asarray(curves[row["leg"]][side])
                label = (f"{row['mode']} {side}: {row['classification']}"
                         f" (M={float(row['m_final']):.4f})")
                ax.plot(points[:, 0], points[:, 1], line, lw=1.0, label=label)
        ax.invert_yaxis()
        ax.set_title(f"medium r_far={radius}")
        ax.set_xlabel("x/c")
        ax.grid(alpha=0.3)
    axes.flat[0].set_ylabel("Cp")
    by_label = {}
    for ax in axes.flat:
        handles, labels = ax.get_legend_handles_labels()
        by_label.update(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), fontsize=7,
               loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(HERE / "cp_medium.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.analyze_only:
        plot_results()
        return

    self_test()
    write_environment()
    generate_matrix()
    if args.generate_only:
        return

    previous = read_rows("summary.csv")
    solve_seconds = max((float(r["cumulative_solve_s"]) for r in previous),
                        default=0.0)
    completed = {r["leg"] for r in previous}
    top_legs = len(completed)
    invoked = 0
    for level in LEVELS:
        for radius in RADII:
            mode = "baseline"
            leg = f"{mode}_{level}_r{radius}"
            if leg in completed:
                continue
            if solve_seconds >= TIME_BUDGET_S or top_legs >= 10:
                break
            solve_seconds += solve_leg(mode, level, radius, solve_seconds)
            top_legs += 1
            invoked += 1

    if not args.baseline_only:
        baseline = {(r["level"], int(r["r_far"])): r
                    for r in read_rows("summary.csv") if r["mode"] == "baseline"}
        for radius in (20, 30):
            base = baseline.get(("medium", radius))
            if base is None or base["converged"] == "1":
                continue
            for mode in ("continuation", "entropy_off"):
                leg = f"{mode}_medium_r{radius}"
                if leg in completed:
                    continue
                if solve_seconds >= TIME_BUDGET_S or top_legs >= 10:
                    break
                solve_seconds += solve_leg(mode, "medium", radius, solve_seconds)
                top_legs += 1
                invoked += 1

    with open(HERE / "run_state.json", "w") as fh:
        json.dump(dict(top_legs_total=top_legs,
                       top_legs_this_invocation=invoked,
                       cumulative_solver_seconds=solve_seconds,
                       time_budget_exhausted=solve_seconds >= TIME_BUDGET_S,
                       baseline_only=args.baseline_only), fh, indent=2)
        fh.write("\n")
    plot_results()


if __name__ == "__main__":
    main()
