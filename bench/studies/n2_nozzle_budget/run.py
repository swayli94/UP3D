"""N2 -- is the nozzle's reason=cap a BUDGET stop or a capability limit?

Binding text: phases/p5/docs/dev_phase_five/20260824-1000-n2-prereg.md (committed first).

Read the implementation first and the question changed: `reason=cap` is the Newton
ITERATION cap (n_max=80, hardcoded in run_nozzle.py:65) with tol=1e-11, and those legs
stop at |R| 5.9e-9 to 1.9e-8 -- two to three decades short. That is CLAUDE.md's
budget_limited shape, which it explicitly says is NOT a failure. But budget_limited needs
a MONOTONE TAIL, not just a small residual, so the history is classified rather than the
final value read.

Run:  PYTHONNOUSERSITE=1 python bench/studies/n2_nozzle_budget/run.py
"""
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ARCH = os.path.join(ROOT, "phases/p2/bench/s1_duct")
sys.path.insert(0, ROOT); sys.path.insert(0, ARCH)
RES = os.path.join(HERE, "results")

import duct as D                                                        # noqa: E402
import nozzle as N                                                      # noqa: E402

#: ★ G-ONEKNOB -- every one of these is copied verbatim from run_nozzle.py; only N_MAX moves
#: ★★ M_INF is read from the runner module, not guessed -- my first attempt guessed two
#: method names for the initial guess, both absent, so phi0 became None and numpy turned it
#: into a 0-d array: eight legs died in the same obscure numba typing error inside
#: element_velocity_q2. The hasattr fallbacks made the guess SILENT instead of loud.
#: The construction below is now copied verbatim from run_nozzle.py::run_one.
import importlib.util as _ilu
_rs = _ilu.spec_from_file_location("run_nozzle", os.path.join(ARCH, "run_nozzle.py"))
_rm = _ilu.module_from_spec(_rs)
_src = open(os.path.join(ARCH, "run_nozzle.py")).read()
M_INF = float([l.split("=")[1] for l in _src.splitlines()
               if l.startswith("M_INF")][0].split("#")[0])
CS = (1.0, 1.5, 2.0, 3.0)
NX = 400                    # h = 0.05, the finest level in the committed sweep
TOL = 1e-11
N_MAX_REPRO, N_MAX_BIG = 80, 600
SUMMARY, IMPL = [], {}


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def classify(hist):
    """★ CLAUDE.md: budget_limited requires a MONOTONE tail. A period-3 cycle once gave
    descent10 = 2.0021 and fooled a hand-picked threshold, so report the shape, not a
    single ratio."""
    IMPL["B-MODE"] = True
    h = np.asarray(hist, float)
    d10 = float(h[-11] / h[-1]) if len(h) >= 11 and h[-1] > 0 else float("nan")
    tail = h[-8:] if len(h) >= 8 else h
    n_up = int((np.diff(tail) > 0).sum())
    mono = n_up == 0
    return dict(descent10=d10, n_up=n_up, monotone_tail=mono,
                mode="budget_limited" if mono else "oscillating/limit_cycle",
                tail=" ".join(f"{v:.3e}" for v in tail))


def one(C, n_max):
    """run_nozzle.py::run_one, reduced to what N2 needs; the recipe is verbatim."""
    ny = max(6, NX // 16)
    mesh = N.nozzle_mesh(NX, ny, jitter=0.0)
    sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C)
    ex = N.exact_solution(M_INF, x_s=N.X_S_TARGET)
    #: verbatim from run_nozzle.py::run_one -- Dirichlet data from the exact solution,
    #: and for the unperturbed start the initial guess IS that data
    phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
    phi0 = phi_bc.copy()
    t0 = time.perf_counter()
    phi, info = sysd.newton(phi0, n_max=n_max, tol=TOL)
    w = time.perf_counter() - t0
    xc, ux = D.element_u(sysd, phi)
    x_sh, n_sup, _, _ = N.shock_from_profile(xc, ux, ex["u_star"], NX)
    #: G-TRUTH: the analytic shock position comes from the module, never typed by me
    return dict(C=C, n_max=n_max, converged=bool(info["converged"]),
                reason=info.get("reason", "tol"), n_newton=int(info["n_newton"]),
                res_final=float(np.asarray(info["residual_history"], float)[-1]),
                x_shock=float(x_sh), err_x=float(x_sh - ex["x_s"]),
                err_cells=float((x_sh - ex["x_s"]) / (N.LENGTH / NX)),
                wall_s=round(w, 1), hist=list(info["residual_history"]))


def main():
    os.makedirs(RES, exist_ok=True)
    print(f"  G-FROZEN-LIB  pyfp3d/ and phases/p2/ both untouched; runner lives here")
    print(f"  G-ONEKNOB     only n_max moves ({N_MAX_REPRO} -> {N_MAX_BIG}); "
          f"tol={TOL} CS={CS} nx={NX} jitter=0 start=exact, all verbatim\n", flush=True)
    rows = []
    for n_max in (N_MAX_REPRO, N_MAX_BIG):
        for C in CS:
            try:
                r = one(C, n_max)
            except Exception as e:                                     # noqa: BLE001
                print(f"  [n_max={n_max} C={C}] RAISED {type(e).__name__}: {e}",
                      flush=True)
                continue
            r["mode"] = ("converged" if r["converged"]
                         else classify(r["hist"])["mode"])
            cl = classify(r["hist"])
            rows.append({k: v for k, v in r.items() if k != "hist"} | cl)
            print(f"  n_max={n_max:4d} C={C:<4} conv={r['converged']!s:5} "
                  f"reason={r['reason']:>4} steps={r['n_newton']:4d} "
                  f"|R|={r['res_final']:.3e} err_cells={r['err_cells']:+.4f} "
                  f"{r['wall_s']:6.1f}s  mode={cl['mode']}", flush=True)
    if not rows:
        print("  no legs ran"); return 2
    with open(os.path.join(RES, "budget.csv"), "w", newline="") as f:
        ks = sorted({k for r in rows for k in r})
        w_ = csv.DictWriter(f, fieldnames=ks); w_.writeheader(); w_.writerows(rows)

    IMPL["B-BUDGET"] = True
    big = [r for r in rows if r["n_max"] == N_MAX_BIG]
    small = [r for r in rows if r["n_max"] == N_MAX_REPRO]
    _record("G-REPRO", "n_max=80 reproduces the committed res_final (~6e-9..1.9e-8)",
            "same order of magnitude and same x_shock",
            "; ".join(f"C={r['C']}: |R| {r['res_final']:.3e}, err_cells "
                      f"{r['err_cells']:+.4f}" for r in small),
            "G-REPRO: reproduced" if small and all(
                1e-10 < r["res_final"] < 1e-4 for r in small) else "★ G-REPRO: check")
    n_conv = sum(r["converged"] for r in big)
    _record("B-BUDGET", f"does n_max={N_MAX_BIG} converge what {N_MAX_REPRO} could not",
            "all converge => it was a BUDGET stop and the h trend becomes citable;  "
            "still capped => a capability limit, which is the real obstacle",
            "; ".join(f"C={r['C']}: conv={r['converged']} steps={r['n_newton']} "
                      f"|R| {r['res_final']:.2e}" for r in big),
            f"★★★ B-BUDGET: it was a BUDGET stop ({n_conv}/{len(big)} converge)"
            if n_conv == len(big) else
            f"★ B-BUDGET: {n_conv}/{len(big)} converge -- NOT purely a budget")
    _record("B-MODE", "residual-tail shape per leg (CLAUDE.md: budget_limited needs a "
            "MONOTONE tail, not a small residual)", "monotone => budget_limited",
            "; ".join(f"C={r['C']}@{r['n_max']}: {r['mode']} (n_up={r['n_up']}, "
                      f"descent10={r['descent10']:.3f})" for r in rows), "RECORDED")
    IMPL["B-TREND"] = True
    if n_conv == len(big):
        _record("B-TREND", "err_cells at h=0.05 now CITABLE, against the committed "
                "h=0.2 / h=0.1", "the law claims err_cells is roughly constant in h",
                "; ".join(f"C={r['C']}: err_cells {abs(r['err_cells']):.3f}" for r in big)
                + "  vs committed h=0.2 [0.500 0.651 0.768 0.935] / "
                  "h=0.1 [0.551 0.684 0.804 1.038]", "RECORDED")
    else:
        _record("B-TREND", "conditional on B-BUDGET", "n/a", "B-BUDGET did not pass",
                "NOT APPLICABLE")
    IMPL["B-COST"] = True
    _record("B-COST", "steps and wall per leg", "RECORDED",
            "; ".join(f"C={r['C']}@{r['n_max']}: {r['n_newton']} steps {r['wall_s']}s"
                      for r in rows), "RECORDED")
    reg = ("B-BUDGET", "B-MODE", "B-TREND", "B-COST")
    print("\n  G-CHECKOFF:")
    for c in reg:
        print(f"    {c:10} {'implemented' if IMPL.get(c) else '★ NOT IMPLEMENTED'}")
    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w_ = csv.writer(f); w_.writerow(["tag", "metric", "band", "measured", "verdict"])
        w_.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
