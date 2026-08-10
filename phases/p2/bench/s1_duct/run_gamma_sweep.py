"""GS1.7: sweep the circulation on each mesh level and separate the FIELD map
from the CLOSURE root.

Protocol pre-registered in docs/dev_phase_two/20260729-0010-s1-gamma-sweep.md.

For every mesh level and every Gamma in a per-level window, the FIELD is solved
self-consistently with the Kutta row removed (so no compatibility assumption is
made -- the failure mode of GS1.6's fixed-Gamma legs, where a Gamma taken from
another level's equilibrium drove the trailing-edge flow into the limiter). Each
state then yields:

  * the field map        Gamma -> x_shock, cl_p, M_max
  * BOTH closure residuals at that state:
        F_probe(Gamma) = T(phi) - Gamma          (T = TE potential jump)
        F_press(Gamma) = |q_u|^2 - |q_l|^2       (raw, sigma-free)
  * the clamp counters (a clamped point is recorded and DISCARDED from the fits)

Readouts:
  b = dT/dGamma from the sweep slope, compared with the matrix value from GS1.3b
  (0.933 / 0.955 / 0.969); the root of each closure; and whether the root shift
  between levels equals dF/(dF/dGamma) for BOTH renderings -- which would mean
  the root drift is a property of the trailing-edge flow state, not of the
  closure algebra.

Outputs: results/gs1_7_gamma_sweep.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.constraints.wake import kutta_targets              # noqa: E402
from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.physics.isentropic import mach_squared_field       # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients        # noqa: E402
from pyfp3d.solve.newton import (NewtonWorkspace,              # noqa: E402
                                 solve_newton_lifting)

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
UPWIND_C, M_CRIT, M_CAP, RHO_FLOOR = 1.5, 0.95, 3.0, 0.05
LEVELS = ("coarse", "medium", "fine")

#: per-condition Gamma windows. Chosen to OVERLAP across levels (criterion 0)
#: while staying near the levels' own equilibria: the equilibria are
#: 0.1892/0.2643/0.3575 at M0.7875 and 0.1218/0.1267/0.1278 at M0.72.
SWEEPS = {
    0.7875: np.array([0.19, 0.22, 0.25, 0.28, 0.31, 0.34]),
    0.72: np.array([0.110, 0.118, 0.1267, 0.135, 0.143]),
}
#: the fine level is ~70 s per point, so it samples a subset
FINE_SUBSET = {0.7875: (0.22, 0.28, 0.34), 0.72: (0.118, 0.1267, 0.135)}


def field_newton(ws, phi_free, gamma, n_max=60, tol=1e-9):
    R, F, state = ws.eval_residual(phi_free, gamma, UPWIND_C, M_CRIT, M_CAP,
                                   RHO_FLOOR)
    hist = [float(np.max(np.abs(R)))]
    for _ in range(n_max):
        J_ff, _ = ws.assemble_coupled(state, UPWIND_C, M_CRIT, RHO_FLOOR)
        try:
            d = spla.spsolve(J_ff.tocsc(), -R)
        except Exception:                                      # noqa: BLE001
            return phi_free, state, hist, "linear solve failed"
        if not np.all(np.isfinite(d)):
            return phi_free, state, hist, "non-finite step"
        lam, best, r0 = 1.0, None, hist[-1]
        for _ in range(10):
            trial = phi_free + lam * d
            Rt, Ft, st = ws.eval_residual(trial, gamma, UPWIND_C, M_CRIT,
                                          M_CAP, RHO_FLOOR)
            rt = float(np.max(np.abs(Rt)))
            if np.isfinite(rt) and rt < r0:
                break
            if best is None or (np.isfinite(rt) and rt < best[0]):
                best = (rt, trial, Rt, Ft, st)
            lam *= 0.5
        else:
            if best is None:
                return phi_free, state, hist, "line search"
            rt, trial, Rt, Ft, st = best
        phi_free, R, F, state = trial, Rt, Ft, st
        hist.append(rt)
        if rt < tol:
            return phi_free, state, hist, "tol"
    return phi_free, state, hist, "cap"


def main():
    rows = []
    for m_inf, gammas in SWEEPS.items():
        print(f"\n########## M {m_inf} ##########", flush=True)
        for level in LEVELS:
            path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
            if not path.exists():
                continue
            mc, wc = cut_wake(read_mesh(path))
            dz = float(np.ptp(mc.nodes[:, 2]))
            ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA,
                                 kutta_estimator="pressure")
            ws.set_mach(m_inf)
            # Seed from this level's OWN coupled equilibrium and walk Gamma
            # outward from it. A Gamma taken from another level's equilibrium
            # drives the TE flow into the limiter (GS1.6 4.3), and a 5-step
            # Picard seed is not good enough at M0.7875 on the finer meshes
            # (measured: every medium point capped at M_max 3.0).
            eq = solve_newton_lifting(
                mc, wc, m_inf=m_inf, alpha_deg=ALPHA, upwind_c=UPWIND_C,
                m_crit=M_CRIT, freeze_tol=1e-6, freeze_refresh_max=8,
                precond="direct", direct_refactor_every=4, n_newton_max=80)
            g_eq = float(np.atleast_1d(eq["gamma"])[0])
            phi_seed = np.asarray(eq["phi"],
                                  dtype=np.float64)[:ws.n_red][ws.free].copy()
            todo = (gammas if level != "fine"
                    else np.array(FINE_SUBSET[m_inf]))
            todo = sorted(todo, key=lambda x: abs(x - g_eq))
            print(f"--- {level} ({len(mc.nodes)} nodes), {len(todo)} points, "
                  f"equilibrium Gamma {g_eq:.6f} (conv={eq['converged']}) ---",
                  flush=True)
            phi_free = phi_seed.copy()
            for g in todo:
                gamma = np.array([g], dtype=np.float64)
                t0 = time.perf_counter()
                # continue from the previous Gamma's field: cheaper and keeps
                # the sweep on one branch
                phi_free, state, hist, why = field_newton(ws, phi_free, gamma)
                wall = time.perf_counter() - t0
                phi_cut = state["phi_cut"]
                # both closure residuals at this state
                f_probe = float(kutta_targets(phi_cut, wc)[0] - g)
                f_press = float(np.atleast_1d(
                    ws.cvs.residual_stations(phi_cut))[0])
                rep = shock_report(wall_cp_curve(mc, phi_cut, z=0.5 * dz,
                                                m_inf=m_inf), m_inf)
                fo = wall_force_coefficients(
                    mc.nodes, mc.elements, mc.boundary_faces["wall"], phi_cut,
                    alpha_deg=ALPHA, s_ref=dz, m_inf=m_inf)
                usable = (state["n_limited"] == 0 and state["n_floored"] == 0
                          and hist[-1] < 1e-5)
                row = dict(m_inf=m_inf, level=level, n_dof=len(mc.nodes),
                           gamma=round(float(g), 6), reason=why,
                           res_final=hist[-1], usable=usable,
                           T_probe=round(f_probe + g, 8),
                           F_probe=round(f_probe, 8),
                           F_press=round(f_press, 8),
                           cl_p=round(fo["cl"], 6),
                           x_shock=rep["upper"].get("x_shock"),
                           m_max=round(float(np.sqrt(np.max(
                               mach_squared_field(state["q2l"], m_inf)))), 5),
                           n_limited=int(state["n_limited"]),
                           n_floored=int(state["n_floored"]),
                           wall_s=round(wall, 1))
                print(f"   G={g:.4f} {why:8s} |R|={hist[-1]:.1e} "
                      f"use={str(usable):5s} T={row['T_probe']:.6f} "
                      f"F_pr={f_probe:+.6f} F_ps={f_press:+.6f} "
                      f"x_sh={row['x_shock']} M_max={row['m_max']} "
                      f"({wall:.0f}s)", flush=True)
                rows.append(row)

    with open(OUT / "gs1_7_gamma_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_7_gamma_sweep.csv")

    # ---------------- readouts ----------------
    print("\n=== b = dT/dGamma from the sweep (matrix values 0.933/0.955/0.969"
          " at M0.7875) ===")
    for m_inf in SWEEPS:
        for level in LEVELS:
            sub = [r for r in rows if r["m_inf"] == m_inf
                   and r["level"] == level and r["usable"]]
            if len(sub) < 2:
                print(f"  M{m_inf} {level:7s}: {len(sub)} usable point(s)")
                continue
            g = np.array([r["gamma"] for r in sub])
            T = np.array([r["T_probe"] for r in sub])
            fp = np.array([r["F_press"] for r in sub])
            b = np.polyfit(g, T, 1)[0]
            dfp = np.polyfit(g, fp, 1)[0]
            root_probe = np.polyfit(g, np.array([r["F_probe"] for r in sub]),
                                    1)
            root_press = np.polyfit(g, fp, 1)
            print(f"  M{m_inf} {level:7s}: n={len(sub)}  b={b:.4f}  "
                  f"1-b={1 - b:+.4f}  dF_press/dG={dfp:.4f}  "
                  f"root_probe={-root_probe[1] / root_probe[0]:.5f}  "
                  f"root_press={-root_press[1] / root_press[0]:.5f}")

    print("\n=== field map: x_shock at common Gamma ===")
    for m_inf in SWEEPS:
        common = sorted(set.intersection(*[
            {r["gamma"] for r in rows if r["m_inf"] == m_inf
             and r["level"] == lv and r["usable"]}
            for lv in LEVELS
            if any(r["m_inf"] == m_inf and r["level"] == lv for r in rows)]
        )) if rows else []
        print(f"  M{m_inf}: common usable Gamma = {common}")
        for g in common:
            vals = []
            for lv in LEVELS:
                m = [r for r in rows if r["m_inf"] == m_inf
                     and r["level"] == lv and r["gamma"] == g and r["usable"]]
                vals.append(m[0]["x_shock"] if m else None)
            spread = (max(v for v in vals if v) - min(v for v in vals if v)
                      if all(v for v in vals) else float("nan"))
            print(f"    G={g:.4f}  x_shock " + " / ".join(
                f"{v:.4f}" if v else "--" for v in vals)
                + f"   spread {spread:.4f}")


if __name__ == "__main__":
    main()
