"""GS1.2: quantify what controls the shock position, and whether LIFT converges
under mesh refinement.

Two parts (pre-registered in
docs/dev_phase_two/20260728-1825-s1-resolution.md):

  Q1  nozzle bench (exact reference): position error and smear width expressed
      in CELLS, vs h and C. If the error is a roughly constant number of cells,
      h is the only lever the present operator offers.
  Q2  airfoil at the M1 condition: the same case on h_wall = 0.02 / 0.01 /
      0.005 chord (coarse / medium / fine). Does cl converge? The audit already
      measured coarse -> medium = +40 % at M0.7875; if medium -> fine moves by
      the same amount, refinement is not a route to metric M1.

The fine leg uses precond="amg" (audit exp2 measured amg and direct give the
identical solution on M6 medium, only the wall clock differs), everything else
is the recipe used by the audit and by tests/test_p8_newton.py.

Outputs: results/gs1_2_nozzle_cells.csv, results/gs1_2_airfoil_h.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

import duct as D                                              # noqa: E402
import nozzle as N                                            # noqa: E402
from pyfp3d.mesh.reader import read_mesh                      # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                     # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve             # noqa: E402
from pyfp3d.post.shock import shock_report                    # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients       # noqa: E402
from pyfp3d.solve.newton import solve_newton_transonic        # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

# ---- Q1: nozzle ----------------------------------------------------------
M_INF_N = 0.80
NXS = (100, 200, 400, 800)
CS_N = (1.0, 1.5, 3.0)

# ---- Q2: airfoil ---------------------------------------------------------
ALPHA = 1.25
CASES_A = [(0.7875, 1.5), (0.80, 1.5), (0.7875, 3.0)]
LEVELS = ("coarse", "medium", "fine")
H_WALL = {"coarse": 0.02, "medium": 0.01, "fine": 0.005}


def smear_width(sysd, phi, ex, n_bins, hx):
    """Transition width in CELLS: bins whose u lies strictly between the two
    exact branch values (5 % .. 95 % of the jump)."""
    xc, ux = D.element_u(sysd, phi)
    edges = np.linspace(0.0, N.LENGTH, n_bins + 1)
    idx = np.clip(np.digitize(xc, edges) - 1, 0, n_bins - 1)
    ub = np.full(n_bins, np.nan)
    for b in range(n_bins):
        m = idx == b
        if m.any():
            ub[b] = ux[m].mean()
    ok = ~np.isnan(ub)
    xb = (0.5 * (edges[:-1] + edges[1:]))[ok]
    ub = ub[ok]
    # local branch values at each bin, from the exact solution
    u_sup = np.interp(xb, ex["x"], ex["u"])          # exact profile (has jump)
    x_sh, _, _, _ = N.shock_from_profile(xc, ux, ex["u_star"], n_bins)
    if not np.isfinite(x_sh):
        return -1
    # use the exact jump size local to the shock
    i = int(np.argmin(np.abs(ex["x"] - ex["x_s"])))
    hi = float(np.max(ex["u"][max(0, i - 200):i]))
    lo = float(np.min(ex["u"][i:i + 200]))
    band = (ub < hi - 0.05 * (hi - lo)) & (ub > lo + 0.05 * (hi - lo))
    near = np.abs(xb - x_sh) < 10 * hx
    return int(np.count_nonzero(band & near))


def q1_nozzle(rows):
    ex = N.exact_solution(M_INF_N)
    print(f"\n=== Q1 nozzle (exact x_s = {ex['x_s']}, "
          f"M_shock {ex['m_shock_up']:.3f}) ===", flush=True)
    for nx in NXS:
        ny = max(6, nx // 16)
        hx = N.LENGTH / nx
        mesh = N.nozzle_mesh(nx, ny, jitter=0.0)
        phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
        for C in CS_N:
            sysd = D.DuctSystem(mesh, m_inf=M_INF_N, upwind_c=C)
            t0 = time.perf_counter()
            phi, info = sysd.newton(phi_bc.copy(), n_max=80, tol=1e-11)
            wall = time.perf_counter() - t0
            x_sh, _, _, _ = N.shock_from_profile(
                *D.element_u(sysd, phi), ex["u_star"], nx)
            w = smear_width(sysd, phi, ex, nx, hx)
            err = x_sh - ex["x_s"] if np.isfinite(x_sh) else None
            row = dict(part="nozzle", nx=nx, h=round(hx, 5), C=C,
                       n_dof=len(mesh.nodes),
                       converged=info["converged"], reason=info["reason"],
                       res_final=info["residual_history"][-1],
                       wall_s=round(wall, 1),
                       x_shock=round(x_sh, 5) if err is not None else None,
                       err_x=round(err, 5) if err is not None else None,
                       err_cells=round(err / hx, 3) if err is not None else None,
                       smear_cells=w)
            print("  ", {k: row[k] for k in ("nx", "h", "C", "converged",
                                             "err_x", "err_cells",
                                             "smear_cells")}, flush=True)
            rows.append(row)


def q2_airfoil(rows):
    print("\n=== Q2 airfoil: does cl converge with h? ===", flush=True)
    for m_inf, C in CASES_A:
        for level in LEVELS:
            path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
            if not path.exists():
                print(f"  skip {level}: mesh missing")
                continue
            mc, wc = cut_wake(read_mesh(path))
            dz = float(np.ptp(mc.nodes[:, 2]))
            precond = "amg" if level == "fine" else "direct"
            recipe = dict(m_start=0.70, dm=0.025, dm_min=0.003,
                          freeze_tol=1e-6, upwind_c=C,
                          newton_kw=dict(freeze_refresh_max=8,
                                         precond=precond, n_newton_max=60))
            t0 = time.perf_counter()
            try:
                r = solve_newton_transonic(mc, wc, m_inf=m_inf,
                                           alpha_deg=ALPHA, **recipe)
                err = ""
            except Exception as exc:                          # noqa: BLE001
                rows.append(dict(part="airfoil", m_inf=m_inf, C=C, level=level,
                                 n_dof=len(mc.nodes), h_wall=H_WALL[level],
                                 converged=False, error=type(exc).__name__))
                print("  ", rows[-1], flush=True)
                continue
            wall = time.perf_counter() - t0
            rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                            m_inf=m_inf), m_inf)
            f = wall_force_coefficients(mc.nodes, mc.elements,
                                        mc.boundary_faces["wall"], r["phi"],
                                        alpha_deg=ALPHA, s_ref=dz, m_inf=m_inf)
            done = [lr["m"] for lr in r["level_results"] if lr["converged"]]
            row = dict(part="airfoil", m_inf=m_inf, C=C, level=level,
                       n_dof=len(mc.nodes), h_wall=H_WALL[level],
                       precond=precond, converged=r["converged"],
                       m_reached=max(done) if done else None,
                       wall_s=round(wall, 1),
                       cl_p=round(f["cl"], 6),
                       cl_kj=round(2 * float(r["gamma"][0]), 6),
                       x_shock_up=rep["upper"].get("x_shock"),
                       n_cells_up=rep["upper"].get("n_cells"),
                       m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                       n_limited=r["n_limited"], n_floored=r["n_floored"],
                       res_final=r["residual_history"][-1], error=err)
            print("  ", {k: row[k] for k in
                         ("m_inf", "C", "level", "converged", "m_reached",
                          "cl_p", "x_shock_up", "n_cells_up", "m_max",
                          "n_floored", "wall_s")}, flush=True)
            rows.append(row)


def report(rows):
    print("\n=== Q1 headline: error and smear in CELLS ===")
    print(f"  {'nx':>4s} {'h':>7s} " + " ".join(
        f"{'C=' + str(c) + ' err/smear':>18s}" for c in CS_N))
    for nx in NXS:
        cells = []
        for C in CS_N:
            m = [r for r in rows if r.get("part") == "nozzle"
                 and r["nx"] == nx and r["C"] == C]
            if m and m[0]["err_cells"] is not None:
                cells.append(f"{m[0]['err_cells']:10.2f} /{m[0]['smear_cells']:6d}")
            else:
                cells.append(f"{'--':>18s}")
        h = N.LENGTH / nx
        print(f"  {nx:4d} {h:7.4f} " + " ".join(cells))

    print("\n=== Q2 headline: cl vs h_wall ===")
    for m_inf, C in CASES_A:
        sub = [r for r in rows if r.get("part") == "airfoil"
               and r.get("m_inf") == m_inf and r.get("C") == C]
        if not sub:
            continue
        print(f"  M {m_inf}, C {C}")
        prev = None
        for level in LEVELS:
            m = [r for r in sub if r["level"] == level]
            if not m:
                continue
            r0 = m[0]
            cl = r0.get("cl_p")
            d = (f"{100 * (cl - prev) / abs(prev):+7.1f} %"
                 if (prev is not None and cl is not None) else "     --  ")
            print(f"    {level:7s} h_wall {r0['h_wall']:.3f}  "
                  f"conv={str(r0.get('converged')):5s} "
                  f"cl_p {cl if cl is not None else float('nan'):.4f}  "
                  f"d(cl) {d}   x_shock "
                  f"{r0.get('x_shock_up')}  cells {r0.get('n_cells_up')}  "
                  f"M_max {r0.get('m_max')}")
            if cl is not None:
                prev = cl


def main():
    rows = []
    q1_nozzle(rows)
    q2_airfoil(rows)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(OUT / "gs1_2_resolution.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_2_resolution.csv")
    report(rows)


if __name__ == "__main__":
    main()
