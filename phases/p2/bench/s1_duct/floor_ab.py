"""GS1.1 decisive A/B: is the spurious second solution MANUFACTURED by the
artificial-density floor?

diagnose_two_states.py measured, on one nozzle BVP with identical Dirichlet
data, two states with ||R|| ~ 1e-12:

    start @ 12 (exact)  -> x_shock 11.93, n_floored  0, rho_tilde_min 0.670
    start @ 10          -> x_shock  7.98, n_floored 72, rho_tilde_min 0.050  <- == rho_floor

`rho_floor` is a hard clamp: where it binds, the element no longer discretises
rho(q^2) -- it discretises the constant 0.05, a DIFFERENT (degenerate) equation.
A clamped region can therefore host a state that is a genuine root of the
discrete residual while not being a solution of the full-potential equation.

This script sweeps `rho_floor` over four decades from the same two starts. If
the floor is the cause, the spurious branch must disappear (or fail to
converge) as the floor is lowered, while the correct branch -- which never
touches the floor -- must be unaffected.

Outputs: results/gs1_1_floor_ab.csv
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

import duct as D                                              # noqa: E402
import nozzle as N                                            # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
M_INF, C, NX, NY = 0.80, 1.5, 200, 12
FLOORS = (0.05, 0.02, 5e-3, 1e-4, 1e-8)
STARTS = (("exact@12", 12.0), ("perturbed@10", 10.0))


def main():
    ex = N.exact_solution(M_INF)
    mesh = N.nozzle_mesh(NX, NY, jitter=0.0)
    phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
    rows = []
    print(f"exact x_s = {ex['x_s']}, mdot {ex['mdot']:.6f}, "
          f"shipped default rho_floor = 0.05")
    for tag, x0 in STARTS:
        for floor in FLOORS:
            sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C,
                                rho_floor=floor)
            if x0 == 12.0:
                phi0 = phi_bc.copy()
            else:
                e2 = N.exact_solution(M_INF, x_s=x0)
                phi0 = e2["phi_of_x"](mesh.nodes[:, 0])
                phi0[sysd.dir_nodes] = phi_bc[sysd.dir_nodes]
            phi, info = sysd.newton(phi0, n_max=80, tol=1e-11)
            x_sh, n_sup, _, _ = N.shock_from_profile(
                *D.element_u(sysd, phi), ex["u_star"], NX)
            _, _, _, rho_t = sysd.state(phi)
            row = dict(start=tag, rho_floor=floor,
                       converged=info["converged"], reason=info["reason"],
                       n_newton=info["n_newton"],
                       res_final=info["residual_history"][-1],
                       x_shock=round(x_sh, 4) if np.isfinite(x_sh) else None,
                       err_x=round(x_sh - ex["x_s"], 4)
                       if np.isfinite(x_sh) else None,
                       n_floored=int(sysd.upw.n_floored),
                       rho_tilde_min=round(float(rho_t.min()), 6),
                       n_supersonic_bins=n_sup)
            rows.append(row)
            print(f"  {tag:13s} floor {floor:<8g} conv={str(row['converged']):5s} "
                  f"({row['reason']:14s}) |R| {row['res_final']:.2e}  "
                  f"x_shock {row['x_shock']}  n_floored {row['n_floored']:4d}  "
                  f"rho_t_min {row['rho_tilde_min']}", flush=True)

    with open(OUT / "gs1_1_floor_ab.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_1_floor_ab.csv")

    print("\n=== reading ===")
    for tag, _ in STARTS:
        sub = [r for r in rows if r["start"] == tag]
        xs = [r["x_shock"] for r in sub]
        print(f"  {tag:13s} x_shock over floors {FLOORS}: {xs}")


if __name__ == "__main__":
    main()
