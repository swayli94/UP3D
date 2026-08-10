"""GS1.6 criterion 0 (the ACTION PREMISE): does the exact-solution nozzle bench
also lose h-convergence when the shock is strong?

The entropy-fix hypothesis (GS1.3b 4.3) says the lift divergence at M_max ~ 1.4
comes from the ISENTROPIC shock over-predicting as the artificial dissipation
vanishes. If that is the mechanism, then the SAME degradation must be visible in
a case with no lift and no Kutta coupling -- the Laval nozzle, where the exact
shock position is known -- once the shock is pushed to the same strength.

If instead the nozzle stays O(h)-convergent at M_shock 1.4-1.5, then isentropic
over-prediction is NOT what breaks the airfoil, and per the pre-registered kill
criterion this round stops before any entropy term is written.

Sweeps the shock strength (M_shock = 1.15, 1.30, 1.45, 1.60) over three mesh
levels at fixed C, and reports the position error and its convergence ratio --
the same statistic (Delta2/Delta1) used for the airfoil envelope.

Outputs: results/gs1_6_strong_shock.csv
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

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

M_INF = 0.80
M_SHOCKS = (1.15, 1.30, 1.45, 1.60)
NXS = (100, 200, 400)
C = 1.5


def area_ratio_for(m_shock, gamma=1.4):
    """A/A* that puts the supersonic branch at `m_shock` -- used to retune the
    nozzle so the shock at x_s = 12 has the requested strength."""
    m2 = m_shock * m_shock
    return ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0))) \
        * (1.0 + 0.5 * (gamma - 1.0) * m2) ** ((gamma + 1.0)
                                               / (2.0 * (gamma - 1.0))) / m_shock


def main():
    rows = []
    a_out_default = N.A_OUT
    for m_shock in M_SHOCKS:
        # retune the diverging wall so that H(x_s)/H_t equals the area ratio
        # of the requested supersonic Mach (the throat stays sonic, so the
        # shock at x_s = 12 then sits at m_shock)
        ar = area_ratio_for(m_shock)
        frac = ((N.X_S_TARGET - N.X_T) / (N.LENGTH - N.X_T)) ** 2
        N.A_OUT = (ar - 1.0) / frac
        ex = N.exact_solution(M_INF)
        dphi_s, dphi_sub, uniq = N.verify_uniqueness(M_INF)
        print(f"\n=== M_shock target {m_shock} -> a_out {N.A_OUT:.4f}, "
              f"exact M just upstream {ex['m_shock_up']:.4f}, "
              f"unique={uniq} ===", flush=True)
        if not uniq:
            print("   (boundary data not uniquely shocked -- skipping)")
            continue
        errs = {}
        for nx in NXS:
            ny = max(6, nx // 16)
            hx = N.LENGTH / nx
            mesh = N.nozzle_mesh(nx, ny, jitter=0.0)
            sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C)
            phi0 = ex["phi_of_x"](mesh.nodes[:, 0])
            t0 = time.perf_counter()
            phi, info = sysd.newton(phi0, n_max=80, tol=1e-11)
            wall = time.perf_counter() - t0
            x_sh, _, _, _ = N.shock_from_profile(
                *D.element_u(sysd, phi), ex["u_star"], nx)
            err = (x_sh - ex["x_s"]) if np.isfinite(x_sh) else None
            errs[nx] = err
            row = dict(m_shock_target=m_shock,
                       m_shock_exact=round(ex["m_shock_up"], 4),
                       a_out=round(N.A_OUT, 4), nx=nx, h=round(hx, 5), C=C,
                       n_dof=len(mesh.nodes), converged=info["converged"],
                       reason=info["reason"],
                       res_final=info["residual_history"][-1],
                       x_shock=round(x_sh, 5) if err is not None else None,
                       err_x=round(err, 5) if err is not None else None,
                       err_cells=round(err / hx, 3) if err is not None else None,
                       wall_s=round(wall, 1))
            print(f"   nx={nx:4d} h={hx:.4f} conv={str(info['converged']):5s} "
                  f"x_shock={row['x_shock']} err={row['err_x']} "
                  f"err_cells={row['err_cells']}", flush=True)
            rows.append(row)
        if all(errs.get(nx) is not None for nx in NXS):
            d1 = abs(errs[NXS[1]]) - abs(errs[NXS[0]])
            d2 = abs(errs[NXS[2]]) - abs(errs[NXS[1]])
            ratio = d2 / d1 if d1 != 0 else float("nan")
            print(f"   |err| {abs(errs[NXS[0]]):.5f} -> "
                  f"{abs(errs[NXS[1]]):.5f} -> {abs(errs[NXS[2]]):.5f}   "
                  f"ratio {ratio:.3f}", flush=True)
    N.A_OUT = a_out_default

    with open(OUT / "gs1_6_strong_shock.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_6_strong_shock.csv")

    print("\n=== headline: position error (in cells) vs shock strength ===")
    print(f"  {'M_shock':>8s} " + " ".join(f"{'nx=' + str(n):>10s}"
                                           for n in NXS) + "   |err| ratio")
    for m_shock in M_SHOCKS:
        sub = {r["nx"]: r for r in rows if r["m_shock_target"] == m_shock}
        if not sub:
            continue
        cells = [f"{sub[n]['err_cells']:10.3f}"
                 if n in sub and sub[n]["err_cells"] is not None
                 else f"{'--':>10s}" for n in NXS]
        if all(n in sub and sub[n]["err_x"] is not None for n in NXS):
            a = [abs(sub[n]["err_x"]) for n in NXS]
            d1, d2 = a[1] - a[0], a[2] - a[1]
            ratio = f"{d2 / d1:.3f}" if d1 != 0 else "--"
        else:
            ratio = "--"
        print(f"  {m_shock:8.2f} " + " ".join(cells) + f"   {ratio}")


if __name__ == "__main__":
    main()
