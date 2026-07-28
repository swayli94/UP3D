"""GS1.1 follow-up: the two machine-zero states of the SAME nozzle BVP.

run_nozzle.py found that starting the Newton from a shock at x = 10 (with the
Dirichlet data of the x_s = 12 exact solution) converges to ||R|| ~ 1e-15 with
the shock near x = 8, while starting from the exact solution converges to
x ~ 11.93. The continuous BVP is provably unique (verify_uniqueness). Before
calling that a spurious discrete solution, rule out the alternatives:

  (a) a clamp is active -- rho_tilde hitting `rho_floor` turns the operator
      into a non-physical algebraic clamp that can manufacture solutions;
  (b) the second state is not choked / not quasi-1-D (e.g. two supersonic
      zones, or a transverse structure the 1-D uniqueness argument misses);
  (c) the two states carry DIFFERENT discrete mass flux, which the
      Dirichlet-Dirichlet boundary data does not pin.

Prints, for both states: the conserved discrete flux rho_tilde*u*H(x)/mdot
(the correct normalisation -- the constant-area version used by
duct.mass_flux_profile is wrong for a nozzle), the Mach profile, the number of
supersonic bins, the throat Mach, floor/limiter counters, and the imposed
Delta_phi.

Outputs: results/gs1_1_two_states.csv
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
from pyfp3d.physics.isentropic import mach_number_squared     # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
M_INF, C, NX, NY = 0.80, 1.5, 200, 12


def profile(sysd, phi, ex, n_bins):
    """Bin-averaged Mach and CONSERVED flux rho_tilde*u_x*H(x)/mdot."""
    grad, q2, rho, rho_t = sysd.state(phi)
    xc = sysd.nodes[sysd.elements].mean(axis=1)[:, 0]
    m2 = mach_number_squared(q2, M_INF)
    flux = rho_t * grad[:, 0] * N.height(xc)
    edges = np.linspace(0.0, N.LENGTH, n_bins + 1)
    idx = np.clip(np.digitize(xc, edges) - 1, 0, n_bins - 1)
    mb = np.full(n_bins, np.nan)
    fb = np.full(n_bins, np.nan)
    for b in range(n_bins):
        s = idx == b
        if s.any():
            mb[b] = np.sqrt(m2[s].mean())
            fb[b] = flux[s].mean() / ex["mdot"]
    xb = 0.5 * (edges[:-1] + edges[1:])
    return xb, mb, fb, rho_t


def main():
    ex = N.exact_solution(M_INF)
    mesh = N.nozzle_mesh(NX, NY, jitter=0.0)
    phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
    rows = []
    print(f"exact: mdot {ex['mdot']:.6f}  delta_phi {ex['delta_phi']:.6f}  "
          f"u* {ex['u_star']:.6f}")

    for tag, x_init in (("start@12(exact)", 12.0), ("start@10", 10.0)):
        sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C)
        if x_init == 12.0:
            phi0 = phi_bc.copy()
        else:
            e2 = N.exact_solution(M_INF, x_s=x_init)
            phi0 = e2["phi_of_x"](mesh.nodes[:, 0])
            phi0[sysd.dir_nodes] = phi_bc[sysd.dir_nodes]
        phi, info = sysd.newton(phi0, n_max=80, tol=1e-11)
        xb, mb, fb, rho_t = profile(sysd, phi, ex, NX)
        ok = ~np.isnan(mb)
        sup = mb[ok] > 1.0
        # count the supersonic RUNS (a second pocket would show up here)
        runs = int(np.count_nonzero(sup[1:] & ~sup[:-1])) + int(sup[0])
        i_t = int(np.argmin(np.abs(xb - N.X_T)))
        x_sh, n_sup, _, _ = N.shock_from_profile(
            *D.element_u(sysd, phi), ex["u_star"], NX)
        dphi = float(phi[np.argmax(mesh.nodes[:, 0])]
                     - phi[np.argmin(mesh.nodes[:, 0])])
        row = dict(
            state=tag, converged=info["converged"],
            n_newton=info["n_newton"], res_final=info["residual_history"][-1],
            x_shock=round(x_sh, 4), n_supersonic_bins=n_sup,
            n_supersonic_runs=runs,
            mach_max=round(float(np.nanmax(mb)), 4),
            mach_at_throat=round(float(mb[i_t]), 4),
            mach_inlet=round(float(mb[ok][0]), 4),
            mach_exit=round(float(mb[ok][-1]), 4),
            flux_min=round(float(np.nanmin(fb)), 5),
            flux_max=round(float(np.nanmax(fb)), 5),
            flux_mean_smooth=round(float(np.nanmean(
                fb[ok][(mb[ok] < 0.97) | (mb[ok] > 1.03)])), 5),
            n_floored=int(sysd.upw.n_floored),
            nu_max=round(float(sysd.upw.nu_max), 4),
            rho_tilde_min=round(float(rho_t.min()), 5),
            delta_phi_imposed=round(dphi, 5))
        rows.append(row)
        print(f"\n--- {tag} ---")
        for k, v in row.items():
            print(f"   {k:22s} {v}")
        # a compact profile dump around the shock
        s = slice(max(0, int(NX * 0.30)), int(NX * 0.70))
        print("   x        M       flux/mdot")
        for x, m, f in list(zip(xb[s], mb[s], fb[s]))[::8]:
            print(f"   {x:6.2f}  {m:6.3f}   {f:8.5f}")

    with open(OUT / "gs1_1_two_states.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1_1_two_states.csv")


if __name__ == "__main__":
    main()
