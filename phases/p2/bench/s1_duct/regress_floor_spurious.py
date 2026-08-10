"""GS1.4 regression: the rho_floor spurious solution must stay reproducible.

GS1.1 measured that the artificial-density floor manufactures a machine-zero
solution: on the Laval nozzle BVP (unique exact answer, shock at x = 12) the
Newton started from a shock at x = 10 converges to |R| ~ 9e-13 with the shock at
x ~ 8 and 72 elements pinned at exactly rho_floor, while the physical branch is
insensitive to the floor over four decades.

This script is the standing regression for that finding -- not a pytest (it is a
~1 min bench run), but an assert-carrying script so a silent change shows up:

  * the PHYSICAL branch must stay floor-independent (same shock to 1e-4 over
    floor 0.05 .. 1e-8, n_floored == 0 throughout);
  * the SPURIOUS branch must still exist at the shipped floor 0.05 (converged,
    clamped, far from the exact position) -- if it ever stops existing, the
    operator changed and GS1.4's motivation must be re-measured.

Exit code 0 = both hold.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

import duct as D                                              # noqa: E402
import nozzle as N                                            # noqa: E402

M_INF, C, NX, NY = 0.80, 1.5, 200, 12
FLOORS = (0.05, 0.02, 5e-3, 1e-4, 1e-8)


def main():
    ex = N.exact_solution(M_INF)
    mesh = N.nozzle_mesh(NX, NY, jitter=0.0)
    phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
    ok = True

    print("physical branch (start = exact):")
    xs = []
    for floor in FLOORS:
        sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C, rho_floor=floor)
        phi, info = sysd.newton(phi_bc.copy(), n_max=80, tol=1e-11)
        x, _, _, _ = N.shock_from_profile(*D.element_u(sysd, phi),
                                          ex["u_star"], NX)
        _, _, _, rho_t = sysd.state(phi)
        n_fl = int(sysd.upw.n_floored)
        print(f"  floor {floor:<8g} conv={str(info['converged']):5s} "
              f"x_shock={x:.5f} n_floored={n_fl}")
        xs.append(x)
        if n_fl != 0:
            print("   !! the physical branch touched the floor")
            ok = False
    if max(xs) - min(xs) > 1e-4:
        print(f"  !! physical branch moved with the floor "
              f"({max(xs) - min(xs):.2e} > 1e-4)")
        ok = False
    else:
        print(f"  physical branch floor-independent "
              f"(spread {max(xs) - min(xs):.2e})")

    print("spurious branch (start = shock at x = 10, shipped floor 0.05):")
    e2 = N.exact_solution(M_INF, x_s=10.0)
    sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C, rho_floor=0.05)
    phi0 = e2["phi_of_x"](mesh.nodes[:, 0])
    phi0[sysd.dir_nodes] = phi_bc[sysd.dir_nodes]
    phi, info = sysd.newton(phi0, n_max=80, tol=1e-11)
    x, _, _, _ = N.shock_from_profile(*D.element_u(sysd, phi), ex["u_star"], NX)
    n_fl = int(sysd.upw.n_floored)
    print(f"  conv={info['converged']} |R|={info['residual_history'][-1]:.2e} "
          f"x_shock={x:.5f} (exact {ex['x_s']}) n_floored={n_fl}")
    if not (info["converged"] and n_fl > 0 and abs(x - ex["x_s"]) > 1.0):
        print("   !! the spurious clamped branch no longer reproduces -- "
              "re-measure GS1.4's motivation")
        ok = False
    else:
        print("  spurious clamped branch reproduces (this is why GS1.4 exists)")

    print("\nOK" if ok else "\nFAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
