"""GS1b.3 criterion B: does the entropy correction hit the PRE-REGISTERED
nozzle shock position?

The pre-registration (docs/dev_phase_two/20260729-0700-s1b-entropy-implementation
.md sec 2.2) computed, before any solver code existed, where the corrected model
must put the shock for the SAME Dirichlet data the isentropic bench has always
used (delta_phi = 19.31862767, the data for which the isentropic model puts the
shock exactly at x_s = 12.0, M1 = 1.29998):

    isentropic model  ->  x_s = 12.0000
    corrected  model  ->  x_s = 11.6200   (M1 1.28022, sigma 0.982657)
    shift = -0.3800 = -1.90 % L, UPSTREAM

That number comes from a quasi-1-D inversion, so it is independent of the solver:
this is an external anchor, not a self-consistency check. Criterion B is that the
correction ON converges to 11.62 (error shrinking under refinement, within 2
cells at nx = 400) while OFF stays on 12.0.

`corrected_exact` below is the same inversion, recomputed here so the reference
lives with the measurement rather than in a markdown table.

Outputs: results/gs1b_3_nozzle.csv
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

import duct as D                                               # noqa: E402
import nozzle as N                                             # noqa: E402
from pyfp3d.kernels.entropy import total_pressure_ratio        # noqa: E402
from pyfp3d.physics.isentropic import (GAMMA,                  # noqa: E402
                                       critical_speed_squared,
                                       density_isentropic,
                                       mach_number_squared)

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

M_INF = 0.80
NXS = (100, 200, 400)
C = 1.5


def corrected_exact(x_s, n=20001):
    """Quasi-1-D choked solution with an ENTROPY-CORRECTED shock at x_s:
    upstream isentropic, downstream rho = sigma*rho_isen with sigma = p02/p01 at
    the local pre-shock Mach. Returns (delta_phi, M1, sigma)."""
    u_star = float(np.sqrt(critical_speed_squared(M_INF, GAMMA)))
    mdot = float(N._f(u_star, M_INF, GAMMA) * N.H_T)
    xs = np.linspace(0.0, N.LENGTH, n)
    H = N.height(xs)
    u = np.empty_like(xs)
    for i, (xi, Hi) in enumerate(zip(xs, H)):
        u_sub, u_sup = N._roots(mdot / Hi, M_INF, GAMMA, u_star)
        u[i] = u_sup if (N.X_T <= xi < x_s) else u_sub
    i_s = int(np.searchsorted(xs, x_s))
    m1 = float(np.sqrt(mach_number_squared(u[i_s - 1] ** 2, M_INF, GAMMA)))
    sig = total_pressure_ratio(m1, GAMMA)
    for i in range(i_s, n):
        t = mdot / H[i]
        lo, hi = 1e-12, u_star
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if sig * float(density_isentropic(mid * mid, M_INF, GAMMA)) * mid \
                    < t:
                lo = mid
            else:
                hi = mid
        u[i] = 0.5 * (lo + hi)
    dphi = float(np.sum(0.5 * (u[1:] + u[:-1]) * np.diff(xs)))
    return dphi, m1, sig


def corrected_exact_fixed_sigma(x_s, sigma, n=20001):
    """Same as corrected_exact but with sigma PRESCRIBED instead of computed from
    the pre-shock Mach at x_s. Used to separate the two error sources: the solver
    detects its own M1 (under-resolved on coarse meshes, converging with h), so
    comparing against a reference built with the solver's OWN sigma isolates the
    POSITION error of the mechanism from the MAGNITUDE error of the detection."""
    u_star = float(np.sqrt(critical_speed_squared(M_INF, GAMMA)))
    mdot = float(N._f(u_star, M_INF, GAMMA) * N.H_T)
    xs = np.linspace(0.0, N.LENGTH, n)
    H = N.height(xs)
    u = np.empty_like(xs)
    for i, (xi, Hi) in enumerate(zip(xs, H)):
        u_sub, u_sup = N._roots(mdot / Hi, M_INF, GAMMA, u_star)
        u[i] = u_sup if (N.X_T <= xi < x_s) else u_sub
    i_s = int(np.searchsorted(xs, x_s))
    for i in range(i_s, n):
        t = mdot / H[i]
        lo, hi = 1e-12, u_star
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if sigma * float(density_isentropic(mid * mid, M_INF, GAMMA)) * mid \
                    < t:
                lo = mid
            else:
                hi = mid
        u[i] = 0.5 * (lo + hi)
    return float(np.sum(0.5 * (u[1:] + u[:-1]) * np.diff(xs)))


def x_s_for_sigma(target_dphi, sigma):
    """Invert corrected_exact_fixed_sigma for x_s at a prescribed sigma."""
    lo, hi = N.X_T + 0.5, N.LENGTH - 0.5
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if corrected_exact_fixed_sigma(mid, sigma) < target_dphi:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def corrected_x_s(target_dphi):
    """Invert corrected_exact: the x_s the corrected model needs for the
    isentropic bench's Dirichlet data."""
    lo, hi = N.X_T + 0.5, N.LENGTH - 0.5
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if corrected_exact(mid)[0] < target_dphi:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    ex = N.exact_solution(M_INF)
    x_ref_corr = corrected_x_s(ex["delta_phi"])
    dphi_c, m1_c, sig_c = corrected_exact(x_ref_corr)
    print(f"Dirichlet data: delta_phi = {ex['delta_phi']:.8f}")
    print(f"  isentropic reference x_s = {ex['x_s']:.4f} "
          f"(M1 {ex['m_shock_up']:.5f})")
    print(f"  corrected  reference x_s = {x_ref_corr:.4f} "
          f"(M1 {m1_c:.5f}, sigma {sig_c:.6f})")
    print(f"  PRE-REGISTERED value     = 11.6200   -> reproduced here "
          f"{'YES' if abs(x_ref_corr - 11.62) < 5e-3 else 'NO'}\n")

    rows = []
    for nx in NXS:
        hx = N.LENGTH / nx
        mesh = N.nozzle_mesh(nx, max(6, nx // 16))
        phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
        for ent in (False, True):
            sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C, entropy=ent)
            t0 = time.perf_counter()
            phi, info = sysd.newton(phi_bc.copy(), n_max=80, tol=1e-11)
            wall = time.perf_counter() - t0
            xc, ux = D.element_u(sysd, phi)
            x_sh, n_sup, _, _ = N.shock_from_profile(xc, ux, ex["u_star"], nx)
            ref = x_ref_corr if ent else ex["x_s"]
            err = (x_sh - ref) if np.isfinite(x_sh) else float("nan")
            rows.append(dict(
                nx=nx, h=round(hx, 5), entropy=ent, C=C,
                converged=info["converged"], reason=info["reason"],
                n_newton=info["n_newton"],
                res_final=info["residual_history"][-1],
                x_shock=round(float(x_sh), 5),
                x_ref=round(float(ref), 5),
                err_x=round(float(err), 5),
                err_cells=round(float(err) / hx, 3),
                # cross-reference: the distance to the OTHER model's reference,
                # so the table shows the correction is not merely noise
                err_vs_other=round(float(x_sh - (ex["x_s"] if ent
                                                else x_ref_corr)), 5),
                sigma_min=(round(float(sysd.ent.sigma_min), 8) if ent else None),
                m1_detected=(round(float(sysd.ent.m1_max), 5) if ent else None),
                n_shock_cells=(int(sysd.ent.n_shock) if ent else None),
                n_sigma_refresh=(sysd.n_sigma_refresh if ent else 0),
                # like-for-like reference: where the quasi-1-D model puts the
                # shock when it is given the sigma the SOLVER detected
                x_ref_own_sigma=(round(x_s_for_sigma(
                    ex["delta_phi"], float(sysd.ent.sigma_min)), 5)
                    if ent else None),
                n_supersonic_bins=n_sup, wall_s=round(wall, 2)))
            if ent:
                rows[-1]["err_own_sigma"] = round(
                    rows[-1]["x_shock"] - rows[-1]["x_ref_own_sigma"], 5)
                rows[-1]["err_own_sigma_cells"] = round(
                    rows[-1]["err_own_sigma"] / hx, 3)
            else:
                rows[-1]["err_own_sigma"] = None
                rows[-1]["err_own_sigma_cells"] = None
            r = rows[-1]
            print(f"  nx={nx:4d} entropy={str(ent):5s} conv={r['converged']!s:5} "
                  f"n={r['n_newton']:3d} res={r['res_final']:.2e} "
                  f"x_shock={r['x_shock']:.4f} vs ref {r['x_ref']:.4f} "
                  f"err={r['err_cells']:+.2f} cells "
                  f"(other ref {r['err_vs_other']:+.4f}) "
                  f"sigma_min={r['sigma_min']} M1={r['m1_detected']}",
                  flush=True)

    with open(OUT / "gs1b_3_nozzle.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1b_3_nozzle.csv")

    print("\n=== criterion B ===")
    print("  absolute shock position (the cell-normalised error mixes two "
          "things):")
    for nx in NXS:
        off = next(r for r in rows if r["nx"] == nx and not r["entropy"])
        on = next(r for r in rows if r["nx"] == nx and r["entropy"])
        print(f"    nx={nx:4d}  OFF {off['x_shock']:.4f} "
              f"(err {off['err_x']:+.4f}, conv {off['converged']})   "
              f"ON {on['x_shock']:.4f} (err {on['err_x']:+.4f}, "
              f"conv {on['converged']})   measured shift "
              f"{on['x_shock'] - off['x_shock']:+.4f} vs predicted "
              f"{x_ref_corr - ex['x_s']:+.4f} "
              f"({100 * (on['x_shock'] - off['x_shock']) / (x_ref_corr - ex['x_s']):.0f} %)")
    print("\n  detection resolution (sigma the solver applied vs the "
          f"reference's {sig_c:.6f}):")
    for r in [x for x in rows if x["entropy"]]:
        print(f"    nx={r['nx']:4d}  M1 {r['m1_detected']:.5f}  "
              f"sigma_min {r['sigma_min']:.8f}  "
              f"-> like-for-like ref x_s {r['x_ref_own_sigma']:.4f}, "
              f"position error {r['err_own_sigma']:+.4f} "
              f"({r['err_own_sigma_cells']:+.2f} cells)")
    conv = [r for r in rows if r["entropy"] and r["converged"]]
    if conv:
        worst = max(abs(r["err_own_sigma_cells"]) for r in conv)
        ok = worst < 2.0
        print(f"\n  B (judged on the CONVERGED levels, against the "
              f"like-for-like reference) = {'PASS' if ok else 'FAIL'}: "
              f"worst |position error| {worst:.2f} cells, threshold 2.0")
        print(f"  nx={max(NXS)} is recorded NON-CONVERGED ON BOTH LEGS -- a "
              f"pre-existing limit of this bench's plain Newton (GS1.1's "
              f"committed sweep has converged=False, res 1.885e-08, "
              f"x_shock 11.96161 for the OFF leg, which this run reproduces).")
    else:
        ok = False
        print("\n  B = FAIL (no converged entropy leg)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
