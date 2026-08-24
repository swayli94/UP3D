"""GS1.1 Case B sweep: where does the SHIPPED shock operator put a shock whose
exact position is uniquely determined?

Protocol (pre-registered, phases/p2/docs/dev_phase_two/20260728-1640-s1-shock-bench.md):

  * Laval nozzle, choked, exact quasi-1-D shock at x_s = 12.0 of a duct of
    length 20, shock upstream Mach 1.30; uniqueness of the boundary-value
    problem verified numerically (nozzle.verify_uniqueness);
  * Dirichlet phi at inlet and outlet from the exact solution; walls natural;
  * solve the shipped discrete system (PicardOperator + UpwindOperator walk
    flux) with a damped Newton;
  * sweep the artificial-dissipation constant C, the mesh size h, and mesh
    regularity (structured vs interior-jittered);
  * a second start (shock initialised at x = 10 instead of the exact 12) tests
    whether the DISCRETE problem has spurious multiple solutions.

Outputs: results/gs1_1_nozzle_sweep.csv, results/gs1_1_nozzle.png,
         results/gs1_1_nozzle_profiles.csv
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
CS = (1.0, 1.5, 2.0, 3.0)
NXS = (100, 200, 400)
LEGS = (("regular", 0.0), ("irregular", 0.35))
PERTURBED_NX = 200          # the multiple-solution probe runs on this level
X_S_PERTURBED = 10.0


def run_one(nx, jitter, C, ex, x_s_init, rows, prof_rows, tag):
    ny = max(6, nx // 16)
    hx = N.LENGTH / nx
    mesh = N.nozzle_mesh(nx, ny, jitter=jitter)
    sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C)

    # Dirichlet data ALWAYS from the exact (x_s = 12) solution; only the
    # initial guess changes for the perturbed probe.
    phi_bc = ex["phi_of_x"](mesh.nodes[:, 0])
    if x_s_init == ex["x_s"]:
        phi0 = phi_bc.copy()
    else:
        ex2 = N.exact_solution(M_INF, x_s=x_s_init)
        phi0 = ex2["phi_of_x"](mesh.nodes[:, 0])
        phi0[sysd.dir_nodes] = phi_bc[sysd.dir_nodes]     # keep the true BCs

    R0, _ = sysd.residual(phi0)
    t0 = time.perf_counter()
    phi, info = sysd.newton(phi0, n_max=80, tol=1e-11)
    wall = time.perf_counter() - t0

    xc, ux = D.element_u(sysd, phi)
    x_sh, n_sup, xb, ub = N.shock_from_profile(xc, ux, ex["u_star"], nx)
    _, _, _, ub_ex = N.shock_from_profile(
        *D.element_u(sysd, phi_bc), ex["u_star"], nx)
    _, fb = D.mass_flux_profile(sysd, phi, ex, n_bins=nx)
    ok = ~np.isnan(fb)
    defect = float(np.sum(np.abs(1.0 - fb[ok])) * (N.LENGTH / nx))

    err = (x_sh - ex["x_s"]) if np.isfinite(x_sh) else None
    row = dict(
        leg="regular" if jitter == 0.0 else "irregular", start=tag, nx=nx,
        h=round(hx, 5), n_dof=len(mesh.nodes), n_tets=len(mesh.elements), C=C,
        converged=info["converged"], reason=info["reason"],
        n_newton=info["n_newton"], wall_s=round(wall, 2),
        res_exact=float(np.max(np.abs(R0[sysd.free]))),
        res_final=info["residual_history"][-1],
        x_shock=round(x_sh, 5) if np.isfinite(x_sh) else None,
        err_x=round(err, 5) if err is not None else None,
        err_cells=round(err / hx, 3) if err is not None else None,
        n_supersonic_bins=n_sup,
        massflux_min=round(float(np.nanmin(fb)), 5),
        massflux_max=round(float(np.nanmax(fb)), 5),
        massflux_defect_L1=round(defect, 6))
    print("  ", {k: row[k] for k in ("leg", "start", "nx", "C", "converged",
                                     "n_newton", "x_shock", "err_cells",
                                     "massflux_defect_L1", "res_final")},
          flush=True)
    rows.append(row)
    if C == 1.5 and jitter == 0.0 and tag == "exact":
        for x, u, u_ex in zip(xb, ub, ub_ex):
            prof_rows.append(dict(nx=nx, x=round(float(x), 4),
                                  u=round(float(u), 6),
                                  u_exact_on_mesh=round(float(u_ex), 6)))
    return row


def figure(rows, ex):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:                                          # noqa: BLE001
        print("(matplotlib unavailable)")
        return
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5))
    for leg, mk in (("regular", "o-"), ("irregular", "s--")):
        for C in CS:
            sub = sorted([r for r in rows
                          if r["leg"] == leg and r["C"] == C
                          and r["start"] == "exact" and r["err_x"] is not None],
                         key=lambda r: r["h"])
            if not sub:
                continue
            hs = [r["h"] for r in sub]
            axes[0].plot(hs, [abs(r["err_x"]) for r in sub], mk,
                         label=f"{leg} C={C}", ms=4)
            axes[1].plot(hs, [abs(r["err_cells"]) for r in sub], mk, ms=4)
            axes[2].plot(hs, [r["massflux_defect_L1"] for r in sub], mk, ms=4)
    h_ref = np.array([0.05, 0.2])
    axes[0].plot(h_ref, 0.5 * h_ref / h_ref[1] * 1.0, "k:", lw=1,
                 label="O(h) slope")
    for ax, lab, ttl in zip(
            axes,
            ["|x_shock - 12.0|", "|error| / cell size",
             "L1 mass-flux defect"],
            ["shock-position error vs h", "error in cells",
             "conservation error of the operator"]):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("h_x = L/nx")
        ax.set_ylabel(lab)
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.3)
    axes[1].axhline(1.0, color="k", lw=0.8, ls=":")
    axes[0].legend(fontsize=6, ncol=2)
    fig.suptitle("GS1.1 Case B: choked nozzle, exact shock at x=12 "
                 f"(M_shock {ex['m_shock_up']:.2f}); Dirichlet phi from the "
                 "exact solution", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "gs1_1_nozzle.png", dpi=130)
    print("wrote", OUT / "gs1_1_nozzle.png")


def main():
    ex = N.exact_solution(M_INF)
    dphi_s, dphi_sub, uniq = N.verify_uniqueness(M_INF)
    print(f"exact: u* {ex['u_star']:.6f} mdot {ex['mdot']:.6f} "
          f"delta_phi {ex['delta_phi']:.6f} M_inlet {ex['m_inlet']:.4f} "
          f"M_shock {ex['m_shock_up']:.4f}")
    print(f"uniqueness: shocked {dphi_s:.6f} > unchoked-max {dphi_sub:.6f} "
          f"-> {uniq}", flush=True)
    assert uniq, "boundary data is not uniquely satisfied by the shocked branch"

    rows, prof_rows = [], []
    for leg, jitter in LEGS:
        print(f"\n=== {leg} mesh, start = exact ===", flush=True)
        for nx in NXS:
            for C in CS:
                run_one(nx, jitter, C, ex, ex["x_s"], rows, prof_rows, "exact")
    print(f"\n=== multiple-solution probe: start with the shock at "
          f"x = {X_S_PERTURBED} (nx = {PERTURBED_NX}) ===", flush=True)
    for leg, jitter in LEGS:
        for C in CS:
            run_one(PERTURBED_NX, jitter, C, ex, X_S_PERTURBED, rows,
                    prof_rows, "perturbed")

    with open(OUT / "gs1_1_nozzle_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    if prof_rows:
        with open(OUT / "gs1_1_nozzle_profiles.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(prof_rows[0]))
            w.writeheader()
            w.writerows(prof_rows)
    print("\nwrote", OUT / "gs1_1_nozzle_sweep.csv")
    figure(rows, ex)

    print("\n=== headline: x_shock (exact 12.0), start = exact ===")
    print(f"{'leg':10s} {'nx':>4s} " + " ".join(f"{'C=' + str(c):>10s}"
                                               for c in CS) + "   spread")
    for leg, _ in LEGS:
        for nx in NXS:
            vals, cells = [], []
            for C in CS:
                m = [r for r in rows if r["leg"] == leg and r["nx"] == nx
                     and r["C"] == C and r["start"] == "exact"]
                v = m[0]["x_shock"] if m else None
                vals.append(v)
                cells.append(f"{v:10.4f}" if v is not None else f"{'--':>10s}")
            good = [v for v in vals if v is not None]
            spread = (max(good) - min(good)) if len(good) > 1 else float("nan")
            print(f"{leg:10s} {nx:4d} " + " ".join(cells)
                  + f"   {spread:.4f}")
    print("\n=== multiple-solution probe: x_shock from a start at x=10 ===")
    for leg, _ in LEGS:
        for C in CS:
            m = [r for r in rows if r["leg"] == leg and r["C"] == C
                 and r["start"] == "perturbed"]
            e = [r for r in rows if r["leg"] == leg and r["C"] == C
                 and r["start"] == "exact" and r["nx"] == PERTURBED_NX]
            if m and e:
                print(f"  {leg:10s} C={C:<4} start@10 -> "
                      f"{m[0]['x_shock']}   (start@12 -> {e[0]['x_shock']})")


if __name__ == "__main__":
    main()
