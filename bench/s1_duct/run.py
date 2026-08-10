"""GS1.1: does the shipped shock operator hold a shock where mass conservation
puts it? Sweep over the artificial-dissipation constant C, mesh size h and mesh
regularity.

Protocol (pre-registered in
docs/dev_phase_two/20260728-1640-s1-shock-bench.md):

  * constant-area duct, exact piecewise-linear stationary shock at x_s = 2.0
    of a duct of length 4 -- an EXACT weak solution of the full-potential
    equation for any x_s, with zero background artificial dissipation
    (rho_e == rho_upstream in both uniform states);
  * Dirichlet phi at inlet and outlet taken from the exact solution, walls
    natural; initial guess = the exact solution;
  * solve the SHIPPED discrete system (PicardOperator + UpwindOperator walk
    flux) with a damped Newton to ||R||_inf ~ 1e-12;
  * measure where the shock ends up. The exact problem has a one-parameter
    family (any x_s), so a consistent operator must leave it within ~1 cell:
    every unit of drift is spurious.

Outputs: results/gs1_1_duct_sweep.csv, results/gs1_1_drift.png,
         results/gs1_1_massflux.csv

Run:
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
        python bench/s1_duct/run.py
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

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

M_INF, M_SUP = 0.80, 1.30
LENGTH, HEIGHT, DZ = 4.0, 1.0, 0.1
X_S = 2.0
CS = (1.0, 1.5, 2.0, 3.0)
NXS = (40, 80, 160)
LEGS = (("regular", 0.0), ("irregular", 0.35))


def one(nx, jitter, C, st, rows, flux_rows):
    ny = max(4, nx // 8)
    h = LENGTH / nx
    mesh = D.duct_mesh(nx, ny, LENGTH, HEIGHT, DZ, jitter=jitter, seed=7)
    sysd = D.DuctSystem(mesh, m_inf=M_INF, upwind_c=C)
    phi0 = D.phi_exact(mesh.nodes[:, 0], X_S, st)
    R0, _ = sysd.residual(phi0)
    r_exact = float(np.max(np.abs(R0[sysd.free])))
    x0, w0 = D.shock_position(sysd, phi0, st, n_bins=nx)

    t0 = time.perf_counter()
    phi, info = sysd.newton(phi0, n_max=60, tol=1e-12)
    wall = time.perf_counter() - t0
    x1, w1 = D.shock_position(sysd, phi, st, n_bins=nx)
    xb, fb = D.mass_flux_profile(sysd, phi, st, n_bins=nx)
    ok = ~np.isnan(fb)
    defect = float(np.sum((1.0 - fb[ok]) * (xb[1] - xb[0])))   # integral defect

    row = dict(
        leg="regular" if jitter == 0.0 else "irregular", nx=nx, h=round(h, 5),
        n_dof=len(mesh.nodes), n_tets=len(mesh.elements), C=C,
        converged=info["converged"], reason=info["reason"],
        n_newton=info["n_newton"], wall_s=round(wall, 2),
        res_exact=r_exact, res_final=info["residual_history"][-1],
        x_shock_exact=round(x0, 5) if np.isfinite(x0) else None,
        x_shock_final=round(x1, 5) if np.isfinite(x1) else None,
        drift=round(x1 - X_S, 5) if np.isfinite(x1) else None,
        drift_cells=round((x1 - X_S) / h, 3) if np.isfinite(x1) else None,
        width_bins_exact=w0, width_bins_final=w1,
        massflux_min=round(float(np.nanmin(fb)), 5),
        massflux_defect_int=round(defect, 6))
    print("  ", {k: v for k, v in row.items()
                 if k in ("leg", "nx", "C", "converged", "n_newton",
                          "x_shock_final", "drift_cells", "massflux_min",
                          "massflux_defect_int", "res_final")}, flush=True)
    rows.append(row)
    for x, f in zip(xb, fb):
        flux_rows.append(dict(leg=row["leg"], nx=nx, C=C, x=round(float(x), 5),
                              flux_over_mdot=round(float(f), 6)
                              if np.isfinite(f) else None))
    return row


def figure(rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:                                          # noqa: BLE001
        print("(matplotlib unavailable; skipping figure)")
        return
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for leg, mk in (("regular", "o-"), ("irregular", "s--")):
        for C in CS:
            sub = [r for r in rows if r["leg"] == leg and r["C"] == C
                   and r["drift"] is not None]
            if not sub:
                continue
            hs = [r["h"] for r in sub]
            axes[0].plot(hs, [abs(r["drift"]) for r in sub], mk,
                         label=f"{leg} C={C}")
            axes[1].plot(hs, [abs(r["drift_cells"]) for r in sub], mk,
                         label=f"{leg} C={C}")
            axes[2].plot(hs, [r["massflux_defect_int"] for r in sub], mk,
                         label=f"{leg} C={C}")
    for ax, lab, ttl in zip(
            axes,
            ["|drift| (chord-equivalent length)", "|drift| / cell size",
             "integral mass-flux defect"],
            ["shock drift vs h", "drift in cells (should be <= ~1)",
             "conservation error of the shock operator"]):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("h = L/nx")
        ax.set_ylabel(lab)
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.3)
    axes[1].axhline(1.0, color="k", lw=0.8, ls=":")
    axes[0].legend(fontsize=6, ncol=2)
    fig.suptitle("GS1.1 constant-area duct, exact stationary shock at x=2.0 "
                 f"(M_inf {M_INF}, M_shock {M_SUP})", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "gs1_1_drift.png", dpi=130)
    print("wrote", OUT / "gs1_1_drift.png")


def main():
    st = D.duct_states(M_INF, M_SUP)
    print("exact states:", {k: round(v, 6) for k, v in st.items()}, flush=True)
    rows, flux_rows = [], []
    for leg, jitter in LEGS:
        print(f"\n=== {leg} mesh (jitter {jitter}) ===", flush=True)
        for nx in NXS:
            for C in CS:
                one(nx, jitter, C, st, rows, flux_rows)

    with open(OUT / "gs1_1_duct_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(OUT / "gs1_1_massflux.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(flux_rows[0]))
        w.writeheader()
        w.writerows(flux_rows)
    print("\nwrote", OUT / "gs1_1_duct_sweep.csv")
    figure(rows)

    print("\n=== headline: drift in cells ===")
    print(f"{'leg':10s} {'nx':>4s} " + " ".join(f"{'C=' + str(c):>9s}"
                                               for c in CS))
    for leg, _ in LEGS:
        for nx in NXS:
            cells = []
            for C in CS:
                m = [r for r in rows if r["leg"] == leg and r["nx"] == nx
                     and r["C"] == C]
                cells.append(f"{m[0]['drift_cells']:9.2f}"
                             if m and m[0]["drift_cells"] is not None
                             else f"{'--':>9s}")
            print(f"{leg:10s} {nx:4d} " + " ".join(cells))


if __name__ == "__main__":
    main()
