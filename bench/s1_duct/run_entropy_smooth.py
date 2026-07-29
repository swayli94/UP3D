"""GS1b.4 criterion S1 + S5: is sigma a STABLE function of the state, and how
sensitive is that to the smoothing widths?

The defect being fixed (measured in GS1b.3b sec 7.3): walking the fine mesh up in
0.005 Mach steps, adjacent converged states read M1 = 1.305, 1.306, 1.380, 1.317,
1.468, 1.492, 1.259 with sigma_min swinging 0.932-0.986 -- 5.4 percentage points of
density correction jumping between neighbouring conditions, which is not physical.
S1 requires max|d sigma_min| between adjacent Mach steps <= 1 percentage point.

Protocol: branch-continue with a FIXED 0.005 step (no halving -- a halving would
change the step between points and confound the comparison), record sigma_min and
M1 at each converged state, and report the largest adjacent jump. Repeated for a
sweep of eps_m so the knob's sensitivity is measured rather than asserted (S5).

Outputs: results/gs1b_4_smooth.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.kernels.entropy import EntropyOperator                # noqa: E402
from pyfp3d.mesh.reader import read_mesh                          # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                         # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting              # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA, C = 1.25, 1.5
LEVEL = "fine"
M_FROM, M_TO, DM = 0.7300, 0.7700, 0.005
EPS_SWEEP = (0.02, 0.05, 0.10)


def walk(mc, wc, eps_m):
    """Fixed-step continuation from M_FROM to M_TO, recording sigma per state."""
    ent_kw = dict(eps_m=eps_m)
    out = []

    def one(m, phi=None, gam=None):
        kw = dict(m_inf=m, alpha_deg=ALPHA, upwind_c=C, m_crit=0.95,
                  freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
                  direct_refactor_every=4, n_newton_max=80,
                  entropy_correction=True, entropy_kwargs=ent_kw)
        if phi is not None:
            kw.update(phi_init=phi, gamma_init=gam, n_picard_seed=0)
        r = solve_newton_lifting(mc, wc, **kw)
        return r, bool(r["converged"]) and not r.get("clamped", False)

    # ★ fine cannot be cold-started at M_FROM (measured: it does not converge),
    # so reach it by the same adaptive walk the other GS1b scripts use, and only
    # then switch to the FIXED step the stability reading needs.
    m, phi, gam, halv = 0.72, None, None, 0
    while halv <= 4 and m < M_FROM - 1e-12:
        m_next = min(m + 0.02 / (2 ** halv), M_FROM)
        r, ok = one(m_next, phi, gam)
        if ok:
            phi, gam, m = r["phi"], r["gamma"], m_next
        else:
            halv += 1
    if abs(m - M_FROM) > 1e-12:
        print(f"    seed walk stalled at M={m:.5f} (target {M_FROM})")
        return out
    out.append((M_FROM, r["sigma_min"], r["m1_max"], r["n_shock_cells"]))
    print(f"    M={M_FROM:.5f} OK  sigma_min={r['sigma_min']:.6f} "
          f"M1={r['m1_max']:.4f} n_shock={r['n_shock_cells']} (seed)",
          flush=True)
    m = M_FROM
    while m < M_TO - 1e-12:
        m = round(m + DM, 6)
        t0 = time.perf_counter()
        r, ok = one(m, phi, gam)
        print(f"    M={m:.5f} {'OK ' if ok else 'BAD'} "
              f"sigma_min={r['sigma_min']:.6f} M1={r['m1_max']:.4f} "
              f"n_shock={r['n_shock_cells']} ({time.perf_counter()-t0:.0f}s)",
              flush=True)
        if not ok:
            break
        phi, gam = r["phi"], r["gamma"]
        out.append((m, r["sigma_min"], r["m1_max"], r["n_shock_cells"]))
    return out


def main():
    mc, wc = cut_wake(read_mesh(REPO / f"cases/meshes/naca0012_2.5d/{LEVEL}.msh"))
    rows = []
    print(f"{LEVEL}: {len(mc.nodes)} nodes, fixed step {DM} from {M_FROM} to "
          f"{M_TO}\n")
    for eps in EPS_SWEEP:
        print(f"  eps_m = {eps}", flush=True)
        pts = walk(mc, wc, eps)
        if len(pts) < 2:
            print(f"    (only {len(pts)} state(s) -- cannot form a jump)")
            continue
        sig = np.array([p[1] for p in pts])
        m1 = np.array([p[2] for p in pts])
        dsig = np.abs(np.diff(sig))
        for (m, s_, mm, ns) in pts:
            rows.append(dict(level=LEVEL, eps_m=eps, m_inf=m,
                             sigma_min=round(s_, 8), m1=round(mm, 5),
                             n_shock_cells=ns))
        print(f"    -> sigma_min range {sig.min():.6f}..{sig.max():.6f}, "
              f"max adjacent |d sigma| = {100*dsig.max():.3f} pp, "
              f"M1 range {m1.min():.4f}..{m1.max():.4f} "
              f"({len(pts)} states)\n", flush=True)

    if not rows:
        print("no states recorded -- nothing to write")
        return
    with open(OUT / "gs1b_4_smooth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "gs1b_4_smooth.csv")

    print("\n=== criterion S1 (<= 1.0 pp) and S5 (eps_m sensitivity) ===")
    print("  reference, HARD version (GS1b.3b sec 7.3): sigma_min 0.932..0.986,"
          " max adjacent jump 4.6 pp, M1 1.259..1.492")
    for eps in EPS_SWEEP:
        sub = [r for r in rows if r["eps_m"] == eps]
        if len(sub) < 2:
            continue
        sig = np.array([r["sigma_min"] for r in sub])
        m1 = np.array([r["m1"] for r in sub])
        jump = 100.0 * np.abs(np.diff(sig)).max()
        print(f"  eps_m={eps:<5} max adjacent |d sigma| = {jump:6.3f} pp  "
              f"{'PASS' if jump <= 1.0 else 'FAIL'}   "
              f"sigma {sig.min():.6f}..{sig.max():.6f}  "
              f"M1 {m1.min():.4f}..{m1.max():.4f}  ({len(sub)} states)")


if __name__ == "__main__":
    main()
