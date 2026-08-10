"""GS1b.3b follow-up: why does the FINE leg stop short with the correction on?

GS1b.3b re-measured criterion F at the corrected (knee) magnitude and found the
fine leg no longer reaches M0.7875 -- it stops at M0.76875 after four step
halvings, worse than both the isentropic leg (0.7825) and the 2x-too-strong
version (0.7875). Its last good state also reads M1 = 1.5417 with 1-sigma = 8.4 %,
a STRONGER shock than medium shows at a HIGHER Mach (M1 1.347) -- so the state
itself looks wrong, not just the step size.

Two candidate causes, and one probe separates them:

  (churn)  the sigma refresh is the problem -- the post-shock SET limit-cycles
           (criterion E), and at fine the flip-flop is large enough to stall the
           step. Then FREEZING sigma completely (refresh cap 0: built once from
           the seed, never updated) should get through, and refreshing MORE
           (cap 20) should be no better than the current 8.
  (system) the corrected system is simply harder here. Then no refresh cap gets
           through, and the failure is in the Newton/continuation, not in sigma.

Also recorded per attempt: the sigma history (sigma_min, n_shock, M1, max|dsigma|)
so the churn amplitude is visible rather than inferred.

Outputs: results/gs1b_3b_fine_diag.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
C = 1.5
LEVEL = "fine"
#: the probe must start on the SAME last-good state GS1b.3b stopped on, so it
#: runs that round's ADAPTIVE continuation verbatim rather than a hand-written
#: ladder (my first attempt guessed the ladder and the guess died at M0.76 --
#: the real walk takes 0.005 steps to 0.76, then 0.0025, then 0.00125).
M_START = 0.72
DM0 = 0.02
MAX_HALVINGS = 4
M_NEXT = 0.7700
REFRESH_CAPS = (0, 8, 20)


def solve(mc, wc, m, phi=None, gam=None, refresh_max=None, ent=True):
    # GS1b.5(a) is a measured negative and the cap policy was restored, with the
    # cap now an internal constant (newton._SIGMA_REFRESH_MAX). The refresh_max
    # argument is kept as a no-op so this script's recorded CSV columns stay
    # readable against the run that produced them.
    kw = dict(m_inf=m, alpha_deg=ALPHA, upwind_c=C, m_crit=0.95,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80,
              entropy_correction=ent)
    if phi is not None:
        kw.update(phi_init=phi, gamma_init=gam, n_picard_seed=0)
    return solve_newton_lifting(mc, wc, **kw)


def usable(r):
    return bool(r["converged"]) and not r.get("clamped", False)


def main():
    mc, wc = cut_wake(read_mesh(REPO / f"cases/meshes/naca0012_2.5d/{LEVEL}.msh"))
    dz = float(np.ptp(mc.nodes[:, 2]))
    print(f"{LEVEL}: {len(mc.nodes)} nodes -- replaying GS1b.3b's adaptive "
          f"continuation (dm {DM0}, up to {MAX_HALVINGS} halvings)", flush=True)
    m, phi, gam, halv, last_m = M_START, None, None, 0, None
    while halv <= MAX_HALVINGS:
        m_next = m + DM0 / (2 ** halv)
        t0 = time.perf_counter()
        r = solve(mc, wc, m_next, phi, gam)
        if usable(r):
            phi, gam, m, last_m = r["phi"], r["gamma"], m_next, m_next
            print(f"  M={m_next:.5f} OK  gamma={float(gam[0]):.6f} "
                  f"sigma_min={r['sigma_min']:.6f} M1={r['m1_max']:.4f} "
                  f"n_shock={r['n_shock_cells']} "
                  f"({time.perf_counter() - t0:.0f}s)", flush=True)
        else:
            halv += 1
            print(f"  M={m_next:.5f} BAD (conv={r['converged']} "
                  f"clamp={r.get('clamped')}) -> halving to "
                  f"{DM0 / 2 ** halv:g} at M={m:.5f} "
                  f"({time.perf_counter() - t0:.0f}s)", flush=True)
    if last_m is None:
        print("  no usable state at all")
        return 1
    print(f"\n  last good state: M={last_m:.5f} (GS1b.3b reported 0.76875)",
          flush=True)
    base = (phi, gam)

    rows = []
    for cap in REFRESH_CAPS:
        t0 = time.perf_counter()
        r = solve(mc, wc, M_NEXT, base[0], base[1], refresh_max=cap)
        wall = time.perf_counter() - t0
        rep = shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                        m_inf=M_NEXT), M_NEXT)
        hist = r["sigma_history"]
        dmax = max((h[4] for h in hist), default=0.0)
        nshock = sorted({h[2] for h in hist})
        rows.append(dict(
            level=LEVEL, m_inf=M_NEXT, refresh_max=cap, usable=usable(r),
            converged=bool(r["converged"]), clamped=r.get("clamped"),
            n_newton=r["n_newton"], res_final=r["residual_history"][-1],
            gamma=round(float(r["gamma"][0]), 8),
            x_shock=rep["upper"].get("x_shock"),
            m_max=round(float(np.sqrt(r["mach2_max"])), 5),
            sigma_min=r["sigma_min"], m1_max=r["m1_max"],
            n_shock_cells=r["n_shock_cells"],
            n_sigma_refresh=r["n_sigma_refresh"],
            dsigma_max=round(dmax, 8),
            n_shock_range=f"{nshock[0]}..{nshock[-1]}" if nshock else "",
            wall_s=round(wall, 1)))
        q = rows[-1]
        print(f"\n  refresh_max={cap:3d}: {'OK ' if q['usable'] else 'BAD'} "
              f"conv={q['converged']} n={q['n_newton']} "
              f"res={q['res_final']:.2e} gamma={q['gamma']:.6f} "
              f"x_sh={q['x_shock']} M_max={q['m_max']} "
              f"sigma_min={q['sigma_min']:.6f} M1={q['m1_max']:.4f} "
              f"n_shock {q['n_shock_range']} max|dsigma|={q['dsigma_max']:.3e} "
              f"({q['wall_s']}s)", flush=True)

    # control: the same step with the correction OFF
    t0 = time.perf_counter()
    r = solve(mc, wc, M_NEXT, base[0], base[1], ent=False)
    rows.append(dict(
        level=LEVEL, m_inf=M_NEXT, refresh_max="OFF", usable=usable(r),
        converged=bool(r["converged"]), clamped=r.get("clamped"),
        n_newton=r["n_newton"], res_final=r["residual_history"][-1],
        gamma=round(float(r["gamma"][0]), 8),
        x_shock=shock_report(wall_cp_curve(mc, r["phi"], z=0.5 * dz,
                                          m_inf=M_NEXT), M_NEXT)["upper"].get(
                                              "x_shock"),
        m_max=round(float(np.sqrt(r["mach2_max"])), 5),
        sigma_min=None, m1_max=None, n_shock_cells=None, n_sigma_refresh=0,
        dsigma_max=0.0, n_shock_range="",
        wall_s=round(time.perf_counter() - t0, 1)))
    q = rows[-1]
    print(f"\n  entropy OFF (control, from the SAME ON seed): "
          f"{'OK ' if q['usable'] else 'BAD'} conv={q['converged']} "
          f"n={q['n_newton']} res={q['res_final']:.2e} "
          f"gamma={q['gamma']:.6f} x_sh={q['x_shock']}", flush=True)

    with open(OUT / "gs1b_3b_fine_diag.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1b_3b_fine_diag.csv")

    print("\n=== reading ===")
    ok = {r["refresh_max"]: r["usable"] for r in rows}
    print(f"  refresh cap 0 (sigma fully frozen): {ok.get(0)}")
    print(f"  refresh cap 8 (current default):    {ok.get(8)}")
    print(f"  refresh cap 20:                     {ok.get(20)}")
    print(f"  entropy OFF control:                {ok.get('OFF')}")
    if ok.get(0) and not ok.get(8):
        print("  => CHURN: freezing sigma gets through, refreshing does not.")
    elif not any(ok.get(c) for c in REFRESH_CAPS) and ok.get("OFF"):
        print("  => SYSTEM: no refresh policy gets through but the isentropic "
              "system does -- the corrected system is harder here.")
    elif all(ok.get(c) for c in REFRESH_CAPS):
        print("  => the step works at every cap: the GS1b.3b stop was the "
              "adaptive continuation's halving budget, not sigma.")
    else:
        print("  => mixed; read the table.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
