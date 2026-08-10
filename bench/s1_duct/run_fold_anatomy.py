"""GS1b.2 Q4: is the branch termination a LIMIT POINT, or is it my continuation
control / the field solve failing?

Why this run exists (it questions my own GS1b.1 verdict). Q2's C = 3.0 fine leg
died at M_inf 0.75 -- EARLIER than the C = 1.5 fine leg (0.77) -- which breaks
the monotone-in-C story that holds at medium, and its last good point carried
dGamma/dM = 1.14 (medium's terminal value was 13.6) while M_max JUMPED
discontinuously 1.258 -> 1.525 on the failing step. Two possible readings:

  (F) a genuine limit point: the coupled system has no solution past M_last;
  (S) my walker: `run_fold_branch.py` refines the step ONCE (one-shot flag), and
      at C = 3.0 fine that single refinement was consumed early by a 0.02-step
      failure at M 0.74, so the walk ended at the FIRST 0.005-step failure.

Two probes separate them, both seeded from the same converged last-good state:

  P1  step refinement -- retry M_last + delta for delta in {2.5, 1.0, 0.5}e-3.
      If a smaller step converges, the termination was (S), not (F). The walk
      then continues on the smallest delta to see how far the branch really goes.
  P2  fixed-Gamma field probe -- at M_last + 2.5e-3 freeze Gamma at its last-good
      value and solve the FIELD block alone (the GS1.6 / A2 discriminator).
        field converges + coupled fails  => the singularity is in the GAMMA
                                            direction: a fold (F)
        field also fails                 => the failure is in the field/shock,
                                            NOT in the lift closure
        both fine at a smaller step      => (S), continuation control

P2 is what makes this implementation-independent: it does not read any Schur
complement or condition number (GS1b.1 sec 1.1 showed those are
rendering-dependent), only "does a solution exist when Gamma is held".

Outputs: results/gs1b_2_q4_anatomy.csv
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                      # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve              # noqa: E402
from pyfp3d.post.shock import shock_report                     # noqa: E402
from pyfp3d.solve.newton import NewtonWorkspace                # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting           # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 1.25
M_CRIT, M_CAP, RHO_FLOOR = 0.95, 3.0, 0.05

#: (level, C, mach ladder up to the measured last-good point). The ladders
#: reproduce the GS1b.1 / Q2 walks exactly, so the probe starts from the SAME
#: state those runs terminated on.
CASES = (
    ("medium", 1.5, [0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.7925]),
    ("fine",   1.5, [0.72, 0.74, 0.76, 0.765, 0.77]),
    ("fine",   3.0, [0.72, 0.725, 0.73, 0.735, 0.74, 0.745, 0.75]),
)
DELTAS = (0.0025, 0.001, 0.0005)
N_CONTINUE = 8          # extra steps to take once a delta works


def coupled(mc, wc, m_inf, C, phi=None, gam=None):
    kw = dict(m_inf=m_inf, alpha_deg=ALPHA, upwind_c=C, m_crit=M_CRIT,
              freeze_tol=1e-6, freeze_refresh_max=8, precond="direct",
              direct_refactor_every=4, n_newton_max=80)
    if phi is not None:
        kw.update(phi_init=phi, gamma_init=gam, n_picard_seed=0)
    return solve_newton_lifting(mc, wc, **kw)


def usable(r):
    return bool(r["converged"]) and not r.get("clamped", False)


def field_only(mc, wc, m_inf, C, phi_seed, gamma_fix, n_max=60, tol=1e-9):
    """Newton on the field block with Gamma FROZEN (no Kutta row)."""
    ws = NewtonWorkspace(mc, wc, alpha_deg=ALPHA)
    ws.set_mach(m_inf)
    phi_free = np.asarray(phi_seed, dtype=np.float64)[:ws.n_red][ws.free].copy()
    gamma = np.asarray([gamma_fix], dtype=np.float64)
    R, _, state = ws.eval_residual(phi_free, gamma, C, M_CRIT, M_CAP, RHO_FLOOR)
    r0 = float(np.max(np.abs(R)))
    for _ in range(n_max):
        J_ff, _ = ws.assemble_coupled(state, C, M_CRIT, RHO_FLOOR)
        try:
            d = spla.spsolve(J_ff.tocsc(), -R)
        except Exception:                                      # noqa: BLE001
            return False, r0, "linear solve failed", state
        if not np.all(np.isfinite(d)):
            return False, r0, "non-finite step", state
        lam, best = 1.0, None
        for _ in range(12):
            trial = phi_free + lam * d
            Rt, _, st = ws.eval_residual(trial, gamma, C, M_CRIT, M_CAP,
                                         RHO_FLOOR)
            rt = float(np.max(np.abs(Rt)))
            if np.isfinite(rt) and rt < r0:
                break
            if best is None or (np.isfinite(rt) and rt < best[0]):
                best = (rt, trial, Rt, st)
            lam *= 0.5
        else:
            if best is None:
                return False, r0, "line search", state
            rt, trial, Rt, st = best
        phi_free, R, state, r0 = trial, Rt, st, rt
        if rt < tol:
            return True, rt, "tol", state
    return False, r0, "cap", state


def main():
    rows = []
    for level, C, ladder in CASES:
        path = REPO / f"cases/meshes/naca0012_2.5d/{level}.msh"
        if not path.exists():
            print(f"skip {level}: mesh missing")
            continue
        mc, wc = cut_wake(read_mesh(path))
        dz = float(np.ptp(mc.nodes[:, 2]))
        print(f"\n=== {level}, C = {C}: walking to M {ladder[-1]} ===",
              flush=True)
        phi = gam = None
        for m in ladder:
            r = coupled(mc, wc, m, C, phi, gam)
            if not usable(r):
                print(f"   ! ladder point M={m} not usable "
                      f"(conv={r['converged']} clamp={r.get('clamped')}) "
                      f"-- anatomy invalid for this case", flush=True)
                phi = None
                break
            phi, gam = r["phi"], r["gamma"]
        if phi is None:
            continue
        m_last = ladder[-1]
        g_last = float(gam[0])
        mmax_last = float(np.sqrt(r["mach2_max"]))
        print(f"   arrived: gamma={g_last:.6f} M_max={mmax_last:.4f}",
              flush=True)

        # ---- P1: does a smaller step get past the termination? ----
        worked = None
        for d in DELTAS:
            t0 = time.perf_counter()
            rr = coupled(mc, wc, m_last + d, C, phi, gam)
            ok = usable(rr)
            mm = float(np.sqrt(rr["mach2_max"]))
            rows.append(dict(
                level=level, upwind_c=C, probe="P1_step", m_last=m_last,
                delta=d, m_inf=round(m_last + d, 6), usable=ok,
                converged=bool(rr["converged"]), clamped=rr.get("clamped"),
                gamma=round(float(rr["gamma"][0]), 8),
                m_max=round(mm, 5), dm_max_jump=round(mm - mmax_last, 5),
                res_final=rr["residual_history"][-1], n_newton=rr["n_newton"],
                note="", wall_s=round(time.perf_counter() - t0, 1)))
            print(f"   P1 delta={d:<7} M={m_last + d:.5f} "
                  f"{'OK ' if ok else 'BAD'} gamma={rows[-1]['gamma']:.6f} "
                  f"M_max={mm:.4f} (jump {mm - mmax_last:+.4f}) "
                  f"res={rows[-1]['res_final']:.2e}", flush=True)
            if ok and worked is None:
                worked = (d, rr)

        # ---- P1b: if a delta worked, keep walking on it ----
        if worked is not None:
            d, rr = worked
            p, g = rr["phi"], rr["gamma"]
            m = m_last + d
            for _ in range(N_CONTINUE):
                m += d
                r2 = coupled(mc, wc, m, C, p, g)
                ok = usable(r2)
                mm = float(np.sqrt(r2["mach2_max"]))
                rows.append(dict(
                    level=level, upwind_c=C, probe="P1b_continue",
                    m_last=m_last, delta=d, m_inf=round(m, 6), usable=ok,
                    converged=bool(r2["converged"]),
                    clamped=r2.get("clamped"),
                    gamma=round(float(r2["gamma"][0]), 8), m_max=round(mm, 5),
                    dm_max_jump="", res_final=r2["residual_history"][-1],
                    n_newton=r2["n_newton"], note="", wall_s=""))
                print(f"   P1b M={m:.5f} {'OK ' if ok else 'BAD'} "
                      f"gamma={rows[-1]['gamma']:.6f} M_max={mm:.4f}",
                      flush=True)
                if not ok:
                    break
                p, g = r2["phi"], r2["gamma"]

        # ---- P2: fixed-Gamma field probe at the first failing step ----
        t0 = time.perf_counter()
        ok_f, res_f, why_f, st = field_only(mc, wc, m_last + DELTAS[0], C,
                                            phi, g_last)
        rep = shock_report(wall_cp_curve(mc, st["phi_cut"], z=0.5 * dz,
                                        m_inf=m_last + DELTAS[0]),
                           m_last + DELTAS[0])
        rows.append(dict(
            level=level, upwind_c=C, probe="P2_fixed_gamma", m_last=m_last,
            delta=DELTAS[0], m_inf=round(m_last + DELTAS[0], 6), usable=ok_f,
            converged=ok_f, clamped="", gamma=round(g_last, 8),
            m_max=round(float(np.sqrt(st["mach2_max"]))
                        if "mach2_max" in st else float("nan"), 5),
            dm_max_jump="", res_final=res_f, n_newton="",
            note=f"{why_f}; x_shock={rep['upper'].get('x_shock')}",
            wall_s=round(time.perf_counter() - t0, 1)))
        print(f"   P2 fixed-Gamma at M={m_last + DELTAS[0]:.5f}: "
              f"{'CONVERGED' if ok_f else 'FAILED'} |R|={res_f:.3e} ({why_f})",
              flush=True)

    with open(OUT / "gs1b_2_q4_anatomy.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT / "gs1b_2_q4_anatomy.csv")

    print("\n=== reading ===")
    for level, C, ladder in CASES:
        sub = [r for r in rows if r["level"] == level and r["upwind_c"] == C]
        if not sub:
            continue
        p1 = [r for r in sub if r["probe"] == "P1_step"]
        p1b = [r for r in sub if r["probe"] == "P1b_continue"]
        p2 = [r for r in sub if r["probe"] == "P2_fixed_gamma"]
        got = [r for r in p1 if r["usable"]]
        reach = max([r["m_inf"] for r in p1b + p1 if r["usable"]],
                    default=None)
        print(f"  {level:7s} C={C}: smaller step works = "
              f"{'YES @ delta ' + str(got[0]['delta']) if got else 'NO'};"
              f" branch reaches {reach}; fixed-Gamma field at the failing step ="
              f" {'CONVERGES' if p2 and p2[0]['usable'] else 'FAILS'}")


if __name__ == "__main__":
    main()
