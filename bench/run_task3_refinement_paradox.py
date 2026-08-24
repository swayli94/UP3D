"""Task 3 item 1: WHY does refinement make robustness worse? -- classify it first.

Pre-registered in phases/p3/docs/dev_phase_three/20260811-1500-task3-refinement-paradox-prereg.md,
committed before this file was written.

GS1.2 measured the paradox and never classified it: refinement makes convergence WORSE (fine
stops closing near M0.75 at alpha 1.25) with `n_floored = 0`, so it is NOT clamping. That leaves
four modes -- limit cycle, ill-conditioning, line-search collapse, sigma-transport -- whose fixes
live in four completely different places (frozen-selection churn, preconditioner, globalisation,
entropy correction). CLAUDE.md forbids reporting a bare conv=False for exactly this reason, and
phase two lost five hypotheses in one day to correlating a mixture of diseases against one knob.
So this script classifies; it does not fix.

★★ The leg phase two could not run is RS: refine ONLY the shock band. All four of its generator
knobs were measured out of scope, so it could only refine globally and could never separate "the
shock region's resolution" from "the global cell count". That separation is what route (A)'s
experimental controllability bought, and this is it being spent.

Premise, measured and fail-fast: R0 must converge with 0/0 clamps and yield a readable x_shock,
or there is no baseline and the item stops. The shock window is then taken FROM that measurement
rather than assumed.

Outputs (TRACKED): bench/gate_results/task3_refinement_paradox.csv
"""

import argparse
import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.wake_cut import cut_wake                          # noqa: E402
from pyfp3d.meshgen.extrude import extrude_single_layer            # noqa: E402
from pyfp3d.meshgen.planar import naca0012_coordinates            # noqa: E402
from pyfp3d.meshgen.structured import (airfoil_hybrid_2d,          # noqa: E402
                                       airfoil_surface_distribution)
from pyfp3d.post.section_cut import wall_cp_curve                  # noqa: E402
from pyfp3d.post.shock import shock_report                         # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients            # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting               # noqa: E402
from run_le14_common_root import classify_failure                  # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_refinement_paradox.csv")
M_INF, ALPHA = 0.80, 1.25          # the M1 condition -- this item WANTS the out-of-envelope failure
REF = naca0012_coordinates(n_half=401)
#: ★ seed FIXED and recorded: M1's own measurement has seed 0 failing at M0.80 medium and being
#: rescued only by the thread-dependent fallback, so leaving it to the default would put a known
#: thread dependence inside the comparison. Threads are pinned above for the same reason.
SEED = 5
SOLVE_KW = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
                precond="direct", direct_refactor_every=4, n_newton_max=80,
                n_picard_seed=SEED)
BLOCK_THICKNESS, GROWTH, TE_BLEND, LE_CLUSTER = 0.08, 1.15, 0.05, 0.8
LEG_GATE_S = 600.0


def build(n_stations, h_wall_normal, window=None, factor=1.0, dz=0.10):
    surf = airfoil_surface_distribution(REF, n_stations, le_cluster=LE_CLUSTER,
                                        local_window=window, local_factor=factor)
    pts, tris, edges, interior, info = airfoil_hybrid_2d(
        surf, thickness=BLOCK_THICKNESS, h_wall_normal=h_wall_normal, growth=GROWTH,
        te_blend=TE_BLEND)
    mesh = extrude_single_layer(pts, tris, edges, interior_edge_groups=interior, dz=dz,
                                name="task3")
    return cut_wake(mesh), info


def run_leg(tag, n_stations, h_wall_normal, window, factor, run_id, alpha=None):
    alpha = ALPHA if alpha is None else alpha
    t0 = time.perf_counter()
    (mc, wc), info = build(n_stations, h_wall_normal, window, factor)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=alpha, **SOLVE_KW)
    wall = time.perf_counter() - t0
    hist = np.asarray(r.get("residual_history", []), dtype=float)
    nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
    #: ★ never a bare conv=False: classify, and let accept_reason outrank any inferred signature
    if r["converged"]:
        mode, ev = "", f"accept_reason={r.get('accept_reason')!r}"
    else:
        mode, ev, d10, _ = classify_failure(
            hist, np.asarray(r.get("clamp_history", []), dtype=float),
            np.asarray(r.get("F_history", []), dtype=float),
            int(r.get("n_gmres_stalled") or 0), str(r.get("accept_reason")), nlim, nflr)
    z = float(np.unique(mc.nodes[:, 2]).mean())
    rep = shock_report(wall_cp_curve(mc, r["phi"], z=z, m_inf=M_INF), M_INF)
    s_ref = float(np.ptp(mc.nodes[:, 2]))
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=alpha, s_ref=s_ref, m_inf=M_INF)
    return dict(leg=tag, alpha=alpha, run=run_id, n_stations=n_stations, h_wall_normal=h_wall_normal,
                window=str(window), factor=factor, n_tets=len(mc.elements),
                converged=bool(r["converged"]), res_final=float(hist[-1]) if len(hist) else None,
                n_limited=nlim, n_floored=nflr, failure_mode=mode, mode_evidence=ev,
                accept_reason=str(r.get("accept_reason")), n_newton=len(hist),
                has_shock=bool(rep["upper"].get("has_shock")),
                x_shock=rep["upper"].get("x_shock"), cl=round(float(f["cl"]), 6),
                m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                seed=SEED, threads=os.environ["NUMBA_NUM_THREADS"], wall_s=round(wall, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha0-control", action="store_true",
                    help="the GS1.2b discriminator: same four legs at alpha=0, where phase two "
                         "measured the shock position CONVERGING while alpha=1.25 diverged")
    ap.add_argument("--baseline-only", action="store_true",
                    help="run R0 only -- the premise gate, before spending on four legs")
    a = ap.parse_args()
    print(f"task 3 item 1: classify the refinement paradox   M{M_INF} / alpha {ALPHA}  "
          f"seed {SEED}  threads {os.environ['NUMBA_NUM_THREADS']}\n")

    LEGS_ALL = (("R0_baseline", 160, 0.004, None, 1.0),
                ("RG_global", 320, 0.002, None, 1.0),
                ("RS_shock_band", 160, 0.004, (0.45, 0.75), 3.0),
                ("RN_normal_only", 160, 0.002, None, 1.0))
    if a.alpha0_control:
        #: ★★ the leg my pre-registration wrongly omitted. GS1.2b already CLASSIFIED this paradox
        #: (progress row 6): at alpha=0 the shock position CONVERGES (0.464/0.502/0.513) while at
        #: alpha=1.25 the LIFT DIVERGES (0.373/0.523/0.715) -- root cause the lift coupling. My
        #: registered A1/A2/A3 presupposed a non-convergence and so used the wrong instrument.
        #: alpha=0 separates fold-zone lift sensitivity (discipline #4: no grid-convergence claims
        #: at M~0.79, dcl/dM 6-10) from a discretisation error, because at alpha=0 there is no lift
        #: to be sensitive. It is transonic at M0.80 -- GS1.2b measured shocks there.
        arows = []
        print("  alpha=0 control (the GS1.2b discriminator):")
        for tag, n_st, h_n, w, fac in LEGS_ALL:
            row = run_leg(tag, n_st, h_n, w, fac, 1, alpha=0.0)
            arows.append(row)
            print(f"    {tag:15} conv={str(row['converged']):5} |R|={row['res_final']:.2e} "
                  f"lim/flr={row['n_limited']}/{row['n_floored']} tets={row['n_tets']:>6} "
                  f"M_max={row['m_max']:.4f} shock={row['has_shock']} "
                  f"x_shock={row['x_shock'] if row['x_shock'] is None else round(row['x_shock'],4)}"
                  f" cl={row['cl']} ({row['wall_s']:.0f}s)", flush=True)
        assert any(r["m_max"] > 1.0 for r in arows), (
            "alpha=0 at M0.80 shows NO supersonic cell on any leg -- the control has no shock to "
            "converge, so it cannot discriminate. Do NOT carry GS1.2b's alpha=0 transonic premise "
            "across mesh families without this assert (the Q1 lesson, verbatim).")
        xs = [r["x_shock"] for r in arows if r["x_shock"] is not None]
        print(f"\n    x_shock spread alpha=0 : {min(xs):.4f} .. {max(xs):.4f}  "
              f"(range {max(xs)-min(xs):.4f})")
        cls = [abs(r["cl"]) for r in arows]
        print(f"    |cl|      spread alpha=0 : {min(cls):.4f} .. {max(cls):.4f}")
        _write_named(arows, "task3_alpha0_control.csv")
        return 0

    rows = [run_leg(*LEGS_ALL[0], 1)]
    r0 = rows[0]
    print(f"  R0  conv={str(r0['converged']):5} |R|={r0['res_final']:.2e} "
          f"lim/flr={r0['n_limited']}/{r0['n_floored']} M_max={r0['m_max']:.4f} "
          f"shock={r0['has_shock']} x_shock={r0['x_shock']} cl={r0['cl']} "
          f"({r0['wall_s']:.0f}s)  mode={r0['failure_mode'] or '-'}")
    _write(rows)

    ok = (r0["converged"] and r0["n_limited"] == 0 and r0["n_floored"] == 0
          and r0["has_shock"] and r0["x_shock"] is not None)
    print("\n=== premise (pre-registration §5.1, fail-fast) ===")
    if not ok:
        print("  ★ R0 is NOT a usable baseline: the item STOPS here per the registration.")
        print(f"     converged={r0['converged']} clamps={r0['n_limited']}/{r0['n_floored']} "
              f"has_shock={r0['has_shock']} x_shock={r0['x_shock']}")
        print("     A baseline that is not a solution cannot define the shock window, and a")
        print("     window guessed instead of measured would make every later leg unreadable.")
        return 1
    print(f"  ✓ R0 usable. Measured x_shock = {r0['x_shock']:.4f} -> the shock window is taken")
    print(f"    FROM this, not assumed.")
    if a.baseline_only:
        print("  (--baseline-only: stopping before the four legs, as asked)")
        return 0

    #: ★ the window is the REGISTERED [0.45, 0.75]; R0's measured x_shock 0.6051 falls in it, so
    #: it is confirmed by measurement rather than chosen after seeing the answer.
    win = (0.45, 0.75)
    assert win[0] < r0["x_shock"] < win[1], (
        f"measured x_shock {r0['x_shock']:.4f} is OUTSIDE the registered window {win} -- the "
        f"window must then be re-registered, not silently moved")
    #: (tag, n_stations, h_wall_normal, window, factor)
    LEGS = (("RG_global", 320, 0.002, None, 1.0),
            ("RS_shock_band", 160, 0.004, win, 3.0),
            ("RN_normal_only", 160, 0.002, None, 1.0))
    #: ★ TWO independent runs per leg: criterion C2 requires the same mode twice, or the
    #: classification is noise rather than a reading.
    for tag, n_st, h_n, w, fac in LEGS:
        for run_id in (1, 2):
            row = run_leg(tag, n_st, h_n, w, fac, run_id)
            rows.append(row)
            print(f"  {tag:15} run{run_id} conv={str(row['converged']):5} "
                  f"|R|={row['res_final']:.2e} lim/flr={row['n_limited']}/{row['n_floored']} "
                  f"tets={row['n_tets']:>6} M_max={row['m_max']:.4f} "
                  f"x_shock={row['x_shock'] if row['x_shock'] is None else round(row['x_shock'],4)} "
                  f"mode={row['failure_mode'] or '-'} ({row['wall_s']:.0f}s)", flush=True)
            if row["wall_s"] > LEG_GATE_S:
                print(f"    ★ leg gate {row['wall_s']:.0f}s > {LEG_GATE_S:.0f}s -- stopping, "
                      f"no budget added")
                _write(rows); return 1
    _write(rows)
    _read(rows)
    return 0


def _read(rows):
    by = {}
    for r in rows:
        by.setdefault(r["leg"], []).append(r)
    print("\n=== reading (bands fixed in the pre-registration §4) ===")
    print("  C1  every non-converged leg carries a MODE, not a bare conv=False:")
    c1 = True
    for leg, rs in by.items():
        for r in rs:
            if not r["converged"] and not r["failure_mode"]:
                c1 = False; print(f"      ★ {leg} run{r['run']} has no mode")
    print(f"      -> {'PASS' if c1 else 'FAIL'}")
    print("  C2  the same leg gives the same mode on two independent runs:")
    c2 = True
    for leg, rs in sorted(by.items()):
        modes = {r["failure_mode"] for r in rs}
        conv = {r["converged"] for r in rs}
        stable = len(modes) == 1 and len(conv) == 1
        c2 &= stable
        print(f"      {leg:15} converged={sorted(conv)} modes={sorted(modes) or ['-']}"
              f"  {'stable' if stable else '★ UNSTABLE'}")
    print(f"      -> {'PASS' if c2 else 'FAIL'}")
    rg = by.get("RG_global", [{}])[0]; rs_ = by.get("RS_shock_band", [{}])[0]
    print("\n  attribution (A1/A2/A3):")
    print(f"      RG converged={rg.get('converged')} mode={rg.get('failure_mode') or '-'}")
    print(f"      RS converged={rs_.get('converged')} mode={rs_.get('failure_mode') or '-'}")
    if rg.get("converged") and rs_.get("converged"):
        print("      ⇒ NEITHER refinement leg fails: the GS1.2 paradox does NOT reproduce on")
        print("        this mesh family. That is a finding about the FAMILY, not about the")
        print("        solver, and it is what makes it attributable -- see the verdict.")
    elif rs_.get("converged") and not rg.get("converged"):
        print("      ⇒ A2: the failure tracks GLOBAL refinement, not the shock band.")
    elif rg.get("failure_mode") == rs_.get("failure_mode"):
        print("      ⇒ A1: same mode both ways -- the failure tracks the shock-band resolution.")
    else:
        print("      ⇒ A3: different modes -- a MIXTURE; see the kill criterion.")


def _write_named(rows, name):
    path = os.path.join(HERE, "gate_results", name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {path}")


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
