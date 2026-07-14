"""Track B / B15 demo -- LS Newton transonic ramp + N5 freeze-selection.

The claim: the level-set TRANSONIC PICARD parks its top Mach levels on the
shock-position residual PLATEAU (the P4/B6/N5 soft mode) and burns its whole
outer budget there -- that plateau IS the 24.5/38.4 min of the M6-medium M0.84
workflow solve. Newton has no such soft mode. B15 gives the LS Newton the two
things it lacked: a per-side FROZEN upwind selection (N5) and a Mach-ramp
wrapper.

Parts (all cheap except part 3):
  1. GB15.2 -- the freeze cures a genuine LIMIT CYCLE. With the live per-side
     selection the LS Newton does not converge on NACA coarse M0.75: it parks
     in a period-6 cycle at |R| ~ 3e-7 with 0 limited/floored (a CLEAN stall =
     assignment churn). Arming the freeze converges it, WITHOUT moving gamma.
  2. GB15.3 -- the ramp beats the Picard, and `intermediate_tol` is free.
     NACA coarse M0.80 / alpha 1.25 (the B6 gate condition, whose same-mesh
     CONFORMING-Newton truth is shock 0.658 / cl_p 0.459 / M_max 1.408).
  3. GB15.4 (gated, PYFP3D_TRANSONIC_GATES=1) -- ONERA M6 MEDIUM M0.84, the
     case the committed Picard takes 2304.7 s (wake-free) to leave
     bounded-but-not-converged. The Picard baseline is NOT re-run: it is
     committed evidence (cases/demo/m6_medium_ls_workflow/results/summary.csv).

Run:
    NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
        python cases/demo/b15_ls_newton_ramp/run_demo.py
    (add PYFP3D_TRANSONIC_GATES=1 for part 3, ~10-40 min)
"""
import csv
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.solve.newton_ls import (
    B_NEWTON_M6_DEFAULTS,
    solve_multivalued_newton,
    solve_multivalued_newton_transonic,
)
from pyfp3d.solve.picard_ls import solve_multivalued_transonic
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
GATED = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"

NACA = ROOT / "cases" / "meshes" / "naca0012_2.5d" / "coarse.msh"
M6_WF = ROOT / "cases" / "meshes" / "onera_m6_wakefree" / "medium.msh"

# committed Picard reference for part 3 (m6_medium_ls_workflow/summary.csv)
PICARD_M6 = dict(wall_s=2304.7, cl_kj=0.27648, m_max=2.4549)

checks, rows = [], []


def check(name, ok, detail=""):
    checks.append({"check": name, "status": "PASS" if ok else "FAIL",
                   "detail": detail})
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def naca_mvop(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
                       direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ---------------------------------------------------------------------------
# Part 1 -- GB15.2: the freeze cures the limit cycle.
# ---------------------------------------------------------------------------
print("\n=== Part 1 (GB15.2): freeze vs the live limit cycle "
      "(NACA coarse M0.75) ===")
mesh = read_mesh(NACA)
common = dict(mesh=mesh, m_inf=0.75, alpha_deg=2.0, farfield="neumann",
              upwind_c=1.5, m_crit=0.95, n_seed=25, n_newton_max=60,
              tol_residual=1e-10)
live = solve_multivalued_newton(mvop=naca_mvop(mesh), **common)
frz = solve_multivalued_newton(mvop=naca_mvop(mesh), freeze_tol=1e-6, **common)

check("live Newton does NOT converge (the limit cycle)",
      not live["converged"],
      f"|R|={live['residual_history'][-1]:.2e} after {live['n_newton']} steps")
check("...and it is a CLEAN stall (0 limited/floored) = assignment churn",
      live["n_limited"] == 0 and live["n_floored"] == 0,
      f"lim/flr={live['n_limited']}/{live['n_floored']}")
check("frozen Newton CONVERGES", frz["converged"] and frz["froze"],
      f"|R|={frz['residual_history'][-1]:.2e} in {frz['n_newton']} steps "
      f"(reason={frz['accept_reason']})")
check("...in fewer steps than the live budget",
      frz["n_newton"] < live["n_newton"],
      f"{frz['n_newton']} vs {live['n_newton']}")
check("...and the freeze does NOT move the solution",
      abs(frz["gamma"] - live["gamma"]) < 1e-4,
      f"gamma {frz['gamma']:.6f} vs live {live['gamma']:.6f}")
check("...with no freeze reverts", frz["n_freeze_reverts"] == 0, "")
rows.append(dict(part="GB15.2", case="NACA coarse M0.75", method="live",
                 wall_s="", n_iter=live["n_newton"],
                 residual=f"{live['residual_history'][-1]:.3e}",
                 converged=int(live["converged"]), gamma=f"{live['gamma']:.6f}",
                 m_max=f"{np.sqrt(live['mach2_max']):.4f}"))
rows.append(dict(part="GB15.2", case="NACA coarse M0.75", method="frozen",
                 wall_s="", n_iter=frz["n_newton"],
                 residual=f"{frz['residual_history'][-1]:.3e}",
                 converged=int(frz["converged"]), gamma=f"{frz['gamma']:.6f}",
                 m_max=f"{np.sqrt(frz['mach2_max']):.4f}"))

fig, ax = plt.subplots(figsize=(7, 4.4))
ax.semilogy(live["residual_history"], "o-", ms=3, lw=1.1,
            label=f"live selection ({live['n_newton']} steps, NOT converged)")
ax.semilogy(frz["residual_history"], "s-", ms=3, lw=1.1,
            label=f"frozen selection ({frz['n_newton']} steps, converged)")
ax.axhline(1e-10, color="k", ls=":", lw=0.9, label="tol_residual")
ax.set_xlabel("Newton step")
ax.set_ylabel(r"$\|R\|_\infty$")
ax.set_title("B15/GB15.2 — the frozen selection removes the assignment\n"
             "limit cycle (NACA coarse, M∞0.75)")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "gb152_freeze_limit_cycle.png", dpi=130)
plt.close(fig)

# ---------------------------------------------------------------------------
# Part 2 -- GB15.3: the ramp vs the Picard, and intermediate_tol.
# ---------------------------------------------------------------------------
print("\n=== Part 2 (GB15.3): Newton ramp vs Picard "
      "(NACA coarse M0.80, alpha 1.25 = the B6 gate) ===")
RAMP = dict(mesh=mesh, m_target=0.80, alpha_deg=1.25, farfield="neumann",
            m_start=0.60, dm=0.05, dm_min=0.005, freeze_tol=1e-6,
            n_seed=30, n_newton_max=60, tol_residual=1e-10,
            direct_refactor_every=1000)

t = time.perf_counter()
pic = solve_multivalued_transonic(
    mvop=naca_mvop(mesh), mesh=mesh, m_target=0.80, alpha_deg=1.25,
    farfield="neumann", m_start=0.60, dm=0.05, n_outer_seed=120,
    n_outer_level=400, tol_residual=1e-7, direct_refactor_every=1000)
t_pic = time.perf_counter() - t
pic_outers = sum(l["n_outer"] for l in pic["levels"])
pic_conv = sum(1 for l in pic["levels"] if l["converged"])

t = time.perf_counter()
nwt = solve_multivalued_newton_transonic(mvop=naca_mvop(mesh),
                                         intermediate_tol=None, **RAMP)
t_nwt = time.perf_counter() - t

t = time.perf_counter()
nwi = solve_multivalued_newton_transonic(mvop=naca_mvop(mesh),
                                         intermediate_tol=1e-5, **RAMP)
t_nwi = time.perf_counter() - t

check("PICARD leaves top Mach levels UNCONVERGED (the plateau)",
      pic_conv < len(pic["levels"]),
      f"{pic_conv}/{len(pic['levels'])} levels, |R|={pic['residual_norm']:.2e}, "
      f"{pic_outers} outers, {t_pic:.1f}s")
check("NEWTON ramp reaches the target STRICTLY",
      nwt["target_reached"] and nwt["converged"]
      and nwt["residual_history"][-1] < 1e-10,
      f"|R|={nwt['residual_history'][-1]:.2e}, {t_nwt:.1f}s")
check("...0 limited / 0 floored",
      nwt["n_limited"] == 0 and nwt["n_floored"] == 0, "")
check("...and it is FASTER than the (stalled) Picard",
      t_nwt < t_pic, f"{t_nwt:.1f}s vs {t_pic:.1f}s = {t_pic/t_nwt:.1f}x")
check("M_max consistent with the same-mesh CONFORMING-Newton truth (1.408)",
      abs(np.sqrt(nwt["mach2_max"]) - 1.408) < 0.05,
      f"M_max={np.sqrt(nwt['mach2_max']):.4f} vs 1.408")
check("intermediate_tol does NOT move the final solution",
      abs(nwi["gamma"] - nwt["gamma"]) < 1e-6,
      f"gamma {nwi['gamma']:.6f} vs strict {nwt['gamma']:.6f}")
check("...and it is faster still",
      t_nwi < t_nwt,
      f"{t_nwi:.1f}s vs {t_nwt:.1f}s; newton steps "
      f"{sum(l['n_newton'] for l in nwi['levels'])} vs "
      f"{sum(l['n_newton'] for l in nwt['levels'])}")

for tag, res, wall in (("picard", pic, t_pic), ("newton", nwt, t_nwt),
                       ("newton+interm_tol", nwi, t_nwi)):
    n_it = (sum(l["n_outer"] for l in res["levels"]) if tag == "picard"
            else sum(l["n_newton"] for l in res["levels"]))
    rnorm = (res["residual_norm"] if tag == "picard"
             else res["residual_history"][-1])
    rows.append(dict(part="GB15.3", case="NACA coarse M0.80 a1.25", method=tag,
                     wall_s=f"{wall:.1f}", n_iter=n_it,
                     residual=f"{rnorm:.3e}",
                     converged=int(res.get("target_reached",
                                           res["converged"])),
                     gamma=f"{res['gamma']:.6f}",
                     m_max=f"{np.sqrt(res['mach2_max']):.4f}"))

fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.bar(["Picard\n(stalled)", "Newton\nramp", "Newton\n+interm_tol"],
       [t_pic, t_nwt, t_nwi],
       color=["#c44", "#48a", "#4a8"])
for i, (v, r) in enumerate(((t_pic, pic["residual_norm"]),
                            (t_nwt, nwt["residual_history"][-1]),
                            (t_nwi, nwi["residual_history"][-1]))):
    ax.text(i, v + 0.6, f"{v:.1f}s\n|R|={r:.0e}", ha="center", fontsize=8)
ax.set_ylabel("wall clock (s)")
ax.set_title("B15/GB15.3 — the Newton ramp removes the Picard shock plateau\n"
             "(NACA coarse, M∞0.80, α1.25 — the B6 gate condition)")
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(OUT / "gb153_ramp_vs_picard.png", dpi=130)
plt.close(fig)

# ---------------------------------------------------------------------------
# Part 3 -- GB15.4: ONERA M6 medium M0.84 (gated).
# ---------------------------------------------------------------------------
print("\n=== Part 3 (GB15.4): ONERA M6 MEDIUM M0.84 "
      f"({'RUNNING' if GATED else 'SKIPPED — set PYFP3D_TRANSONIC_GATES=1'}) ===")
if GATED and M6_WF.exists():
    # The TE polyline MUST come from the authoritative geometry: CutElementMap
    # finds TE nodes by matching the polyline onto WALL NODES, and the M6 TE
    # endpoints are exact wall nodes (M2: tip corner 0.0, junction 1.5e-9). A
    # hand-rolled x_te that is off by even ~2e-4 in x matches NOTHING -> 0 TE
    # nodes -> no Kutta -> the solve blows up (measured: 340k limited, NaN).
    from pyfp3d.meshgen.wing3d import B_SEMI, x_te
    ALPHA = 3.06

    m6 = read_mesh(M6_WF)
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(m6.nodes, m6.elements, wls,
                       wall_nodes=np.unique(m6.boundary_faces["wall"]))
    mv6 = MultivaluedOperator(m6.nodes, m6.elements, cm, levelset=wls)

    # B_NEWTON_M6_DEFAULTS -- the MEASURED recipe. freeze_tol=1e-3 sits ABOVE
    # the churn floor (which RISES with Mach: <1e-6 / 8.6e-6 / 2.7e-4 at M0.60 /
    # 0.65 / 0.70) and freeze_max_clamped=8 tolerates the lone floored cell that
    # would otherwise block the freeze at ANY freeze_tol (the P9/G9.1 wall).
    t = time.perf_counter()
    r6 = solve_multivalued_newton_transonic(
        mvop=mv6, mesh=m6, m_target=0.84, alpha_deg=ALPHA, farfield="neumann",
        n_seed=40, n_newton_max=80, tol_residual=1e-10,
        **B_NEWTON_M6_DEFAULTS)
    t6 = time.perf_counter() - t
    lv = r6["levels"]

    check("M6 medium M0.84: the Newton ramp REACHES the target",
          r6["target_reached"],
          f"m_final={r6['m_final']}, {t6:.0f}s ({t6/60:.1f} min)")
    check("...every ramp level converged",
          all(l["converged"] for l in lv),
          f"{sum(l['converged'] for l in lv)}/{len(lv)} levels")
    check("...the freeze armed at every level, with no reverts",
          all(l["froze"] for l in lv)
          and sum(l["n_freeze_reverts"] for l in lv) == 0,
          f"froze={[int(l['froze']) for l in lv]}")
    check("...FASTER than the committed Picard (2304.7 s)",
          t6 < PICARD_M6["wall_s"],
          f"{t6:.0f}s vs {PICARD_M6['wall_s']}s = "
          f"{PICARD_M6['wall_s']/t6:.2f}x")
    check("...and the PLATEAU is gone: |R| ~1e-11, not the Picard 1e-5..1e-4",
          r6["residual_history"][-1] < 1e-9,
          f"|R|={r6['residual_history'][-1]:.2e} "
          f"(Picard parks at 1e-5..1e-4)")
    check("...M_max agrees with the committed Picard (2.4549) within 5%",
          abs(np.sqrt(r6["mach2_max"]) - PICARD_M6["m_max"])
          / PICARD_M6["m_max"] < 0.05,
          f"M_max={np.sqrt(r6['mach2_max']):.4f} vs {PICARD_M6['m_max']:.4f}")
    # HONESTY: most levels accept at the assignment-discontinuity floor
    # ("assignment_cycle"), i.e. the FROZEN system is converged to ~1e-11 and
    # the LIVE residual has stopped improving across refreshes. That is the N5
    # semantics the conforming path also uses -- it is NOT a claim that the live
    # residual is below 1e-10. Report it rather than paper over it.
    print(f"    [note] accept reasons: {[l['accept_reason'] for l in lv]}")
    print(f"    [note] final clamped cells (live): "
          f"{r6['n_limited']}/{r6['n_floored']} of {len(m6.elements)} "
          f"(committed Picard: <=3)")
    # ★ RECORD the clamped cells as a committed artifact. `freeze_max_clamped>0`
    # RELAXES the convergence semantics: the assignment_cycle / refresh_budget
    # accept routes do not re-check the clamp count, so this converged=True state
    # CARRIES clamped cells. They do NOT "clear themselves" -- they persist at
    # every level from M0.70 up. Anyone quoting the M6 number must see this.
    per_level_clamped = ", ".join(
        "M{:.2f}:{}/{}".format(l["m_inf"], l["n_limited"], l["n_floored"])
        for l in lv)
    accept_routes = sorted({l["accept_reason"] for l in lv})
    check("...HONEST: the converged state CARRIES clamped cells (semantics "
          "relaxed by freeze_max_clamped>0, NOT a 0-clamped solution)",
          True,
          f"live lim/flr={r6['n_limited']}/{r6['n_floored']} of "
          f"{len(m6.elements)} tets (Picard <=3); per level "
          f"[{per_level_clamped}]; accept routes {accept_routes}")
    for l in lv:
        rows.append(dict(part="GB15.4", case=f"  level M{l['m_inf']:.4f}",
                         method=f"newton ({l['accept_reason']})",
                         wall_s="", n_iter=l["n_newton"],
                         residual=f"{l['residual_norm']:.3e}",
                         converged=int(l["converged"]),
                         gamma=f"{l['gamma']:.6f}",
                         m_max=f"{l['mach_max']:.4f}",
                         clamped=f"{l['n_limited']}/{l['n_floored']}"))
    rows.append(dict(part="GB15.4", case="M6 medium M0.84 wake-free",
                     method="picard (committed)",
                     wall_s=f"{PICARD_M6['wall_s']:.1f}", n_iter="",
                     residual="1e-5..1e-4 (plateau)", converged=0,
                     gamma="", m_max=f"{PICARD_M6['m_max']:.4f}"))
    rows.append(dict(part="GB15.4", case="M6 medium M0.84 wake-free",
                     method="newton ramp",
                     wall_s=f"{t6:.1f}",
                     n_iter=sum(l["n_newton"] for l in lv),
                     residual=f"{r6['residual_history'][-1]:.3e}",
                     converged=int(r6["target_reached"]),
                     gamma=f"{r6['gamma']:.6f}",
                     m_max=f"{np.sqrt(r6['mach2_max']):.4f}"))
    np.savez(OUT / "m6_medium_newton.npz", phi_ext=r6["phi_ext"], wall_s=t6)

# ---------------------------------------------------------------------------
with open(OUT / "summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["part", "case", "method", "wall_s",
                                      "n_iter", "residual", "converged",
                                      "gamma", "m_max", "clamped"],
                       restval="")
    w.writeheader()
    w.writerows(rows)
with open(OUT / "checks.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["check", "status", "detail"])
    w.writeheader()
    w.writerows(checks)

n_pass = sum(1 for c in checks if c["status"] == "PASS")
print(f"\n=== {n_pass}/{len(checks)} checks PASS ===")
print(f"artifacts -> {OUT}")
