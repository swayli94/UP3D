"""Outstanding item 2: the coarse 2.0 % -- is it the flat->round CAP switch, or a real pipeline gap?

Pre-registered in phases/p3/docs/dev_phase_three/20260816-1800-coarse-pipeline-gap-prereg.md. Phase 3's last debt.

★★★ The entry check relocated the debt before any code was written. "2.0 % against the P14 coarse anchor"
had no citable referent in the round file; g82_anchor_check.csv is MEDIUM (run_g82_anchor_check.py:51 calls
_m6_case("medium")). The real referent is in run_m3_budget.py itself:
P14_ANCHOR["coarse"] = (0.262778, 0.268813), "GV5.3's P14_ANCHOR verbatim". Against the triage's HEAD coarse
reading of 0.268115 that is +2.03 % -- the 2.0 %.

And the two "pipelines" turn out to be ONE recipe: GV5.3 used the P14 recipe verbatim (M0.70 probe seed, the
NEWTON_M6_RECIPE ramp, pressure Kutta, NO taper), which is exactly the triage's production_solve. What
differs is that coarse.msh was REGENERATED on 2026-08-04 when the level names flipped flat -> round. So this
round runs the same leg on coarse_flat.msh and asks whether it lands on the anchor.

★ Prediction, written in the registration before this file existed: F-MESH. It fails if the flat leg does NOT
land near 0.262778, which would mean a genuine cross-pipeline inconsistency remains.

Outputs (TRACKED): bench/gate_results/task3_coarse_pipeline_gap.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.metrics import precompute_element_geometry            # noqa: E402
from pyfp3d.mesh.reader import read_mesh                               # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                              # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared              # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,              # noqa: E402
                                 wall_force_coefficients)
import run_m3_budget as MB                                             # noqa: E402
import run_task3_m6_triage as TR                                       # noqa: E402
from failure_modes import classify_failure                      # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_coarse_pipeline_gap.csv")
#: the referent, read from run_m3_budget (GV5.3's P14_ANCHOR verbatim), NOT recalled
ANCHOR_CL_P, ANCHOR_CL_KJ = MB.P14_ANCHOR["coarse"]
#: G-R: the triage's own HEAD coarse reading, from its round file
TRIAGE_COARSE_CL_P = 0.268115
BAND = 0.005                       #: +-0.5 %, the triage's own G-R band -- no new number
CELLS_MAX_DIFF = 0.20
LEG_GATE_S, TOTAL_GATE_S = 900.0, 1800.0
LEGS = (("R_round", "coarse.msh"), ("F_flat", "coarse_flat.msh"))


def leg(mesh_name):
    """The triage's production_solve(seed=0, ramp=True) VERBATIM -- imported, not re-typed."""
    p = os.path.join(REPO, "cases", "meshes", "onera_m6", mesh_name)
    if not os.path.exists(p):
        return None, dict(leg=mesh_name, error="mesh missing")
    mc, wc = cut_wake(read_mesh(p))
    wall = mc.boundary_faces["wall"]
    #: G-M's registered metric is the raw cell count. The tip/non-tip split is REPORTED ALONGSIDE
    #: it as extra evidence -- it does NOT rescue the threshold, which stays where it was registered.
    zc = mc.nodes[wall].mean(axis=1)[:, 2]
    tip = zc > 0.95 * B_SEMI
    info = dict(mesh=mesh_name, n_tet=int(len(mc.elements)), n_node=int(len(mc.nodes)),
                n_wall=int(len(wall)), n_wall_tip=int(tip.sum()),
                n_wall_nontip=int((~tip).sum()))
    t0 = time.perf_counter()
    r, _r0 = TR.production_solve(mc, wc, 0, True)
    wall = time.perf_counter() - t0
    hist = list(r.get("residual_history", []))
    conv = bool(r.get("converged"))
    try:
        mode, ev, _d, _v = classify_failure(
            np.asarray(hist, dtype=float),
            np.asarray(r.get("clamp_history", []), dtype=float),
            np.asarray(r.get("F_history", []), dtype=float),
            int(r.get("n_gmres_stalled", 0) or 0), str(r.get("accept_reason")),
            int(r.get("n_limited", 0) or 0), int(r.get("n_floored", 0) or 0))
    except Exception as exc:                                           # noqa: BLE001
        mode, ev = "classifier_failed", f"{type(exc).__name__}: {exc}"
    row = dict(info, converged=conv, res_final=(hist[-1] if hist else None),
               n_newton=len(hist), n_limited=r.get("n_limited"), n_floored=r.get("n_floored"),
               failure_mode=(None if conv else mode), failure_evidence=(None if conv else ev),
               solve_s=round(wall, 1))
    if conv:
        phi = np.asarray(r["phi"])
        gam = np.atleast_1d(np.asarray(r["gamma"]))
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                    alpha_deg=MB.ALPHA, s_ref=s_ref, m_inf=MB.M_INF)
        Bg, _V = precompute_element_geometry(mc.nodes, mc.elements)
        gg = np.einsum("eaj,ea->ej", Bg, phi[mc.elements])
        m2 = mach_number_squared(np.einsum("ej,ej->e", gg, gg), MB.M_INF)
        o = np.argsort(wc.station_z)
        row.update(cl_p=float(f["cl"]),
                   cl_kj=float(cl_kj_3d(gam[o], wc.station_z[o], s_ref, B_SEMI)),
                   m_max=float(np.sqrt(m2.max())))
    return row, None


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}\n")
    print(f"★ referent (read from run_m3_budget.P14_ANCHOR, GV5.3 verbatim): "
          f"coarse cl_p {ANCHOR_CL_P:.6f} / cl_KJ {ANCHOR_CL_KJ:.6f}")
    print(f"★ G-R target (the triage's HEAD coarse reading): cl_p {TRIAGE_COARSE_CL_P:.6f}")
    print(f"★ the recipe is run_task3_m6_triage.production_solve(seed=0, ramp=True), IMPORTED "
          f"verbatim -- no taper on either side, as on the anchor's side\n")

    rows, t_all = [], time.perf_counter()
    for tag, mesh_name in LEGS:
        if time.perf_counter() - t_all > TOTAL_GATE_S:
            print(f"★ total gate exceeded -- {tag} NOT run (kill clause 4)")
            break
        print(f"=== {tag} ({mesh_name}) ===")
        row, err = leg(mesh_name)
        if row is None:
            print(f"  ★ {err}")
            rows.append(err)
            continue
        row["leg"] = tag
        rows.append(row)
        _write(rows)
        print(f"  {row['n_tet']} tets / {row['n_node']} nodes / {row['n_wall']} wall tris "
              f"(tip {row['n_wall_tip']} / non-tip {row['n_wall_nontip']})   "
              f"conv={row['converged']} |R|={(row['res_final'] or float('nan')):.3e} "
              f"lim={row['n_limited']} flr={row['n_floored']} ({row['solve_s']:.0f}s)")
        if row["converged"]:
            print(f"  cl_p {row['cl_p']:.6f}   cl_KJ {row['cl_kj']:.6f}   "
                  f"M_max {row['m_max']:.4f}")
        else:
            print(f"  ★ mode={row['failure_mode']}   {row['failure_evidence']}")
        if row["solve_s"] > LEG_GATE_S:
            print(f"  ★ leg gate exceeded -- stopping here (kill clause 3)")
            break
    _write(rows)
    return report(rows)


def report(rows):
    by = {r.get("leg"): r for r in rows if r.get("leg")}
    R, F = by.get("R_round"), by.get("F_flat")

    print("\n=== G-M: are the two meshes different only in the cap? ===")
    if R and F and R.get("n_tet") and F.get("n_tet"):
        d = abs(R["n_tet"] - F["n_tet"]) / max(R["n_tet"], F["n_tet"])
        print(f"  round {R['n_tet']} tets / flat {F['n_tet']} tets -> {100 * d:.1f} % apart "
              f"(REGISTERED limit {100 * CELLS_MAX_DIFF:.0f} %) -> "
              f"{'PASS' if d <= CELLS_MAX_DIFF else '★ FAIL'}")
        nt = abs(R["n_wall_nontip"] - F["n_wall_nontip"]) / max(R["n_wall_nontip"],
                                                                F["n_wall_nontip"])
        print(f"  ★ where the difference lives: wall tris tip {R['n_wall_tip']} vs "
              f"{F['n_wall_tip']}, NON-tip {R['n_wall_nontip']} vs {F['n_wall_nontip']} "
              f"= {100 * nt:.2f} % apart")
        print("    (reported as extra evidence; it does NOT move the registered threshold)")
        if d > CELLS_MAX_DIFF:
            print("  ★ beyond the limit: they differ by MORE than the cap ⇒ the verdict below is")
            print("    degraded to RECORDED (kill clause 2).")
    else:
        print("  (one leg missing)")

    print("\n=== G-R: does the round leg reproduce the triage's HEAD coarse reading? ===")
    if not (R and R.get("converged")):
        print(f"  ★ round leg not converged/missing -> G-R UNDEFINED "
              f"(mode={R.get('failure_mode') if R else 'missing'})")
    else:
        rel = abs(R["cl_p"] - TRIAGE_COARSE_CL_P) / TRIAGE_COARSE_CL_P
        print(f"  cl_p {R['cl_p']:.6f} vs {TRIAGE_COARSE_CL_P:.6f} = {100 * rel:+.3f} % -> "
              f"{'PASS' if rel <= BAND else '★ FAIL'}")
        if rel > BAND:
            print("  -> STOP: not the same leg (kill clause 1).")
            return 1

    print("\n=== the verdict (binding = the FLAT leg against the anchor) ===")
    if not (F and F.get("converged")):
        print(f"  -> ★ F-UNDEF  the flat leg did not converge "
              f"(mode={F.get('failure_mode') if F else 'missing'}) ⇒ UNDEFINED.")
        print("     ★ Non-convergence is NOT evidence for either side.")
        if F and F.get("failure_evidence"):
            print(f"     evidence: {F['failure_evidence']}")
        return 0
    rel_f = (F["cl_p"] - ANCHOR_CL_P) / ANCHOR_CL_P
    rel_r = ((R["cl_p"] - ANCHOR_CL_P) / ANCHOR_CL_P) if (R and R.get("converged")) else None
    print(f"  FLAT  cl_p {F['cl_p']:.6f} vs anchor {ANCHOR_CL_P:.6f} = {100 * rel_f:+.3f} %")
    if rel_r is not None:
        print(f"  ROUND cl_p {R['cl_p']:.6f} vs anchor {ANCHOR_CL_P:.6f} = {100 * rel_r:+.3f} %")
    print(f"  (cl_KJ: flat {F.get('cl_kj', float('nan')):.6f} vs anchor {ANCHOR_CL_KJ:.6f} = "
          f"{100 * (F.get('cl_kj', float('nan')) - ANCHOR_CL_KJ) / ANCHOR_CL_KJ:+.3f} %)")
    if abs(rel_f) <= BAND:
        print(f"  -> ★★ F-MESH  the flat leg lands on the anchor within {100 * BAND:.1f} % ⇒ the")
        print("     2.0 % IS the flat->round cap switch of coarse.msh, not a pipeline inconsistency.")
        print("     ⇒ P14_ANCHOR is a FLAT-CAP-ERA number and must be labelled as such.")
    else:
        print("  -> ★★★ F-OTHER  the flat leg does NOT land on the anchor ⇒ MY PREDICTION IS WRONG:")
        print("     the cap does not explain the gap, so a genuine cross-pipeline inconsistency")
        print("     remains and has to be chased (next round: diff recipe and post-processing).")
    return 0


def _write(rows):
    if not rows:
        return
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
