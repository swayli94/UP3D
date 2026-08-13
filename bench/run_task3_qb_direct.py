"""Outstanding item 1: the Q-B DIRECT arm -- does the SYSTEM have several solutions without the ramp?

Pre-registered in docs/dev_phase_three/20260816-1400-qb-direct-arm-prereg.md. Ruling D7 keeps this in
phase 3.

Arm A (the ramp, three seeds) is answered: 0.0011 % cl_p spread on the ramp path. Arm B asks whether the
RAMP is what pins that -- because the project's own note is that a ramp and a seed are two implementations of
one function, so "reproducible through the ramp" is not "the system is unique".

★★★ The triage script's own arm B does NOT apply the production tip taper, and its committed CSVs are the
configuration known to die (cl_p 0.2493, all legs non-converged). Running it verbatim would give a
non-converged direct arm because the TAPER is missing rather than the RAMP -- two causes in one bare
conv=False. So this arm is built on run_m3_budget.solve's production leg (which HAS the taper) with the ramp
step removed and nothing else changed, and G-S asserts by source inspection that the replication is faithful.

★★ QB-UNDEF (fewer than two legs converging) is the most likely outcome and is a first-class result; its
value is the per-leg failure classification. The forbidden sentence is carried over verbatim from the triage
registration: NON-CONVERGENCE IS NOT EVIDENCE OF UNIQUENESS.

Outputs (TRACKED): bench/gate_results/task3_qb_direct.csv
"""

import csv
import inspect
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
from pyfp3d.post.section_cut import section_cp_curve                   # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,              # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import solve_newton_lifting                   # noqa: E402
import run_m3_budget as MB                                             # noqa: E402
from run_le14_common_root import classify_failure                      # noqa: E402
from tests.test_p8_newton import NEWTON_M6_RECIPE                      # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_qb_direct.csv")
#: G-P: arm A's reference, one file per seed (the triage-redo round's own artifacts)
REF = {0: "m3_budget_head_medium.csv", 5: "m3_budget_head_medium_seed5.csv",
       12: "m3_budget_head_medium_seed12.csv"}
SEEDS = (0, 5, 12)
MESH = os.path.join(REPO, "cases", "meshes", "onera_m6", "medium.msh")
PINS, UNIQUE = 10.0, 2.0
LEG_GATE_S, TOTAL_GATE_S = 1200.0, 3600.0


def gs_source_guard():
    """★★ G-S: assert the direct arm replicates run_m3_budget.solve's pressure branch.

    Two arms are only comparable if they share everything but the ramp, so the lines this file
    replicates are checked against that function's source rather than trusted from memory.
    """
    src = inspect.getsource(MB.solve)
    need = [
        'tip_taper=tip_taper_factors(wc.station_z, B_SEMI,',
        'kw["newton_kw"] = dict(kw["newton_kw"], entropy_correction=entropy)',
        'r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,',
        'kw["newton_kw"].update(kutta_estimator="pressure", phi_init=r0["phi"],',
        'return solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)',
    ]
    missing = [n for n in need if n not in src]
    assert not missing, f"G-S: run_m3_budget.solve no longer contains {missing}"
    print("★ G-S: the pressure branch of run_m3_budget.solve is unchanged -- the direct arm below "
          "replicates it and removes ONLY the ramp (the last line above).")


def direct_leg(mc, wc, seed):
    """run_m3_budget.solve's production pressure leg with the RAMP removed."""
    from pyfp3d.constraints.wake import tip_taper_factors
    kw = dict(NEWTON_M6_RECIPE)
    kw["newton_kw"] = dict(kw["newton_kw"],
                           tip_taper=tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                                                       0.05 * B_SEMI))
    for k, v in MB.M6_NEWTON_KW.items():                   # the recipe drift guard, verbatim
        assert kw["newton_kw"][k] == v, f"recipe drift: newton_kw[{k}]"
    kw["newton_kw"] = dict(kw["newton_kw"], entropy_correction=True)
    r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=MB.ALPHA, entropy_correction=True,
                              **dict(MB.M6_NEWTON_KW, n_picard_seed=seed))
    kw["newton_kw"].update(kutta_estimator="pressure", phi_init=r0["phi"],
                           gamma_init=r0["gamma"], n_picard_seed=0)
    #: ★ the ONLY difference from run_m3_budget.solve: lifting at M_INF instead of the ramp
    r = solve_newton_lifting(mc, wc, m_inf=MB.M_INF, alpha_deg=MB.ALPHA, **dict(kw["newton_kw"]))
    return r, r0


def metrics(mc, wc, r, s_ref, exp):
    phi = np.asarray(r["phi"])
    gam = np.atleast_1d(np.asarray(r["gamma"]))
    Bg, _V = precompute_element_geometry(mc.nodes, mc.elements)
    gg = np.einsum("eaj,ea->ej", Bg, phi[mc.elements])
    m2 = mach_number_squared(np.einsum("ej,ej->e", gg, gg), MB.M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                alpha_deg=MB.ALPHA, s_ref=s_ref, m_inf=MB.M_INF)
    o = np.argsort(wc.station_z)
    ss = nn = 0.0
    for eta in MB.ETAS[:MB.N_UNMASKED]:
        try:
            c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=MB.M_INF)
        except Exception:                                              # noqa: BLE001
            continue
        got = MB.band_rms({eta: c}, exp, eta).get("LE_upper")
        if got:
            ss += got[0]
            nn += got[1]
    return dict(cl_p=float(f["cl"]),
                cl_kj=float(cl_kj_3d(gam[o], wc.station_z[o], s_ref, B_SEMI)),
                m_max=float(np.sqrt(m2.max())), m1_max=r.get("m1_max"),
                sigma_min=r.get("sigma_min"),
                band_LE_upper=(float(np.sqrt(ss / nn)) if nn else None))


def spread(vals):
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return None
    return 100.0 * (max(vals) - min(vals)) / abs(np.mean(vals))


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}  (RECORDED only -- wall clock is not a criterion here)\n")
    gs_source_guard()

    #: --- G-P: arm A's reference, read from its own committed files -------------------------------
    print("\n★ G-P: arm A (RAMP) reference, one file per seed:")
    ref = {}
    for sd, fn in REF.items():
        p = os.path.join(HERE, "gate_results", fn)
        if not os.path.exists(p):
            print(f"  ★ missing {fn} -- STOP (kill clause 5; no substitute allowed)")
            return 1
        row = next(iter(csv.DictReader(open(p))))
        ref[sd] = dict(cl_p=float(row["cl_p"]), band_LE_upper=float(row["band_LE_upper"]),
                       converged=row["converged"])
        print(f"  seed {sd:2}: cl_p {ref[sd]['cl_p']:.6f}  LE_upper {ref[sd]['band_LE_upper']:.6f}"
              f"  conv {ref[sd]['converged']}  [{fn}]")
    ref_cl = spread([ref[s]["cl_p"] for s in SEEDS])
    ref_le = spread([ref[s]["band_LE_upper"] for s in SEEDS])
    print(f"  ⇒ RAMP-arm spread: cl_p {ref_cl:.4f} %   LE_upper {ref_le:.4f} %")

    if not os.path.exists(MESH):
        print(f"★ mesh missing: {MESH} -- regenerate via generate_onera_m6.py")
        return 1
    exp = MB.parse_experiment()
    mc, wc = cut_wake(read_mesh(MESH))
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    print(f"\n=== the DIRECT arm (no Mach ramp), medium {len(mc.elements)} tets, "
          f"M{MB.M_INF}/alpha {MB.ALPHA} ===")

    rows, t_all = [], time.perf_counter()
    for sd in SEEDS:
        if time.perf_counter() - t_all > TOTAL_GATE_S:
            print(f"  ★ total gate exceeded -- seed {sd} NOT run (kill clause 4)")
            break
        t0 = time.perf_counter()
        try:
            r, r0 = direct_leg(mc, wc, sd)
        except Exception as exc:                                       # noqa: BLE001
            wall = time.perf_counter() - t0
            print(f"  seed {sd}: RAISED {type(exc).__name__}: {exc}  ({wall:.0f}s)")
            rows.append(dict(arm="B_direct", seed=sd, converged=False,
                             error=f"{type(exc).__name__}: {exc}", solve_s=round(wall, 1)))
            continue
        wall = time.perf_counter() - t0
        hist = list(r.get("residual_history", []))
        conv = bool(r.get("converged"))
        #: ★★★ TWO defects fixed here, and both are the SAME family the project keeps recording:
        #: "read signatures and import paths, do not recall them", and "guard REPORT code the way
        #: solver code is guarded". The first call passed PYTHON LISTS (classify_failure does
        #: `tail > 0`, so it needs arrays -- the caller at run_le14_common_root.py:226-228 builds
        #: them with np.asarray) and unpacked TWO returns where the function returns FOUR
        #: (mode, ev, d10, revisits, per its own caller at line 260). It raised AFTER seed 0's
        #: direct solve had finished and BEFORE rows.append, so the solve was lost -- the same
        #: hazard as the 40-minute solve thrown away by a float(None) in a reporting line.
        #: ⇒ the classification is now (a) called with the right types and arity and (b) wrapped,
        #: so a reporting-layer error can never again destroy a completed solve.
        try:
            mode, ev, _d10, _rev = classify_failure(
                np.asarray(hist, dtype=float),
                np.asarray(r.get("clamp_history", []), dtype=float),
                np.asarray(r.get("F_history", []), dtype=float),
                int(r.get("n_gmres_stalled", 0) or 0), str(r.get("accept_reason")),
                int(r.get("n_limited", 0) or 0), int(r.get("n_floored", 0) or 0))
        except Exception as exc:                                       # noqa: BLE001
            mode, ev = "classifier_failed", f"{type(exc).__name__}: {exc}"
            print(f"    ★ classifier raised ({exc}) -- recorded, the leg's numbers are kept")
        row = dict(arm="B_direct", seed=sd, converged=conv,
                   res_final=(hist[-1] if hist else None), n_newton=len(hist),
                   n_limited=r.get("n_limited"), n_floored=r.get("n_floored"),
                   accept_reason=r.get("accept_reason"), failure_mode=(None if conv else mode),
                   failure_evidence=(None if conv else ev), solve_s=round(wall, 1),
                   over_leg_gate=bool(wall > LEG_GATE_S),
                   probe_cl=float(np.atleast_1d(np.asarray(r0["gamma"])).sum()))
        rows.append(row)          #: ★ CACHE BEFORE YOU REPORT -- the row exists before metrics run
        _write(rows)
        if conv:
            try:
                row.update(metrics(mc, wc, r, s_ref, exp))
                _write(rows)
            except Exception as exc:                                   # noqa: BLE001
                row["metrics_error"] = f"{type(exc).__name__}: {exc}"
                print(f"    ★ metrics raised ({exc}) -- the leg's convergence data is kept")
        print(f"  seed {sd:2}: conv={conv} |R|={(hist[-1] if hist else float('nan')):.3e} "
              f"n={len(hist)} lim={r.get('n_limited')} flr={r.get('n_floored')} "
              f"reason={r.get('accept_reason')} ({wall:.0f}s)"
              + ("" if conv else f"   ★ mode={mode}"))
        if conv:
            print(f"           cl_p {row['cl_p']:.6f}  cl_kj {row['cl_kj']:.6f}  "
                  f"LE_upper {row['band_LE_upper']:.6f}  M_max {row['m_max']:.4f}")
        _write(rows)
        if wall > LEG_GATE_S:
            print(f"  ★ leg gate {wall:.0f}s > {LEG_GATE_S:.0f}s -- stopping the arm here "
                  f"(kill clause 3; the ramp exists because direct solves are slow)")
            break

    _write(rows)
    return report(rows, ref_cl, ref_le)


def report(rows, ref_cl, ref_le):
    print("\n=== G-D: did the seed reach the solve? (probe circulation per leg) ===")
    pc = {r["seed"]: r.get("probe_cl") for r in rows if r.get("probe_cl") is not None}
    for sd, v in pc.items():
        print(f"  seed {sd:2}: probe gamma sum {v:.9f}")
    if len(pc) >= 2 and len({round(v, 12) for v in pc.values()}) == 1:
        print("  -> ★ G-D FAIL: every leg's probe state is identical ⇒ the seed did NOT reach the")
        print("     solve, so any spread below is meaningless. STOP (kill clause 2).")
        return 1
    dupes = [s for s in pc if list(pc.values()).count(pc[s]) > 1]
    if dupes:
        print(f"  ★ seeds {sorted(set(dupes))} share a probe state -- recorded as 'this seed may not "
              f"have entered the solve' (the redo round saw the same for seed 12), NOT counted as "
              f"independent samples")
    else:
        print("  -> PASS (all probe states distinct)")

    ok = [r for r in rows if r.get("converged")]
    print(f"\n=== QB verdict (binding = cl_p spread; reference = the RAMP arm's {ref_cl:.4f} %) ===")
    print(f"  converged direct legs: {len(ok)}/{len(rows)}")
    if len(ok) < 2:
        print("  -> ★★★ QB-UNDEF  fewer than two converged legs ⇒ UNDEFINED.")
        print("     ★★ FORBIDDEN SENTENCE (carried verbatim from the triage registration): a")
        print("     non-converging direct arm is NOT evidence of uniqueness. It says the direct")
        print("     path is hard, which is why the ramp exists -- nothing about how many solutions")
        print("     the system has.")
        print("\n  per-leg failure classification (this is the round's product):")
        for r in rows:
            if r.get("converged"):
                continue
            print(f"    seed {r['seed']:2}: mode={r.get('failure_mode')}  "
                  f"reason={r.get('accept_reason')}  |R|={r.get('res_final')}  "
                  f"lim/flr={r.get('n_limited')}/{r.get('n_floored')}  n={r.get('n_newton')}")
            if r.get("failure_evidence"):
                print(f"             evidence: {r['failure_evidence']}")
        return 0
    s_cl = spread([r["cl_p"] for r in ok])
    s_le = spread([r.get("band_LE_upper") for r in ok])
    ratio = s_cl / ref_cl if ref_cl else float("nan")
    print(f"  DIRECT-arm spread: cl_p {s_cl:.4f} %   LE_upper "
          f"{(s_le if s_le is not None else float('nan')):.4f} %")
    print(f"  ratio to the ramp arm: {ratio:.2f}x   (>= {PINS} -> ramp pins; <= {UNIQUE} -> unique)")
    if ratio >= PINS:
        print("  -> ★★ QB-RAMP-PINS  the ramp IS pinning the selection ⇒ arm A's reproducibility is")
        print("     NOT 'the system is unique'; there are solutions the ramp hides. This connects to")
        print("     the 2.5-D non-uniqueness result.")
    elif ratio <= UNIQUE:
        print("  -> ★ QB-UNIQUE  as far as seeds probe, the system itself is unique at this")
        print("     condition even without continuation ⇒ the 2.5-D 'M1 (b)(c) unmeasurable'")
        print("     finding does NOT carry over to 3-D here.")
    else:
        print("  -> QB-MIX  RECORDED, no direction claimed.")
    for r in ok:
        print(f"    seed {r['seed']:2}: cl_p {r['cl_p']:.6f}  LE_upper "
              f"{(r.get('band_LE_upper') or float('nan')):.6f}  M_max {r['m_max']:.4f}  "
              f"|R| {r['res_final']:.2e}")
    return 0


def _write(rows):
    if not rows:
        rows = [dict(note="no legs run")]
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
