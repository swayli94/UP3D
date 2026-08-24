"""S4's two founding premises, checked on HEAD. Near-zero cost.

Pre-registered in phases/p3/docs/dev_phase_three/20260815-2200-s4-premise-check-prereg.md.

S4 ("viscous restart by subtraction") rests on two PHASE-ONE measurements:
  P1  the boundary layer costs 9x the inviscid solve (1032 s against 115 s);
  P2  the crossflow it pays for is unused (max|B|/|A| at or below 0.072).
This session has three times shown that a stale premise puts a whole line on sand, so both are
re-measured on HEAD before anything is built on them.

★ Cheap by construction: phi is already cached, so the FP driver is a STUB that returns it instantly,
and n_outer_max = 1 gives exactly ONE IBL solve. No coupling loop is run -- both premises are
single-solve quantities (P1 compares single solves, P2 is a single-solve diagnostic), so the loop adds
cost without changing what is measured.

★★ Declared: with a stub driver the IBL sees the INVISCID u_e -- precisely the state at outer k = 1 of a
real loop. Nothing here is a converged coupled state and nothing is read as one.

★ G-U: the crossflow ratio's definition is READ from the V5 bridge (bench/studies/v5_m6_bridge/run.py
lines 235-247) rather than recalled: live = ~tip, a_mag = max|U[live, 1]|, b_mag = max|U[live, 2]|, and
the tip mask is z > z_tip * (1 - 0.05).

Outputs (TRACKED): bench/gate_results/task3_s4_premise_check.csv
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

from pyfp3d.mesh.reader import read_mesh                              # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                             # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, chord_at, x_le         # noqa: E402
from pyfp3d.viscous.coupling import (CouplingConfig, build_wing_case,  # noqa: E402
                                     run_loose_coupling)
import run_task3_fixed_budget as FB                                   # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_s4_premise_check.csv")
INVISCID_CSV = os.path.join(HERE, "gate_results", "task3_selfconvergence.csv")
#: the same leg as this phase's rounds
M_INF, ALPHA = 0.8395, 3.06
RE_MAC, X_TR, TIP_FRAC = 11.72e6, 0.05, 0.05
#: ★ 20260816-0200 round: S2 then S4 in ONE process so the two wall clocks share load and threads
#: (wall clock is a calibration -- the project has recorded a 5.4x spread on identical physics).
#: S0 is carried from the previous round's committed CSV and labelled cross-run.
LEVELS = tuple(
    x for x in (("S2", 0.0300, "rnd_0.03_0.015_0.015_3.6.msh"),
                ("S4", 0.0150, "rnd_0.015_0.0075_0.0075_1.8.msh"),
                ("S0", 0.0600, "rnd_0.06_0.03_0.03_7.2.msh"),
                #: ★ 20260816-0600: the REPEAT leg. G-L2 replaces the defective load proxy by
                #: measuring the property it was for -- repeatability -- with S4 sandwiched
                #: between two S2 legs in ONE process.
                ("S2r", 0.0300, "rnd_0.03_0.015_0.015_3.6.msh"))
    if x[0] in os.environ.get("PYFP3D_S4_LEVELS", "S2,S4,S2r").split(","))
#: phase-one values, quoted as ORDER-OF-MAGNITUDE references only
P1_REF, P2_REF = 9.0, 0.072
P1_BAND, P2_MAX = (4.0, 20.0), 0.15
#: ★ DEFECT FIXED 20260816: the registration for the repeatability round set this to 2400 s, but the
#: constant was left at 900 from the first premise round -- so the REPEAT leg (S2r) was skipped by my
#: own stale gate and G-L2 could not be evaluated. Implementation-versus-registration mismatch, mine.
TOTAL_GATE_S = 2400.0


def inviscid_wall(tag):
    """G-P: the inviscid wall clock must come from THIS phase's own committed CSV, with the path
    printed. Phase one's 115 s is cross-time and cross-pipeline and may not substitute.

    ★ DEFECT FIXED 20260816 (the SECOND self-inflicted blocker in this round, both in the
    report/instrument layer -- the project's own lesson is to guard report code the way solver code is
    guarded): the repeat leg is S2 BY CONSTRUCTION, so its inviscid reference is S2's row. The phi path
    already stripped the trailing 'r'; this lookup did not, so the leg died on kill clause 3 after 17
    minutes of compute.
    """
    tag = tag.rstrip("r")
    with open(INVISCID_CSV) as fh:
        for r in csv.DictReader(fh):
            if r["leg"] == tag:
                return float(r["solve_s"]), int(float(r["n_tet"])), r["res_final"]
    return None, None, None


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}")
    print("★ G-U: crossflow ratio read from bench/studies/v5_m6_bridge/run.py:235-247 --")
    print("   live = ~tip;  a = max|U[live,1]|;  b = max|U[live,2]|;  tip: z > z_tip*(1-0.05)\n")
    la0 = os.getloadavg()
    #: ★ G-R is built on NON-TIMING quantities: wall clock is not reproducible (5.4x recorded),
    #: so the chain of custody uses S2's surface node count and P2 ratio instead.
    GR_S2 = dict(n_node_surf=3557, p2_ratio=0.1471)
    rows, t_all = [], time.perf_counter()

    for tag, h, mesh_name in LEVELS:
        if time.perf_counter() - t_all > TOTAL_GATE_S:
            print(f"★ total gate exceeded before {tag} (kill clause 4)")
            break
        #: ★ the repeat leg is S2 by construction -- same phi, same mesh, same knobs
        ppath = os.path.join(FB.SCRATCH, f"phi_{tag.rstrip(chr(114))}.npz")
        mpath = os.path.join(FB.SCRATCH, mesh_name)
        if not (os.path.exists(ppath) and os.path.exists(mpath)):
            print(f"★ {tag}: phi or mesh cache missing -- STOP (kill clause 1; no re-solve)")
            return 1
        z = np.load(ppath)
        phi_c, gam_c = np.asarray(z["phi"]), np.asarray(z["gamma"])
        print(f"=== {tag} (h_wall {h}) ===")
        print(f"  G-phi: {ppath.split('/')[-1]}  phi {phi_c.shape}  gamma {gam_c.shape} "
              f"(cached -- no inviscid re-solve)")
        inv_s, n_tet, res_f = inviscid_wall(tag)
        if inv_s is None:
            print(f"  ★ no same-provenance inviscid wall clock for {tag} -- STOP (kill clause 3)")
            return 1
        print(f"  G-P: inviscid wall {inv_s:.1f} s from {INVISCID_CSV.split('/')[-1]} "
              f"(leg {tag}, {n_tet} tets, |R| {res_f})")

        mc, wc = cut_wake(read_mesh(mpath))
        wall = mc.boundary_faces["wall"]
        cfg = CouplingConfig(re_chord=RE_MAC / MAC, m_inf=M_INF, alpha_deg=ALPHA,
                             x_tr_upper=X_TR, x_tr_lower=X_TR, n_outer_max=1)
        case = build_wing_case(mc.nodes, mc.elements, wall, cfg,
                               x_le=x_le, chord_at=chord_at, tip_mask_frac=TIP_FRAC)
        sm = case.sm
        print(f"  IBL surface: {sm.n_node} nodes / {sm.n_tri} tris; "
              f"tip-masked {int(case.outflow_pin_surf.sum())}; "
              f"LE-band {int(case.inflow_candidates.sum())}; "
              f"BL unknowns = 6 x {sm.n_node} = {6 * sm.n_node}")

        calls = [0]

        def stub(rhs, seed):
            """Zero-cost FP driver: returns the cached inviscid state. ★ So the IBL sees the
            INVISCID u_e -- exactly outer k = 1 of a real loop, and nothing here is a converged
            coupled state."""
            calls[0] += 1
            return phi_c, gam_c, dict(converged=True, stub=True)

        t0 = time.perf_counter()
        try:
            res = run_loose_coupling(stub, case, cfg)
        except Exception as exc:                                       # noqa: BLE001
            print(f"  ★ {tag}: one IBL solve FAILED on HEAD: {type(exc).__name__}: {exc}")
            rows.append(dict(level=tag, ibl_ok=False, error=f"{type(exc).__name__}: {exc}"))
            continue
        ibl_s = time.perf_counter() - t0
        print(f"  one IBL solve: {ibl_s:.1f} s   (fp_solve stub calls {calls[0]}, "
              f"n_outer {res.n_outer}, converged {res.converged})")

        #: --- P2: the crossflow ratio, definition read from the bridge ---------------------------
        U = res.U
        zz = sm.xyz[:, 2]
        z_tip = float(np.max(zz))
        tip = zz > z_tip * (1.0 - TIP_FRAC)
        live = ~tip
        a = np.abs(U[live, 1])
        b = np.abs(U[live, 2])
        a_mag, b_mag = float(a.max()), float(b.max())
        ratio = b_mag / max(a_mag, 1e-30)
        pw = np.abs(U[live, 2]) / np.maximum(np.abs(U[live, 1]), 1e-30)
        print(f"  P2 crossflow: max|B| {b_mag:.4e} / max|A| {a_mag:.4e} = {ratio:.4f}   "
              f"(n_live {int(live.sum())})")
        print(f"     pointwise |B|/|A| median {np.median(pw):.4f}  p90 {np.percentile(pw, 90):.4f}"
              f"  max {pw.max():.4f}")

        #: ★ CACHE BEFORE YOU REPORT: last round did not store U, which is exactly why the
        #: B-field question could not be answered without re-running. Columns per the G-U read:
        #: 1 = A (streamwise), 2 = B (crossflow).
        upath = os.path.join(FB.SCRATCH, f"U_{tag}.npz")
        np.savez(upath, U=U, xyz=sm.xyz, tip=tip)
        print(f"  G-cache: U written to {upath.split(chr(47))[-1]}  U {U.shape}")

        #: --- P1: the cost ratio -----------------------------------------------------------------
        p1 = ibl_s / inv_s
        print(f"  P1 cost: one IBL {ibl_s:.1f} s / inviscid {inv_s:.1f} s = {p1:.2f}x   "
              f"(phase-one reference {P1_REF}x)")
        rows.append(dict(level=tag, h_wall=h, n_node_surf=sm.n_node, n_tri=sm.n_tri,
                         bl_unknowns=6 * sm.n_node, n_tet=n_tet, ibl_ok=True,
                         ibl_s=round(ibl_s, 2), inviscid_s=inv_s, p1_ratio=p1,
                         n_outer=res.n_outer, ibl_converged=bool(res.converged),
                         a_mag=a_mag, b_mag=b_mag, p2_ratio=ratio,
                         p2_median=float(np.median(pw)), p2_p90=float(np.percentile(pw, 90)),
                         p2_max=float(pw.max()), n_live=int(live.sum()),
                         n_tip_masked=int(tip.sum())))
        print()

    la1 = os.getloadavg()
    #: ★★ G-L2 REPLACES G-L. The old guard read load average, which INCLUDES this process, so it
    #: measured its own effect and could never pass a compute-heavy round (defect #8, the third
    #: self-reference this session). The new guard measures the property the old one was for --
    #: repeatability -- by comparing the two S2 legs that sandwich S4. It reads only its own two
    #: wall clocks and no global state, so it is not self-referential (G-B, declared in advance).
    print(f"  (load average before {la0[0]:.2f} -> after {la1[0]:.2f}; RECORDED only -- it "
          f"includes this process, which is why it is no longer a guard)")
    t_a = next((r["ibl_s"] for r in rows if r["level"] == "S2" and r.get("ibl_ok")), None)
    t_b = next((r["ibl_s"] for r in rows if r["level"] == "S2r" and r.get("ibl_ok")), None)
    judgeable = False
    if t_a is None or t_b is None:
        print("★ G-L2: no repeat leg -> the trend stays UNJUDGED")
    else:
        rel = abs(t_b - t_a) / t_a
        judgeable = rel <= 0.10
        print(f"★ G-L2 (repeatability): S2 {t_a:.1f} s -> S2r {t_b:.1f} s = {100 * rel:+.1f} % "
              f"-> {chr(9733) + chr(9733) + ' PASS, the trend MAY be judged' if judgeable else chr(9733) + ' FAIL, timings not repeatable to 10 %: trend stays UNJUDGED'}")
    s2 = next((r for r in rows if r["level"] == "S2" and r.get("ibl_ok")), None)
    if s2 is not None:
        okn = s2["n_node_surf"] == GR_S2["n_node_surf"]
        okp = abs(s2["p2_ratio"] - GR_S2["p2_ratio"]) <= 5e-5
        print(f"★ G-R (non-timing): S2 surface nodes {s2['n_node_surf']} vs "
              f"{GR_S2['n_node_surf']} / P2 {s2['p2_ratio']:.4f} vs {GR_S2['p2_ratio']} -> "
              f"{'PASS' if okn and okp else chr(9733) + ' FAIL'}")
        if not (okn and okp):
            print("  -> STOP: not the same leg (kill clause 3).")
            _write(rows)
            return 1
    _write(rows)
    ok = [r for r in rows if r.get("ibl_ok")]
    if not ok:
        print("★ Q UNDEFINED: one IBL solve is not usable on HEAD at any level tried "
              "(kill clause 2). Reported as a reading, not a failure.")
        return 0
    b = ok[0]                       #: the finest usable level is the binding one
    print("=== the verdict (binding = the finest usable level) ===")
    if not judgeable:
        print("  ★ G-L2 did not pass -> the LEVEL TREND is not judged; the Q bands below are")
        print("    read on the binding level only, as in the previous round.")
    print(f"  level {b['level']}   P1 = {b['p1_ratio']:.2f}x (band {P1_BAND})   "
          f"P2 = {b['p2_ratio']:.4f} (band <= {P2_MAX})")
    if b["p1_ratio"] > 100.0 or b["p1_ratio"] < 1.0:
        print("  -> ★ SUSPICION, not a verdict (kill clause 5): a ratio this far out means I divided")
        print("     two different things. Check what each side measures before reading anything.")
        return 1
    p1ok = P1_BAND[0] <= b["p1_ratio"] <= P1_BAND[1]
    p2ok = b["p2_ratio"] <= P2_MAX
    if p1ok and p2ok:
        print("  -> ★★ Q-BOTH  both premises reproduce on HEAD ⇒ S4's subtraction has HEAD evidence:")
        print("     the BL is the cost, and the crossflow it pays for is unused ⇒ GS4.1 (chordwise")
        print("     2-D strips) is the right first gate.")
    elif p1ok:
        print("  -> Q-COST  only the cost premise reproduces ⇒ the crossflow may NOT be dropped;")
        print("     GS4.1's 2-D strip argument needs re-making.")
    elif p2ok:
        print("  -> Q-CROSS  only the crossflow premise reproduces ⇒ the BL is NOT the cost driver,")
        print("     so subtraction buys little speed and S4's payoff must be re-estimated.")
    else:
        print("  -> ★★★ Q-NEITHER  neither premise reproduces ⇒ S4's founding evidence needs")
        print("     REWRITING before anything is built on it. The most valuable outcome.")
    if len(ok) > 1:
        print("\n  other levels (RECORDED): "
              + "; ".join(f"{r['level']} P1 {r['p1_ratio']:.2f}x P2 {r['p2_ratio']:.4f}"
                          for r in ok[1:]))
    print("\n  ★ wall clock is a CALIBRATION, not a guarantee -- P1's band is an order of magnitude,")
    print("    and the thread count and load average are printed above.")
    return 0


def _write(rows):
    if not rows:
        rows = [dict(note="no rows")]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
