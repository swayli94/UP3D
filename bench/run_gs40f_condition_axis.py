"""
Does k*=1 survive a change of FLOW CONDITION?

Pre-registration: docs/dev_phase_four/20260817-1900-condition-axis-prereg.md
Previous round:   docs/dev_phase_four/20260817-1700-level-ladder-verdict.md (K-FIXED)

The previous verdict wrote its own limit as the fourth forbidden sentence: k*=1
across three LEVELS is not k*=1 across CONDITIONS. Here geometry, mesh level and
recipe are frozen and ONLY the condition moves -- a cross design, M-axis at
alpha 3.06 and alpha-axis at M 0.88, sharing the state already on disk.

★ TWO VACUITY TRAPS CLOSED IN ADVANCE (the 8th question):
  (1) a SUBCRITICAL leg has no sonic set, so C is undefined; counting it as "no
      counterexample" would pass the criterion for free. Only legs with >= 100
      crossing faces on the CONTROL arm count.
  (2) "k* unchanged across conditions" is worthless if the conditions did not
      differ. The supersonic-fraction spread must reach 2x or W-FIXED is
      DOWNGRADED to RECORDED.
"""

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "bench"))

import run_capability_matrix as cap                                 # noqa: E402
import run_gs40c_coherence_ladder as L                              # noqa: E402
import run_gs40e_level_ladder as E                                  # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import (build_face_adjacency,              # noqa: E402
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

MSH = REPO / "cases/meshes/onera_m6_wingbody_conforming/coarse.msh"
CACHE = REPO / "bench/gate_results/gs40f_states"
OUT = REPO / "bench/gate_results/gs40f_condition_axis.csv"
ANCHOR_NPZ = REPO / "bench/gate_results/gs40d_levels/coarse.npz"
GAMMA, R_C = 1.4, 0.05

#: legs fixed BEFORE knowing which converge
LEGS = (("M0.78", 0.78, 3.06), ("M0.82", 0.82, 3.06), ("M0.86", 0.86, 3.06),
        ("ANCHOR", 0.88, 3.06),
        ("A1.53", 0.88, 1.53), ("A4.59", 0.88, 4.59))
LADDER = ("P0", "P1", "P2", "P3")
ARMS = LADDER + ("S1",)

#: every threshold carried over verbatim
PEAK_FACTOR, FLAT_RATIO, MIN_CROSS = L.PEAK_FACTOR, L.FLAT_RATIO, L.MIN_CROSS
STRENGTH_SPREAD_MIN = 2.0          # the premise check, registered as a kill condition
ANCHOR_C = dict(P0=0.221, P1=0.479, P2=0.396, P3=0.291)
GATE_S = 25 * 60


def solve_leg(tag, m_inf, alpha, mc, wc):
    npz = CACHE / f"{tag}.npz"
    if npz.exists():
        d = np.load(npz)
        return (np.asarray(d["phi"], float), bool(d["conv"]), int(d["nlim"]),
                int(d["nflr"]), float(d["res"]), str(d["note"]), 0.0)
    t = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", R_C * B_SEMI)
    t0 = time.perf_counter()
    try:
        seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART, alpha_deg=alpha,
                                    **cap.CONF_SEED_KW)
        nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                  phi_init=seed["phi"], gamma_init=seed["gamma"],
                  n_picard_seed=0, tip_taper=t)
        r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=alpha,
                                   m_start=cap.WB_MSTART, dm=cap.DM, dm_min=0.01,
                                   freeze_tol=1e-5, intermediate_tol=1e-4,
                                   newton_kw=nk)
    except Exception as exc:                                       # noqa: BLE001
        return None, False, -1, -1, float("nan"), f"{type(exc).__name__}: {exc}", \
            time.perf_counter() - t0
    wall = time.perf_counter() - t0
    #: classify a failure rather than reporting a bare conv=False (CLAUDE.md's table)
    note = ""
    if not r["converged"]:
        hist = np.asarray(r["residual_history"], float)
        note = (f"accept_reason={r.get('accept_reason')}; "
                f"clamps={int(r['n_limited'])}/{int(r['n_floored'])}; "
                f"gmres_stalled={r.get('n_gmres_stalled')}; "
                f"m_last_converged={r.get('m_last_converged')}; "
                f"target_reached={r.get('target_reached')}; "
                f"tail={hist[-3:] if len(hist) >= 3 else hist}")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz, phi=r["phi"], conv=bool(r["converged"]),
                        nlim=int(r["n_limited"]), nflr=int(r["n_floored"]),
                        res=float(r["residual_history"][-1]), note=note,
                        mesh_sha=mesh_fingerprint(MSH)["sha256"])
    return (np.asarray(r["phi"], float), bool(r["converged"]), int(r["n_limited"]),
            int(r["n_floored"]), float(r["residual_history"][-1]), note, wall)


def main():
    t0 = time.perf_counter()
    E.guard_recipe()
    man = read_manifest(MSH)
    sha = mesh_fingerprint(MSH)["sha256"]
    assert man is not None and sha == man["sha256"], "G-MESH: mesh moved"
    print(f"G-MESH coarse sha {sha[:12]}  PASS")

    mc, wc = cut_wake(read_mesh(MSH))
    nodes, elements = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, elements)
    face_nb, _ = build_face_adjacency(elements)
    fa, fb = L.faces(face_nb)
    P, p1e = {}, None
    for k in (0, 1, 2, 3):
        po, pi = L.khop_patches(elements, face_nb, k)
        P[f"P{k}"] = L.patch_nodes(elements, po, pi)
        if k == 1:
            p1e = (po, pi)
    print(f"mesh {len(nodes)} nodes / {len(elements)} tets; patches built "
          f"({time.perf_counter() - t0:.1f} s)\n")

    rows, keep, C_all, strength = [], {}, {}, {}
    for tag, m_inf, alpha in LEGS:
        if tag == "ANCHOR":
            d = np.load(ANCHOR_NPZ)
            phi, conv, nlim, nflr, res, note, wall = (
                np.asarray(d["phi"], float), bool(d["conv"]), int(d["nlim"]),
                int(d["nflr"]), float(d["res"]), "reused from gs40d", 0.0)
        else:
            phi, conv, nlim, nflr, res, note, wall = solve_leg(tag, m_inf, alpha, mc, wc)
        if phi is None or not conv or nlim or nflr:
            print(f"[{tag}] M{m_inf} a{alpha}  NOT USABLE  conv={conv} "
                  f"clamps={nlim}/{nflr}  {note}")
            rows.append(dict(leg=tag, m_inf=m_inf, alpha=alpha, status="NOT_CONVERGED",
                             converged=conv, n_limited=nlim, n_floored=nflr,
                             res_final=res, wall_s=round(wall, 1), note=note))
            continue

        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)
        m2A = mach_number_squared(q2, m_inf, GAMMA)
        sup_frac = float((m2A >= 1.0).mean()); m_max = float(np.sqrt(m2A.max()))

        g = {n: L.refit(nodes, elements, phi, *P[n]) for n in LADDER}
        g["S1"] = L.average_over_patch(nodes, elements, gA, *p1e)
        e_a = float(np.max(np.abs(g["P0"] - gA)))
        assert e_a < 1e-10, f"G-A {tag}: {e_a:.2e}"

        C, res_c = {}, {}
        for n in ARMS:
            m2 = mach_number_squared((g[n] * g[n]).sum(axis=1), m_inf, GAMMA)
            r = L.coherence(np.sqrt(np.maximum(m2, 0.0)), fa, fb)
            C[n] = r["C"]; res_c[n] = r
            rows.append(dict(leg=tag, m_inf=m_inf, alpha=alpha, status="OK", arm=n,
                             converged=True, n_limited=0, n_floored=0, res_final=res,
                             sup_frac=round(sup_frac, 5), m_max=round(m_max, 4),
                             n_cross=r["n_cross"], max_comp=r["max_comp"],
                             n_ge50=r["n_usable"], C=round(r["C"], 4),
                             wall_s=round(wall, 1), note=note))
        n_cross0 = res_c["P0"]["n_cross"]
        counted = n_cross0 >= MIN_CROSS
        C_all[tag] = C; strength[tag] = (sup_frac, m_max)
        if counted:
            keep[tag] = C
        print(f"[{tag}] M{m_inf} a{alpha}  {wall:5.1f}s  sup_frac {sup_frac:.4f} "
              f"M_max {m_max:.3f}  n_cross(P0) {n_cross0:>6}  "
              f"{'COUNTED' if counted else '★ EXCLUDED-SUBCRITICAL'}")
        print("        C = " + "  ".join(f"{n} {C[n]:.3f}" for n in ARMS))

        if tag == "ANCHOR":
            bad = {k: (C[k], ANCHOR_C[k]) for k in LADDER
                   if abs(C[k] - ANCHOR_C[k]) > 1e-3}
            print(f"G-ANCHOR vs the published coarse read: "
                  f"{'PASS' if not bad else '** FAIL ** ' + str(bad)}")
            assert not bad, "G-ANCHOR: kill criterion 2 -- instrument, not a finding"

    # ------------------------------------------------------------- verdict --
    print("\n" + "=" * 76)
    if len(keep) < 2:
        verdict = f"UNDEFINED (only {len(keep)} counting leg(s)) -- kill criterion 3"
        per = {}
    else:
        per = {}
        for tag, C in keep.items():
            vals = [C[k] for k in LADDER]
            if max(vals) / max(min(vals), 1e-30) < FLAT_RATIO:
                per[tag] = "NO-PEAK (flat)"
            elif max(C["P1"], C["P2"]) >= PEAK_FACTOR * max(C["P0"], C["P3"]):
                kstar = min(LADDER, key=lambda k: (-round(C[k], 6), LADDER.index(k)))
                per[tag] = f"PEAK k*={kstar}"
            else:
                per[tag] = "NO-PEAK (no interior max)"
            print(f"{tag:8} C = " + " ".join(f"{C[k]:.3f}" for k in LADDER)
                  + f"   -> {per[tag]}")
        codes = list(per.values())
        if any(c.startswith("NO-PEAK") for c in codes):
            verdict = "W-VANISH (a counting leg has no interior maximum) -- scope cut back"
        elif all(c == "PEAK k*=P1" for c in codes):
            verdict = "W-FIXED (k*=1 on every counting leg)"
        else:
            verdict = "W-DRIFT (interior max everywhere, but k* moves with condition)"

    # the registered PREMISE check on strength spread
    sf = [strength[t][0] for t in keep]
    spread = (max(sf) / max(min(sf), 1e-30)) if sf else float("nan")
    print("=" * 76)
    print(f"strength spread across counting legs: sup_frac "
          f"{min(sf):.4f}..{max(sf):.4f} = {spread:.2f}x "
          f"(premise needs >= {STRENGTH_SPREAD_MIN}x)")
    if verdict.startswith("W-FIXED") and spread < STRENGTH_SPREAD_MIN:
        verdict = ("W-FIXED DOWNGRADED to RECORDED -- the conditions did not "
                   f"differ enough ({spread:.2f}x < {STRENGTH_SPREAD_MIN}x); "
                   "no cross-condition claim (kill criterion 4)")
    print(f"VERDICT: {verdict}")
    print("\n★ COHERENT IS NOT CORRECT; zero statement about behaviour inside the\n"
          "  solver; and geometry is still ONE family (wing-body conforming coarse).")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict=verdict, per_leg=per, strength_spread=round(float(spread), 3),
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 5"


if __name__ == "__main__":
    main()
