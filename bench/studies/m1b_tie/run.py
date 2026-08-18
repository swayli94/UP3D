"""M1-b + M1-c -- measure S-TIE at the registered quantity, and ask what sets D's floor.

Binding text: docs/dev_phase_five/20260822-2000-m1b-prereg.md  (+ addendum #1,
both committed before this file existed).

★★ Round 8 implemented S-TIE against the FOOTPRINT VECTOR; the registered quantity
is the argmin of D, and D needs J1 -- the shock position read off the section Cp
curve. So this script computes D the way round 18 formed it, by importing round 18's
own foot_from_field rather than re-deriving it.

★★★ M1-c comes from what building D exposes: it compares a volume field's M = 1 wall
footprint against a wall-Cp shock position. Those are two definitions of a shock, so
D's floor may be a definition difference rather than the extractor's resolution.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/m1b_tie/run.py
"""

import csv
import hashlib
import importlib.util
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "bench"))
RESULTS = os.path.join(HERE, "results")

import run_gs40c_coherence_ladder as L                                  # noqa: E402
import run_gs40i_position_and_blastradius as I                          # noqa: E402
import run_gs40j_position_primary as J                                  # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import build_face_adjacency, precompute_element_geometry  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared               # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                    # noqa: E402
from pyfp3d.post.shock import shock_report                              # noqa: E402

#: ★ the sub-cell rule is round 8's, imported rather than re-typed
_spec = importlib.util.spec_from_file_location(
    "m1a_subcell_run", os.path.join(ROOT, "bench/studies/m1a_subcell/run.py"))
M1A = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(M1A)

GAMMA = 1.4
#: G-REPRO, addendum #1: BOTH scopes, both must reproduce
ANCHOR_TIES = {"LADDER": 4, "ARMS": 5}
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def d_median(feet, j1):
    """Round 18's D, bench/run_gs40j_position_primary.py:172-174, verbatim."""
    dd = [abs(feet[e_] - j1[e_]) for e_ in J.ETAS
          if j1[e_] is not None and np.isfinite(feet[e_])]
    return float(np.median(dd)) if len(dd) >= J.MIN_STATIONS else np.nan


def tie_scope(D, scope):
    dm = min(D[a] for a in scope)
    return [a for a in scope if abs(D[a] - dm) < 1e-9]


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    J.guard_no_solve_executed()
    print("  G-FROZEN-LIB  the extractors live in bench/, pyfp3d/ is untouched")
    print(f"  G-LADDER      ties are taken over LADDER={J.LADDER} (round 18's own "
          f"scope) AND over ARMS={J.ARMS}; anchors {ANCHOR_TIES}")

    rows, drows, ties = [], [], {"LADDER": {"percell": 0, "subcell": 0},
                                 "ARMS": {"percell": 0, "subcell": 0}}
    n_states = 0
    for geom, label, sp, msh, m_inf in J.STATES:
        if not os.path.exists(str(sp)):
            continue
        d = np.load(sp)
        if not bool(d["conv"]) or int(d["nlim"]) or int(d["nflr"]):
            continue
        phi = np.asarray(d["phi"], float)
        sha = hashlib.sha256(phi.tobytes()).hexdigest()[:12]
        mc, _ = cut_wake(read_mesh(str(msh)))
        nodes, elements = mc.nodes, mc.elements
        if len(phi) != len(nodes):
            continue
        wall_faces = mc.boundary_faces["wall"]

        # ---- J1: exogenous, bit-identical across arms, computed ONCE ---------
        j1 = {}
        for e_ in J.ETAS:
            try:
                cur = section_cp_curve(mc, phi, eta=e_, b_semi=J.B_SEMI, m_inf=m_inf)
                x = shock_report(cur, m_inf)["upper"].get("x_shock")
                j1[e_] = None if x is None or not np.isfinite(x) else float(x)
            except Exception:                                          # noqa: BLE001
                j1[e_] = None
        n_ok = sum(v is not None for v in j1.values())
        if n_ok < J.MIN_STATIONS:
            print(f"  [{label:11}] EXCLUDED-NO-WALL-SHOCK  J1 {n_ok}/7")
            continue

        # ---- the arms, exactly as round 18 builds them ------------------------
        B, _ = precompute_element_geometry(nodes, elements)
        fn, _ = build_face_adjacency(elements)
        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)
        P, pel = {}, {}
        for k in (0, 1, 2, 3):
            po, pi = L.khop_patches(elements, fn, k)
            P[f"P{k}"] = L.patch_nodes(elements, po, pi); pel[f"P{k}"] = (po, pi)
        g = {n: L.refit(nodes, elements, phi, *P[n]) for n in J.LADDER}
        for s, p in J.SMOOTH.items():
            g[s] = L.average_over_patch(nodes, elements, gA, *pel[p])
        assert np.max(np.abs(g["P0"] - gA)) < 1e-10, f"G-A {label}"

        cent = nodes[elements].mean(axis=1)
        wn = np.unique(wall_faces.reshape(-1))
        isw = np.zeros(len(nodes), bool); isw[wn] = True
        adj = isw[elements].any(axis=1)

        Dp, Ds = {}, {}
        for n in J.ARMS:
            M = np.sqrt(np.maximum(mach_number_squared(
                (g[n] * g[n]).sum(axis=1), m_inf, GAMMA), 0.0))
            fp = I.foot_from_field(nodes, elements, M, wall_faces, J.ETAS, J.B_SEMI)
            fs, _ = M1A.foot_subcell(cent, adj, M, J.ETAS, J.B_SEMI)
            Dp[n], Ds[n] = d_median(fp, j1), d_median(fs, j1)
            for e_ in J.ETAS:                       # M1-c needs the SIGNED residual
                if j1[e_] is None:
                    continue
                for rule, ff in (("percell", fp), ("subcell", fs)):
                    if np.isfinite(ff[e_]):
                        drows.append({"state": label, "arm": n, "eta": e_,
                                      "rule": rule, "j1": j1[e_], "foot": ff[e_],
                                      "signed": ff[e_] - j1[e_]})
        if any(not np.isfinite(Dp[n]) for n in J.LADDER):
            print(f"  [{label:11}] EXCLUDED  an arm has < {J.MIN_STATIONS} stations")
            continue
        n_states += 1
        for scope, names in (("LADDER", J.LADDER), ("ARMS", J.ARMS)):
            for rule, D in (("percell", Dp), ("subcell", Ds)):
                if len(tie_scope(D, names)) > 1:
                    ties[scope][rule] += 1
        rows.append({"state": label, "geom": geom, "m_inf": m_inf, "phi_sha": sha,
                     "n_j1": n_ok,
                     **{f"D_pc_{n}": round(Dp[n], 6) for n in J.ARMS},
                     **{f"D_sc_{n}": round(Ds[n], 6) for n in J.ARMS},
                     "kpos_pc": min(J.LADDER, key=lambda a: Dp[a]),
                     "kpos_sc": min(J.LADDER, key=lambda a: Ds[a])})
        print(f"  [{label:11}] sha {sha}  J1 {n_ok}/7   D_pc "
              + " ".join(f"{n}:{Dp[n]:.4f}" for n in J.LADDER)
              + f"   ties_pc {len(tie_scope(Dp, J.LADDER))}"
              f"  ties_sc {len(tie_scope(Ds, J.LADDER))}")

    if not rows:
        print("no usable states"); return 1
    for name, rr in (("D.csv", rows), ("signed.csv", drows)):
        with open(os.path.join(RESULTS, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rr[0].keys()))
            w.writeheader(); w.writerows(rr)

    # ---- G-REPRO ----------------------------------------------------------
    got = {s: ties[s]["percell"] for s in ("LADDER", "ARMS")}
    ok = all(got[s] == ANCHOR_TIES[s] for s in got)
    _record("G-REPRO", "per-cell argmin-tie count reproduces round 18, BOTH scopes",
            f"LADDER {ANCHOR_TIES['LADDER']}/13, ARMS {ANCHOR_TIES['ARMS']}/13",
            f"LADDER {got['LADDER']}/{n_states}, ARMS {got['ARMS']}/{n_states}",
            "G-REPRO PASS" if ok else "★ G-REPRO FAIL -- kill criterion 1")
    if not ok:
        print("\n★★ kill criterion 1: the baseline does not reproduce. Stopping "
              "before any conclusion -- the recomputation is what to check first.")
        with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
            w.writerows(SUMMARY)
        return 2

    # ---- T-TIE ------------------------------------------------------------
    tp, ts = ties["LADDER"]["percell"], ties["LADDER"]["subcell"]
    _record("T-TIE", "argmin-tie states, per-cell vs sub-cell (LADDER scope)",
            "<=1 => resolution; still >=4 => attribution DEFINITIVELY excluded",
            f"per-cell {tp}/{n_states} -> sub-cell {ts}/{n_states}  "
            f"(ARMS: {ties['ARMS']['percell']} -> {ties['ARMS']['subcell']})",
            "T-TIE PASS -- ties were resolution" if ts <= 1
            else "★ T-TIE: ties are NOT a resolution artefact")
    dpc = np.array([[r[f"D_pc_{n}"] for n in J.LADDER] for r in rows]).ravel()
    dsc = np.array([[r[f"D_sc_{n}"] for n in J.LADDER] for r in rows]).ravel()
    _record("T-D", "D itself under the two extractors", "RECORDED",
            f"median D per-cell {np.median(dpc):.5f} -> sub-cell {np.median(dsc):.5f}"
            f"  (median |change| {np.median(np.abs(dsc - dpc)):.5f})", "RECORDED")

    # ---- M1-c: is D's floor a BIAS (definition) or SCATTER (resolution)? ----
    for rule in ("percell", "subcell"):
        s = np.array([r["signed"] for r in drows if r["rule"] == rule])
        bias, mad = float(np.median(s)), float(np.median(np.abs(s - np.median(s))))
        pos = float((s > 0).mean())
        _record("T-DEF" if rule == "percell" else "T-DEF(sub)",
                f"bias vs scatter of (foot - J1), {rule}",
                "|bias| >> MAD => D's floor is a DEFINITION difference",
                f"n={len(s)}  bias {bias:+.5f}  MAD {mad:.5f}  ratio {abs(bias)/max(mad,1e-12):.2f}",
                "T-DEF: BIAS-dominated" if abs(bias) > 2.0 * mad else
                "T-DEF: scatter-dominated")
        _record("T-SIGN" if rule == "percell" else "T-SIGN(sub)",
                f"fraction of readings with foot > J1, {rule}",
                "near 0 or 1 => a fixed sign => definition difference",
                f"{pos:.3f} of {len(s)}",
                "T-SIGN: fixed sign" if (pos > 0.9 or pos < 0.1) else "T-SIGN: mixed")

    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n  {n_states} states, {len(drows)} signed readings, "
          f"{time.perf_counter() - t0:.1f} s")
    print("\n★ This does NOT fix M1: its three gates are untouched and the recorded "
          "'not reachable on the current discretisation' verdict stands.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
