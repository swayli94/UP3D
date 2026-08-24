"""R11 -- give the footprint rule failure detection, then re-measure the arm ranking.

Binding text: docs/dev_phase_five/20260823-0000-r11-prereg.md (committed before this
file existed).

Round 10 measured that `foot_from_field` returns the downstream end of its own search
band whenever the supersonic region reaches it -- a "not found" returned as a value,
20.9% of readings, at a rate that is monotone in patch size (P0 36.3% -> P3 12.5%). So
the arm ranking D exists to produce is confounded with the instrument's failure rate.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/r11_saturation/run.py
"""

import collections
import csv
import hashlib
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

GAMMA = 1.4
SAT, SAT_ALT = 0.95, 0.90            # G-SATDEF: both must classify identically
DETA = 0.01
ANCHOR = {"sat": 131, "n": 627, "ties": 4}
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def foot_flagged(nodes, elements, M, wall_faces, etas, b_semi):
    """I.foot_from_field, plus THE FLAG IT NEVER HAD.

    The rule reports the largest x that is still supersonic inside the station band.
    When the supersonic set reaches the band's downstream end there is no shock inside
    the window and the value returned is the window's own edge -- a non-answer. The
    flag says so instead of letting the caller read it as a position.
    """
    cent = nodes[elements].mean(axis=1)
    wn = np.unique(wall_faces.reshape(-1))
    isw = np.zeros(len(nodes), bool); isw[wn] = True
    adj = isw[elements].any(axis=1)
    out, sat = {}, {}
    for e_ in etas:
        z0 = e_ * b_semi
        inb = adj & (np.abs(cent[:, 2] - z0) <= 0.02 * b_semi)
        band = inb & (M >= 1.0)
        if band.sum() < 3:
            out[e_] = np.nan; sat[e_] = False; continue
        x0, x1 = cent[inb, 0].min(), cent[inb, 0].max()
        out[e_] = float((cent[band, 0].max() - x0) / max(x1 - x0, 1e-9))
        sat[e_] = out[e_] > SAT
    return out, sat


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    J.guard_no_solve_executed()
    print("  G-FROZEN-LIB  the flag lives in bench/, pyfp3d/ is untouched")

    rows, j1rows, statedat = [], [], []
    n_sat = n_tot = n_alt = 0
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

        def j1_at(e_, _mc=mc, _phi=phi, _m=m_inf):
            try:
                cur = section_cp_curve(_mc, _phi, eta=e_, b_semi=J.B_SEMI, m_inf=_m)
                x = shock_report(cur, _m)["upper"].get("x_shock")
                return None if x is None or not np.isfinite(x) else float(x)
            except Exception:                                          # noqa: BLE001
                return None
        j1 = {e_: j1_at(e_) for e_ in J.ETAS}
        j1m = {e_: j1_at(e_ - DETA) for e_ in J.ETAS}
        j1p = {e_: j1_at(e_ + DETA) for e_ in J.ETAS}
        if sum(v is not None for v in j1.values()) < J.MIN_STATIONS:
            continue
        #: R-J1 -- the SIGNED one-sided differences round 10 threw away
        for e_ in J.ETAS:
            if all(x[e_] is not None for x in (j1, j1m, j1p)):
                j1rows.append({"state": label, "eta": e_,
                               "fwd": j1p[e_] - j1[e_], "bwd": j1[e_] - j1m[e_]})

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

        feet, sats = {}, {}
        for n in J.ARMS:
            M = np.sqrt(np.maximum(mach_number_squared(
                (g[n] * g[n]).sum(axis=1), m_inf, GAMMA), 0.0))
            f_, s_ = foot_flagged(nodes, elements, M, wall_faces, J.ETAS, J.B_SEMI)
            ref = I.foot_from_field(nodes, elements, M, wall_faces, J.ETAS, J.B_SEMI)
            for e_ in J.ETAS:                    # the flag must not change the value
                a, b = f_[e_], ref[e_]
                assert (np.isnan(a) and np.isnan(b)) or a == b, f"G-VALUE {label} {n}"
            feet[n], sats[n] = f_, s_
            for e_ in J.ETAS:
                if j1[e_] is None or not np.isfinite(f_[e_]):
                    continue
                n_tot += 1; n_sat += bool(s_[e_]); n_alt += bool(f_[e_] > SAT_ALT)
                rows.append({"state": label, "arm": n, "eta": e_, "foot": f_[e_],
                             "j1": j1[e_], "d": f_[e_] - j1[e_], "sat": int(s_[e_])})
        #: G-COMMONSET -- stations unsaturated for EVERY arm, so all arms are
        #: compared on identical stations (the sample set is fixed BEFORE comparing)
        common = [e_ for e_ in J.ETAS
                  if j1[e_] is not None
                  and all(np.isfinite(feet[n][e_]) and not sats[n][e_] for n in J.ARMS)]
        allset = [e_ for e_ in J.ETAS
                  if j1[e_] is not None and all(np.isfinite(feet[n][e_]) for n in J.ARMS)]

        def dmed(n, ES):
            v = [abs(feet[n][e_] - j1[e_]) for e_ in ES]
            return float(np.median(v)) if len(v) >= J.MIN_STATIONS else np.nan
        Dall = {n: dmed(n, allset) for n in J.ARMS}
        Dcln = {n: dmed(n, common) for n in J.ARMS}
        ok = len(common) >= J.MIN_STATIONS
        ka = min(J.LADDER, key=lambda a: Dall[a]) if all(
            np.isfinite(Dall[a]) for a in J.LADDER) else None
        kc = (min(J.LADDER, key=lambda a: Dcln[a]) if ok and all(
            np.isfinite(Dcln[a]) for a in J.LADDER) else None)
        statedat.append({"state": label, "sha": sha, "n_all": len(allset),
                         "n_common": len(common), "usable": int(ok),
                         "kpos_all": ka, "kpos_clean": kc,
                         **{f"Dall_{n}": round(Dall[n], 6) for n in J.LADDER},
                         **{f"Dcln_{n}": (round(Dcln[n], 6)
                                          if np.isfinite(Dcln[n]) else "")
                            for n in J.LADDER}})
        print(f"  [{label:11}] sha {sha}  stations all {len(allset)}/7 -> "
              f"common-clean {len(common)}/7   k_pos {ka} -> {kc}"
              f"{'' if ok else '   ★ UNDEFINED (<3)'}")

    for name, rr in (("readings.csv", rows), ("states.csv", statedat),
                     ("j1_signed.csv", j1rows)):
        with open(os.path.join(RESULTS, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rr[0].keys()))
            w.writeheader(); w.writerows(rr)

    # ---- guards -----------------------------------------------------------
    _record("G-REPRO", "saturation census reproduces round 10",
            f"{ANCHOR['sat']}/{ANCHOR['n']}", f"{n_sat}/{n_tot}",
            "G-REPRO PASS" if (n_sat, n_tot) == (ANCHOR["sat"], ANCHOR["n"])
            else "★ G-REPRO FAIL -- kill criterion 1")
    _record("G-SATDEF", f"foot>{SAT} and foot>{SAT_ALT} classify identically",
            "identical counts (round 10 measured a clean gap)",
            f"{n_sat} vs {n_alt}",
            "G-SATDEF PASS" if n_sat == n_alt else "★ G-SATDEF FAIL -- kill criterion 1")
    if (n_sat, n_tot) != (ANCHOR["sat"], ANCHOR["n"]) or n_sat != n_alt:
        print("\n★★ kill criterion 1 -- stopping before any conclusion.")
        with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
            w.writerows(SUMMARY)
        return 2

    # ---- R-N --------------------------------------------------------------
    usable = [s for s in statedat if s["usable"]]
    _record("R-N", "states retaining >=3 common-clean stations",
            "if most states fall below 3, D cannot be computed cleanly at all",
            f"{len(usable)}/{len(statedat)}; common-clean per state "
            + " ".join(str(s["n_common"]) for s in statedat),
            "R-N: D is computable" if len(usable) >= 7 else
            "★ R-N: D is NOT cleanly computable on most states")

    # ---- R-RANK -----------------------------------------------------------
    chg = [s for s in usable if s["kpos_all"] != s["kpos_clean"]]
    ha = collections.Counter(s["kpos_all"] for s in usable)
    hc = collections.Counter(s["kpos_clean"] for s in usable)
    _record("R-RANK", "states whose D-argmin changes once saturated stations go",
            ">=5/13 => round 18's ranking was contaminated;  <=2/13 => it survives",
            f"{len(chg)}/{len(usable)} changed "
            f"({', '.join(s['state'] for s in chg) if chg else 'none'});  "
            f"k_pos* all {dict(ha)} -> clean {dict(hc)}",
            "★ R-RANK: CONTAMINATED" if len(chg) >= 5 else
            ("R-RANK: ranking survives" if len(chg) <= 2 else "R-RANK: UNDECIDABLE"))

    # ---- R-J1 (folds M1-f) ------------------------------------------------
    fw = np.array([r["fwd"] for r in j1rows]); bw = np.array([r["bwd"] for r in j1rows])
    keep = (np.abs(fw) > 1e-12) & (np.abs(bw) > 1e-12)
    same = float((np.sign(fw[keep]) == np.sign(bw[keep])).mean())
    _record("R-J1", "signed one-sided differences of J1 agree in sign",
            ">=70% same sign => real spanwise slope => round 10's (a) is REFUTED;  "
            "<=30% => noise => (a) stands",
            f"{100*same:.1f}% of {int(keep.sum())} station-slots  "
            f"(median |fwd| {np.median(np.abs(fw)):.5f}, "
            f"|bwd| {np.median(np.abs(bw)):.5f})",
            "R-J1: REAL SLOPE -- round 10's (a) refuted" if same >= 0.70 else
            ("R-J1: noise -- (a) stands" if same <= 0.30 else "R-J1: UNDECIDABLE"))

    # ---- R-SAT (RECORDED) -------------------------------------------------
    per = {a: np.mean([r["sat"] for r in rows if r["arm"] == a]) for a in J.ARMS}
    cl = np.abs([r["d"] for r in rows if not r["sat"]])
    _record("R-SAT", "saturation rate by arm; |d| on unsaturated readings", "RECORDED",
            "  ".join(f"{a} {100*per[a]:.1f}%" for a in J.ARMS)
            + f";  clean |d| median {np.median(cl):.5f} = {np.median(cl)/0.041:.2f} cell",
            "RECORDED")

    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s")
    print("\n★ kill criterion 5: NO arm is adopted here. Ranking on M6 is not passing "
          "on M1 -- adoption needs NACA0012 M0.80/alpha1.25 (R12).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
