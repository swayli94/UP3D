"""M1-e -- what carries D's 0.36-cell station-to-station scatter.

Binding text: phases/p5/docs/dev_phase_five/20260822-2200-m1e-prereg.md (committed before
this file existed).

Three candidates with different consequences: J1's own noise (a), the footprint's
+-2% b_semi station band (b), or a genuine smooth spanwise structure (c). (a) and (b)
mean D can be repaired; (c) means aggregating seven stations into a median is the
wrong construction and the metric's FORM is what needs replacing.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/m1e_scatter/run.py
"""

import collections
import csv
import hashlib
import itertools
import math
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
BANDS = (0.01, 0.02, 0.04)          # E-BAND; 0.02 is the shipped value
DETA = 0.01                          # E-J1
PERM_MIN = 5                         # a p-value on < 5 stations is not worth reading
#: G-REPRO, round 9's committed base readings
ANCHOR = {"bias": 0.00257, "mad": 0.01459, "n": 627}
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def foot_band(cent, adj, M, etas, b_semi, half):
    """I.foot_from_field with the station half-width exposed.

    ★ G-BAND-IDENT asserts this reproduces I.foot_from_field BITWISE at half=0.02,
    so the band arm is a re-parameterisation and not a re-implementation.
    """
    out = {}
    for e_ in etas:
        z0 = e_ * b_semi
        inb = adj & (np.abs(cent[:, 2] - z0) <= half * b_semi)
        band = inb & (M >= 1.0)
        if band.sum() < 3:
            out[e_] = np.nan; continue
        x0, x1 = cent[inb, 0].min(), cent[inb, 0].max()
        out[e_] = float((cent[band, 0].max() - x0) / max(x1 - x0, 1e-9))
    return out


def perm_p(vals):
    """Exact smoothness p-value: T = sum of squared consecutive differences over the
    span-ordered stations, against ALL n! orderings. No random numbers."""
    v = np.asarray(vals, float)
    t_obs = float(((v[1:] - v[:-1]) ** 2).sum())
    n_le = 0
    for p in itertools.permutations(range(len(v))):
        w = v[list(p)]
        if float(((w[1:] - w[:-1]) ** 2).sum()) <= t_obs + 1e-15:
            n_le += 1
    return n_le / float(math.factorial(len(v))), t_obs


def mad(x):
    x = np.asarray(x, float)
    return float(np.median(np.abs(x - np.median(x))))


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    J.guard_no_solve_executed()
    print("  G-FROZEN-LIB  extractors live in bench/, pyfp3d/ is untouched")

    base, bandrows, j1rows, permrows, pairrows = [], [], [], [], []
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

        # ---- J1 at eta and at eta +- DETA (E-J1); footprint untouched ---------
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
            print(f"  [{label:11}] EXCLUDED-NO-WALL-SHOCK"); continue

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
        cent = nodes[elements].mean(axis=1)
        wn = np.unique(wall_faces.reshape(-1))
        isw = np.zeros(len(nodes), bool); isw[wn] = True
        adj = isw[elements].any(axis=1)

        dfun = {}
        for n in J.ARMS:
            M = np.sqrt(np.maximum(mach_number_squared(
                (g[n] * g[n]).sum(axis=1), m_inf, GAMMA), 0.0))
            fb = {h: foot_band(cent, adj, M, J.ETAS, J.B_SEMI, h) for h in BANDS}
            #: G-BAND-IDENT -- the shipped width must be the shipped function
            ref = I.foot_from_field(nodes, elements, M, wall_faces, J.ETAS, J.B_SEMI)
            for e_ in J.ETAS:
                a, b = fb[0.02][e_], ref[e_]
                assert (np.isnan(a) and np.isnan(b)) or a == b, \
                    f"G-BAND-IDENT {label} {n} {e_}: {a} vs {b}"
            dd = {}
            for e_ in J.ETAS:
                if j1[e_] is None or not np.isfinite(fb[0.02][e_]):
                    continue
                dd[e_] = fb[0.02][e_] - j1[e_]
                base.append({"state": label, "arm": n, "eta": e_, "d": dd[e_]})
            dfun[n] = dd
            # E-BAND, on the set valid at ALL widths (G-FIXEDSET)
            keep = [e_ for e_ in J.ETAS if j1[e_] is not None
                    and all(np.isfinite(fb[h][e_]) for h in BANDS)]
            for h in BANDS:
                for e_ in keep:
                    bandrows.append({"state": label, "arm": n, "eta": e_, "half": h,
                                     "d": fb[h][e_] - j1[e_]})
            if n == "P0":
                print(f"  [{label:11}] sha {sha}  band-fixed {len(keep)}/7", end="")
            # E-PERM (binding on P0)
            v = [dd[e_] for e_ in sorted(dd)]
            if len(v) >= PERM_MIN:
                p_, t_ = perm_p(v)
                permrows.append({"state": label, "arm": n, "n_st": len(v),
                                 "p": p_, "T": t_})
                if n == "P0":
                    print(f"   E-PERM n={len(v)} p={p_:.4f}")
            elif n == "P0":
                print(f"   E-PERM SKIPPED (n={len(v)} < {PERM_MIN})")
        # E-J1, on the set valid at all three eta (G-FIXEDSET)
        keep1 = [e_ for e_ in J.ETAS
                 if all(x[e_] is not None for x in (j1, j1m, j1p))]
        for e_ in keep1:
            j1rows.append({"state": label, "eta": e_,
                           "dminus": abs(j1m[e_] - j1[e_]),
                           "dplus": abs(j1p[e_] - j1[e_])})
        # E-PAIRED (RECORDED ONLY -- must not be used to rank)
        Dm = {n: float(np.median([abs(x) for x in dfun[n].values()]))
              for n in J.LADDER if len(dfun[n]) >= J.MIN_STATIONS}
        if len(Dm) == len(J.LADDER):
            srt = sorted(Dm, key=lambda a: Dm[a]); w_, r_ = srt[0], srt[1]
            com = [e_ for e_ in J.ETAS if e_ in dfun[w_] and e_ in dfun[r_]]
            nw = sum(abs(dfun[w_][e_]) < abs(dfun[r_][e_]) for e_ in com)
            pairrows.append({"state": label, "winner": w_, "runner": r_,
                             "n_common": len(com), "n_winner_closer": nw})

    for name, rr in (("base_d.csv", base), ("band.csv", bandrows),
                     ("j1_eta.csv", j1rows), ("perm.csv", permrows),
                     ("paired.csv", pairrows)):
        if rr:
            with open(os.path.join(RESULTS, name), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rr[0].keys()))
                w.writeheader(); w.writerows(rr)

    # ---- G-REPRO ----------------------------------------------------------
    dv = np.array([r["d"] for r in base])
    bias, m_ = float(np.median(dv)), mad(dv)
    ok = (abs(bias - ANCHOR["bias"]) < 5e-5 and abs(m_ - ANCHOR["mad"]) < 5e-5
          and len(dv) == ANCHOR["n"])
    _record("G-REPRO", "base d reproduces round 9",
            f"bias {ANCHOR['bias']}, MAD {ANCHOR['mad']}, n {ANCHOR['n']}",
            f"bias {bias:+.5f}, MAD {m_:.5f}, n {len(dv)}",
            "G-REPRO PASS" if ok else "★ G-REPRO FAIL -- kill criterion 1")
    if not ok:
        print("\n★★ kill criterion 1: check the recomputation before concluding.")
        with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
            w.writerows(SUMMARY)
        return 2

    # ---- E-PERM (binding on P0) -------------------------------------------
    p0 = [r for r in permrows if r["arm"] == "P0"]
    sig = sum(r["p"] < 0.05 for r in p0)
    allp = [r for r in permrows if r["arm"] in J.LADDER]
    _record("E-PERM", "d(eta) smoother than a random reordering (exact, all n!)",
            ">=5/13 => (c) structure;  <=2/13 => not (c);  else UNDECIDABLE",
            f"P0 {sig}/{len(p0)} states with p<0.05 (null expects "
            f"{0.05*len(p0):.2f});  pooled LADDER "
            f"{sum(r['p'] < 0.05 for r in allp)}/{len(allp)}",
            "E-PERM: (c) STRUCTURE" if sig >= 5 else
            ("E-PERM: not (c)" if sig <= 2 else "E-PERM: UNDECIDABLE"))

    # ---- E-BAND -----------------------------------------------------------
    stats = {}
    for h in BANDS:
        v = np.array([r["d"] for r in bandrows if r["half"] == h and r["arm"] == "P0"])
        stats[h] = (float(np.median(v)), mad(v), len(v))
    rel = abs(stats[0.01][1] - stats[0.02][1]) / max(stats[0.02][1], 1e-12)
    _record("E-BAND", "MAD of d vs station half-width (P0, set fixed across widths)",
            ">=30% drop at 0.01 => (b);  <10% change => not (b)",
            "  ".join(f"h={h}: bias {stats[h][0]:+.5f} MAD {stats[h][1]:.5f} "
                      f"n={stats[h][2]}" for h in BANDS)
            + f"   |dMAD|/MAD = {rel:.1%}",
            "E-BAND: (b) CARRIES IT" if rel >= 0.30 else
            ("E-BAND: not (b)" if rel < 0.10 else "E-BAND: UNDECIDABLE"))

    # ---- E-J1 -------------------------------------------------------------
    dj = np.array([r["dminus"] for r in j1rows] + [r["dplus"] for r in j1rows])
    md = float(np.median(dj))
    _record("E-J1", f"|J1(eta+-{DETA}) - J1(eta)| (set fixed across the three eta)",
            ">=0.007 => (a);  <=0.002 => not (a)",
            f"median {md:.5f}, mean {dj.mean():.5f}, n={len(dj)} "
            f"(of {14*len({r['state'] for r in j1rows})} station-slots)",
            "E-J1: (a) CARRIES IT" if md >= 0.007 else
            ("E-J1: not (a)" if md <= 0.002 else "E-J1: UNDECIDABLE"))

    # ---- E-VAR (RECORDED, nearly a deduction -- context only) --------------
    byse = collections.defaultdict(dict)
    for r in base:
        byse[(r["state"], r["eta"])][r["arm"]] = r["d"]
    across_arm = [np.std(list(v.values())) for v in byse.values() if len(v) >= 4]
    byar = collections.defaultdict(list)
    for r in base:
        byar[(r["state"], r["arm"])].append(r["d"])
    across_eta = [np.std(v) for v in byar.values() if len(v) >= 4]
    _record("E-VAR", "spread across eta vs across arms (declared in advance to be "
            "nearly a deduction; context only)", "RECORDED",
            f"median sigma across eta {np.median(across_eta):.5f} vs across arms "
            f"{np.median(across_arm):.5f}", "RECORDED")

    # ---- E-PAIRED (RECORDED, ranking forbidden) ---------------------------
    tot = sum(r["n_winner_closer"] for r in pairrows)
    ntot = sum(r["n_common"] for r in pairrows)
    _record("E-PAIRED", "stations where the D-argmin arm is closer than the "
            "runner-up (RECORDED -- no significance machinery registered)",
            "RECORDED; ranking claims are FORBIDDEN this round",
            f"{tot}/{ntot} station-comparisons over {len(pairrows)} states "
            f"(chance {ntot/2:.0f})", "RECORDED")

    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s")
    print("\n★ This does NOT fix M1: its three gates are untouched and the recorded "
          "'not reachable on the current discretisation' verdict stands.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
