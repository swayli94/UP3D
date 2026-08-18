"""M1-a -- make the shock-position instrument SUB-CELL resolvable.

Binding text: docs/dev_phase_five/20260822-1800-m1a-prereg.md (committed before
this file existed).

★★★ This does NOT fix M1. M1's three gates are untouched and its recorded verdict
-- "not reachable on the current discretisation" -- stands. What this round asks is
narrower: round 18 measured the extractor resolving 0.29 of a cell while being asked
to distinguish a quarter of one, so until the instrument is sub-cell, no fix can be
shown to work.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/m1a_subcell/run.py
"""

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

import run_gs40j_position_primary as J                                  # noqa: E402
# ★ import paths READ from bench/run_gs40j_position_primary.py, not recalled --
# my first attempt guessed pyfp3d.wake_cut and pyfp3d.assembly.geometry, both wrong.
from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared               # noqa: E402

GAMMA = 1.4
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def foot_percell(cent, adj, M, etas, b_semi):
    """Round 18's rule, IMPORTED IN SPIRIT and re-implemented here verbatim from
    bench/run_gs40i_position_and_blastradius.py:114 so both readings run on the
    same field in one process (G-SAMEDATA)."""
    out = {}
    for e_ in etas:
        z0 = e_ * b_semi
        inb = adj & (np.abs(cent[:, 2] - z0) <= 0.02 * b_semi)
        band = inb & (M >= 1.0)
        if band.sum() < 3:
            out[e_] = np.nan
            continue
        x0, x1 = cent[inb, 0].min(), cent[inb, 0].max()
        out[e_] = float((cent[band, 0].max() - x0) / max(x1 - x0, 1e-9))
    return out


def foot_subcell(cent, adj, M, etas, b_semi):
    """★ The sub-cell rule: sort the band by centroid x, take the LAST adjacent
    pair straddling M = 1, and interpolate. The landing point lies BETWEEN two
    centroids, which is what sub-cell means."""
    out, diag = {}, {}
    for e_ in etas:
        z0 = e_ * b_semi
        inb = adj & (np.abs(cent[:, 2] - z0) <= 0.02 * b_semi)
        if (inb & (M >= 1.0)).sum() < 3:
            out[e_] = np.nan
            continue
        idx = np.where(inb)[0]
        o = idx[np.argsort(cent[idx, 0])]
        x, m = cent[o, 0], M[o]
        cross = np.where((m[:-1] >= 1.0) & (m[1:] < 1.0))[0]
        if cross.size == 0:
            out[e_] = np.nan
            continue
        i = int(cross[-1])
        xs = x[i] + (1.0 - m[i]) * (x[i + 1] - x[i]) / (m[i + 1] - m[i])
        x0, x1 = x.min(), x.max()
        out[e_] = float((xs - x0) / max(x1 - x0, 1e-9))
        diag[e_] = (float(x[i]), float(xs), float(x[i + 1]))   # S-MONO
    return out, diag


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    J.guard_no_solve_executed()
    print("  G-FROZEN-LIB  the extractor lives here, pyfp3d/ is untouched")

    rows, arms_rows, mono_bad = [], [], 0
    for geom, label, sp, msh, m_inf in J.STATES:
        if not sp.exists():
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
        # ★ the SAME computation round 18 used (its lines 159-161), not my own
        B, _ = precompute_element_geometry(nodes, elements)
        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)
        M = np.sqrt(np.maximum(mach_number_squared(q2, m_inf, GAMMA), 0.0))
        cent = nodes[elements].mean(axis=1)
        wall_nodes = np.unique(mc.boundary_faces["wall"].reshape(-1))
        isw = np.zeros(len(nodes), bool); isw[wall_nodes] = True
        adj = isw[elements].any(axis=1)
        b_semi = J.B_SEMI      # ★ imported, not re-derived
        etas = J.ETAS          # ★ round 18's own station list, imported

        # ★★ S-TIE is the REGISTERED primary criterion and needs the per-ARM
        # readings, not just the base field. My first implementation computed only
        # S-MONO and S-SHIFT on the unrefit field, which does not answer it.
        import run_gs40c_coherence_ladder as L
        from pyfp3d.mesh.metrics import build_face_adjacency
        fn, _ = build_face_adjacency(elements)
        Pp, pel = {}, {}
        for k in (0, 1, 2, 3):
            po, pi = L.khop_patches(elements, fn, k)
            Pp[f"P{k}"] = L.patch_nodes(elements, po, pi); pel[f"P{k}"] = (po, pi)
        garm = {n: L.refit(nodes, elements, phi, *Pp[n]) for n in J.LADDER}
        for sm, p in J.SMOOTH.items():
            garm[sm] = L.average_over_patch(nodes, elements, gA, *pel[p])
        for arm in J.ARMS:
            Ma = np.sqrt(np.maximum(mach_number_squared(
                (garm[arm] * garm[arm]).sum(axis=1), m_inf, GAMMA), 0.0))
            a_pc = foot_percell(cent, adj, Ma, etas, b_semi)
            a_sc, _ = foot_subcell(cent, adj, Ma, etas, b_semi)
            for e_ in etas:
                if np.isfinite(a_pc.get(e_, np.nan)) and np.isfinite(a_sc.get(e_, np.nan)):
                    arms_rows.append({"state": label, "arm": arm, "eta": e_,
                                      "x_percell": a_pc[e_], "x_subcell": a_sc[e_]})

        pc = foot_percell(cent, adj, M, etas, b_semi)
        sc, diag = foot_subcell(cent, adj, M, etas, b_semi)
        # a cell, in the same normalised units, for S-SHIFT
        h_cell = float(np.median(np.diff(np.sort(cent[adj, 0])))) if adj.sum() else np.nan
        for e_ in etas:
            if np.isnan(pc.get(e_, np.nan)) or np.isnan(sc.get(e_, np.nan)):
                continue
            if e_ in diag:
                lo, xs, hi = diag[e_]
                if not (min(lo, hi) - 1e-12 <= xs <= max(lo, hi) + 1e-12):
                    mono_bad += 1
            rows.append({"state": label, "geom": geom, "eta": e_, "phi_sha": sha,
                         "x_percell": pc[e_], "x_subcell": sc[e_],
                         "shift": sc[e_] - pc[e_]})
        print(f"  [{label}] phi sha {sha}  stations {len(etas)}  "
              f"median |shift| {np.nanmedian([abs(sc[e]-pc[e]) for e in etas if e in sc and e in pc and np.isfinite(sc[e]) and np.isfinite(pc[e])]):.5f}")

    if not rows:
        print("no usable states"); return 1
    with open(os.path.join(RESULTS, "subcell.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # ---- S-TIE: how many states have a TIED set of arm readings? ------------
    import collections
    def tie_count(key):
        by = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in arms_rows:
            by[r["state"]][r["arm"]].append(r[key])
        n = 0
        for st, a in by.items():
            sig = {arm: tuple(round(v, 9) for v in vs) for arm, vs in a.items()}
            if len(set(sig.values())) < len(sig):
                n += 1
        return n, len(by)
    tp, ns = tie_count("x_percell"); ts, _ = tie_count("x_subcell")
    _record("S-TIE", "states where two arms give an IDENTICAL footprint vector",
            "per-cell rule's count is the target to beat",
            f"per-cell {tp}/{ns}  ->  sub-cell {ts}/{ns}",
            "S-TIE PASS" if ts <= 1 else f"S-FAIL -- ties are NOT a resolution artefact")
    with open(os.path.join(RESULTS, "arms.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(arms_rows[0].keys()))
        w.writeheader(); w.writerows(arms_rows)

    sh = np.array([abs(r["shift"]) for r in rows])
    _record("S-MONO", "the interpolated x lies between the two centroids",
            "every station", f"{mono_bad} violations of {len(rows)}",
            "S-MONO PASS" if mono_bad == 0 else "S-FAIL")
    _record("S-SHIFT", "|sub-cell - per-cell|, in the normalised chord the "
            "per-cell rule uses", "RECORDED (expected ~0.3 cell)",
            f"median {np.median(sh):.5f}, max {sh.max():.5f}", "RECORDED")
    print(f"\n  {len(rows)} station-readings over "
          f"{len({r['state'] for r in rows})} states, {time.perf_counter()-t0:.1f} s")
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print("\n★ This does NOT fix M1: its three gates are untouched and the recorded "
          "'not reachable on the current discretisation' verdict stands.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
