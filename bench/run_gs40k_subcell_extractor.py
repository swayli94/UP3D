"""
Sub-cell footprint extraction: change the RULER, not the thing being measured.

Pre-registration: docs/dev_phase_four/20260818-1700-subcell-extractor-prereg.md
Previous round:   docs/dev_phase_four/20260818-1500-position-primary-verdict.md

Same 13 states, same 7 arms, same J1 judge, same D. The ONLY variable is the
footprint extractor: per-cell (resolution h) versus sub-cell (interpolate the
M = 1 crossing on supersonic/subsonic face pairs).

★ Second meaning, and it is the point: if sub-cell interpolation makes the
position metric discriminate where the per-cell one could not, that is a direct
demonstration that SUB-CELL REPRESENTATION MATTERS -- which is what registered
item 1 is about. If it does not, my explanation for POS-MIXED was wrong.
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

import run_gs40c_coherence_ladder as L                              # noqa: E402
import run_gs40i_position_and_blastradius as I                      # noqa: E402
import run_gs40j_position_primary as J                              # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import build_face_adjacency, precompute_element_geometry  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402

OUT = REPO / "bench/gate_results/gs40k_subcell.csv"
GAMMA = 1.4
ETAS, LADDER, ARMS = J.ETAS, J.LADDER, J.ARMS
MIN_STATIONS, MIN_STATES = 3, 4
TIE_MAX, GATE_S = 2, 20 * 60
G_SUB_SLACK = 0.002


def foot_subcell(nodes, elements, M, wall_faces, face_nb, etas, b_semi):
    """Sub-cell footprint: interpolate the M = 1 crossing on face pairs.

    For a face between e (M_e) and e' (M_e') straddling 1, the crossing sits at
    t = (1 - M_e') / (M_e - M_e') along the centroid-to-centroid segment. The
    footprint is the largest x/c among those crossings inside the station band.

    ★ Registered in advance: M_e == M_e' means no crossing on that face and it is
    SKIPPED (a zero denominator is not a crossing).
    """
    cent = nodes[elements].mean(axis=1)
    isw = np.zeros(len(nodes), bool); isw[np.unique(wall_faces.reshape(-1))] = True
    adj = isw[elements].any(axis=1)
    e0 = np.repeat(np.arange(len(elements)), 4)
    e1 = face_nb.reshape(-1)
    m = (e1 >= 0) & (e1 > e0) & adj[e0] & adj[e1]
    a, b = e0[m], e1[m]
    Ma, Mb = M[a], M[b]
    cross = ((Ma >= 1.0) & (Mb < 1.0)) | ((Mb >= 1.0) & (Ma < 1.0))
    a, b, Ma, Mb = a[cross], b[cross], Ma[cross], Mb[cross]
    den = Ma - Mb
    ok = np.abs(den) > 0.0
    a, b, Ma, Mb, den = a[ok], b[ok], Ma[ok], Mb[ok], den[ok]
    t = np.clip((1.0 - Mb) / den, 0.0, 1.0)
    pt = cent[b] + t[:, None] * (cent[a] - cent[b])
    out, nsta = {}, {}
    for e_ in etas:
        z0 = e_ * b_semi
        band = np.abs(pt[:, 2] - z0) <= 0.02 * b_semi
        wb = adj & (np.abs(cent[:, 2] - z0) <= 0.02 * b_semi)
        if band.sum() < 1 or wb.sum() < 2:
            out[e_] = np.nan; nsta[e_] = 0; continue
        x0, x1 = cent[wb, 0].min(), cent[wb, 0].max()
        out[e_] = float((pt[band, 0].max() - x0) / max(x1 - x0, 1e-9))
        nsta[e_] = int(band.sum())
    return out, nsta


def main():
    t0 = time.perf_counter()
    J.guard_no_solve_executed()

    rows, per = [], {}
    n_tie_cell = n_tie_sub = 0
    for geom, label, sp, msh, m_inf in J.STATES:
        if not sp.exists():
            continue
        man = read_manifest(msh)
        assert man is None or mesh_fingerprint(msh)["sha256"] == man["sha256"], \
            f"G-MESH {label}"
        d = np.load(sp)
        if not bool(d["conv"]) or int(d["nlim"]) or int(d["nflr"]):
            continue
        phi = np.asarray(d["phi"], float)
        mc, wc = cut_wake(read_mesh(msh))
        nodes, elements = mc.nodes, mc.elements
        if len(phi) != len(nodes):
            continue
        B, _ = precompute_element_geometry(nodes, elements)
        fn, _ = build_face_adjacency(elements)
        wall_faces = mc.boundary_faces["wall"]

        j1 = {}
        for e_ in ETAS:
            try:
                cur = section_cp_curve(mc, phi, eta=e_, b_semi=B_SEMI, m_inf=m_inf)
                x = shock_report(cur, m_inf)["upper"].get("x_shock")
                j1[e_] = None if x is None or not np.isfinite(x) else float(x)
            except Exception:                                      # noqa: BLE001
                j1[e_] = None
        if sum(v is not None for v in j1.values()) < MIN_STATIONS:
            print(f"[{label}] EXCLUDED-NO-WALL-SHOCK"); continue

        P, pel = {}, {}
        for k in (0, 1, 2, 3):
            po, pi = L.khop_patches(elements, fn, k)
            P[f"P{k}"] = L.patch_nodes(elements, po, pi); pel[f"P{k}"] = (po, pi)
        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)
        g = {n: L.refit(nodes, elements, phi, *P[n]) for n in LADDER}
        for s, p in J.SMOOTH.items():
            g[s] = L.average_over_patch(nodes, elements, gA, *pel[p])
        assert np.max(np.abs(g["P0"] - gA)) < 1e-10, f"G-A {label}"

        Dc, Ds, nst = {}, {}, {}
        for n in ARMS:
            m2 = mach_number_squared((g[n] * g[n]).sum(axis=1), m_inf, GAMMA)
            M = np.sqrt(np.maximum(m2, 0.0))
            fc = I.foot_from_field(nodes, elements, M, wall_faces, ETAS, B_SEMI)
            fs, ns = foot_subcell(nodes, elements, M, wall_faces, fn, ETAS, B_SEMI)
            nst[n] = sum(1 for e_ in ETAS if j1[e_] is not None and np.isfinite(fs[e_]))
            for tag, ff, store in (("cell", fc, Dc), ("sub", fs, Ds)):
                dd = [abs(ff[e_] - j1[e_]) for e_ in ETAS
                      if j1[e_] is not None and np.isfinite(ff[e_])]
                store[n] = float(np.median(dd)) if len(dd) >= MIN_STATIONS else np.nan
        if any(np.isnan(Ds[n]) or np.isnan(Dc[n]) for n in LADDER):
            print(f"[{label}] EXCLUDED: an arm has < {MIN_STATIONS} usable stations")
            continue

        #: exact ties only -- no tolerance, so "tie" is not a number I choose
        tc = len({round(Dc[n], 12) for n in LADDER}) < len(LADDER)
        ts = len({round(Ds[n], 12) for n in LADDER}) < len(LADDER)
        n_tie_cell += tc; n_tie_sub += ts
        kc = min(LADDER, key=lambda a: Dc[a])
        ks = min(LADDER, key=lambda a: Ds[a])
        per[label] = dict(Dc=Dc, Ds=Ds, kc=kc, ks=ks, tie_c=tc, tie_s=ts)
        print(f"[{label:11}] cell D: " + " ".join(f"{Dc[n]:.4f}" for n in LADDER)
              + f"  k={kc}{' TIE' if tc else ''}")
        print(f"{'':14}sub  D: " + " ".join(f"{Ds[n]:.4f}" for n in LADDER)
              + f"  k={ks}{' TIE' if ts else ''}   stations "
              + "/".join(str(nst[n]) for n in LADDER))
        for n in ARMS:
            rows.append(dict(state=label, geom=geom, m_inf=m_inf, arm=n,
                             D_cell=round(Dc[n], 5), D_sub=round(Ds[n], 5),
                             n_stations_sub=nst[n], k_cell=kc, k_sub=ks))

    # ------------------------------------------------------------- guards ----
    print("\n" + "=" * 76)
    n = len(per)
    if n < MIN_STATES:
        print(f"VERDICT: SUB-UNDEF (only {n} states)"); return
    better = sum(1 for v in per.values()
                 if v["Ds"]["P0"] <= v["Dc"]["P0"] + G_SUB_SLACK)
    print(f"G-SUB   sub-cell no worse than per-cell on the CONTROL arm: "
          f"{better}/{n} states  {'PASS' if better > n / 2 else '** FAIL **'}")
    assert better > n / 2, "G-SUB: the interpolation is suspect -- kill criterion 1"

    medc = float(np.median([v["Dc"][a] for v in per.values() for a in LADDER]))
    meds = float(np.median([v["Ds"][a] for v in per.values() for a in LADDER]))
    print(f"        median D: per-cell {medc:.4f} -> sub-cell {meds:.4f} c "
          f"({100*(meds-medc)/medc:+.1f}%)")

    votes = {}
    for lb, v in per.items():
        votes.setdefault("MIXED" if v["tie_s"] else v["ks"], []).append(lb)
    win = max(votes, key=lambda k: len(votes[k]))
    maj = len(votes[win]) > n / 2 and win != "MIXED"
    print(f"\nexact ties across arms: per-cell {n_tie_cell}/{n}  ->  "
          f"sub-cell {n_tie_sub}/{n}   (SUB needs <= {TIE_MAX})")
    print("k_pos* votes (sub-cell): " + "  ".join(
        f"{k}:{len(v)}" for k, v in sorted(votes.items())))
    if n_tie_sub > TIE_MAX:
        verdict = ("SUB-SAME (ties persist) -- my POS-MIXED explanation was WRONG "
                   "-- kill criterion 2")
    elif maj:
        verdict = f"SUB-SHARP (ties gone AND k_pos* has a majority at {win})"
    else:
        verdict = "SUB-NOISY (ties gone, but k_pos* still has no majority)"
    print("=" * 76)
    print(f"VERDICT: {verdict}")

    agree = sum(1 for v in per.values() if v["ks"] == v["kc"])
    print(f"\nRECORDED: sub-cell and per-cell pick the same arm on {agree}/{n} states")
    allsub = [v["Ds"][a] for v in per.values() for a in LADDER]
    if max(allsub) < 0.002:
        print("RECORDED: ★ every arm is within 0.2% chord -- differences are tiny")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict=verdict, ties_cell=n_tie_cell, ties_sub=n_tie_sub,
             median_D_cell=medc, median_D_sub=meds, n_states=n,
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 4"


if __name__ == "__main__":
    main()
