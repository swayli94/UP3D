"""F: does the level-to-level change align with the FOOT displacement or the PEAK location?
G: which spanwise quantity also TURNS OVER near eta 0.90?  Both ZERO-SOLVE.

Pre-registered in docs/dev_phase_three/20260815-1400-alignment-and-spanwise-prereg.md.

F ★★ The sampling problem is faced first, not discovered afterwards: the foot-swept interval is 0.002 to
0.020 chord wide while each station holds only 5-10 experimental points across [0, 0.15), so it would
contain 0-2 of them -- the trap that returned n = 0 two rounds ago. So the PRIMARY form evaluates the
difference on a DENSE GRID of the computed curves (the mechanism question is about the SOLUTION, not about
the error metric), and the experimental-abscissae form is a consistency check only.
★★★ D(x) on the dense grid and d_self on the experimental points are DIFFERENT POINT SETS: only position
conclusions are compared, never the numbers.
★★ The two candidate windows are forced EQUAL WIDTH, and a station whose windows overlap by more than 50 %
is declared to have no discriminating power and is excluded from the hit count -- registered in advance,
because unequal generosity is the defect the previous round caught in itself.

G ★ The discriminator is "which candidate also TURNS OVER": d_self rises then falls, peaking at eta 0.90,
while most spanwise quantities climb monotonically to the tip. No correlation coefficients are computed --
over seven stations they have no discriminating power, and this session has twice been misled by exactly
that.

★ G-Z: no solver is called. `cut_wake` IS used, for station_z only -- a mesh operation, declared in the
registration section 2, not a solve.

Outputs (TRACKED): bench/gate_results/task3_alignment_spanwise.csv
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                              # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                             # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at                    # noqa: E402
from run_m3_budget import ETAS, M_INF, parse_experiment               # noqa: E402
import run_task3_fixed_budget as FB                                   # noqa: E402
import run_task3_nonconv_discriminator as ND                          # noqa: E402
import run_task3_shape_and_taper as ST                                # noqa: E402
import run_task3_tipstations_mechanism as TM                          # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_alignment_spanwise.csv")
PAIR = ("S2", "S4")
#: the seven SEEN d_self values live in the shape/taper round (TM carries only five)
GR7 = ST.GR7
UNMASKED = (0.20, 0.44, 0.65, 0.80, 0.90)
N_GRID, BAND_HI = 301, 0.15
OVERLAP_MAX, F_HITS = 0.50, 3
ETA_TOL = 0.10                     #: one (coarse-side) station spacing near eta 0.90
S4_MESH = "rnd_0.015_0.0075_0.0075_1.8.msh"


def dense_D(c2, c4, eta, side="upper"):
    """D(x) on a dense grid inside the COMMON x-range of both curves (G-GRID)."""
    x2, p2 = np.asarray(c2[eta][f"x_{side}"]), np.asarray(c2[eta][f"cp_{side}"])
    x4, p4 = np.asarray(c4[eta][f"x_{side}"]), np.asarray(c4[eta][f"cp_{side}"])
    o2, o4 = np.argsort(x2), np.argsort(x4)
    x2, p2, x4, p4 = x2[o2], p2[o2], x4[o4], p4[o4]
    lo = max(x2.min(), x4.min(), 0.0)
    hi = min(x2.max(), x4.max(), BAND_HI)
    if not (hi > lo):
        return None
    g = np.linspace(lo, hi, N_GRID)
    return dict(x=g, D=np.interp(g, x2, p2) - np.interp(g, x4, p4), lo=lo, hi=hi)


def overlap_frac(a, b):
    """|A cap B| / |A| for two intervals of equal width."""
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    inter = max(0.0, hi - lo)
    return inter / (a[1] - a[0]) if a[1] > a[0] else 0.0


def main():
    src = open(__file__).read().split(chr(34) * 3, 2)[2]
    for forbidden in ("solve" + "_newton", "run_m3_budget." + "solve"):
        assert forbidden not in src, f"G-Z: no solver ({forbidden})"
    print("★ G-Z: zero solves (cut_wake is used for station_z only -- declared in the registration)\n")
    exp = parse_experiment()
    cp_star = ND.cp_critical(M_INF)
    curves = {}
    for t in ("S0", "S2", "S4"):
        c, _r = TM.solve_all_stations(t, exp)
        if c is None:
            print(f"★ {t} curves unavailable -- STOP")
            return 1
        curves[t] = c

    #: --- G-R ------------------------------------------------------------------------------------
    print("=== G-R: the seven d_self values must still reproduce bit-for-bit ===")
    st = {}
    bad = []
    for eta in ETAS:
        got = TM.band_points(curves, exp, eta)
        if got[0] is None:
            bad.append(eta)
            continue
        d = got[2][PAIR[0]] - got[2][PAIR[1]]
        st[eta] = dict(x=got[0], d=d, n=len(d))
        r = float(np.sqrt(np.mean(d ** 2)))
        ok = abs(r - GR7[eta]) <= 5e-7
        print(f"  eta {eta:.2f}: {r:.6f} vs {GR7[eta]:.6f}  {'PASS' if ok else '★ FAIL'}")
        if not ok:
            bad.append(eta)
    if bad:
        print(f"  -> ★ G-R FAIL at {bad}. STOP (kill clause 1).")
        return 1
    print("  -> PASS")

    #: --- F --------------------------------------------------------------------------------------
    print("\n=== F: alignment -- foot-swept window A vs EQUAL-WIDTH peak window B ===")
    print(f"  {'eta':>6}{'x_argmax':>10}{'x_cen':>9}{'w':>9}{'window A':>18}{'window B':>18}"
          f"{'ovl':>6}{'in A':>6}{'in B':>6}{'expN':>6}   verdict")
    rows, hitA, hitB, usable = [], 0, 0, 0
    for eta in ETAS:
        if eta not in st:
            continue
        dd = dense_D(curves[PAIR[0]], curves[PAIR[1]], eta)
        g2 = TM.geometry(curves[PAIR[0]][eta], cp_star)
        g4 = TM.geometry(curves[PAIR[1]][eta], cp_star)
        if dd is None or g2[0] is None or g4[0] is None:
            print(f"  {eta:>6.2f}   (grid or geometry unavailable -- skipped)")
            continue
        A = (min(g2[0], g4[0]), max(g2[0], g4[0]))
        w = A[1] - A[0]
        B = (g4[1] - 0.5 * w, g4[1] + 0.5 * w)
        m = np.abs(dd["D"])
        xam = float(dd["x"][int(np.argmax(m))])
        mass = dd["D"] ** 2
        xcen = float(np.sum(dd["x"] * mass) / np.sum(mass))
        ovl = overlap_frac(A, B)
        inA, inB = A[0] <= xam <= A[1], B[0] <= xam <= B[1]
        expN = int(((st[eta]["x"] >= A[0]) & (st[eta]["x"] <= A[1])).sum())
        disc = ovl <= OVERLAP_MAX
        if eta in UNMASKED and disc:
            usable += 1
            hitA += int(inA)
            hitB += int(inB)
        verdict = ("★ no discriminating power (overlap > 50 %)" if not disc
                   else "in A" if inA and not inB else "in B" if inB and not inA
                   else "in both" if inA and inB else "in neither")
        mark = "  " if eta in UNMASKED else " ★"
        print(f"{mark}{eta:>6.2f}{xam:>10.5f}{xcen:>9.5f}{w:>9.5f}"
              f"{f'[{A[0]:.4f},{A[1]:.4f}]':>18}{f'[{B[0]:.4f},{B[1]:.4f}]':>18}"
              f"{ovl:>6.2f}{str(inA):>6}{str(inB):>6}{expN:>6}   {verdict}")
        rows.append(dict(kind="F", eta=eta, x_argmax=xam, x_cen=xcen, w=w, A_lo=A[0], A_hi=A[1],
                         B_lo=B[0], B_hi=B[1], overlap=ovl, in_A=inA, in_B=inB,
                         exp_pts_in_A=expN, grid_lo=dd["lo"], grid_hi=dd["hi"],
                         x_peak_S4=g4[1], x_foot_S2=g2[0], x_foot_S4=g4[0]))

    print(f"\n  usable (discriminating) unmasked stations: {usable}/5   in A: {hitA}   in B: {hitB}")
    print(f"  ★ experimental points inside window A: "
          f"{[r['exp_pts_in_A'] for r in rows if r['eta'] in UNMASKED]} "
          f"-- the auxiliary form is {'unusable (as predicted)' if all(r['exp_pts_in_A'] <= 2 for r in rows if r['eta'] in UNMASKED) else 'usable'}")
    if usable < F_HITS:
        print(f"  -> ★ F UNDEFINED: fewer than {F_HITS} discriminating stations (kill clause 3).")
    elif hitA >= F_HITS and hitA > hitB:
        print("  -> ★★ F-FOOT  the change aligns with the FOOT-SWEPT interval ⇒ the synthesis's")
        print("     prediction is confirmed; the target is the sub-cell sonic interface.")
    elif hitB >= F_HITS and hitB > hitA:
        print("  -> ★★ F-PEAK  the change aligns with the PEAK LOCATION ⇒ ★★★ MY PREVIOUS ROUND'S")
        print("     SYNTHESIS IS WRONG and the verdict must say so and reorganise the readings.")
    elif hitA >= F_HITS and hitB >= F_HITS:
        print("  -> F-BOTH  both windows hit ⇒ RECORDED; the test cannot separate them here.")
    else:
        print("  -> ★ F-NEITHER  the change sits in a THIRD place ⇒ RECORDED, and informative.")

    #: --- G --------------------------------------------------------------------------------------
    print("\n=== G: which spanwise quantity also TURNS OVER near eta 0.90? ===")
    dprof = {e: float(np.sqrt(np.mean(st[e]["d"] ** 2))) for e in st}
    cands = {}
    cands["x_foot(S4)"] = {e: TM.geometry(curves["S4"][e], cp_star)[0] for e in st}
    cands["|Cp_min|(S4)"] = {e: abs(TM.geometry(curves["S4"][e], cp_star)[2]) for e in st}
    cands["|dCp_min|"] = {e: abs(TM.geometry(curves["S2"][e], cp_star)[2]
                                 - TM.geometry(curves["S4"][e], cp_star)[2]) for e in st}
    cands["chord c(eta)"] = {e: float(chord_at(e * B_SEMI)) for e in st}
    gpath = os.path.join(FB.SCRATCH, "phi_S4.npz")
    mpath = os.path.join(FB.SCRATCH, S4_MESH)
    if os.path.exists(gpath) and os.path.exists(mpath):
        gam = np.asarray(np.load(gpath)["gamma"])
        _mc, wc = cut_wake(read_mesh(mpath))
        z = np.asarray(wc.station_z)
        o = np.argsort(z)
        z, gam = z[o], gam[o]
        eta_z = z / B_SEMI
        cands["Gamma(z)"] = {e: float(np.interp(e, eta_z, gam)) for e in st}
        dg = np.gradient(gam, z)
        cands["|dGamma/dz|"] = {e: abs(float(np.interp(e, eta_z, dg))) for e in st}
        cands["Gamma/c"] = {e: float(np.interp(e, eta_z, gam)) / float(chord_at(e * B_SEMI))
                            for e in st}
        print(f"  [gamma loaded: {len(gam)} stations, eta {eta_z.min():.3f}-{eta_z.max():.3f}]")
    else:
        print("  ★ phi/gamma or mesh cache missing -> the Gamma candidates are UNDEFINED "
              "(kill clause 4; no re-solve)")

    print(f"\n  {'candidate':16}{'argmax eta':>12}{'turns over?':>13}{'|argmax-0.90|':>15}"
          f"{'matches?':>10}   normalised profile")
    matches = []
    for name, prof in [("d_self (target)", dprof)] + sorted(cands.items()):
        vals = [prof[e] for e in ETAS if e in prof]
        ets = [e for e in ETAS if e in prof]
        am = ets[int(np.argmax(vals))]
        turns = am != ets[-1]
        dist = abs(am - 0.90)
        ok = turns and dist <= ETA_TOL
        if name != "d_self (target)" and ok:
            matches.append(name)
        mx = max(vals)
        prof_s = " ".join(f"{v / mx:.2f}" for v in vals)
        print(f"  {name:16}{am:>12.2f}{str(turns):>13}{dist:>15.2f}"
              f"{('★ YES' if ok else 'no'):>10}   {prof_s}")
        rows.append(dict(kind="G", candidate=name, argmax_eta=am, turns_over=turns,
                         dist_to_090=dist, matches=ok,
                         **{f"prof_eta{e:g}": prof[e] for e in ets}))
    print(f"\n  matching candidates: {matches if matches else 'NONE'}")
    if not matches:
        print("  -> ★★ G-NONE  no available spanwise quantity turns over near eta 0.90 ⇒ the")
        print("     turnover is NOT explained by loading, chord, foot position or peak depth.")
        print("     That is a NEW question, registered as such.")
    elif len(matches) == 1:
        print(f"  -> ★★ G-MATCH  {matches[0]} ⇒ the peak gets a loading/geometry explanation.")
    else:
        print(f"  -> G-MULTI  {matches} ⇒ RECORDED; these are mutually related, so this test")
        print("     cannot say which.")
    _write(rows)
    return 0


def _write(rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
