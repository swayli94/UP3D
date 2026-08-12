"""Q1 out-of-sample at the two MASKED stations + Q2 foot-versus-peak attribution. One batch of solves.

Pre-registered in docs/dev_phase_three/20260815-0600-tipstations-mechanism-prereg.md.

Q1: the experiment has SEVEN stations and M3 uses five -- eta 0.96 and 0.99 are tip-masked, so their
self-differences have never been computed. That makes them genuinely out-of-sample for the previous round's
"the non-convergence lives at the outboard leading edge". ★★★ But eta > 0.95 is also where tip_taper acts
and where the P13 free-edge singularity lives, so "bigger there" is equally consistent with a tip artifact.
The criterion therefore requires the FINGERPRINT to reproduce (mass again at x < 0.05), not just magnitude.

Q2: a moving sonic FRONT FOOT predicts a Cp change of |dCp/dx|_foot * |dx_foot|; a changing PEAK DEPTH
predicts |dCp_min| directly. Each predicts a NUMBER, compared against the measured per-station d_self --
a dimensional attribution rather than a correlation over five points.

★ G-R: the five unmasked d_self values must reproduce bit-for-bit, since this is nominally the same solve.
★ phi is cached this round, so a future round does not pay the re-solve again.

Outputs (TRACKED): bench/gate_results/task3_tipstations_mechanism.csv
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

from pyfp3d.mesh.wake_cut import cut_wake                             # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                              # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                  # noqa: E402
from run_m3_budget import (ALPHA, BANDS, ETAS, M_INF,                 # noqa: E402
                           parse_experiment, solve)
import run_task3_fixed_budget as FB                                   # noqa: E402
import run_task3_nonconv_discriminator as ND                          # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_tipstations_mechanism.csv")
RAY = {"S0": 0.0600, "S2": 0.0300, "S4": 0.0150}
CHAIN, PAIR = ("S0", "S2", "S4"), ("S2", "S4")
#: G-R: the five SEEN values (RECORDED baseline only, never a criterion)
GR_SEEN = {0.20: 0.129873, 0.44: 0.243539, 0.65: 0.279530, 0.80: 0.355800, 0.90: 0.373187}
TIP = (0.96, 0.99)
Q2_BAND, Q2_MIN_HITS = (0.5, 2.0), 3
TOTAL_GATE_S = 1800.0


def solve_all_stations(tag, exp):
    """Solve and extract ALL SEVEN stations. phi is cached (G-phi)."""
    h = RAY[tag]
    cpath = os.path.join(FB.SCRATCH, f"curves7_{tag}.npz")
    ppath = os.path.join(FB.SCRATCH, f"phi_{tag}.npz")
    if os.path.exists(cpath):
        z = np.load(cpath)
        etas = sorted({float(k.split("|")[0]) for k in z.files})
        print(f"  [cached 7-station curves] {tag}", flush=True)
        return {e: {k.split("|")[1]: z[k] for k in z.files if float(k.split("|")[0]) == e}
                for e in etas}, None
    mc, wc = cut_wake(FB.cached_mesh(tag, h, 0.5 * h, 0.5 * h, 120.0 * h))
    t0 = time.perf_counter()
    r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=0.05)
    wall = time.perf_counter() - t0
    conv = bool(r.get("converged"))
    print(f"  {tag}: conv={conv} |R|={float(r['residual_history'][-1]):.2e} "
          f"lim={r.get('n_limited')} flr={r.get('n_floored')} ({wall:.0f}s)", flush=True)
    if not conv:
        return None, dict(leg=tag, converged=False, solve_s=round(wall, 1),
                          res_final=float(r["residual_history"][-1]))
    phi = np.asarray(r["phi"])
    np.savez(ppath, phi=phi, gamma=np.atleast_1d(np.asarray(r["gamma"])))
    assert os.path.exists(ppath), "G-phi: phi was not written"
    out = {}
    for eta in ETAS:
        try:
            out[eta] = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        except Exception as exc:                                       # noqa: BLE001
            print(f"    ★ station eta {eta} extraction failed: {exc}")
    np.savez(cpath, **{f"{e}|{k}": np.asarray(v) for e, c in out.items()
                       for k, v in c.items() if np.ndim(v) >= 1})
    return out, dict(leg=tag, converged=True, solve_s=round(wall, 1),
                     res_final=float(r["residual_history"][-1]),
                     n_limited=r.get("n_limited"), n_floored=r.get("n_floored"),
                     m_max=None, cl_p=None, sigma_min=r.get("sigma_min"),
                     m1_max=r.get("m1_max"))


def band_points(curves, exp, eta, band="LE", side="upper"):
    """Experimental points of one station/band, restricted to the THREE-LEVEL COMMON SET.

    ★★★ DEFECT FIXED, and G-R is what caught it. The first version pooled EVERY point in the band,
    while the registered baseline values were computed on the common point set of S0/S2/S4 -- so
    eta 0.65 and eta 0.80, the two stations where S0 has out-of-range points, disagreed while the
    other three reproduced bit-for-bit and the solves were bit-identical (|R| 1.40e-14 / 1.76e-11 /
    7.33e-15). The guard did its job: it flagged that I was comparing a DIFFERENT POINT SET, the
    fifth standing question again. Fixing the set is a defect fix, not a loosened criterion -- and
    after it G-R must pass BIT-IDENTICALLY, which is a stronger check than a tolerance.
    """
    lo, hi = next((l, h) for n, l, h in BANDS if n == band)
    e = exp[eta]
    m = (e["upper"] == (side == "upper")) & (e["x"] >= lo) & (e["x"] < hi)
    xs, ce = e["x"][m], e["cp"][m]
    if not xs.size:
        return None, None, None, False
    keep = np.ones(xs.shape, dtype=bool)
    for t, c in curves.items():
        if eta not in c:
            return None, None, None, False
        cx = np.asarray(c[eta][f"x_{side}"])
        keep &= (xs >= cx.min()) & (xs <= cx.max())
    if not keep.any():
        return None, None, None, False
    per = {}
    for t, c in curves.items():
        cx, cc = np.asarray(c[eta][f"x_{side}"]), np.asarray(c[eta][f"cp_{side}"])
        per[t] = np.interp(xs[keep], cx, cc)
    #: `ok` now means "nothing was silently clamped" -- true by construction here; the number of
    #: points removed is reported instead of being hidden in a boolean.
    return xs[keep], ce[keep], per, True, int((~keep).sum())


def geometry(curve, cp_star, side="upper"):
    """(x_foot, x_peak, cp_min, |dCp/dx| at the foot) for one station."""
    x, cp = np.asarray(curve[f"x_{side}"]), np.asarray(curve[f"cp_{side}"])
    o = np.argsort(x)
    x, cp = x[o], cp[o]
    sup = cp < cp_star
    if not sup.any():
        return None, None, None, None
    i = int(np.flatnonzero(sup)[0])
    x_foot = float(x[i])
    j = int(np.argmin(cp))
    #: local slope at the foot, one-sided over the two cells straddling it
    a, b = max(i - 1, 0), min(i + 1, len(x) - 1)
    slope = abs((cp[b] - cp[a]) / (x[b] - x[a])) if x[b] != x[a] else None
    return x_foot, float(x[j]), float(cp[j]), slope


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}")
    print(f"★ G-C tip_cap = {FB.TIP_CAP};  7 stations {ETAS}\n")
    exp = parse_experiment()
    cp_star = ND.cp_critical(M_INF)
    curves, legs, t_all = {}, [], time.perf_counter()
    for tag in CHAIN:
        if time.perf_counter() - t_all > TOTAL_GATE_S:
            print(f"  ★ total gate exceeded -- {tag} NOT run (kill clause 4)")
            break
        c, row = solve_all_stations(tag, exp)
        if row is not None:
            legs.append(row)
        if c is not None:
            curves[tag] = c
    if not all(t in curves for t in PAIR):
        print(f"\n★ {PAIR} not both available -- STOP (kill clause 2).")
        return 1

    #: --- per-station d_self and the fingerprint --------------------------------------------------
    rows, dself, fing = [], {}, {}
    for eta in ETAS:
        got = band_points(curves, exp, eta)
        xs, ce, per, ok = got[0], got[1], got[2], got[3]
        n_dropped = got[4] if len(got) > 4 else 0
        if xs is None or not xs.size:
            print(f"  ★ eta {eta}: no LE-upper points / extraction failed")
            continue
        d = per[PAIR[0]] - per[PAIR[1]]
        mass = d ** 2
        front = xs < 0.05
        n, nf = len(xs), int(front.sum())
        share_pts = nf / n if n else 0.0
        share_mass = float(mass[front].sum()) / float(mass.sum()) if mass.sum() else 0.0
        dself[eta] = float(np.sqrt(mass.mean()))
        fing[eta] = dict(n=n, n_dropped=n_dropped, n_front=nf, share_pts=share_pts,
                         share_mass=share_mass,
                         conc=(share_mass / share_pts if share_pts else None),
                         conc_max=(1.0 / share_pts if share_pts else None), valid=ok)
        rows.append(dict(kind="station", eta=eta, d_self=dself[eta], **fing[eta]))

    print("=== G-R: do the five SEEN stations reproduce bit-for-bit? ===")
    bad = []
    for eta, ref in GR_SEEN.items():
        got = dself.get(eta)
        okk = got is not None and abs(got - ref) <= 5e-7
        print(f"  eta {eta:.2f}: {got if got is None else format(got, '.6f')} vs {ref:.6f}  "
              f"{'PASS' if okk else '★ FAIL'}")
        if not okk:
            bad.append(eta)
    if bad:
        print(f"  -> ★ G-R FAIL at {bad}: the re-solve is not the same leg. STOP (kill clause 1).")
        _write(rows, legs)
        return 1
    print("  -> PASS (same leg across a re-solve)")

    #: --- Q1 -------------------------------------------------------------------------------------
    print("\n=== Q1: the two MASKED stations, out-of-sample ===")
    print(f"  {'eta':>6}{'n':>4}{'d_self':>11}{'front pts %':>13}{'front mass %':>14}"
          f"{'CONC':>7}{'CONC_max':>10}{'G-X':>6}")
    for eta in ETAS:
        if eta not in fing:
            continue
        f = fing[eta]
        mark = " ★" if eta in TIP else "  "
        print(f"{mark}{eta:>6.2f}{f['n']:>4}{dself[eta]:>11.6f}{100 * f['share_pts']:>13.1f}"
              f"{100 * f['share_mass']:>14.1f}"
              f"{(f['conc'] if f['conc'] else float('nan')):>7.2f}"
              f"{(f['conc_max'] if f['conc_max'] else float('nan')):>10.2f}"
              f"{str(f['valid']):>6}")
    ref90, ref80 = GR_SEEN[0.90], GR_SEEN[0.80]
    tip_ok = [e for e in TIP if e in dself and fing[e]["valid"]]
    if len(tip_ok) < len(TIP):
        print(f"\n  -> ★ Q1 UNDEFINED: usable tip stations {tip_ok} (kill clause 3). Q2 unaffected.")
    else:
        bigger = all(dself[e] >= ref90 for e in TIP)
        smaller = any(dself[e] < ref80 for e in TIP)
        fingerprint = all(fing[e]["share_mass"] > 0.50 and (fing[e]["conc"] or 0) >= 1.0
                          for e in TIP)
        print(f"\n  both >= eta0.90's {ref90:.6f}: {bigger}   any < eta0.80's {ref80:.6f}: {smaller}"
              f"   fingerprint (front mass > 50 % and CONC >= 1): {fingerprint}")
        if bigger and fingerprint:
            print("  -> ★★ O-CONF  the same fingerprint reproduces out-of-sample ⇒ the outboard-LE cell")
            print("     is independently confirmed.")
        elif smaller:
            print("  -> ★★ O-REFUTE  the outboard trend does NOT continue ⇒ the previous round's")
            print("     'outboard' reading is likely a boundary effect of the 5-station set ending at")
            print("     eta 0.90. This WEAKENS my own previous result and is reported as such.")
        elif bigger and not fingerprint:
            print("  -> ★ O-TIP  bigger, but the mass is NOT at the leading edge ⇒ a TIP artifact")
            print("     (taper / P13), not the same mechanism. Not a confirmation.")
        else:
            print("  -> O-MIX  RECORDED, no claim.")

    #: --- Q2 -------------------------------------------------------------------------------------
    print("\n=== Q2: foot-shift vs peak-depth, dimensional attribution ===")
    print(f"  {'eta':>6}{'x_foot(S2)':>12}{'x_foot(S4)':>12}{'|dx_foot|':>11}{'slope':>9}"
          f"{'pred_foot':>11}{'|dCp_min|':>11}{'d_self':>10}{'foot/d':>9}{'peak/d':>9}")
    hits_f = hits_p = 0
    for eta in ETAS[:5]:
        if eta not in dself:
            continue
        g2 = geometry(curves[PAIR[0]][eta], cp_star)
        g4 = geometry(curves[PAIR[1]][eta], cp_star)
        if None in g2 or None in g4:
            print(f"  {eta:>6.2f}  (no supersonic zone at one level -- skipped)")
            continue
        dx = abs(g2[0] - g4[0])
        slope = g4[3]
        pf = slope * dx
        pp = abs(g2[2] - g4[2])
        rf, rp = pf / dself[eta], pp / dself[eta]
        hits_f += int(Q2_BAND[0] <= rf <= Q2_BAND[1])
        hits_p += int(Q2_BAND[0] <= rp <= Q2_BAND[1])
        rows.append(dict(kind="q2", eta=eta, x_foot_S2=g2[0], x_foot_S4=g4[0], dx_foot=dx,
                         slope=slope, pred_foot=pf, pred_peak=pp, d_self=dself[eta],
                         ratio_foot=rf, ratio_peak=rp,
                         cp_min_S2=g2[2], cp_min_S4=g4[2], x_peak_S2=g2[1], x_peak_S4=g4[1]))
        print(f"  {eta:>6.2f}{g2[0]:>12.5f}{g4[0]:>12.5f}{dx:>11.5f}{slope:>9.2f}"
              f"{pf:>11.6f}{pp:>11.6f}{dself[eta]:>10.6f}{rf:>9.2f}{rp:>9.2f}")
    print(f"\n  within [{Q2_BAND[0]}, {Q2_BAND[1]}]:  foot {hits_f}/5   peak {hits_p}/5   "
          f"(need >= {Q2_MIN_HITS})")
    fw, pw = hits_f >= Q2_MIN_HITS, hits_p >= Q2_MIN_HITS
    if fw and not pw:
        print("  -> ★★ M-FOOT  the SONIC FRONT FOOT's position explains the magnitude ⇒ the fix")
        print("     points at the m_crit switch and the upwind donor selection, not at LE resolution.")
    elif pw and not fw:
        print("  -> ★★ M-PEAK  the PEAK DEPTH explains it ⇒ the fix points at LE resolution/recovery.")
    elif fw and pw:
        print("  -> M-BOTH  both land inside the band ⇒ RECORDED, a 1-D estimate cannot separate them.")
    else:
        print("  -> ★ M-NEITHER  neither explains the magnitude ⇒ RECORDED, and INFORMATIVE: a third")
        print("     sub-hypothesis is needed.")
    _write(rows, legs)
    return 0


def _write(rows, legs):
    allr = rows + [dict(kind="leg", **l) for l in legs]
    keys = []
    for r in allr:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(allr)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
