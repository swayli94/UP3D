"""A: what does the peak-depth estimate MISS (zero-solve)?  B: is the eta-0.90 peak taper-made?

Pre-registered in docs/dev_phase_three/20260815-1000-shape-and-taper-prereg.md.

A decomposes the self-difference ITSELF with the identity RMS^2 = bias^2 + var. If f_bias is high and
|bias| matches |dCp_min|, the level-to-level change is a band-wide LEVEL SHIFT that the peak already
measures, and the 24-88 % "gap" is not a new mechanism -- just RMS >= |mean|. If f_bias is low, the curve
is DEFORMING: a third sub-hypothesis. ★ f_bias is compared against 1/n, never against zero (n = 1 makes it
identically 1 -- the vacuity this session already nearly published).

B reuses W2's validated structure: sweep the taper RADIUS and read a self-normalising spanwise profile
ratio. W2 measured cl_p moving 3.13 % while the LE-RMS ratio moved 0.15 %, which is also what proved the
knob reaches the solve. Here the swept quantity is the self-difference profile.

Outputs (TRACKED): bench/gate_results/task3_shape_and_taper.csv
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
from pyfp3d.post.surface import planform_area                         # noqa: E402
from run_m3_budget import (ALPHA, ETAS, M_INF, parse_experiment,       # noqa: E402
                           solve)
import run_task3_fixed_budget as FB                                   # noqa: E402
import run_task3_nonconv_discriminator as ND                          # noqa: E402
import run_task3_tipstations_mechanism as TM                          # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_shape_and_taper.csv")
CHAIN, PAIR = ("S0", "S2", "S4"), ("S2", "S4")
#: G-R(A): the seven SEEN d_self values -- pure re-analysis must reproduce them bit-for-bit
GR7 = {0.20: 0.129873, 0.44: 0.243539, 0.65: 0.279530, 0.80: 0.355800,
       0.90: 0.373187, 0.96: 0.355123, 0.99: 0.268246}
UNMASKED = (0.20, 0.44, 0.65, 0.80, 0.90)
A_LEVEL_MIN, A_SHAPE_MAX, A_HITS = 0.7, 0.3, 3
RATIO_BAND = (0.5, 2.0)
RC_SWEEP = (0.025, 0.05, 0.10)
B_TAPER_MIN, B_GEOM_MAX = 0.10, 0.03
RTOL_IDENTITY = 1e-12
TOTAL_GATE_S = 1800.0


def dself_per_station(curves, exp):
    """Per-station d = Cp_S2 - Cp_S4 on the three-level common point set (TM.band_points)."""
    out = {}
    for eta in ETAS:
        got = TM.band_points(curves, exp, eta)
        if got[0] is None:
            continue
        xs, ce, per = got[0], got[1], got[2]
        d = per[PAIR[0]] - per[PAIR[1]]
        out[eta] = dict(x=xs, cp_exp=ce, d=d, n=len(d))
    return out


def decompose(d):
    bias = float(np.mean(d))
    var = float(np.mean((d - bias) ** 2))
    rms = float(np.sqrt(np.mean(d ** 2)))
    return dict(bias=bias, var=var, sd=float(np.sqrt(var)), rms=rms,
                f_bias=(bias * bias / (rms * rms) if rms > 0 else None),
                identity_rel=(abs(bias * bias + var - rms * rms) / max(rms * rms, 1e-300)))


def solve_pair(rc, exp):
    """B: S2 and S4 at one taper radius. Curves cached per (level, rc)."""
    out, rows = {}, []
    for tag in PAIR:
        h = TM.RAY[tag]
        cpath = os.path.join(FB.SCRATCH, f"curves7_{tag}_rc{rc:g}.npz")
        if os.path.exists(cpath):
            z = np.load(cpath)
            etas = sorted({float(k.split("|")[0]) for k in z.files})
            print(f"    [cached] {tag} rc {rc:g}", flush=True)
            out[tag] = {e: {k.split("|")[1]: z[k] for k in z.files
                            if float(k.split("|")[0]) == e} for e in etas}
            continue
        mc, wc = cut_wake(FB.cached_mesh(tag, h, 0.5 * h, 0.5 * h, 120.0 * h))
        t0 = time.perf_counter()
        r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=rc)
        wall = time.perf_counter() - t0
        conv = bool(r.get("converged"))
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        row = dict(kind="B_leg", leg=tag, taper_rc=rc, converged=conv, solve_s=round(wall, 1),
                   res_final=float(r["residual_history"][-1]),
                   n_limited=r.get("n_limited"), n_floored=r.get("n_floored"))
        if conv:
            row.update({k: v for k, v in FB.flow_row(mc, wc, r, s_ref, exp).items()
                        if k in ("cl_p", "cl_kj", "m_max", "sigma_min")})
            phi = np.asarray(r["phi"])
            cs = {}
            for eta in ETAS:
                try:
                    cs[eta] = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
                except Exception:                                      # noqa: BLE001
                    pass
            np.savez(cpath, **{f"{e}|{k}": np.asarray(v) for e, c in cs.items()
                               for k, v in c.items() if np.ndim(v) >= 1})
            out[tag] = cs
        print(f"    {tag} rc {rc:g}: conv={conv} |R|={float(r['residual_history'][-1]):.2e} "
              f"cl_p={row.get('cl_p')} ({wall:.0f}s)", flush=True)
        rows.append(row)
    return out, rows


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}\n")
    exp = parse_experiment()
    cp_star = ND.cp_critical(M_INF)
    curves = {}
    for tag in CHAIN:
        c, _row = TM.solve_all_stations(tag, exp)
        if c is None:
            print(f"★ {tag} unavailable -- STOP")
            return 1
        curves[tag] = c
    st = dself_per_station(curves, exp)
    rows = []

    #: --- G-R(A) + G-I ---------------------------------------------------------------------------
    print("=== G-R(A): pure re-analysis must reproduce the seven seen d_self values ===")
    bad = []
    for eta, ref in GR7.items():
        got = float(np.sqrt(np.mean(st[eta]["d"] ** 2))) if eta in st else None
        okk = got is not None and abs(got - ref) <= 5e-7
        print(f"  eta {eta:.2f}: {'-' if got is None else format(got, '.6f')} vs {ref:.6f}  "
              f"{'PASS' if okk else '★ FAIL'}")
        if not okk:
            bad.append(eta)
    if bad:
        print(f"  -> ★ G-R FAIL at {bad}. STOP (kill clause 1).")
        return 1
    print("  -> PASS (bit-identical)")

    print("\n=== A: decompose d_self itself (bias vs scatter), with 1/n beside f_bias ===")
    print(f"  {'eta':>6}{'n':>4}{'1/n':>7}{'bias':>11}{'sd':>11}{'rms':>11}{'f_bias':>9}"
          f"{'|dCp_min|':>11}{'|bias|/dCpm':>13}   note")
    worst_id, hits_level, hits_ratio, hits_shape, le_1n = 0.0, 0, 0, 0, 0
    for eta in ETAS:
        if eta not in st:
            continue
        D = decompose(st[eta]["d"])
        worst_id = max(worst_id, D["identity_rel"])
        g2 = TM.geometry(curves[PAIR[0]][eta], cp_star)
        g4 = TM.geometry(curves[PAIR[1]][eta], cp_star)
        dcpm = abs(g2[2] - g4[2]) if (g2[2] is not None and g4[2] is not None) else None
        ratio = abs(D["bias"]) / dcpm if dcpm else None
        inv = 1.0 / st[eta]["n"]
        note = ("★ f_bias <= 1/n (not even a level shift is measurable)" if D["f_bias"] <= inv
                else "LEVEL-dominated" if D["f_bias"] >= A_LEVEL_MIN
                else "SHAPE-dominated" if D["f_bias"] <= A_SHAPE_MAX else "mixed")
        if eta in UNMASKED:
            hits_level += int(D["f_bias"] >= A_LEVEL_MIN)
            hits_shape += int(D["f_bias"] <= A_SHAPE_MAX and D["f_bias"] > inv)
            hits_ratio += int(ratio is not None and RATIO_BAND[0] <= ratio <= RATIO_BAND[1])
            le_1n += int(D["f_bias"] <= inv)
        mark = "  " if eta in UNMASKED else " ★"
        print(f"{mark}{eta:>6.2f}{st[eta]['n']:>4}{inv:>7.3f}{D['bias']:>11.6f}{D['sd']:>11.6f}"
              f"{D['rms']:>11.6f}{D['f_bias']:>9.3f}"
              f"{(dcpm if dcpm else float('nan')):>11.6f}"
              f"{(ratio if ratio else float('nan')):>13.2f}   {note}")
        rows.append(dict(kind="A", eta=eta, n=st[eta]["n"], inv_n=inv, dcp_min=dcpm,
                         ratio_bias_dcpm=ratio, **D))
    print(f"\n  G-I: worst |bias^2+var-RMS^2|/RMS^2 = {worst_id:.2e} (rtol {RTOL_IDENTITY:.0e}) -> "
          f"{'PASS' if worst_id <= RTOL_IDENTITY else '★ FAIL'}")
    if worst_id > RTOL_IDENTITY:
        return 1
    signs = {int(np.sign(decompose(st[e]['d'])['bias'])) for e in UNMASKED if e in st}
    print(f"  bias signs across the 5 unmasked stations: {sorted(signs)}"
          + ("   ★ CONSISTENT (a directional level shift)" if len(signs) == 1
             else "   ★ SIGNS FLIP (not a single shift)"))
    print(f"\n  hits: f_bias >= {A_LEVEL_MIN}: {hits_level}/5   |bias|/|dCp_min| in "
          f"{RATIO_BAND}: {hits_ratio}/5   f_bias <= {A_SHAPE_MAX} (and > 1/n): {hits_shape}/5   "
          f"f_bias <= 1/n: {le_1n}/5")
    if le_1n == 5:
        print("  -> ★★★ KILL CLAUSE 5: every station's f_bias is at or below its own 1/n ⇒ not even a")
        print("     LEVEL SHIFT is measurable; the level-to-level change is PURE SCATTER, which")
        print("     CONTRADICTS 'the peak gets deeper'. Reported as such.")
    elif hits_level >= A_HITS and hits_ratio >= A_HITS:
        print("  -> ★ A-LEVEL  the change is a band-wide LEVEL SHIFT whose size the peak measures ⇒")
        print("     the 24-88 % gap is NOT a new mechanism, only RMS >= |mean|.")
    elif hits_shape >= A_HITS:
        print("  -> ★★ A-SHAPE  the curve DEFORMS ⇒ a third sub-hypothesis, distinct from both")
        print("     'the peak deepens' and 'the front foot moves'.")
    else:
        print("  -> A-MIX  RECORDED, no direction claimed.")

    #: --- B --------------------------------------------------------------------------------------
    print("\n=== B: taper radius dose-response on the d_self spanwise profile (W2's structure) ===")
    t0 = time.perf_counter()
    prof, blegs = {}, []
    for rc in RC_SWEEP:
        if time.perf_counter() - t0 > TOTAL_GATE_S:
            print(f"  ★ gate exceeded -- rc {rc:g} NOT run")
            break
        print(f"  rc = {rc:g}")
        if rc == 0.05:
            prof[rc] = {e: float(np.sqrt(np.mean(st[e]["d"] ** 2))) for e in st}
            print("    [production leg -- reusing the cached S2/S4 curves, declared]")
            continue
        cs, br = solve_pair(rc, exp)
        blegs += br
        if not all(t in cs for t in PAIR):
            print(f"    ★ rc {rc:g}: pair incomplete -- dropped")
            continue
        s2 = dself_per_station({**{t: curves[t] for t in ("S0",)}, **cs}, exp)
        prof[rc] = {e: float(np.sqrt(np.mean(s2[e]["d"] ** 2))) for e in s2}
    rows += [dict(kind="B_profile", taper_rc=rc, eta=e, d_self=v)
             for rc, p in prof.items() for e, v in p.items()]
    rows += [dict(**b) for b in blegs]

    usable = [rc for rc in RC_SWEEP if rc in prof and 0.20 in prof[rc] and 0.90 in prof[rc]]
    print(f"\n  {'rc':>7}{'eta0.20':>10}{'eta0.90':>10}{'R_d':>9}{'argmax':>9}")
    Rd = {}
    for rc in usable:
        Rd[rc] = prof[rc][0.90] / prof[rc][0.20]
        am = max(prof[rc], key=lambda e: prof[rc][e])
        print(f"  {rc:>7.3f}{prof[rc][0.20]:>10.6f}{prof[rc][0.90]:>10.6f}{Rd[rc]:>9.4f}"
              f"{am:>9.2f}")
    if len(usable) < 2:
        print(f"\n  -> ★ B UNDEFINED: usable rc legs {usable} (kill clause 3). A unaffected.")
    else:
        lo, hi = min(usable), max(usable)
        rel = abs(Rd[hi] - Rd[lo]) / Rd[lo]
        mono = (len(usable) < 3 or
                (Rd[usable[0]] < Rd[usable[1]] < Rd[usable[2]]) or
                (Rd[usable[0]] > Rd[usable[1]] > Rd[usable[2]]))
        print(f"\n  R_d over rc {lo:g}->{hi:g}: {Rd[lo]:.4f} -> {Rd[hi]:.4f} = "
              f"{100 * rel:+.1f} % relative   monotone: {mono}")
        if rel >= B_TAPER_MIN and mono:
            print("  -> ★★ B-TAPER  the spanwise peak MOVES with the taper radius ⇒ the eta-0.90 peak")
            print("     is an imprint of a known MODEL BIAS, not a flow structure.")
        elif rel <= B_GEOM_MAX:
            print("  -> ★★ B-GEOM  the profile barely moves ⇒ the peak is the WING'S OWN spanwise")
            print("     structure (taper/sweep/loading), so the next look is geometry, not model bias.")
        else:
            print(f"  -> B-MIX  {100 * rel:.1f} % relative"
                  + ("" if mono else ", NON-MONOTONE") + " ⇒ RECORDED, no claim.")
        cls = [(b["taper_rc"], b["leg"], b.get("cl_p")) for b in blegs if b.get("cl_p")]
        if cls:
            print(f"  ★ cl_p alongside (the taper's model bias moves with r_c -- W2's rule): {cls}")
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
