"""Why does the LE upper surface not converge? A DISCRIMINATOR over four candidates.

Pre-registered in phases/p3/docs/dev_phase_three/20260815-0200-nonconvergence-discriminator-prereg.md.

d_self has so far been read only as an RMS, but it is a POINTWISE quantity, so the distribution of its
squared mass in chord and in span is the fingerprint that separates the candidates:

  (a) LE suction peak      -> mass at the very front (x/c < 0.05), NOT at the sonic line
  (b) sub-cell shock       -> mass AT the per-station sonic crossing
  (d) tip residual         -> mass at outboard stations only
  (c) sigma freeze phase   -> footprint OVERLAPS (b) by construction (sigma acts only where m1 > 1)
                              => not separable by footprint; needs the entropy on/off A/B (arm 2)

★★ Every region is scored by CONCENTRATION = (share of d^2 mass) / (share of points), whose null
hypothesis is 1.0 -- mass in proportion to points. That is the same device as the 1/n comparison in the
bias-versus-scatter round: the criterion compares against a uniform distribution, not against a
hand-picked threshold.

★ Arm 1 is zero-solve (cached curves). Arm 2 solves three levels with entropy_correction=False.

Outputs (TRACKED): bench/gate_results/task3_nonconv_discriminator.csv
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

from run_m3_budget import ALPHA, M_INF, parse_experiment              # noqa: E402
import run_task3_common_pointset as CP                                # noqa: E402
import run_task3_fixed_budget as FB                                   # noqa: E402

CSV = os.path.join(HERE, "gate_results", "task3_nonconv_discriminator.csv")
CHAIN = ("S0", "S2", "S4")
PAIR = ("S2", "S4")                                                   # the binding (growing) pair
#: the ray, needed only by arm 2 (same knobs as the self-convergence round)
RAY = {"S0": 0.0600, "S2": 0.0300, "S4": 0.0150}
GR_REF = 0.242221
CONC_MIN, TIE_FRAC, RTOL_MASS = 2.0, 0.10, 1e-12
SONIC_HALFWIDTH = 0.02
TOTAL_GATE_S = 1800.0


def cp_critical(m_inf, gamma=1.4):
    """Cp* -- the pressure coefficient at which the local flow is sonic."""
    t = (1.0 + 0.5 * (gamma - 1.0) * m_inf ** 2) / (1.0 + 0.5 * (gamma - 1.0))
    return 2.0 / (gamma * m_inf ** 2) * (t ** (gamma / (gamma - 1.0)) - 1.0)


#: ★★ DEFECT FIXED: x_hi defaulted to 0.30, so every station returned ~0.29 -- the WINDOW EDGE,
#: not the crossing. The supersonic zone extends far past 0.30 (measured 0.37-0.65), so the first
#: reading was a clipping artifact that happened to look like a tidy constant. Same family as the
#: project's recorded "a knob's name is not its scope": a window's edge is not a measurement.
def sonic_x(curve, cp_star, side="upper", x_hi=1.0):
    """First x where the computed Cp crosses Cp* from below (more negative than Cp*), forward part."""
    x, cp = np.asarray(curve[f"x_{side}"]), np.asarray(curve[f"cp_{side}"])
    o = np.argsort(x)
    x, cp = x[o], cp[o]
    m = x <= x_hi
    x, cp = x[m], cp[m]
    #: supersonic means Cp < Cp*; find the LAST crossing back to subsonic inside the window
    sup = cp < cp_star
    if not sup.any():
        return None
    idx = np.flatnonzero(sup)
    j = idx[-1]
    if j + 1 >= len(x):
        return float(x[j])
    #: linear crossing between j and j+1
    c0, c1 = cp[j] - cp_star, cp[j + 1] - cp_star
    if c1 == c0:
        return float(x[j])
    return float(x[j] + (x[j + 1] - x[j]) * (0.0 - c0) / (c1 - c0))


def anatomy(curves, exp, tag_a, tag_b, label):
    """Arm 1: where does the d^2 mass live? Returns (rows, regions, notes)."""
    cp_star = cp_critical(M_INF)
    tbl = CP.point_table(curves, exp, "LE", "upper")
    common = [r for r in tbl if all(r[3][t][1] for t in curves)]
    d = np.array([r[3][tag_a][0] - r[3][tag_b][0] for r in common])
    xs = np.array([r[1] for r in common])
    et = np.array([r[0] for r in common])
    mass = d ** 2
    tot = float(mass.sum())

    sx, sup_ext = {}, {}
    for eta in sorted({float(e) for e in et}):
        s = {}
        for t in curves:
            s[t] = sonic_x(curves[t][eta], cp_star, "upper") if eta in curves[t] else None
        sx[eta] = s
        #: ★ the band question is not "where is the shock" but "is the BAND supersonic" -- record the
        #: extent, because the LE band turns out to sit INSIDE the supersonic pocket.
        c = curves[tag_b].get(eta)
        if c is not None:
            x, cp = np.asarray(c["x_upper"]), np.asarray(c["cp_upper"])
            m = cp < cp_star
            sup_ext[eta] = (float(x[m].min()), float(x[m].max())) if m.any() else None
    #: the sonic region uses the FINER level of the pair (the reference state)
    near_sonic = np.zeros(len(common), dtype=bool)
    for i, (eta, x) in enumerate(zip(et, xs)):
        s = sx[float(eta)].get(tag_b)
        if s is not None and abs(x - s) <= SONIC_HALFWIDTH:
            near_sonic[i] = True

    regions = {
        "x<0.05 (LE front)": xs < 0.05,
        "0.05<=x<0.15": (xs >= 0.05) & (xs < 0.15),
        f"sonic +-{SONIC_HALFWIDTH}c": near_sonic,
        "eta<=0.65 (inboard)": et <= 0.65,
        "eta>=0.80 (outboard)": et >= 0.80,
    }
    out = {}
    for name, m in regions.items():
        n = int(m.sum())
        share_pts = n / len(common) if len(common) else 0.0
        share_mass = float(mass[m].sum()) / tot if tot else 0.0
        out[name] = dict(n=n, share_pts=share_pts, share_mass=share_mass,
                         conc=(share_mass / share_pts if share_pts > 0 else None),
                         rms=float(np.sqrt(mass[m].mean())) if n else None)
    return dict(label=label, n_common=len(common), d_rms=float(np.sqrt(mass.mean())),
                total_mass=tot, regions=out, sonic=sx, sup_ext=sup_ext, cp_star=cp_star,
                per_station={float(e): float(np.sqrt(mass[et == e].mean()))
                             for e in sorted({float(x) for x in et})})


def arm2(exp, elapsed):
    """Arm 2: the same three levels with entropy_correction=False -- the ONLY thing that separates
    sigma's footprint from the shock's."""
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.surface import planform_area
    from run_m3_budget import solve as _solve
    out, rows = {}, []
    for tag in CHAIN:
        if elapsed[0] > TOTAL_GATE_S:
            print(f"    ★ total gate exceeded -- {tag} (sigma OFF) NOT run")
            break
        h = RAY[tag]
        path = os.path.join(FB.SCRATCH, f"curves_{tag}_noent.npz")
        mc, wc = cut_wake(FB.cached_mesh(tag, h, 0.5 * h, 0.5 * h, 120.0 * h))
        t0 = time.perf_counter()
        r = _solve(mc, wc, entropy=False, kutta="pressure", taper=True, probe_seed=0, taper_rc=0.05)
        wall = time.perf_counter() - t0
        elapsed[0] += wall
        conv = bool(r.get("converged"))
        print(f"    {tag} sigma OFF: conv={conv} |R|={float(r['residual_history'][-1]):.2e} "
              f"lim={r.get('n_limited')} flr={r.get('n_floored')} ({wall:.0f}s)", flush=True)
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        row = dict(leg=f"{tag}_noent", h_wall=h, converged=conv,
                   res_final=float(r["residual_history"][-1]), solve_s=round(wall, 1),
                   n_limited=r.get("n_limited"), n_floored=r.get("n_floored"))
        if conv:
            row.update(FB.flow_row(mc, wc, r, s_ref, exp))
            phi = np.asarray(r["phi"])
            cs = {}
            from pyfp3d.post.section_cut import section_cp_curve
            from pyfp3d.meshgen.wing3d import B_SEMI
            from run_m3_budget import ETAS, N_UNMASKED
            for eta in ETAS[:N_UNMASKED]:
                try:
                    cs[eta] = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
                except Exception:                                      # noqa: BLE001
                    pass
            np.savez(path, **{f"{e}|{k}": np.asarray(v) for e, c in cs.items()
                              for k, v in c.items() if np.ndim(v) >= 1})
            out[tag] = cs
        rows.append(row)
    return out, rows


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}\n")
    exp = parse_experiment()
    curves = {t: CP.load_curves(t) for t in CHAIN}
    if any(v is None for v in curves.values()):
        print("★ cached curves missing -- STOP")
        return 1

    print("=== arm 1 (ZERO SOLVE): where does the d^2 mass live? ===")
    A = anatomy(curves, exp, PAIR[0], PAIR[1], "sigma ON")
    print(f"  Cp* (M {M_INF}) = {A['cp_star']:.5f}   common points {A['n_common']}   "
          f"d_self({PAIR[0]}->{PAIR[1]}) = {A['d_rms']:.6f}")

    #: --- G-R + G-I --------------------------------------------------------------------------------
    tbl = CP.point_table(curves, exp, "LE", "upper")
    com = [r for r in tbl if all(r[3][t][1] for t in CHAIN)]
    e4 = CP.rms([r[3]["S4"][0] - r[2] for r in tbl if r[3]["S4"][0] is not None])
    print(f"\n  G-R: S4 full-set e(LE upper) {e4:.6f} vs {GR_REF} -> "
          f"{'PASS' if abs(e4 - GR_REF) <= 5e-7 else '★ FAIL'}")
    if abs(e4 - GR_REF) > 5e-7:
        return 1
    m1 = A["regions"]["x<0.05 (LE front)"]["share_mass"]
    m2 = A["regions"]["0.05<=x<0.15"]["share_mass"]
    print(f"  G-I: chordwise partition mass {m1:.12f} + {m2:.12f} = {m1 + m2:.12f} -> "
          f"{'PASS' if abs(m1 + m2 - 1.0) <= RTOL_MASS else '★ FAIL'}")
    if abs(m1 + m2 - 1.0) > RTOL_MASS:
        print("  -> ★ G-I FAIL: partition arithmetic wrong. STOP (kill clause 2).")
        return 1
    print("  ★ the SPANWISE partition has a gap at eta 0.80 -- shares only, NOT exhaustive "
          "(pre-registered)")

    print(f"\n  {'region':24}{'n':>5}{'pts %':>8}{'mass %':>9}{'CONC':>8}{'CONC_max':>10}   reading")
    for name, v in A["regions"].items():
        c = v["conc"]
        note = ("-" if c is None else
                "★ CONCENTRATED" if c >= CONC_MIN else
                "above uniform" if c > 1.1 else
                "≈ uniform" if c > 0.9 else "depleted")
        cmax = (1.0 / v["share_pts"]) if v["share_pts"] > 0 else None
        #: ★★ CONC_max = 1/(point share). A region holding most of the points CANNOT reach a high
        #: concentration -- printed next to every value because the registered 2.0 threshold turned
        #: out to be unreachable BY CONSTRUCTION for the LE-front region (criterion defect #5).
        print(f"  {name:24}{v['n']:>5}{100 * v['share_pts']:>8.1f}{100 * v['share_mass']:>9.1f}"
              f"{(c if c is not None else float('nan')):>8.2f}"
              f"{(cmax if cmax else float('nan')):>10.2f}   {note}"
              + ("  ★ 2.0 UNREACHABLE BY CONSTRUCTION" if cmax and cmax < CONC_MIN else ""))

    print("\n  per-station d_self and the sonic crossing (is it even INSIDE the band?)")
    print(f"  {'eta':>6}{'d_self':>11}{'x_sonic(S2)':>13}{'x_sonic(S4)':>13}{'supersonic ext':>17}   in [0,0.15)?")
    for eta, d in A["per_station"].items():
        s2, s4 = A["sonic"][eta].get("S2"), A["sonic"][eta].get("S4")
        inb = (s4 is not None and s4 < 0.15)
        ext = A["sup_ext"].get(eta)
        print(f"  {eta:>6.2f}{d:>11.6f}"
              f"{(s2 if s2 is not None else float('nan')):>13.5f}"
              f"{(s4 if s4 is not None else float('nan')):>13.5f}"
              f"{(f'[{ext[0]:.3f},{ext[1]:.3f}]' if ext else '-'):>17}   "
              f"{'YES' if inb else '★ NO (footprint test N/A in-band)'}")

    print("\n=== RECORDED cross-tab: are 'LE front' and 'outboard' the SAME cell? ===")
    tblx = CP.point_table(curves, exp, "LE", "upper")
    comx = [r for r in tblx if all(r[3][t][1] for t in CHAIN)]
    dx = np.array([r[3][PAIR[0]][0] - r[3][PAIR[1]][0] for r in comx])
    xx = np.array([r[1] for r in comx]); ee = np.array([r[0] for r in comx])
    mm = dx ** 2
    print(f"  {'cell':26}{'n':>5}{'pts %':>8}{'mass %':>9}{'CONC':>8}{'CONC_max':>10}")
    for lab, sel in (("x<0.05 & eta>=0.80", (xx < 0.05) & (ee >= 0.80)),
                     ("x<0.05 & eta<=0.65", (xx < 0.05) & (ee <= 0.65)),
                     ("x>=0.05 & eta>=0.80", (xx >= 0.05) & (ee >= 0.80)),
                     ("x>=0.05 & eta<=0.65", (xx >= 0.05) & (ee <= 0.65))):
        n = int(sel.sum())
        sp = n / len(comx)
        sm = float(mm[sel].sum()) / float(mm.sum())
        print(f"  {lab:26}{n:>5}{100 * sp:>8.1f}{100 * sm:>9.1f}"
              f"{(sm / sp if sp else float('nan')):>8.2f}{(1 / sp if sp else float('nan')):>10.2f}")
    print("  ★ RECORDED, post-hoc: the registered criterion scored the MARGINALS, so a joint cell")
    print("    exceeding 2.0 is NOT a pre-registered pass -- it needs its own confirmation round.")

    print("\n=== arm 2: entropy OFF on the same three levels (separates sigma from the shock) ===")
    elapsed = [0.0]
    noent, rows2 = arm2(exp, elapsed)
    d_off = None
    if all(t in noent for t in PAIR):
        B = anatomy({t: noent[t] for t in noent}, exp, PAIR[0], PAIR[1], "sigma OFF")
        d_off = B["d_rms"]
        print(f"\n  d_self({PAIR[0]}->{PAIR[1]}): sigma ON {A['d_rms']:.6f}   "
              f"sigma OFF {d_off:.6f}   ratio {d_off / A['d_rms']:.3f}")
    else:
        print(f"\n  ★ arm 2 UNDEFINED: converged sigma-OFF levels = "
              f"{sorted(noent)} (need both of {PAIR}) -- arm 1's verdict is unaffected "
              f"(kill clause 3, pre-registered as separate)")

    #: --- the verdict -----------------------------------------------------------------------------
    print("\n=== C-A / C-B / C-C / C-D / C-MIX ===")
    cand = {"C-A": A["regions"]["x<0.05 (LE front)"]["conc"],
            "C-B": A["regions"][f"sonic +-{SONIC_HALFWIDTH}c"]["conc"],
            "C-D": A["regions"]["eta>=0.80 (outboard)"]["conc"]}
    ranked = sorted(((k, v) for k, v in cand.items() if v is not None),
                    key=lambda kv: -kv[1])
    for k, v in ranked:
        print(f"  {k}: concentration {v:.2f}")
    sigma_wins = d_off is not None and d_off <= 0.5 * A["d_rms"]
    if sigma_wins:
        print(f"  -> ★★ C-C  sigma OFF halves d_self ({d_off:.6f} <= "
              f"{0.5 * A['d_rms']:.6f}) ⇒ the sigma freeze phase is the driver.")
    elif not ranked or ranked[0][1] < CONC_MIN:
        top = ranked[0] if ranked else ("-", float("nan"))
        allunif = all(abs((v or 0) - 1.0) < 0.2 for v in cand.values())
        if allunif:
            print(f"  -> ★★★ KILL CLAUSE 5: every concentration is ~1.0 (top {top[0]} "
                  f"{top[1]:.2f}) ⇒ the d^2 mass is UNIFORM in the band.")
            print("     That is NOT C-MIX: all four candidates predicted concentration, so all four")
            print("     fail together and a FIFTH hypothesis is required.")
        else:
            print(f"  -> C-MIX  no region reaches {CONC_MIN} (top {top[0]} {top[1]:.2f}) ⇒ RECORDED.")
    elif len(ranked) > 1 and (ranked[0][1] - ranked[1][1]) / ranked[0][1] < TIE_FRAC:
        print(f"  -> C-MIX  {ranked[0][0]} and {ranked[1][0]} within {100 * TIE_FRAC:.0f} % "
              f"({ranked[0][1]:.2f} vs {ranked[1][1]:.2f}) ⇒ tie, RECORDED.")
    else:
        k, v = ranked[0]
        print(f"  -> ★★ {k}  concentration {v:.2f} >= {CONC_MIN} and highest")
        if k == "C-B":
            print("     ⇒ this MERGES with the boundary's first unexplained deficit (sub-cell shock")
            print("       positioning) -- two lines become one.")
    _write(A, rows2, d_off)
    return 0


def _write(A, rows2, d_off):
    rows = [dict(kind="region", region=n, **{k: v for k, v in val.items()})
            for n, val in A["regions"].items()]
    rows += [dict(kind="station", region=f"eta{e:g}", rms=d) for e, d in A["per_station"].items()]
    rows += [dict(kind="sonic", region=f"eta{e:g}",
                  **{f"x_sonic_{t}": s.get(t) for t in CHAIN}) for e, s in A["sonic"].items()]
    rows += [dict(kind="summary", region="d_self_sigma_on", rms=A["d_rms"]),
             dict(kind="summary", region="d_self_sigma_off", rms=d_off),
             dict(kind="summary", region="cp_star", rms=A["cp_star"]),
             dict(kind="summary", region="n_common", n=A["n_common"])]
    rows += [dict(kind="arm2", **r) for r in rows2]
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
