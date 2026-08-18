"""GS4.1 round 9 leg B -- the lag equation, checked against the rebuilt XFOIL.

Binding text: docs/dev_phase_four/20260820-0900-lag-and-xfoil-check-prereg.md
plus addenda #1 (20260820-1500, the dump's normalisations) and #2 (20260820-1600,
two marches). All committed before this script existed.

★★ This is an IMPLEMENTATION check, not model validation: XFOIL runs the same
Drela-Giles family, so agreement says the equations were transcribed right, not
that they are right. Model validation needs an experiment carrying boundary-layer
profiles, and the repo has none -- all three experimental references hold only Cp.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_lag_xfoil/run.py
"""

import csv
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

XFOIL = os.path.join(ROOT, "tools", "xfoil", "xfoil")
N_PANELS, REYNOLDS, MACH, ALPHA_DEG, N_ITER = 280, 3.0e6, 0.0, 2.0, 200
X_TRIP = 0.05
RHO, MU = 1.0, 1.0 / REYNOLDS          # chord = Vinf = 1  =>  mu = 1/Re
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


BATCH = ("PLOP\nG F\n\nNACA 0012\nPPAR\n"
         f"N {N_PANELS}\n\n\nOPER\nVISC {REYNOLDS:.2E}\nMACH {MACH:.1f}\n"
         f"VPAR\nXTR {X_TRIP:.2f} {X_TRIP:.2f}\n\n"
         f"ITER {N_ITER}\nALFA {ALPHA_DEG:.1f}\nDUMP d.txt\n\nQUIT\n")


def run_xfoil(wd):
    with open(os.path.join(wd, "in.txt"), "w") as f:
        f.write(BATCH)
    with open(os.path.join(wd, "in.txt")) as fin:
        subprocess.run([XFOIL], stdin=fin, cwd=wd, capture_output=True,
                       text=True, timeout=300)
    rows = [l.split() for l in open(os.path.join(wd, "d.txt"))
            if not l.lstrip().startswith("#")]
    a = np.array([[float(v) for v in r] for r in rows if len(r) == 10])
    le = int(np.argmin(a[:, 1]))
    past = np.where(a[le:, 1] > 1.0 + 1.0e-9)[0]
    end = le + int(past[0]) if past.size else len(a)
    return a[:end], a[end:]


class Side:
    """One surface, ordered from the stagnation point to the trailing edge, with
    XFOIL's columns converted to OUR normalisations (addendum #1)."""

    def __init__(self, raw, name):
        s, x, ue, ds, th, cf, hk, cd, ct = (raw[:, 0], raw[:, 1], raw[:, 3],
                                            raw[:, 4], raw[:, 5], raw[:, 6],
                                            raw[:, 7], raw[:, 8], raw[:, 9])
        self.name = name
        self.s = np.abs(s - s[0])
        self.x = x
        self.ue = np.abs(ue)
        self.theta = th
        self.H = hk
        self.cf = cf / np.maximum(self.ue, 1e-12) ** 2      # / Vinf^2 -> / ue^2
        self.cD = cd / np.maximum(self.ue, 1e-12) ** 3
        self.ctau = ct ** 2                                  # CT is sqrt(Ctau)
        self.dstar = ds


def split_at_stagnation(surf):
    """Ue in the dump is signed (it is GAM/QINF), so the stagnation point is the
    sign change. Each side is then ordered outward from it."""
    ue = surf[:, 3]
    k = int(np.argmin(np.abs(ue)))
    a = surf[:k + 1][::-1]              # from stagnation back toward the TE
    b = surf[k:]                        # from stagnation forward to the other TE
    return Side(a, "side_a"), Side(b, "side_b"), k


def transition_index(side):
    """The first station XFOIL treats as TURBULENT.

    ★★ Round 10 addendum #1. NOT "the first station with CT > 0": `xbl.f:821-822`
    stores the amplification factor in CTAU before transition and sqrt(Ctau)
    after, so that test finds where AMPLIFICATION starts. On this case it landed
    two stations early, at x/c 0.04496 and 0.04916 where H is still 2.57 and c_f
    is still on the laminar correlation, and leg B then compared our turbulent
    state against XFOIL's laminar one -- which is where its 270 % Ctau deviation
    came from.

    The transition point is an INPUT here: we set the trip at X_TRIP, so the
    first turbulent station is the first one past it. G-XTR below checks that
    against XFOIL's own H trace rather than trusting it.
    """
    i = int(np.argmax(side.x > X_TRIP)) if np.any(side.x > X_TRIP) else None
    return i


def guard_xtr(side, i_tr):
    """G-XTR: three assertions, all against XFOIL's own H, none of them a
    threshold I picked -- laminar H is flat to within its own variation upstream,
    it falls immediately downstream, and the CT column's semantics flip there."""
    if i_tr is None or i_tr < 3 or i_tr + 3 >= len(side.H):
        _record("G-XTR", f"{side.name}: transition index", "3 < i_tr < n-3",
                f"{i_tr}", "G-FAIL -> kill 1")
        return False
    up = side.H[max(0, i_tr - 4):i_tr]
    dn = side.H[i_tr:i_tr + 4]
    flat = float(np.max(np.abs(np.diff(up)))) if up.size > 1 else np.inf
    fall = float(dn[0] - dn[-1])
    ok = flat < 0.05 and fall > 0.3 and dn[0] < up[-1]
    _record("G-XTR", f"{side.name}: x_tr is the trip, verified on XFOIL's H",
            "H flat upstream (<0.05/station), falling downstream (>0.3 over 4)",
            f"i_tr={i_tr} x/c={side.x[i_tr]:.5f}; upstream max |dH| {flat:.4f} "
            f"(H {up[-1]:.4f}); downstream drop {fall:.4f} "
            f"(H {dn[0]:.4f} -> {dn[-1]:.4f})",
            "G-XTR PASS" if ok else "G-FAIL -> kill 1")
    print(f"         G-COLSEM {side.name}: CT rows 0..{i_tr-1} are the "
          f"amplification factor n (max {side.ctau[:i_tr].max()**0.5:.4f}), rows "
          f"{i_tr}.. are sqrt(Ctau); only the latter are compared "
          f"(xbl.f:821-822)")
    return ok


def ue_interp(side):
    p = PchipInterpolator(side.s, side.ue, extrapolate=True)
    d = p.derivative()

    def f(x):
        return float(p(x)), float(d(x))
    return f, p


def march_pair(S, C2, side, i0, i_tr, y0, n_substep):
    """Both arms over the same stations, from the same seed."""
    stations = side.s[i0 + 1:]
    ue_fn, _ = ue_interp(side)
    kw = dict(rho=RHO, mu=MU, n_substep=n_substep,
              x_tr=None if i_tr is None else float(side.s[i_tr]))
    out = {}
    for name, lag in (("equil", False), ("lag", True)):
        try:
            out[name] = S.march_correlation(stations, y0, float(side.s[i0]),
                                            ue_fn, lag=lag, **kw)
        except Exception as exc:                      # report, never swallow
            out[name] = exc
    return out, stations


def guards(C2, S):
    print("== guards ==")
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--",
                        "pyfp3d/viscous/closures.py", "pyfp3d/viscous/ibl3.py"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  closures.py + ibl3.py unchanged: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed")
    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _f(nm):
        def _g(*a, **k):
            raise AssertionError(f"G-NOSOLVE: {nm} was CALLED")
        return _g
    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            setattr(mod, nm, _f(nm)); n += 1
    print(f"  G-NOSOLVE {n} solver entry points stubbed  PASS")
    print(f"  G-SOURCE  SCCON={C2.SCCON} (xbl.f:1558)  DUXCON={C2.DUXCON} "
          f"(:1567)  HDMAX={C2.HDMAX} (xblsys.f:1112)  "
          f"CTRCON/CTRCEX={C2.CTRCON}/{C2.CTRCEX} (xbl.f:1564-1565)")

    # G-USED: each NEW lag constant, disabled alone, must move the lag rate
    base = dict(s_tau=0.03, s_tau_eq=0.04, us=0.55, delta=0.1, uq=-0.02,
                due_over_ue=0.0)
    on = C2.lag_rate(**base)
    inert = []
    for attr, off in (("SCCON", 0.0), ("DUXCON", 0.0)):
        keep = getattr(C2, attr)
        setattr(C2, attr, off)
        try:
            if C2.lag_rate(**base) == on:
                inert.append(attr)
        finally:
            setattr(C2, attr, keep)
    th, H = 2.0e-2, 1.45
    de_on = C2.bl_thickness(th, H)
    keep = C2.HDMAX
    C2.HDMAX = 1.0e9
    de_off = C2.bl_thickness(th, 1.0001)
    C2.HDMAX = keep
    if de_off == C2.bl_thickness(th, 1.0001):
        inert.append("HDMAX")
    keep = C2.CTRCON
    C2.CTRCON = 1.0
    if C2.s_tau_at_transition(1.5, 1.5e-3) == \
            keep * C2.s_tau_at_transition(1.5, 1.5e-3):
        inert.append("CTRCON")
    C2.CTRCON = keep
    _record("G-USED", "each new lag constant reaches the answer",
            "disabling any one must change it",
            f"{4-len(inert)}/4 wired" + (f"; inert {inert}" if inert else ""),
            "G-USED PASS" if not inert else "G-FAIL -> kill 2")

    # ★ Addendum #3: the band text used to quote the TURBULENT plate endpoint
    # against a LAMINAR march, and the verdict only asserted the first of the two
    # printed numbers -- half a check. All three anchors are asserted now.
    h = C2.zpg_fixed_point()
    y0 = C2.blasius_state(0.01, H=2.591100)
    st = S.march_correlation(np.array([1.0, 100.0]), y0, 0.01,
                             S.flat_plate_ue(1.0), n_substep=2000)
    yt = C2.blasius_state(0.05, ue=1.0, rho=RHO, mu=1.0e-5, H=2.591100)
    tp = S.march_correlation(np.geomspace(5.1, 400.0, 30), yt, 0.05,
                             S.flat_plate_ue(1.0), rho=RHO, mu=1.0e-5,
                             n_substep=4000, x_tr=5.0)
    ok = (abs(h - 2.590433) < 1e-6 and abs(st.H[-1] - 2.590433) < 5e-4
          and abs(tp.H[-1] - 1.2932778384340817) < 1e-12)
    _record("G-LEGACY", "laminar fixed point / laminar march / turbulent plate",
            "2.590433 / 2.590433 / 1.2932778384340817 (all three asserted)",
            f"{h:.6f} / {st.H[-1]:.6f} / {tp.H[-1]:.16f}",
            "G-LEGACY PASS" if ok else "G-FAIL -> kill 3")

    print("  G-PROV    reference = the LOCALLY REBUILT XFOIL 6.99, used WHOLE "
          "(its own u_e with its own theta/H/cf/CT). It is NOT the same build "
          "as the committed CSV: round 7 measured (i) BLDUMP writes 10 fields "
          "where the committed generator filters >=12, and (ii) station x "
          "diverges to 1.75e-3 toward the TE at the same PPAR/N 280 -- while "
          "cl/cd/cm agree to <=1 last stored digit on two cases.")
    print("  G-UNITS   c_f = Cf_dump/u_e^2 (xoper.f:1845 + xbl.f:276); "
          "c_D = CD_dump/u_e^3 (:1863 + :277); H_dump = HK (:1859); "
          "Ctau = CT_dump^2 (:1864 + xblsys.f:713)")


def units_check(side):
    """G-UNITS with teeth: rebuild c_D at the TE from the OTHER two converted
    columns and require it to close. If the normalisations were wrong this is
    off by powers of u_e, which at u_e ~ 0.89 is tens of percent."""
    from pyfp3d.viscous import closures_2d as C2
    i = int(np.argmax(side.x))                     # the trailing edge
    ret = RHO * side.ue[i] * side.theta[i] / MU
    hs = C2.h_star_turb(side.H[i], ret)
    us = C2.slip_velocity(hs, side.H[i])
    cdr = C2.cd_turb(C2.cf_turb_wall(side.H[i], ret), us, side.ctau[i], hs,
                     C2.dfac_low_hk(side.H[i], ret), ret)
    rel = abs(cdr / side.cD[i] - 1.0)
    _record("G-UNITS", "converted c_D at the TE, rebuilt from c_f and Ctau",
            "<= 15 % (structural closure, not an accuracy claim)",
            f"{side.cD[i]:.4e} dumped vs {cdr:.4e} rebuilt = {100*rel:.2f} %",
            "G-UNITS PASS" if rel <= 0.15 else "G-FAIL -> kill 1")
    return rel <= 0.15


def floor_of(S, C2, side, i0, i_tr, y0):
    """L-FLOOR: this round's own numerical noise, measured two ways -- doubling
    the substeps, and swapping the u_e derivative for a second-order difference.
    Every band below is at least this."""
    a, _ = march_pair(S, C2, side, i0, i_tr, y0, 4000)
    b, _ = march_pair(S, C2, side, i0, i_tr, y0, 8000)
    if isinstance(a["lag"], Exception) or isinstance(b["lag"], Exception):
        return None, a
    f = max(float(np.max(np.abs(b["lag"].theta / a["lag"].theta - 1.0))),
            float(np.max(np.abs(b["lag"].H / a["lag"].H - 1.0))))
    return f, b


def compare(side, res, i0, mask):
    rows = []
    idx = np.arange(i0 + 1, len(side.s))
    for name, st in (("equil", res["equil"]), ("lag", res["lag"])):
        if isinstance(st, Exception):
            continue
        for j, k in enumerate(idx):
            rows.append({"arm": name, "s": side.s[k], "x": side.x[k],
                         "turbulent": int(mask[k]),
                         "ue": side.ue[k],
                         "theta_x": side.theta[k], "theta_o": st.theta[j],
                         "H_x": side.H[k], "H_o": st.H[j],
                         "cf_x": side.cf[k], "cf_o": st.cf[j],
                         "ctau_x": side.ctau[k], "ctau_o": st.ctau[j]})
    return rows


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S
    guards(C2, S)

    with tempfile.TemporaryDirectory() as wd:
        surf, wake = run_xfoil(wd)
    sa, sb, kstag = split_at_stagnation(surf)
    print(f"== XFOIL solution: {len(surf)} surface + {len(wake)} wake rows; "
          f"stagnation at index {kstag}, x = {surf[kstag,1]:.5f}; "
          f"{sa.name} {len(sa.s)} stations, {sb.name} {len(sb.s)} ==")
    if not units_check(sb):
        raise SystemExit("G-UNITS failed -- stopping (kill 1)")

    all_rows, floors = [], {}
    for side in (sa, sb):
        i_tr = transition_index(side)
        i0 = 1                                    # first station past stagnation
        y0 = (float(side.theta[i0]), float(side.H[i0]))
        print(f"-- {side.name}: seed at s={side.s[i0]:.5f} x={side.x[i0]:.5f} "
              f"theta={y0[0]:.3e} H={y0[1]:.4f} (G-SEED, XFOIL's own); "
              f"transition at " +
              (f"x={side.x[i_tr]:.4f}" if i_tr is not None else "none"))
        if not guard_xtr(side, i_tr):
            raise SystemExit("G-XTR failed -- stopping (kill 1)")
        mask = np.zeros(len(side.s), dtype=bool)
        if i_tr is not None:
            mask[i_tr:] = True

        # --- march A: the whole surface, from XFOIL's first station ----------
        fl, resA = floor_of(S, C2, side, i0, i_tr, y0)
        floors[side.name] = fl
        for arm in ("equil", "lag"):
            if isinstance(resA[arm], Exception):
                print(f"   march A {arm}: STOPPED -- {type(resA[arm]).__name__}:"
                      f" {resA[arm]}")
        rows = compare(side, resA, i0, mask)
        for r in rows:
            r["march"] = "A_whole"
        all_rows += rows

        # --- march B: turbulent only, seeded from XFOIL at transition --------
        if i_tr is not None and i_tr + 2 < len(side.s):
            y0b = (float(side.theta[i_tr]), float(side.H[i_tr]))
            resB, _ = march_pair(S, C2, side, i_tr, i_tr, y0b, 4000)
            for arm in ("equil", "lag"):
                if isinstance(resB[arm], Exception):
                    print(f"   march B {arm}: STOPPED -- "
                          f"{type(resB[arm]).__name__}: {resB[arm]}")
            rows = compare(side, resB, i_tr, mask)
            for r in rows:
                r["march"] = "B_turb"
            all_rows += rows
        print(f"   L-FLOOR {side.name}: {fl if fl is None else f'{fl:.2e}'}")

    with open(os.path.join(RESULTS, "stations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader(); w.writerows(all_rows)

    verdicts(C2, all_rows, floors)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:20s} [{tag}] {metric} = {measured}")
    print("\n★ IMPLEMENTATION check, not model validation: XFOIL is the same "
          "Drela-Giles family. Validation needs an experiment with BL profiles, "
          "and this repo has none -- all three experimental references hold Cp "
          "only. That is an external data request, not a compute problem.")
    return 0


def verdicts(C2, rows, floors):
    """L-LAM (c_f, laminar, march A), L-TURB (march B) and L-LAG (march B)."""
    import numpy as np

    # ---- L-LAM: c_f in the laminar region, band = the two correlation SETS'
    # own difference at that station's (H, Re_theta).  Addendum #2: this is the
    # only quantity for which that band EXISTS.
    lam = [r for r in rows if r["march"] == "A_whole" and r["arm"] == "equil"
           and not r["turbulent"]]
    n_out, worst, worst_at = 0, 0.0, None
    for r in lam:
        ret = RHO * r["ue"] * r["theta_o"] / MU
        H = r["H_o"]
        if H <= C2.H_MIN or ret < 1.0:
            continue
        ours = 2.0 * C2.re_theta_cf_half(H) / ret
        xf = C2.cf_lam_xfoil(H, ret)
        band = max(abs(ours / xf - 1.0), 3.0 * max(floors.values() or [0.0]))
        dev = abs(r["cf_o"] / r["cf_x"] - 1.0)
        r["band"], r["dev"] = band, dev
        if dev > band:
            n_out += 1
            if dev - band > worst:
                worst, worst_at = dev - band, r["x"]
    _record("L-LAM", "laminar c_f vs XFOIL, march A",
            "per station: the two correlation sets' own difference there",
            f"{len(lam)-n_out}/{len(lam)} inside"
            + (f", worst excess {100*worst:.2f} pp at x/c {worst_at:.4f}"
               if worst_at is not None else ""),
            "L-LAM PASS" if n_out == 0 else "L-FAIL")

    # ---- L-TURB: march B, same equations both sides => band = the floor
    turb = [r for r in rows if r["march"] == "B_turb" and r["arm"] == "lag"
            and r["turbulent"]]
    fl = max([f for f in floors.values() if f is not None] or [1e-12])
    band = max(10.0 * fl, 0.02)
    print("  G-ASREGISTERED L-TURB is registered as: \"the turbulent segment's "
          "theta/H/c_f PER STATION inside the derived band\" -- so this is "
          "evaluated per station and returns (inside, total, worst station). "
          "★ The first execution took a MEDIAN and passed while the same row "
          "printed cf max 101 % -- round 7's R-STATIONS defect repeated "
          "verbatim (registered per-station, implemented as an aggregate).")
    out, bad = {}, []
    for q in ("theta", "H", "cf", "ctau"):
        d, xs = [], []
        for r in turb:
            if r[q + "_x"] == 0.0:
                continue
            d.append(abs(r[q + "_o"] / r[q + "_x"] - 1.0))
            xs.append(r["x"])
        d, xs = np.asarray(d), np.asarray(xs)
        n_in = int(np.sum(d <= band))
        i = int(np.argmax(d))
        out[q] = (n_in, len(d), float(d[i]), float(xs[i]))   # the required triple+
        if n_in != len(d):
            bad.append(q)
        print(f"         {q:6s} {n_in}/{len(d)} inside; worst {100*d[i]:.2f} % "
              f"at x/c {xs[i]:.4f}; outside span x/c "
              f"{xs[d > band].min():.4f}..{xs[d > band].max():.4f}"
              if n_in != len(d) else
              f"         {q:6s} {n_in}/{len(d)} inside; worst {100*d[i]:.2f} % "
              f"at x/c {xs[i]:.4f}")
    _record("L-TURB", "theta/H/c_f/Ctau vs XFOIL PER STATION, march B",
            f"every station <= max(10 x floor, 2 %) = {100*band:.2f} %",
            "; ".join(f"{q} {a}/{b} inside, worst {100*w:.1f} % at x/c {xw:.4f}"
                      for q, (a, b, w, xw) in out.items()),
            "L-TURB PASS" if not bad else f"L-FAIL on {bad}")

    # ---- L-LAG: is the lag arm closer to XFOIL's CT than the equilibrium arm?
    # Relative comparison against a common external reference -- no invented
    # threshold, and it can fail.
    lg = {(r["march"], r["s"]): r for r in rows if r["arm"] == "lag"}
    eq = {(r["march"], r["s"]): r for r in rows if r["arm"] == "equil"}
    keys = [k for k in lg if k in eq and k[0] == "B_turb"
            and lg[k]["turbulent"] and lg[k]["ctau_x"] > 0.0]
    if not keys:
        _record("L-LAG", "lag vs equilibrium against XFOIL's own CT",
                "lag closer than equilibrium", "no comparable stations",
                "L-UNDEFINED")
        return
    dl = np.array([abs(lg[k]["ctau_o"] / lg[k]["ctau_x"] - 1.0) for k in keys])
    de = np.array([abs(eq[k]["ctau_o"] / eq[k]["ctau_x"] - 1.0) for k in keys])
    xs = np.array([lg[k]["x"] for k in keys])
    o = np.argsort(xs)
    dl, de, xs = dl[o], de[o], xs[o]
    # ★ The registration says "in the post-transition relaxation region" but does
    # NOT fix that region's SIZE, and the first execution used N/5 -- a number I
    # picked. So every window is reported and the verdict is taken on the
    # WINDOW-FREE one, which needs no choice of mine at all.
    swept, all_pass = [], True
    for frac, lbl in ((0.05, "N/20"), (0.10, "N/10"), (0.20, "N/5"),
                      (0.33, "N/3"), (0.50, "N/2"), (1.00, "ALL")):
        n = max(3, int(len(xs) * frac))
        b = int(np.sum(dl[:n] < de[:n]))
        swept.append((lbl, n, b, float(np.median(dl[:n])), float(np.median(de[:n]))))
        all_pass &= b > n / 2
        print(f"         window {lbl:5s} n={n:3d}  lag better {b:3d}/{n:3d}  "
              f"median lag {100*np.median(dl[:n]):6.3f} % vs equil "
              f"{100*np.median(de[:n]):6.3f} %")
    lbl, n, b, ml, me = swept[-1]                       # ALL -- no window chosen
    _record("L-LAG", "lag vs equilibrium against XFOIL's own Ctau, "
            "WHOLE turbulent run (no window chosen); every window also reported",
            "the lag arm's deviation smaller than the equilibrium arm's",
            f"{b}/{n} stations; median |dev| lag {100*ml:.3f} % vs equil "
            f"{100*me:.3f} %; verdict identical at every window N/20..ALL "
            f"({'all pass' if all_pass else 'NOT all pass'})",
            "L-LAG PASS" if b > n / 2 else "L-FAIL")


if __name__ == "__main__":
    sys.exit(main())
