"""GS4.1 round 9 leg A -- the five missing turbulent terms, and what they move.

Binding text: docs/dev_phase_four/20260820-0900-lag-and-xfoil-check-prereg.md
plus addendum #1, 20260820-1000-round9-addendum1-five-more-terms.md (both
committed before this script existed).

★★ Pure re-baseline. No lag equation, no XFOIL comparison -- those are leg B,
whose pre-registration is committed and unchanged. The user's standing rule from
round 8 is that a defect repair and a new feature do not ride in one round.

★ The repair is justified by the SOURCE, not by the outcome: no criterion here
reads "did E-CF improve". The direction is RECORDED.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_five_terms/run.py
"""

import csv
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

#: the XFOIL source is copyrighted and gitignored; A-WHOLE needs it, and reports
#: UNDEFINED rather than PASS when it is absent. The extracted INDEX (names,
#: line numbers, classification) is committed, so the check stays auditable.
SRC_DEFAULT = ("/tmp/claude-1000/-home-lrz-codes-UP3D/"
               "3cd995ca-fd78-494f-80eb-3e7fd575323a/scratchpad/xfsrc/"
               "Xfoil699src/src")
XSRC = os.environ.get("GS41_XFOIL_SRC", SRC_DEFAULT)

RHO, MU, U_INF = 1.0, 1.0e-5, 1.0
X0, X_TR, X1 = 0.05, 5.0, 400.0
RE_HI, ATTRACT_TOL, SEED_PERTURB = 1.0e4, 0.01, 1.15
SUMMARY = []

#: A-WHOLE's table. Every plain (non-derivative) assignment in BLVAR must be in
#: exactly one class, keyed by (name, source line). An unclassified pair FAILS --
#: that is the whole point: the five defects were terms with no counterpart, and
#: no guard over the constants already written could ever have seen them.
BLVAR_TABLE = {
    ("US2", 825): ("HERE", "slip_velocity"),
    ("US2", 838): ("HERE", "US_CLAMP_VAL  (defect 5)"),
    ("US2", 848): ("N/A-WAKE", "ITYP=3 wake clamp 0.99995"),
    ("GCC", 858): ("HERE", "GCCON=0 on the laminar branch"),
    ("GCC", 863): ("HERE", "GCCON in ctau_eq (round 8)"),
    ("HKC", 859): ("HERE", "ctau_eq hkc"),
    ("HKC", 864): ("HERE", "ctau_eq hkc"),
    ("HKC", 868): ("HERE", "ctau_eq hkc 0.01 floor"),
    ("HKB", 874): ("HERE", "ctau_eq hkb"),
    ("USB", 875): ("HERE", "ctau_eq (1-us)"),
    ("CQ2", 876): ("HERE", "ctau_eq (we return its square)"),
    ("CF2", 904): ("N/A-WAKE", "ITYP=3 cf=0"),
    ("CF2", 919): ("HERE", "cf_turb = max(CFT, CFL)  (defect 1)"),
    ("GRT", 972): ("HERE", "dfac_low_hk  (defect 2)"),
    ("HMIN", 973): ("HERE", "dfac_low_hk  (defect 2)"),
    ("FL", 977): ("HERE", "dfac_low_hk  (defect 2)"),
    ("TFL", 981): ("HERE", "dfac_low_hk  (defect 2)"),
    ("DFAC", 982): ("HERE", "dfac_low_hk  (defect 2)"),
    ("DI2", 958): ("HERE", "cd_turb wall term"),
    ("DI2", 994): ("HERE", "cd_turb wall x DFAC  (defect 2)"),
    ("DI2", 999): ("N/A-WAKE", "ITYP=3 zero wall contribution"),
    ("DI2", 1018): ("HERE", "cd_turb outer turbulent  (defect 3)"),
    ("DI2", 1033): ("HERE", "cd_turb outer laminar stress  (defect 4)"),
    ("DI2", 1080): ("N/A-WAKE", "ITYP=3 laminar wake CD"),
    ("DI2", 1093): ("N/A-WAKE", "ITYP=3 doubled for two wake halves"),
    ("DD", 1013): ("HERE", "cd_turb outer turbulent  (defect 3)"),
    ("DD", 1028): ("HERE", "cd_turb outer laminar stress  (defect 4)"),
    ("DE2", 1103): ("DEFERRED", "BL thickness Delta -- LEG B needs it for the "
                                "lag equation's relaxation length"),
    ("HDMAX", 1112): ("DEFERRED", "the Delta cap -- leg B"),
    ("DE2", 1115): ("DEFERRED", "the Delta cap -- leg B"),
}


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def cf_coles_fernholz(ret):
    return 2.0 * ((1.0 / 0.384) * np.log(ret) + 4.127) ** -2


def cf_power_law(ret):
    return 0.024 * ret ** -0.25


def band_at(ret):
    a, b = cf_coles_fernholz(ret), cf_power_law(ret)
    s = abs(a - b) / (0.5 * (a + b))
    return max(0.03, 2.0 * s), s


# ---------------------------------------------------------------------------
# A-WHOLE
# ---------------------------------------------------------------------------

def a_whole():
    path = os.path.join(XSRC, "xblsys.f")
    if not os.path.exists(path):
        _record("A-WHOLE", "every plain assignment in BLVAR is classified",
                "no unclassified (name, line) pair",
                f"source absent at {path}", "A-UNDEFINED")
        return None
    src = open(path).read().splitlines()
    lo = hi = None
    for i, l in enumerate(src, 1):
        if re.match(r"^      SUBROUTINE BLVAR", l):
            lo = i
        elif lo and hi is None and re.match(r"^      END\b\s*(!.*)?$", l) and i > lo:
            hi = i
    pat = re.compile(r"^\s{6,9}([A-Z][A-Z0-9]*)\s*=")
    rows, unclassified = [], []
    for i in range(lo, hi):
        l = src[i - 1]
        if l[:1] in ("C", "c", "#", "*"):
            continue
        m = pat.match(l)
        if not m:
            continue
        key = (m.group(1), i)
        cls, note = BLVAR_TABLE.get(key, ("UNCLASSIFIED", ""))
        rows.append({"var": key[0], "line": key[1], "class": cls, "maps_to": note})
        if cls == "UNCLASSIFIED":
            unclassified.append(key)
    n = {c: sum(r["class"] == c for r in rows) for c in
         ("HERE", "N/A-WAKE", "DEFERRED", "UNCLASSIFIED")}
    print(f"  BLVAR lines {lo}-{hi}: {len(rows)} plain assignments  "
          f"HERE={n['HERE']} N/A-WAKE={n['N/A-WAKE']} "
          f"DEFERRED={n['DEFERRED']} UNCLASSIFIED={n['UNCLASSIFIED']}")
    for r in rows:
        if r["class"] != "HERE":
            print(f"         {r['class']:12s} {r['var']:6s} :{r['line']:5d}  "
                  f"{r['maps_to']}")
    _record("A-WHOLE", "every plain assignment in BLVAR is classified",
            "no unclassified (name, line) pair",
            f"{n['UNCLASSIFIED']} unclassified"
            + (f" {unclassified}" if unclassified else "")
            + f"; {n['DEFERRED']} deferred to leg B",
            "A-WHOLE PASS" if not unclassified else "A-FAIL")
    return rows


# ---------------------------------------------------------------------------
# A-USED: each of the five, disabled on its own, must move the answer
# ---------------------------------------------------------------------------

def a_used(C2):
    """Three classes, per addendum #2 -- registered BEFORE this ever ran.

    ★★ As first written, A-USED read "disabled alone, it must change the
    answer", which is one-sided: it maps "transcribed faithfully but
    unreachable" and "never wired in at all" onto the same FAIL, and those two
    call for opposite actions. The dry run measured that defect 5 is exactly the
    first case -- max raw Us over the whole turbulent H band is 0.922, under the
    0.95 trigger -- so the criterion had to grow a third class before execution.
    """
    def cd(theta, H):
        """★ Addendum #3: the metric is the RHS the march actually consumes,
        not c_D. Defect 1 enters the MOMENTUM equation only -- c_D's wall term
        uses the raw CFT by design -- so reading c_D reported it as unreachable
        when the dry run had already measured it firing. Both components, so a
        change in either shows."""
        r = C2.rhs_turb(theta, H, U_INF, 0.0, rho=RHO, mu=MU)
        return abs(r[0]) + abs(r[1])

    zpg = (2.0e-2, 1.40)                  # Re_theta 2000, on the ZPG plate
    apg = (6.0e-3, 2.90)                  # Re_theta 600, near separation
    legs = [
        ("1 max(CFT, CFL)", apg, "cf_lam_xfoil", lambda: -1.0),
        ("2 DFAC low-Hk", zpg, "dfac_low_hk", lambda: 1.0),
        ("3 0.995 outer", zpg, "CD_OUT_US", 1.0),
        ("4 lam stress", zpg, "CD_LAMSTRESS", 0.0),
        ("5 Us clamp", apg, "US_CLAMP_TRIG", 1.0e9),
    ]
    cls, not_wired = {}, []
    for name, state, attr, off in legs:
        on_z, on_a = cd(*zpg), cd(*state)
        keep = getattr(C2, attr)
        setattr(C2, attr, (lambda *a, **k: off()) if callable(off) else off)
        try:
            off_z, off_a = cd(*zpg), cd(*state)
        finally:
            setattr(C2, attr, keep)
        rz, ra = abs(on_z / off_z - 1.0), abs(on_a / off_a - 1.0)
        if max(rz, ra) > 0.0:
            cls[name] = "WIRED+REACHABLE"
        else:
            # is it wired at the callee level, at an input packet_turb cannot
            # produce?  That distinction is the whole point of addendum #2.
            probe = C2.slip_velocity(2.4, 1.05) == C2.US_CLAMP_VAL
            cls[name] = "WIRED+UNREACHABLE" if probe else "NOT-WIRED"
            if not probe:
                not_wired.append(name)
        print(f"         defect {name:16s} {cls[name]:18s} |RHS| moves "
              f"{100*ra:7.3f} % at its probe state, {100*rz:6.3f} % on ZPG")

    # the unreachability bound, printed as addendum #2 requires
    hi = max(0.5 * C2.h_star_turb(H, ret)
             * (1.0 - (H - 1.0) / (C2.GBCON * H))
             for H in np.linspace(C2.H_TURB_LO, C2.H_TURB_HI, 60)
             for ret in np.geomspace(200.0, 1.0e8, 40))
    print(f"         ★ sup(raw Us) over H in [{C2.H_TURB_LO}, {C2.H_TURB_HI}] "
          f"= {hi:.6f} < trigger {C2.US_CLAMP_TRIG} -- the Us clamp NEVER fires "
          f"in this module. First place it could: leg B's pressure-gradient "
          f"check, or lowering H_TURB_LO.")
    n = {c: sum(v == c for v in cls.values()) for c in
         ("WIRED+REACHABLE", "WIRED+UNREACHABLE", "NOT-WIRED")}
    _record("A-USED", "each correction is wired in (three classes, #2; "
            "metric = the marched RHS, #3)",
            "no NOT-WIRED",
            f"{n['WIRED+REACHABLE']} reachable, {n['WIRED+UNREACHABLE']} wired "
            f"but unreachable (sup Us {hi:.4f} < {C2.US_CLAMP_TRIG}), "
            f"{n['NOT-WIRED']} not wired",
            "A-USED PASS" if not not_wired else "A-FAIL -> kill 2")
    return not not_wired


# ---------------------------------------------------------------------------
# A-SIZE
# ---------------------------------------------------------------------------

def a_size(C2):
    rows = []
    for ret in (300.0, 600.0, 1000.0, 3000.0, 1.0e4, 3.0e4):
        for H in (1.35, 1.45, 1.60, 2.90):
            th = ret * MU / (RHO * U_INF)
            p = C2.packet_turb(th, H, U_INF, rho=RHO, mu=MU)
            us, hs = p["Us"], p["H_star"]
            wall = 0.5 * p["cf_wall"] * us * p["DFAC"]
            outer = p["Ctau_eq"] * (C2.CD_OUT_US - us)
            lam = C2.CD_LAMSTRESS * (C2.CD_OUT_US - us) ** 2 / ret
            tot = wall + outer + lam
            # the pre-round-9 form, for the size of the whole repair
            old = (0.5 * p["cf"] * us + p["Ctau_eq"] * (1.0 - us))
            rows.append({"re_theta": ret, "H": H, "H_star": hs, "Us": us,
                         "DFAC": p["DFAC"], "cf_wall": p["cf_wall"],
                         "cf_max": p["cf"],
                         "cf_max_fires": int(p["cf"] > p["cf_wall"]),
                         "us_clamped": int(us == C2.US_CLAMP_VAL),
                         "cD_wall": wall, "cD_outer": outer, "cD_lam": lam,
                         "cD": tot, "cD_pre_round9": old,
                         "rel_move": tot / old - 1.0})
    for r in rows:
        if r["H"] in (1.45, 2.90):
            print(f"         Re_th={r['re_theta']:7.0f} H={r['H']:.2f} "
                  f"DFAC={r['DFAC']:.4f} wall/outer/lam = "
                  f"{100*r['cD_wall']/r['cD']:5.1f}/"
                  f"{100*r['cD_outer']/r['cD']:5.1f}/"
                  f"{100*r['cD_lam']/r['cD']:4.1f} %   "
                  f"c_D moved {100*r['rel_move']:+6.2f} %"
                  f"{'  cf_max FIRES' if r['cf_max_fires'] else ''}"
                  f"{'  Us CLAMPED' if r['us_clamped'] else ''}")
    mv = [abs(r["rel_move"]) for r in rows]
    _record("A-SIZE", "how far the five move c_D", "RECORDED",
            f"|move| {100*min(mv):.2f}%..{100*max(mv):.2f}%", "RECORDED")
    return rows


# ---------------------------------------------------------------------------
# A-REBASE: round 6's measurement, same pipeline, post-repair
# ---------------------------------------------------------------------------

def a_rebase(C2, S):
    y0 = C2.blasius_state(X0, ue=U_INF, rho=RHO, mu=MU, H=2.591100)
    lead = S.march_correlation(np.array([X_TR * 1.02]), y0, X0,
                               S.flat_plate_ue(U_INF), rho=RHO, mu=MU,
                               n_substep=8000, x_tr=X_TR)
    x_s, th_s, H_s = lead.x[0], lead.theta[0], lead.H[0]
    stations = np.geomspace(x_s * 1.02, X1, 120)
    kw = dict(rho=RHO, mu=MU, n_substep=8000, x_tr=X_TR)
    A = S.march_correlation(stations, (th_s, H_s), x_s, S.flat_plate_ue(U_INF), **kw)
    B = S.march_correlation(stations, (th_s, H_s * SEED_PERTURB), x_s,
                            S.flat_plate_ue(U_INF), **kw)
    d = np.abs(B.H - A.H) / A.H
    start = None
    for i in np.where(d <= ATTRACT_TOL)[0]:
        if np.all(d[i:] <= ATTRACT_TOL):
            start = i
            break
    if start is None:
        _record("E-ATTRACT", "post-repair", "<= 1 % and staying", "never",
                "E-UNDEF")
        return None, None
    _record("E-ATTRACT", "post-repair", "<= 1 % and staying",
            f"Re_theta {A.re_theta[start]:.0f}, final sep {d[-1]:.2e}",
            "E-ATTRACT PASS")

    m = (np.arange(A.x.size) >= start) & (A.re_theta <= RE_HI)
    phys = A.H.min() > 1.05 and A.H.max() < 4.0
    _record("E-PHYS", "post-repair", "1.05 < H < 4 (RECORDED, addendum #3)",
            f"[{A.H.min():.4f}, {A.H.max():.4f}]"
            f" -> standing gate {'holds' if phys else 'BROKEN'}", "RECORDED")
    Hw = A.H[m]
    h_ok = Hw.min() >= 1.25 and Hw.max() <= 1.50
    h_note = "holds" if h_ok else \
        "BROKEN (round 8 already crossed 1.50 at 1.5213)"
    _record("E-H", "post-repair", "[1.25, 1.50] (RECORDED, addendum #3)",
            f"[{Hw.min():.4f}, {Hw.max():.4f}] -> standing gate {h_note}",
            "RECORDED")

    rows, n_out = [], 0
    for i in np.where(m)[0]:
        ret, cf = A.re_theta[i], A.cf[i]
        cfa, cfb = cf_coles_fernholz(ret), cf_power_law(ret)
        b, s = band_at(ret)
        if s > 0.10:
            continue
        dev = max(abs(cf / cfa - 1.0), abs(cf / cfb - 1.0))
        inside = dev <= b
        n_out += (not inside)
        rows.append({"re_theta": ret, "x": A.x[i], "H": A.H[i], "cf": cf,
                     "cf_coles_fernholz": cfa, "cf_power_law": cfb,
                     "band": b, "dev": dev, "inside": int(inside),
                     "excess_pp": max(0.0, 100.0 * (dev - b))})
    _record("E-CF", "post-repair",
            "inside the derived band at every station (RECORDED, addendum #3)",
            f"{len(rows)-n_out}/{len(rows)} inside, worst excess "
            f"{max(r['excess_pp'] for r in rows):.3f} pp -> standing gate "
            + ("holds" if n_out == 0 else "still FAILS, as since round 5"),
            "RECORDED")
    _record("A-STANDING", "the standing round-6 gates, read after the repair",
            "stated separately from the re-baseline rows, not instead of them",
            f"E-PHYS {'holds' if phys else 'BROKEN'}; "
            f"E-H {'holds' if h_ok else 'BROKEN'}; "
            f"E-CF {'holds' if n_out == 0 else 'FAILS'}",
            "A-STANDING RECORDED")
    out = [r for r in rows if not r["inside"]]
    if out:
        print(f"         {len(out)} station(s) outside, Re_theta "
              f"{out[0]['re_theta']:.0f}..{out[-1]['re_theta']:.0f}:")
        for r in out[:12]:
            print(f"           Re_th={r['re_theta']:8.0f} H={r['H']:.4f} "
                  f"dev={r['dev']:.3%} band={r['band']:.3%} "
                  f"excess {r['excess_pp']:.3f} pp")
    return A, rows


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S

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
    k = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            setattr(mod, nm, _f(nm)); k += 1
    print(f"  G-NOSOLVE {k} solver entry points stubbed  PASS")
    print(f"  A-SOURCE  CD_OUT_US={C2.CD_OUT_US} (xblsys.f:1014,1029)  "
          f"CD_LAMSTRESS={C2.CD_LAMSTRESS} (:1029)  "
          f"US_CLAMP={C2.US_CLAMP_TRIG}->{C2.US_CLAMP_VAL} (:836-838)  "
          f"DFAC_C={C2.DFAC_C} (:970)")
    h = C2.zpg_fixed_point()
    y0 = C2.blasius_state(0.01, H=2.591100)
    st = S.march_correlation(np.array([1.0, 100.0]), y0, 0.01,
                             S.flat_plate_ue(U_INF), n_substep=2000)
    ok = abs(h - 2.590433) < 1e-6 and abs(st.H[-1] - 2.590433) < 5e-4
    _record("A-LEGACY", "laminar branch untouched by the repair",
            "ZPG fixed point 2.590433 and x_tr=None march unmoved",
            f"fixed point {h:.6f}, march endpoint {st.H[-1]:.6f}",
            "A-LEGACY PASS" if ok else "A-FAIL -> kill 3")

    print("== A-WHOLE (the source block, read whole) ==")
    whole = a_whole()
    print("== A-USED ==")
    a_used(C2)
    print("== A-SIZE ==")
    size = a_size(C2)
    print("== A-REBASE (round 6's measurement, post-repair) ==")
    A, rows = a_rebase(C2, S)

    if whole:
        with open(os.path.join(RESULTS, "blvar_index.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(whole[0].keys()))
            w.writeheader(); w.writerows(whole)
    with open(os.path.join(RESULTS, "term_sizes.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(size[0].keys()))
        w.writeheader(); w.writerows(size)
    if rows:
        with open(os.path.join(RESULTS, "window_fiveterms.csv"), "w",
                  newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)

    print(f"\n== summary ==  {time.perf_counter()-t0:.2f} s")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:18s} [{tag}] {metric} = {measured}")
    print("\n★ Leg A is a RE-BASELINE. It does not validate the turbulent "
          "closure: being the same equation as XFOIL is not being right, and a "
          "ZPG plate still cannot test pressure gradient or the lag.")
    print("★ Two of the five (max(CFT,CFL), the Us clamp) do NOT fire on a ZPG "
          "plate at all -- rounds 5/6 cannot re-baseline them. Leg B's "
          "pressure-gradient check is where they first act.")
    fails = [r for r in SUMMARY if "FAIL" in r[4]]
    print(f"  {len(fails)} FAIL row(s)")
    return 0
if __name__ == "__main__":
    sys.exit(main())
