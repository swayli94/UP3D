"""GS4.1 round 5 -- the turbulent closure at local equilibrium, on a ZPG plate.

Binding text: phases/p4/docs/dev_phase_four/20260819-1300-turbulent-closure-prereg.md
(committed before this script existed). Codes T-CF / T-H / T-EQUIL / T-CROSS and
the guards are quoted from it and are NOT re-specified here.

★★ SCOPE, stated wherever this is read: local equilibrium (Ctau = CtauEQ). The
LAG equation is absent on purpose -- a zero-pressure-gradient plate cannot test
it, since equilibrium IS Ctau = CtauEQ there. This round does NOT claim to
satisfy the gate's "two-equation + lag"; rounds 5 and 6 do that together.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_turbulent_closure/run.py
"""

import csv
import os
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R1 = os.path.join(ROOT, "bench", "studies", "gs41_strip_core")
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

RHO, MU, U_INF = 1.0, 1.0e-5, 1.0
X0, X_TR, X1 = 0.05, 5.0, 400.0
RE_LO, RE_HI = 500.0, 1.0e4          # the band the criteria are read over
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


# ---------------------------------------------------------------------------
# Established ZPG correlations -- two INDEPENDENT ones, so the band can be
# derived from their own mutual disagreement instead of picked.
# ---------------------------------------------------------------------------

def cf_coles_fernholz(ret):
    """Coles-Fernholz logarithmic form, c_f = 2[(1/kappa) ln Re_theta + C]^-2
    with kappa = 0.384, C = 4.127 (Nagib-Chauhan-Monkewitz)."""
    return 2.0 * ((1.0 / 0.384) * np.log(ret) + 4.127) ** -2


def cf_power_law(ret):
    """1/4 power law, c_f = 0.024 Re_theta^(-1/4) (Schlichting).

    addendum #1: this replaces a misapplied Karman-Schoenherr. KS relates the
    LENGTH-AVERAGED C_F to a LENGTH Reynolds number, not the local c_f to
    Re_theta -- the wrong variable pair, which is question 5 applied to the
    reference rather than to the measurement.

    ★ Ludwieg-Tillmann was considered as the second reference and rejected: it
    depends on H, and H is one of the quantities under test, so the reference
    would partly absorb an H error. Both references here depend on Re_theta
    alone (question 7 -- a guard must not measure itself).
    """
    return 0.024 * ret ** -0.25


def band_at(ret):
    """max(3 %, 2 x the two correlations' own disagreement) -- pre-reg 4.1."""
    a, b = cf_coles_fernholz(ret), cf_power_law(ret)
    s = abs(a - b) / (0.5 * (a + b))
    return max(0.03, 2.0 * s), s


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------

def guard_frozen():
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--",
                        "pyfp3d/viscous/closures.py", "pyfp3d/viscous/ibl3.py"],
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  closures.py + ibl3.py unchanged vs HEAD: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed -- kill criterion 3")


def guard_no_solve():
    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _f(name):
        def _g(*a, **k):
            raise AssertionError(f"G-NOSOLVE: {name} was CALLED")
        return _g
    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            if not hasattr(mod, nm):
                raise SystemExit(f"G-NOSOLVE: {nm} missing")
            setattr(mod, nm, _f(nm)); n += 1
    print(f"  G-NOSOLVE {n} solver entry points stubbed  PASS")


def guard_source(C2):
    """G-SOURCE: every constant printed with its citation, so the transcription
    is auditable from the output rather than by reading the module."""
    rows = [("GACON", C2.GACON, "xbl.f:1559"),
            ("GBCON", C2.GBCON, "xbl.f:1560"),
            ("GCCON", C2.GCCON, "xbl.f:1561"),
            ("CTCON", C2.CTCON, "xbl.f:1569 = 0.5/(GACON^2 GBCON)"),
            ("HSMIN", C2.HSMIN, "xblsys.f:2394 DATA"),
            ("DHSINF", C2.DHSINF, "xblsys.f:2394 DATA")]
    for n, v, src in rows:
        print(f"  G-SOURCE  {n:7s} = {v!r:<22} {src}")
    assert abs(C2.CTCON - 0.5/(C2.GACON**2 * C2.GBCON)) < 1e-15
    print("  G-SOURCE  CTCON is derived from GACON/GBCON, not typed  PASS")


def guard_di(C2):
    """G-DI: the identity the memory attempt got wrong (DI = 2 c_D / H*)."""
    worst = 0.0
    for H, ret in ((1.3, 800.0), (1.4, 2000.0), (1.8, 6000.0), (2.5, 1.0e4)):
        p = C2.packet_turb(1.0e-3, H, 1.0, rho=RHO, mu=MU)
        # rebuild theta so re_theta is the one asked for
        p = C2.packet_turb(ret * MU / (RHO * U_INF), H, U_INF, rho=RHO, mu=MU)
        di = C2.dissipation_identity(p["cD"], p["H_star"])
        lhs = 2.0 * p["cD"] / p["H_star"]
        worst = max(worst, abs(di - lhs) / abs(lhs))
    print(f"  G-DI      2 c_D / H* == DI to {worst:.2e}; ★ c_D is NOT DI -- "
          "reading DIT's return as c_D is what drove the memory attempt to "
          "H = 0.60  PASS")
    if worst > 1e-14:
        raise SystemExit("G-DI failed -- kill criterion 3")


def guard_legacy(S, C2):
    """G-LEGACY: x_tr=None must reproduce round 3 exactly."""
    h = C2.zpg_fixed_point()
    y0 = C2.blasius_state(0.01, H=2.591100)
    st = S.march_correlation(np.array([1.0, 100.0]), y0, 0.01,
                             S.flat_plate_ue(U_INF), n_substep=2000)
    ok = abs(h - 2.590433) < 1e-6 and abs(st.H[-1] - 2.590433) < 5e-4
    print(f"  G-LEGACY  laminar ZPG fixed point {h:.6f} (round 3: 2.590433); "
          f"laminar march H[-1] {st.H[-1]:.6f} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("G-LEGACY failed -- adding turbulence moved a laminar "
                         "reading; kill criterion 2")


# ---------------------------------------------------------------------------
# The ZPG plate
# ---------------------------------------------------------------------------

def plate(arm, n_substep=8000):
    from pyfp3d.viscous import strip2d as S
    from pyfp3d.viscous import closures_2d as C2
    y0 = C2.blasius_state(X0, ue=U_INF, rho=RHO, mu=MU, H=2.591100)
    stations = np.geomspace(X_TR * 1.02, X1, 60)
    t0 = time.perf_counter()
    st = S.march_correlation(stations, y0, X0, S.flat_plate_ue(U_INF),
                             rho=RHO, mu=MU, n_substep=n_substep,
                             x_tr=X_TR, arm=arm)
    return st, time.perf_counter() - t0


def gate(arm, binding):
    print(f"== ZPG turbulent plate, H* arm = '{arm}' "
          f"{'(BINDING)' if binding else '(RECORDED comparison)'} ==")
    st, wall = plate(arm)
    print(f"  G-RES     {st.x.size} stations, {st.n_substep} substeps, "
          f"wall {wall:.3f} s; Re_theta {st.re_theta.min():.0f}"
          f"..{st.re_theta.max():.0f}")

    m = (st.re_theta >= RE_LO) & (st.re_theta <= RE_HI)
    if m.sum() < 5:
        _record("T-CF", f"stations in Re_theta [{RE_LO:.0f}, {RE_HI:.0f}]",
                ">= 5", str(int(m.sum())), "T-UNDEF")
        return st, None

    rows, worst_cf, n_excl = [], 0.0, 0
    for i in np.where(m)[0]:
        ret, cf = st.re_theta[i], st.cf[i]
        cfa, cfb = cf_coles_fernholz(ret), cf_power_law(ret)
        b, s = band_at(ret)
        if s > 0.10:                       # pre-reg section 6: band meaningless
            n_excl += 1
            continue
        d = max(abs(cf/cfa - 1.0), abs(cf/cfb - 1.0))
        worst_cf = max(worst_cf, d - b)    # positive = outside the band
        rows.append({"arm": arm, "x": st.x[i], "re_theta": ret, "H": st.H[i],
                     "cf": cf, "cf_coles_fernholz": cfa,
                     "cf_power_law": cfb, "band": b,
                     "corr_spread": s, "worst_dev": d, "inside": int(d <= b)})
    inside = sum(r["inside"] for r in rows)
    print(f"  band is DERIVED: max(3 %, 2 x the two correlations' own spread); "
          f"their spread here is {np.mean([r['corr_spread'] for r in rows]):.3%}"
          f" -> band {np.mean([r['band'] for r in rows]):.3%}; "
          f"{n_excl} station(s) excluded for spread > 10 %")
    for r in rows[::max(1, len(rows)//6)]:
        print(f"         Re_th={r['re_theta']:8.0f} H={r['H']:.4f} "
              f"cf={r['cf']:.6f} (CF {r['cf_coles_fernholz']:.6f}, "
              f"PL {r['cf_power_law']:.6f}) dev={r['worst_dev']:.3%} "
              f"{'in' if r['inside'] else '★OUT'}")

    tag = "T-CF" if binding else f"T-CF[{arm}]"
    _record(tag, "c_f vs BOTH established ZPG correlations (Coles-Fernholz log + 1/4 power law)",
            "inside the derived band at every station",
            f"{inside}/{len(rows)} stations inside; worst excess "
            f"{100*max(0.0, worst_cf):.2f} pp",
            ("T-CF PASS" if inside == len(rows) else "T-FAIL") if binding
            else "RECORDED")

    Hb = st.H[m]
    tag = "T-H" if binding else f"T-H[{arm}]"
    okH = bool(Hb.min() >= 1.25 and Hb.max() <= 1.50)
    _record(tag, "H over the band", "[1.25, 1.50]",
            f"[{Hb.min():.4f}, {Hb.max():.4f}]",
            ("T-H PASS" if okH else "T-FAIL") if binding else "RECORDED")

    # T-EQUIL: physical throughout + drift over the last half
    half = st.H[st.x >= np.sqrt(st.x[0]*st.x[-1])]
    dec = np.log10(st.x[-1]/np.sqrt(st.x[0]*st.x[-1]))
    drift = abs(half[-1]/half[0] - 1.0)/max(dec, 1e-9)
    phys = bool(np.all(st.H > 1.05) and np.all(st.H < 4.0))
    tag = "T-EQUIL" if binding else f"T-EQUIL[{arm}]"
    _record(tag, "self-consistent equilibrium",
            "1.05 < H < 4 throughout AND drift < 2 %/decade",
            f"H in [{st.H.min():.4f}, {st.H.max():.4f}], "
            f"drift {100*drift:.3f} %/decade",
            ("T-EQUIL PASS" if (phys and drift < 0.02) else "T-FAIL")
            if binding else "RECORDED")
    return st, rows


def cross_check(st):
    """T-CROSS -- ONE-SIDED by registration: agreement supports transcription,
    disagreement is inconclusive. Two different models owe each other nothing."""
    print("== T-CROSS: against the profile family (ONE-SIDED) ==")
    old = list(csv.DictReader(open(os.path.join(
        R1, "results", "turbulent_recorded.csv"))))
    lo = min(float(r["rel_diff"]) for r in old)
    hi = max(float(r["rel_diff"]) for r in old)
    Hlo = min(float(r["H"]) for r in old)
    Hhi = max(float(r["H"]) for r in old)
    m = (st.re_theta >= RE_LO) & (st.re_theta <= RE_HI)
    pw = 0.027 * st.re_x[m] ** (-1.0/7.0)
    rel = (st.cf[m] - pw)/pw
    same_sign_band = (rel.min() > -0.30) and (rel.max() < 0.40)
    Hover = (st.H[m].min() < Hhi) and (st.H[m].max() > Hlo)
    _record("T-CROSS", "vs profile family's committed ZPG reading",
            "ONE-SIDED: supports only",
            f"correlation cf vs 1/7 law {100*rel.min():+.1f}%.."
            f"{100*rel.max():+.1f}% (profile family: {100*lo:+.1f}%.."
            f"{100*hi:+.1f}%); H {st.H[m].min():.3f}-{st.H[m].max():.3f} "
            f"(profile family {Hlo:.3f}-{Hhi:.3f}); overlap="
            f"{bool(same_sign_band and Hover)}", "RECORDED")
    print("  ★ agreement here SUPPORTS correct transcription; disagreement "
          "would be inconclusive, not evidence of an error (pre-reg 4.2)")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t_all = time.perf_counter()
    print("== guards ==")
    guard_frozen()
    guard_no_solve()
    from pyfp3d.viscous import closures_2d as C2
    from pyfp3d.viscous import strip2d as S
    guard_source(C2)
    guard_di(C2)
    guard_legacy(S, C2)

    st_new, rows_new = gate("new", binding=True)
    st_old, rows_old = gate("old", binding=False)
    if rows_new:
        cross_check(st_new)

    if rows_new is not None:
        allrows = (rows_new or []) + (rows_old or [])
        with open(os.path.join(RESULTS, "zpg_plate.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(allrows[0].keys()))
            w.writeheader(); w.writerows(allrows)
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)

    total = time.perf_counter() - t_all
    print(f"\n== summary ==  {total:.2f} s (gate: 20 min)")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:14s} [{tag}] {metric} = {measured}")
    print("\n★ This round does NOT satisfy the gate's 'two-equation + lag': the "
          "lag equation is absent, because a ZPG plate cannot test it.")
    fails = [r for r in SUMMARY if "FAIL" in r[4]]
    print(f"  {len(fails)} FAIL row(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
