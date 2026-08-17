"""GS4.1 round 1 -- verification of the 2-D strip core against analytic
laminar solutions.

Binding text: docs/dev_phase_four/20260818-2100-gs41-strip-core-prereg.md
plus addendum #1 (both committed before this script existed). Criteria
V-BLASIUS / V-ORDER / F-SIMILAR and the guards below are quoted from it and
are NOT re-specified here.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_strip_core/run.py

★ G-ANALYTIC: every reference value is obtained by integrating the Blasius /
Falkner-Skan ODE **in this script**. The textbook numbers (c_f sqrt(Re_x) =
0.664, H = 2.59) are printed as a self-check of that integration; they are
never the criterion, and nothing here is a remembered constant.
"""

import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)

RHO, MU, U_INF = 1.0, 1.0e-5, 1.0          # the V1 study's conventions
SUMMARY = []
WALLS = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def _wall(label, seconds, detail=""):
    WALLS.append((label, seconds, detail))
    print(f"  G-WALL  {label}: {seconds:.4f} s  {detail}")


# ---------------------------------------------------------------------------
# Guards (pre-registration section 3) -- run FIRST, fail fast
# ---------------------------------------------------------------------------

def guard_frozen():
    """G-FROZEN: ibl3.py and closures.py bit-unchanged."""
    paths = ["pyfp3d/viscous/ibl3.py", "pyfp3d/viscous/closures.py"]
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--"] + paths,
                       cwd=ROOT, capture_output=True)
    ok = r.returncode == 0
    print(f"  G-FROZEN  {' '.join(paths)} unchanged vs HEAD: "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("G-FROZEN failed -- kill criterion 3/5 fires")


def guard_reuse():
    """G-REUSE: the new core calls the closure, it does not re-implement it.

    Source assertion: no closure constant appears as a literal in strip2d.py.
    """
    src = open(os.path.join(ROOT, "pyfp3d/viscous/strip2d.py")).read()
    from pyfp3d.viscous import closures as C
    banned = {"KAPPA": C.KAPPA, "B_SPALDING": C.B_SPALDING,
              "A1_BRADSHAW": C.A1_BRADSHAW, "C_L_DEFAULT": C.C_L_DEFAULT,
              "RECOVERY_R": C.RECOVERY_R, "GAMMA_AIR": C.GAMMA_AIR}
    hits = [n for n, v in banned.items() if repr(v) in src]
    calls = sum(src.count(f"C.{n}") for n in
                ("closure_scalar", "stress_source", "blasius_seed"))
    ok = not hits and calls > 0
    print(f"  G-REUSE   closure constant literals in strip2d.py: {hits or 'none'}; "
          f"closure calls: {calls} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("G-REUSE failed -- kill criterion 3 fires")


def guard_no_solve():
    """G-NOSOLVE: every FP solver entry point becomes a raising stub.

    Not a detector -- a solve becomes IMPOSSIBLE (the GS4.0j/k idiom).
    """
    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _forbidden(name):
        def _f(*a, **k):
            raise AssertionError(
                f"G-NOSOLVE: {name} was CALLED -- this round is registered "
                "zero-solve (kill criterion 3)")
        return _f

    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            if not hasattr(mod, nm):
                raise SystemExit(f"G-NOSOLVE: {nm} missing -- library drift, "
                                 "the guard cannot cover what it claims")
            setattr(mod, nm, _forbidden(nm))
            n += 1
    print(f"  G-NOSOLVE {n} solver entry points replaced by raising stubs -- "
          "a solve is IMPOSSIBLE, not merely detected  PASS")


# ---------------------------------------------------------------------------
# G-ANALYTIC: Blasius / Falkner-Skan solved here, from the ODE
# ---------------------------------------------------------------------------

def falkner_skan(m, eta_max=12.0):
    """Solve f''' + f f'' + beta (1 - f'^2) = 0, beta = 2m/(m+1).

    eta = y sqrt((m+1) u_e / (2 nu x)); f(0)=f'(0)=0, f'(inf)=1. Shooting on
    f''(0) by bisection. Returns a dict of the similarity constants.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    beta = 2.0 * m / (m + 1.0)

    def rhs(_, u):
        f, fp, fpp, i1, i2 = u
        return [fp, fpp, -f * fpp - beta * (1.0 - fp * fp),
                1.0 - fp, fp * (1.0 - fp)]

    def shoot(s, full=False):
        sol = solve_ivp(rhs, (0.0, eta_max), [0.0, 0.0, s, 0.0, 0.0],
                        rtol=1e-12, atol=1e-14, dense_output=full,
                        method="DOP853")
        return sol if full else sol.y[1, -1] - 1.0

    s = brentq(shoot, 0.1, 3.0, xtol=1e-13, rtol=1e-15)
    sol = shoot(s, full=True)
    d1_hat, th_hat = sol.y[3, -1], sol.y[4, -1]
    return {
        "m": m, "beta": beta, "fpp0": s,
        "d1_hat": d1_hat, "th_hat": th_hat,
        "H": d1_hat / th_hat,
        # c_f sqrt(Re_x) = 2 f''(0) sqrt((m+1)/2)
        "cf_sqrt_rex": 2.0 * s * np.sqrt((m + 1.0) / 2.0),
        # theta sqrt(Re_x)/x = th_hat sqrt(2/(m+1))
        "theta_coef": th_hat * np.sqrt(2.0 / (m + 1.0)),
    }


def report_analytic(fs):
    print(f"  G-ANALYTIC m={fs['m']:.6f} beta={fs['beta']:.6f}: "
          f"f''(0)={fs['fpp0']:.6f}  H={fs['H']:.6f}  "
          f"cf*sqrt(Re_x)={fs['cf_sqrt_rex']:.6f}  "
          f"theta*sqrt(Re_x)/x={fs['theta_coef']:.6f}")


# ---------------------------------------------------------------------------
# V-BLASIUS + V-ORDER: laminar flat plate
# ---------------------------------------------------------------------------

X_SEED = 0.01                                   # Re_x = 1e3, upstream of the
STATIONS = np.geomspace(0.1, 100.0, 19)         # window Re_x in [1e4, 1e7]:
#                                    3 decades at 6 stations each + endpoint,
#                                    so no decade falls under the >=5-point
#                                    usability floor (pre-registration 4)
N_LADDER = (250, 500, 1000, 2000)               # V-ORDER refinement ladder
N_REF = 32000                                   # self-convergence reference


def _plate_march(n_substep, bl):
    from pyfp3d.viscous import strip2d as S
    theta0 = bl["theta_coef"] * X_SEED / np.sqrt(RHO * U_INF * X_SEED / MU)
    y0 = S.similar_seed(theta0, bl["H"], ue=U_INF, rho=RHO, mu=MU)
    return S.march(STATIONS, y0, X_SEED, S.flat_plate_ue(U_INF),
                   rho=RHO, mu=MU, turbulent=False, n_substep=n_substep)


def gate_v_blasius(bl):
    print("== V-BLASIUS / V-ORDER: laminar flat plate ==")
    st = _plate_march(N_LADDER[-1], bl)
    _wall("V-BLASIUS march", st.wall_time,
          f"({st.n_substep} substeps, {st.x.size} stations)")

    cf_r = st.cf * np.sqrt(st.re_x)
    e_cf = np.abs(cf_r - bl["cf_sqrt_rex"]) / bl["cf_sqrt_rex"]
    e_H = np.abs(st.H - bl["H"]) / bl["H"]

    with open(os.path.join(RESULTS, "v_blasius_stations.csv"), "w") as f:
        f.write("x,re_x,H,H_blasius,rel_err_H,cf_sqrt_rex,"
                "cf_sqrt_rex_blasius,rel_err_cf,A,theta\n")
        for i in range(st.x.size):
            f.write(f"{st.x[i]:.6f},{st.re_x[i]:.6e},{st.H[i]:.8f},"
                    f"{bl['H']:.8f},{e_H[i]:.6e},{cf_r[i]:.8f},"
                    f"{bl['cf_sqrt_rex']:.8f},{e_cf[i]:.6e},{st.A[i]:.8f},"
                    f"{st.theta[i]:.8e}\n")

    decades = np.log10(st.re_x[-1] / st.re_x[0])
    n_per_decade = {}
    for r in st.re_x:
        n_per_decade.setdefault(int(np.floor(np.log10(r))), 0)
        n_per_decade[int(np.floor(np.log10(r)))] += 1
    usable = [d for d, n in n_per_decade.items() if n >= 5]
    print(f"  decades spanned = {decades:.2f}; stations per decade = "
          f"{dict(sorted(n_per_decade.items()))}; usable (>=5 points) = "
          f"{sorted(usable)}")

    if len(usable) < 3:
        _record("V-BLASIUS", "usable Re_x decades", ">= 3",
                f"{len(usable)}", "V-UNDEF")
        return st, bl, None

    ok_cf, ok_H = e_cf.max() <= 0.05, e_H.max() <= 0.05
    _record("V-BLASIUS", "max rel err H vs Blasius", "<= 5 %",
            f"{100*e_H.max():.3f} %", "pass-leg" if ok_H else "fail-leg")
    _record("V-BLASIUS", "max rel err cf*sqrt(Re_x) vs Blasius", "<= 5 %",
            f"{100*e_cf.max():.3f} %", "pass-leg" if ok_cf else "fail-leg")
    _record("V-BLASIUS", "BOTH within +-5 % over >= 3 Re_x decades",
            "both legs pass", f"H {'ok' if ok_H else 'OUT'}, "
            f"cf {'ok' if ok_cf else 'OUT'}",
            "V-BLASIUS PASS" if (ok_cf and ok_H) else "V-FAIL")
    return st, bl, (ok_cf, ok_H)


def gate_v_order(bl):
    print("== V-ORDER: refinement (stations FIXED, only substeps change) ==")
    rows = []
    for n in N_LADDER:
        st = _plate_march(n, bl)
        _wall(f"V-ORDER march n={n}", st.wall_time)
        cf_r = st.cf * np.sqrt(st.re_x)
        rows.append((n, np.abs(st.H - bl["H"]).max() / bl["H"],
                     np.abs(cf_r - bl["cf_sqrt_rex"]).max() / bl["cf_sqrt_rex"],
                     st.H.copy(), cf_r.copy()))
    ref = _plate_march(N_REF, bl)
    _wall(f"V-ORDER reference n={N_REF}", ref.wall_time)
    ref_cf = ref.cf * np.sqrt(ref.re_x)

    with open(os.path.join(RESULTS, "v_order.csv"), "w") as f:
        f.write("n_substep,err_H_vs_blasius,err_cf_vs_blasius,"
                "err_H_self,err_cf_self\n")
        self_rows = []
        for n, eH, ecf, Hs, cfs in rows:
            sH = np.abs(Hs - ref.H).max() / bl["H"]
            scf = np.abs(cfs - ref_cf).max() / bl["cf_sqrt_rex"]
            self_rows.append((n, sH, scf))
            f.write(f"{n},{eH:.8e},{ecf:.8e},{sH:.8e},{scf:.8e}\n")

    def _order(ns, errs):
        ns, errs = np.asarray(ns, float), np.asarray(errs, float)
        good = errs > 0.0
        if good.sum() < 2:
            return float("nan")
        return float(-np.polyfit(np.log(ns[good]), np.log(errs[good]), 1)[0])

    ns = [r[0] for r in rows]
    eH, ecf = [r[1] for r in rows], [r[2] for r in rows]
    monoH = all(eH[i + 1] < eH[i] for i in range(len(eH) - 1))
    monocf = all(ecf[i + 1] < ecf[i] for i in range(len(ecf) - 1))
    pH, pcf = _order(ns, eH), _order(ns, ecf)
    print(f"  err vs Blasius   H: {['%.4e' % e for e in eH]} (order {pH:.3f})")
    print(f"  err vs Blasius  cf: {['%.4e' % e for e in ecf]} (order {pcf:.3f})")

    ok = monoH and monocf and pH >= 0.8 and pcf >= 0.8
    _record("V-ORDER", "err vs Blasius monotone decreasing + order >= 0.8",
            "monotone & >= 0.8",
            f"H mono={monoH} order={pH:.3f}; cf mono={monocf} order={pcf:.3f}",
            "V-ORDER PASS" if ok else "V-FAIL")

    # Pre-registration section 4 row 3 REQUIRES the two explanations be kept
    # apart: this is the discretization error alone, against a fine march.
    sH = [r[1] for r in self_rows]
    scf = [r[2] for r in self_rows]
    qH, qcf = _order(ns, sH), _order(ns, scf)
    print(f"  self-convergence H: {['%.3e' % e for e in sH]} (order {qH:.3f})")
    print(f"  self-convergence cf:{['%.3e' % e for e in scf]} (order {qcf:.3f})")
    _record("V-ORDER", "SEPARATED: self-convergence order (discretization "
            "error alone, vs n=%d)" % N_REF, "RECORDED, not a criterion",
            f"H {qH:.3f}, cf {qcf:.3f}", "RECORDED")
    return {"eH": eH, "ecf": ecf, "sH": sH, "scf": scf, "ns": ns,
            "pH": pH, "pcf": pcf, "qH": qH, "qcf": qcf}


# ---------------------------------------------------------------------------
# F-SIMILAR: Falkner-Skan wedges
# ---------------------------------------------------------------------------

WEDGES = (0.0, 1.0 / 23.0, 1.0)                 # beta = 0, 1/12, 1


def gate_f_similar():
    print("== F-SIMILAR: Falkner-Skan wedges ==")
    from pyfp3d.viscous import strip2d as S

    rows, converged = [], 0
    for m in WEDGES:
        fs = falkner_skan(m)
        report_analytic(fs)
        A_alg, H_alg, cf_alg = S.similarity_fixed_point(m=m, rho=RHO, mu=MU)

        # March from a DISPLACED seed so the reading is the march's own
        # relaxed state, not the seed echoed back.
        ue_fn = S.falkner_skan_ue(m)
        x0, x1 = 0.2, 200.0
        ue0 = ue_fn(x0)[0]
        theta0 = fs["theta_coef"] * x0 / np.sqrt(RHO * ue0 * x0 / MU)
        try:
            y0 = S.similar_seed(theta0, fs["H"], ue=ue0, rho=RHO, mu=MU)
        except ValueError as exc:
            print(f"  m={m:.5f}: seed rejected -- {exc}")
            rows.append((m, fs["H"], float("nan"), H_alg, float("nan"), False))
            continue
        st = S.march(np.geomspace(2.0, x1, 12), y0, x0, ue_fn,
                     rho=RHO, mu=MU, turbulent=False, n_substep=2000)
        _wall(f"F-SIMILAR march m={m:.5f}", st.wall_time)

        H_march = float(st.H[-1])
        settled = abs(st.H[-1] - st.H[-2]) / H_march < 1.0e-4
        # cross-check: the march must land on the closure's own algebraic
        # similar state, else the pressure-gradient terms are wrong
        cross = abs(H_march - H_alg) / H_alg
        good = settled and cross < 1.0e-3
        converged += int(good)
        rel = (H_march - fs["H"]) / fs["H"]
        print(f"  m={m:.5f}: H_march={H_march:.6f} H_algebraic={H_alg:.6f} "
              f"(cross {cross:.2e}) H_FS={fs['H']:.6f} -> {100*rel:+.3f} %"
              f"{'' if good else '  [NOT SETTLED]'}")
        rows.append((m, fs["H"], H_march, H_alg, rel, good))

    with open(os.path.join(RESULTS, "f_similar.csv"), "w") as f:
        f.write("m,beta,H_falkner_skan,H_strip_march,H_algebraic_fixed_point,"
                "rel_err,settled\n")
        for m, hfs, hm, ha, rel, good in rows:
            f.write(f"{m:.8f},{2*m/(m+1):.8f},{hfs:.8f},{hm:.8f},{ha:.8f},"
                    f"{rel:.8e},{int(good)}\n")

    if converged < 2:
        _record("F-SIMILAR", "converged wedges", ">= 2",
                f"{converged}", "F-UNDEF")
        return rows
    worst = max(abs(r[4]) for r in rows if r[5])
    _record("F-SIMILAR", "max |H - H_FS|/H_FS over converged wedges",
            "<= 8 %", f"{100*worst:.3f} % ({converged}/3 wedges)",
            "F-SIMILAR PASS" if worst <= 0.08 else "F-FAIL")
    return rows


# ---------------------------------------------------------------------------
# Turbulent: RECORDED only (pre-registration section 2.3 -- no gate)
# ---------------------------------------------------------------------------

def turbulent_seed(x):
    """Turbulent flat-plate seed (the GV1.1 study's recipe, verbatim).

    The 0.37 / 0.0576 power-law constants are SEEDING correlations, not
    closure constants; they live here rather than in the library so that
    strip2d.py stays free of empirical formulas (G-REUSE).
    """
    from pyfp3d.viscous import closures as C
    re_x = RHO * U_INF * x / MU
    delta = 0.37 * x * re_x ** -0.2
    cf = 0.0576 * re_x ** -0.2
    A = 0.5 * cf * RHO * U_INF * delta / MU
    st = np.array([delta, A, 0.0, 0.0, 1.0e-3, 0.0])
    out, _, _ = C.closure_scalar(st, q=U_INF, rho=RHO, mu=MU, turbulent=True)
    ct = max((C.C_L_DEFAULT * out[C.OUT_SP1] / out[C.OUT_SD]) ** 2, 1.0e-6)
    return np.array([delta, A, ct])


def leg_turbulent():
    print("== turbulent flat plate: RECORDED, no gate (prereg 2.3) ==")
    from pyfp3d.viscous import strip2d as S

    x0 = 0.2
    stations = np.geomspace(0.4, 20.0, 12)
    st = S.march(stations, turbulent_seed(x0), x0, S.flat_plate_ue(U_INF),
                 rho=RHO, mu=MU, turbulent=True, n_substep=2000)
    _wall("turbulent march", st.wall_time)

    cf_pw = 0.027 * st.re_x ** (-1.0 / 7.0)      # Prandtl 1/7 power law
    rel = (st.cf - cf_pw) / cf_pw
    with open(os.path.join(RESULTS, "turbulent_recorded.csv"), "w") as f:
        f.write("x,re_x,re_theta,H,cf,cf_powerlaw_1_7,rel_diff,ctau\n")
        for i in range(st.x.size):
            f.write(f"{st.x[i]:.6f},{st.re_x[i]:.6e},{st.re_theta[i]:.6e},"
                    f"{st.H[i]:.8f},{st.cf[i]:.8e},{cf_pw[i]:.8e},"
                    f"{rel[i]:.6e},{st.ctau[i]:.8e}\n")
    print(f"  H range [{st.H.min():.4f}, {st.H.max():.4f}]; "
          f"cf vs 1/7 power law: {100*rel.min():+.2f} % .. {100*rel.max():+.2f} %")
    _record("turbulent", "cf vs 1/7 power law (Re_x %.1e..%.1e)"
            % (st.re_x[0], st.re_x[-1]), "no gate -- RECORDED",
            f"{100*rel.min():+.2f}%..{100*rel.max():+.2f}%, H in "
            f"[{st.H.min():.3f}, {st.H.max():.3f}]", "RECORDED")
    return st


# ---------------------------------------------------------------------------

def make_figure(plate, order, fs_rows, turb, bl):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    a = ax[0, 0]
    a.semilogx(plate.re_x, plate.H, "o-", label="strip march")
    a.axhline(bl["H"], color="k", ls="--", label=f"Blasius {bl['H']:.4f}")
    a.axhspan(bl["H"] * 0.95, bl["H"] * 1.05, color="g", alpha=0.12,
              label="+-5 % band")
    a.set_xlabel("Re_x"); a.set_ylabel("H"); a.legend(fontsize=8)
    a.set_title("V-BLASIUS: shape factor")

    a = ax[0, 1]
    cf_r = plate.cf * np.sqrt(plate.re_x)
    a.semilogx(plate.re_x, cf_r, "o-", label="strip march")
    a.axhline(bl["cf_sqrt_rex"], color="k", ls="--",
              label=f"Blasius {bl['cf_sqrt_rex']:.4f}")
    a.axhspan(bl["cf_sqrt_rex"] * 0.95, bl["cf_sqrt_rex"] * 1.05, color="g",
              alpha=0.12, label="+-5 % band")
    a.set_xlabel("Re_x"); a.set_ylabel("c_f sqrt(Re_x)"); a.legend(fontsize=8)
    a.set_title("V-BLASIUS: skin friction")

    a = ax[1, 0]
    a.loglog(order["ns"], order["eH"], "o-", label="err H vs Blasius")
    a.loglog(order["ns"], order["ecf"], "s-", label="err cf vs Blasius")
    a.loglog(order["ns"], order["sH"], "o--", label="err H self-conv")
    a.loglog(order["ns"], order["scf"], "s--", label="err cf self-conv")
    a.set_xlabel("substeps"); a.set_ylabel("max rel error")
    a.legend(fontsize=7); a.set_title("V-ORDER: model floor vs discretization")

    a = ax[1, 1]
    ms = [r[0] for r in fs_rows]
    a.plot(ms, [r[1] for r in fs_rows], "ko--", label="Falkner-Skan (exact)")
    a.plot(ms, [r[2] for r in fs_rows], "rs-", label="strip march")
    a.plot(ms, [r[3] for r in fs_rows], "b^:", label="algebraic fixed point")
    a.set_xlabel("wedge exponent m"); a.set_ylabel("H")
    a.legend(fontsize=8); a.set_title("F-SIMILAR: self-similar wedges")

    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "gs41_strip_core.png"), dpi=110)
    plt.close(fig)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.perf_counter()
    print("== guards ==")
    guard_frozen()
    guard_reuse()
    guard_no_solve()

    print("== G-ANALYTIC: references integrated here, not remembered ==")
    bl = falkner_skan(0.0)
    report_analytic(bl)
    print(f"  self-check vs textbook: cf*sqrt(Re_x) {bl['cf_sqrt_rex']:.6f} "
          f"(0.664)  H {bl['H']:.6f} (2.59)  -- printed, NOT the criterion")

    plate, bl, _ = gate_v_blasius(bl)
    order = gate_v_order(bl)
    fs_rows = gate_f_similar()
    turb = leg_turbulent()
    make_figure(plate, order, fs_rows, turb, bl)

    total = time.perf_counter() - t0
    with open(os.path.join(RESULTS, "walls.csv"), "w") as f:
        f.write("label,seconds,detail\n")
        for lab, s, d in WALLS:
            f.write(f"{lab},{s:.6f},{d}\n")
        f.write(f"TOTAL round compute,{total:.6f},"
                f"excludes the close-out lock tier (addendum #1 section 3)\n")
    with open(os.path.join(RESULTS, "summary.csv"), "w") as f:
        f.write("tag,metric,band,measured,verdict\n")
        for row in SUMMARY:
            f.write(",".join('"%s"' % str(c) for c in row) + "\n")

    print(f"\n== summary ==  total round compute {total:.2f} s "
          f"(kill criterion 4 budget: 20 min)")
    fails = [r for r in SUMMARY if "FAIL" in r[4]]
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:18s} [{tag}] {metric} = {measured}")
    print(f"\n  {len(fails)} FAIL row(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
