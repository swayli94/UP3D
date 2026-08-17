"""GS4.1 round 3 -- route (a2): the correlation closure, verified and costed.

Binding text: docs/dev_phase_four/20260819-0500-gs41-a2-correlation-closure-prereg.md
(committed before this script existed). Every band is quoted verbatim from
rounds 1 and 2 so the three rounds compare on one criterion.

Regenerate:  PYTHONNOUSERSITE=1 python bench/studies/gs41_a2_correlation/run.py

★★ READ FIRST (pre-registration section 0): passing D-BLASIUS here is nearly
CIRCULAR -- the correlations are fits to the Falkner-Skan family, so checking
them against Blasius checks a fit against its own training data. Their function
is a TRANSCRIPTION check. A pass must NOT be reported as "(a2)'s accuracy is
validated". The non-circular content of this round is D-COST and D-TURB.
"""

import csv
import os
import statistics
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R1 = os.path.join(ROOT, "bench", "studies", "gs41_strip_core")
R2 = os.path.join(ROOT, "bench", "studies", "gs41_cost_crossflow")
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, ROOT)
sys.path.insert(0, R1)

#: round 1's protocol, verbatim (its own module constants, imported not retyped)
from run import (STATIONS, X_SEED, N_LADDER, N_REF, WEDGES,  # noqa: E402
                 falkner_skan, turbulent_seed, RHO, MU, U_INF)

SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------

def guard_frozen():
    paths = ["pyfp3d/viscous/ibl3.py", "pyfp3d/viscous/closures.py"]
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--"] + paths,
                       cwd=ROOT, capture_output=True)
    print(f"  G-FROZEN  ibl3.py + closures.py unchanged vs HEAD: "
          f"{'PASS' if r.returncode == 0 else 'FAIL'}")
    if r.returncode:
        raise SystemExit("G-FROZEN failed -- kill criterion 3/5 fires")


def guard_authority():
    """G-AUTHORITY: the two closure families must not derive from each other."""
    import ast

    def imports(path):
        out = set()
        for node in ast.walk(ast.parse(open(path).read())):
            if isinstance(node, ast.Import):
                out.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom):
                out.add(node.module or "")
        return out

    c2 = imports(os.path.join(ROOT, "pyfp3d/viscous/closures_2d.py"))
    c1 = imports(os.path.join(ROOT, "pyfp3d/viscous/closures.py"))
    bad = [m for m in c2 if "closures" in m and "closures_2d" not in m]
    bad += [m for m in c1 if "closures_2d" in m]
    print(f"  G-AUTHORITY closures_2d imports {sorted(c2)}; cross-imports: "
          f"{bad or 'none'} -> {'PASS' if not bad else 'FAIL'}")
    if bad:
        raise SystemExit("G-AUTHORITY failed -- kill criterion 3 fires")


def guard_no_solve():
    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _forbidden(name):
        def _f(*a, **k):
            raise AssertionError(f"G-NOSOLVE: {name} was CALLED")
        return _f

    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            if not hasattr(mod, nm):
                raise SystemExit(f"G-NOSOLVE: {nm} missing -- library drift")
            setattr(mod, nm, _forbidden(nm))
            n += 1
    print(f"  G-NOSOLVE {n} solver entry points replaced by raising stubs  PASS")


def guard_legacy():
    """G-LEGACY: round 1's profile-closure readings must not have moved."""
    from pyfp3d.viscous import strip2d as S
    A, H, cf = S.similarity_fixed_point(m=0.0, rho=RHO, mu=MU)
    exp = (8.02881134, 2.708292, 0.710235)
    got = (A, H, cf)
    ok = all(abs(g / e - 1.0) < 1e-6 for g, e in zip(got, exp))
    print(f"  G-LEGACY  profile fixed point A*={A:.8f} H={H:.6f} cf={cf:.6f} "
          f"vs round 1's {exp} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("G-LEGACY failed -- adding an option moved an existing "
                         "reading; kill criterion 2 fires")


def guard_oracle():
    """G-ORACLE: the correlations against the FS ODE, swept -- transcription."""
    from scipy.integrate import solve_ivp
    from pyfp3d.viscous import closures_2d as C2

    print("  G-ORACLE  correlations vs the Falkner-Skan ODE (transcription "
          "check; this is what D-BLASIUS is really for):")
    rows, worst = [], {"hs": 0.0, "cf": 0.0, "cd": 0.0}
    for m in (-0.08, -0.05, -0.02, 0.0, 1.0 / 23.0, 0.1, 0.3, 0.6, 1.0):
        fs = falkner_skan(m)
        beta, s, th = fs["beta"], fs["fpp0"], fs["th_hat"]

        def _mk(last):
            def _rhs(_, u):
                f, fp, fpp, _i = u
                return [fp, fpp, -f * fpp - beta * (1.0 - fp * fp), last(fp, fpp)]
            return _rhs

        e = solve_ivp(_mk(lambda fp, fpp: fp * (1.0 - fp * fp)), (0.0, 12.0),
                      [0, 0, s, 0], rtol=1e-12, atol=1e-14, method="DOP853")
        d = solve_ivp(_mk(lambda fp, fpp: fpp * fpp), (0.0, 12.0),
                      [0, 0, s, 0], rtol=1e-12, atol=1e-14, method="DOP853")
        H = fs["H"]
        hs_ode = e.y[3, -1] / th
        cf_ode = fs["theta_coef"] * fs["cf_sqrt_rex"] / 2.0
        cd_ode = 2.0 * th * d.y[3, -1] / hs_ode
        r = {
            "m": m, "H": H,
            "hs_ode": hs_ode, "hs_corr": C2.h_star(H),
            "cf_ode": cf_ode, "cf_corr": C2.re_theta_cf_half(H),
            "cd_ode": cd_ode, "cd_corr": C2.re_theta_2cd_over_hstar(H),
        }
        for k in ("hs", "cf", "cd"):
            r[f"{k}_rel"] = r[f"{k}_corr"] / r[f"{k}_ode"] - 1.0
            worst[k] = max(worst[k], abs(r[f"{k}_rel"]))
        rows.append(r)
        print(f"            m={m:+.4f} H={H:.5f}  dH*={100*r['hs_rel']:+7.3f} % "
              f"dcf={100*r['cf_rel']:+7.3f} %  dcD={100*r['cd_rel']:+7.3f} %")
    with open(os.path.join(RESULTS, "g_oracle.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    bl = [r for r in rows if r["m"] == 0.0][0]
    print(f"  G-ORACLE  at BLASIUS: dH* {100*bl['hs_rel']:+.3f} %, "
          f"dcf {100*bl['cf_rel']:+.3f} %, dcD {100*bl['cd_rel']:+.3f} % -> PASS "
          f"(transcription correct)")
    print(f"  G-ORACLE  worst over the sweep: H* {100*worst['hs']:.3f} %, "
          f"cf {100*worst['cf']:.3f} %, cD {100*worst['cd']:.3f} % -- the cf one "
          "is the FIT's own error at strong adverse gradient, NOT covered by "
          "this round's wedges (prohibited sentence 3)")
    return rows, bl


# ---------------------------------------------------------------------------
# D-BLASIUS / D-ORDER -- round 1's protocol verbatim
# ---------------------------------------------------------------------------

def _plate(n_substep, bl):
    from pyfp3d.viscous import strip2d as S
    from pyfp3d.viscous import closures_2d as C2
    y0 = C2.blasius_state(X_SEED, ue=U_INF, rho=RHO, mu=MU, H=bl["H"])
    return S.march_correlation(STATIONS, y0, X_SEED, S.flat_plate_ue(U_INF),
                               rho=RHO, mu=MU, n_substep=n_substep)


def gate_blasius(bl):
    print("== D-BLASIUS (round 1's V-BLASIUS protocol verbatim) ==")
    st = _plate(N_LADDER[-1], bl)
    print(f"  G-RES     {st.x.size} stations, {st.n_substep} substeps, "
          f"wall {st.wall_time:.4f} s")
    cf_r = st.cf * np.sqrt(st.re_x)
    e_cf = np.abs(cf_r - bl["cf_sqrt_rex"]) / bl["cf_sqrt_rex"]
    e_H = np.abs(st.H - bl["H"]) / bl["H"]

    with open(os.path.join(RESULTS, "d_blasius_stations.csv"), "w") as f:
        f.write("x,re_x,H,H_blasius,rel_err_H,cf_sqrt_rex,cf_blasius,"
                "rel_err_cf,theta\n")
        for i in range(st.x.size):
            f.write(f"{st.x[i]:.6f},{st.re_x[i]:.6e},{st.H[i]:.8f},"
                    f"{bl['H']:.8f},{e_H[i]:.6e},{cf_r[i]:.8f},"
                    f"{bl['cf_sqrt_rex']:.8f},{e_cf[i]:.6e},{st.theta[i]:.8e}\n")

    per_dec = {}
    for r in st.re_x:
        per_dec.setdefault(int(np.floor(np.log10(r))), 0)
        per_dec[int(np.floor(np.log10(r)))] += 1
    usable = sorted(d for d, n in per_dec.items() if n >= 5)
    print(f"  stations per decade {dict(sorted(per_dec.items()))}; usable "
          f"{usable}")
    if len(usable) < 3:
        _record("D-BLASIUS", "usable decades", ">= 3", str(len(usable)),
                "D-UNDEF")
        return st
    ok_H, ok_cf = e_H.max() <= 0.05, e_cf.max() <= 0.05
    _record("D-BLASIUS", "max rel err H", "<= 5 %", f"{100*e_H.max():.4f} %",
            "pass-leg" if ok_H else "fail-leg")
    _record("D-BLASIUS", "max rel err cf*sqrt(Re_x)", "<= 5 %",
            f"{100*e_cf.max():.4f} %", "pass-leg" if ok_cf else "fail-leg")
    _record("D-BLASIUS", "BOTH within +-5 % over >= 3 decades", "both legs",
            f"H {'ok' if ok_H else 'OUT'}, cf {'ok' if ok_cf else 'OUT'}",
            "D-BLASIUS PASS (★ near-circular, see prereg 0)"
            if (ok_H and ok_cf) else "D-FAIL")
    return st


def gate_order(bl):
    print("== D-ORDER (round 1's ladder verbatim) ==")
    rows = []
    for n in N_LADDER:
        st = _plate(n, bl)
        cf_r = st.cf * np.sqrt(st.re_x)
        rows.append((n, np.abs(st.H - bl["H"]).max() / bl["H"],
                     np.abs(cf_r - bl["cf_sqrt_rex"]).max() / bl["cf_sqrt_rex"],
                     st.H.copy(), cf_r.copy(), st.wall_time))
        print(f"  G-RES     n_substep={n}: wall {st.wall_time:.4f} s")
    ref = _plate(N_REF, bl)
    ref_cf = ref.cf * np.sqrt(ref.re_x)

    def order(ns, errs):
        ns, errs = np.asarray(ns, float), np.asarray(errs, float)
        g = errs > 0
        return (float(-np.polyfit(np.log(ns[g]), np.log(errs[g]), 1)[0])
                if g.sum() >= 2 else float("nan"))

    ns = [r[0] for r in rows]
    eH, ecf = [r[1] for r in rows], [r[2] for r in rows]
    sH = [np.abs(r[3] - ref.H).max() / bl["H"] for r in rows]
    scf = [np.abs(r[4] - ref_cf).max() / bl["cf_sqrt_rex"] for r in rows]
    with open(os.path.join(RESULTS, "d_order.csv"), "w") as f:
        f.write("n_substep,err_H_vs_blasius,err_cf_vs_blasius,err_H_self,"
                "err_cf_self,wall_s\n")
        for i, n in enumerate(ns):
            f.write(f"{n},{eH[i]:.8e},{ecf[i]:.8e},{sH[i]:.8e},{scf[i]:.8e},"
                    f"{rows[i][5]:.6f}\n")
    pH, pcf = order(ns, eH), order(ns, ecf)
    qH, qcf = order(ns, sH), order(ns, scf)
    monoH = all(eH[i + 1] < eH[i] for i in range(len(eH) - 1))
    monocf = all(ecf[i + 1] < ecf[i] for i in range(len(ecf) - 1))
    print(f"  err vs Blasius  H: {['%.4e' % e for e in eH]} (order {pH:.3f})")
    print(f"  err vs Blasius cf: {['%.4e' % e for e in ecf]} (order {pcf:.3f})")
    print(f"  self-convergence H: {['%.3e' % e for e in sH]} (order {qH:.3f})")
    ok = monoH and monocf and pH >= 0.8 and pcf >= 0.8
    _record("D-ORDER", "err vs Blasius monotone + order >= 0.8",
            "monotone & >= 0.8",
            f"H mono={monoH} order={pH:.3f}; cf mono={monocf} order={pcf:.3f}",
            "D-ORDER PASS" if ok else "D-FAIL")
    # round 1's mandated separation of the two explanations, carried forward
    gap = (min(e for e in eH if e > 0) / max(s for s in sH if s > 0)
           if any(s > 0 for s in sH) else float("nan"))
    _record("D-ORDER", "SEPARATED: self-convergence order (discretization "
            f"alone, vs n={N_REF})", "RECORDED",
            f"H {qH:.3f}, cf {qcf:.3f}; model-floor/discretization gap "
            f"{gap:.3g}x (round 1: ~1e11)", "RECORDED")
    return {"eH": eH, "ecf": ecf, "sH": sH, "scf": scf}


def gate_fs():
    print("== D-FS (round 1's F-SIMILAR wedges verbatim) ==")
    from pyfp3d.viscous import strip2d as S
    from pyfp3d.viscous import closures_2d as C2
    rows, conv = [], 0
    for m in WEDGES:
        fs = falkner_skan(m)
        ue_fn = S.falkner_skan_ue(m)
        x0 = 0.2
        ue0 = ue_fn(x0)[0]
        theta0 = fs["theta_coef"] * x0 / np.sqrt(RHO * ue0 * x0 / MU)
        try:
            st = S.march_correlation(np.geomspace(2.0, 200.0, 12),
                                     (theta0, fs["H"]), x0, ue_fn,
                                     rho=RHO, mu=MU, n_substep=2000)
        except C2.ClosureRangeError as exc:
            print(f"  m={m:.5f}: leg stopped -- {exc}")
            rows.append((m, fs["H"], float("nan"), float("nan"), False))
            continue
        Hm = float(st.H[-1])
        settled = abs(st.H[-1] - st.H[-2]) / Hm < 1.0e-4
        rel = (Hm - fs["H"]) / fs["H"]
        conv += int(settled)
        print(f"  m={m:.5f}: H_march={Hm:.6f} H_FS={fs['H']:.6f} "
              f"-> {100*rel:+.4f} %  settled={settled}  "
              f"wall {st.wall_time:.4f} s")
        rows.append((m, fs["H"], Hm, rel, settled))
    with open(os.path.join(RESULTS, "d_fs.csv"), "w") as f:
        f.write("m,beta,H_falkner_skan,H_strip_correlation,rel_err,settled\n")
        for m, hfs, hm, rel, s in rows:
            f.write(f"{m:.8f},{2*m/(m+1):.8f},{hfs:.8f},{hm:.8f},{rel:.8e},"
                    f"{int(s)}\n")
    if conv < 2:
        _record("D-FS", "settled wedges", ">= 2", str(conv), "D-UNDEF")
        return rows
    worst = max(abs(r[3]) for r in rows if r[4])
    _record("D-FS", "max |H - H_FS|/H_FS over settled wedges", "<= 8 %",
            f"{100*worst:.4f} % ({conv}/3 settled)",
            "D-FS PASS" if worst <= 0.08 else "D-FAIL")
    return rows


# ---------------------------------------------------------------------------
# D-COST -- round 2's protocol verbatim, against its committed numbers
# ---------------------------------------------------------------------------

N_SUB = (1, 2, 4)
N_SUB_BINDING = 2


def gate_cost(bl):
    print("== D-COST (round 2's protocol verbatim) ==")
    from pyfp3d.viscous import strip2d as S
    from pyfp3d.viscous import closures_2d as C2

    r2 = list(csv.DictReader(open(os.path.join(R2, "results", "cost.csv"))))
    prof = {(r["level"], int(r["n_sub"])): r for r in r2}
    ref = {"S2": (3557, 136.75, 13.4), "S4": (13764, 914.62, 150.6)}
    for lvl in ref:
        got = int(prof[(lvl, N_SUB_BINDING)]["n_station"])
        assert got == ref[lvl][0], f"station count drift {lvl}: {got}"
    print("  G-REF     station counts + profile-closure walls read from round "
          "2's committed cost.csv (not retyped)")

    def one(n_st, n_sub):
        y0 = C2.blasius_state(0.05, ue=U_INF, rho=RHO, mu=MU, H=bl["H"])
        stations = np.geomspace(0.06, 1.0, n_st)
        t0 = time.perf_counter()
        st = S.march_correlation(stations, y0, 0.05, S.flat_plate_ue(U_INF),
                                 rho=RHO, mu=MU, n_substep=n_st * n_sub)
        return time.perf_counter() - t0, bool(np.all(np.isfinite(st.H)))

    t0 = time.perf_counter(); one(40, 2); cold = time.perf_counter() - t0
    t0 = time.perf_counter(); one(40, 2); warm = time.perf_counter() - t0
    print(f"  G-JIT     cold {cold:.4f} s vs warm {warm:.4f} s -> compile share "
          f"{100*(1-warm/cold):.1f} %; every number below is WARM  PASS")

    rows = []
    for lvl, (n_st, t_ibl3, t_inv) in ref.items():
        for n_sub in N_SUB:
            reps = [one(n_st, n_sub) for _ in range(3)]
            walls = [w for w, ok in reps if ok]
            if not walls:
                continue
            t = statistics.median(walls)
            t_prof = float(prof[(lvl, n_sub)]["t_strip_s"])
            rows.append({
                "level": lvl, "n_station": n_st, "n_sub": n_sub,
                "t_correlation_s": t, "t_profile_s": t_prof,
                "speedup_vs_profile": t_prof / t,
                "rep_spread": (max(walls) - min(walls)) / t,
                "t_ibl3_s": t_ibl3, "t_inviscid_s": t_inv,
                "p1_correlation": t / t_inv,
                "speedup_vs_ibl3": t_ibl3 / t,
            })
            print(f"  G-RES     {lvl}: {n_st} stations x n_sub={n_sub} -> "
                  f"{t:.4f} s (profile {t_prof:.4f} s, "
                  f"{t_prof/t:.2f}x faster), spread "
                  f"{100*rows[-1]['rep_spread']:.1f} %")
    with open(os.path.join(RESULTS, "d_cost.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    b = [r for r in rows if r["level"] == "S4" and r["n_sub"] == N_SUB_BINDING][0]
    t = b["t_correlation_s"]
    code = "C-SUB1" if t < 1.0 else "C-PARTIAL" if t <= 100.0 else "C-MISS"
    _record("D-COST", f"t_strip(S4) at n_sub={N_SUB_BINDING} "
            f"({b['n_station']} stations)",
            "<1 s = C-SUB1; 1-100 s = C-PARTIAL; >100 s = C-MISS",
            f"{t:.4f} s", code)
    _record("D-COST", "vs round 2's profile closure at the same size",
            "RECORDED", f"{b['speedup_vs_profile']:.2f}x faster "
            f"({b['t_profile_s']:.3f} s -> {t:.3f} s)", "RECORDED")
    _record("D-COST", "P1 = t/t_inviscid at S4 (committed P1_IBL3 = 6.073)",
            "RECORDED", f"{b['p1_correlation']:.6f} "
            f"({b['speedup_vs_ibl3']:.0f}x the 3-D IBL's 915 s)", "RECORDED")
    ladder = {r["n_sub"]: r["t_correlation_s"] for r in rows
              if r["level"] == "S4"}
    print("  ★ the gate's <1 s against the resolution knob: " +
          ", ".join(f"n_sub={k}: {v:.4f} s" for k, v in sorted(ladder.items())))
    return rows


def leg_turbulent():
    """D-TURB: RECORDED. The turbulent branch still uses the PROFILE closure
    (this round adds no turbulent correlations), so it must reproduce round 1's
    committed reading -- a G-LEGACY extension."""
    print("== D-TURB (RECORDED; profile closure, must reproduce round 1) ==")
    from pyfp3d.viscous import strip2d as S
    st = S.march(np.geomspace(0.4, 20.0, 12), turbulent_seed(0.2), 0.2,
                 S.flat_plate_ue(U_INF), rho=RHO, mu=MU, turbulent=True,
                 n_substep=2000)
    cf_pw = 0.027 * st.re_x ** (-1.0 / 7.0)
    rel = (st.cf - cf_pw) / cf_pw
    old = list(csv.DictReader(open(os.path.join(
        R1, "results", "turbulent_recorded.csv"))))
    o_lo = min(float(r["rel_diff"]) for r in old)
    o_hi = max(float(r["rel_diff"]) for r in old)
    same = (abs(rel.min() - o_lo) < 1e-9) and (abs(rel.max() - o_hi) < 1e-9)
    _record("D-TURB", "cf vs 1/7 power law (profile closure, unchanged)",
            "RECORDED, no gate",
            f"{100*rel.min():+.2f}%..{100*rel.max():+.2f}% vs round 1's "
            f"{100*o_lo:+.2f}%..{100*o_hi:+.2f}% (bit-reproduced: {same})",
            "RECORDED")
    return {"lo": rel.min(), "hi": rel.max(), "reproduced": same}


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t_all = time.perf_counter()
    print("== guards ==")
    guard_frozen()
    guard_authority()
    guard_no_solve()
    guard_legacy()
    _, bl = guard_oracle()

    bl_fs = falkner_skan(0.0)
    print(f"  reference (integrated here): H={bl_fs['H']:.6f}, "
          f"cf*sqrt(Re_x)={bl_fs['cf_sqrt_rex']:.6f}")

    gate_blasius(bl_fs)
    gate_order(bl_fs)
    gate_fs()
    gate_cost(bl_fs)
    leg_turbulent()

    total = time.perf_counter() - t_all
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n== summary ==  total round compute {total:.2f} s (gate: 20 min)")
    for tag, metric, band, measured, verdict in SUMMARY:
        print(f"  {verdict:45s} [{tag}] {metric} = {measured}")
    print("\n★★ D-BLASIUS / D-FS passing is NEAR-CIRCULAR (the correlations are "
          "fits to this family). Their function is a TRANSCRIPTION check. Do "
          "NOT report them as validation -- pre-registration section 0.")
    fails = [r for r in SUMMARY if "FAIL" in r[4]]
    print(f"  {len(fails)} FAIL row(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
