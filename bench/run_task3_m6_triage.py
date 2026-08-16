"""Triage 1: is M6 coarse at M0.8395/alpha3.06 inside the multi-solution regime?

Pre-registered in docs/dev_phase_three/20260813-0300-m6-triage-prereg.md, with addendum #1 (WHICH of
the two seed sites to vary) committed before this file was written.

The question that matters practically is not "does M6 have multiple solutions" but "are the COMMITTED
M6 anchors a single arbitrary draw". Those are different, because solve_newton_transonic only
cold-seeds LEVEL 0 and every later level warm-starts -- the project's own note is that a ramp and a
seed are two implementations of one function, so the RAMP may be the thing pinning the solution.

  arm A (ramp, production recipe, three seeds) -> Q-A: are the committed anchors reproducible?
  arm B (direct solve at M0.8395, no ramp)     -> Q-B: does the system itself have several solutions?

★ If A is stable and B scatters, the conclusion is "the ramp pins the selection", NOT "M6 has no
non-uniqueness". The registration forbids the second sentence, and so does this docstring.

★★ Addendum #1: the seed is varied at the M0.70 PROBE solve, not at ramp level 0. Level 0 already
receives phi_init, and Picard seeding only acts when there is no warm start -- so varying only that
site would give three BIT-IDENTICAL legs and trip kill clause 4 for the wrong reason (a disconnected
instrument, rather than a seed with no effect). That clause caught this before execution.

Everything else is run_m3_budget's production leg VERBATIM, imported rather than re-typed so the
recipe cannot drift.

Outputs (TRACKED): bench/gate_results/task3_m6_triage.csv
"""

import csv
import os
import sys
import time

#: ★ addendum #3: SIXTEEN, matching run_m3_budget.py's own setdefault and therefore the provenance
#: of the committed reference row. The first execution used 8 and compared its result against a
#: 16-thread committed number -- a cross-provenance comparison, the third instance of that family in
#: this round. The 8-thread arm is kept as RECORDED (its own finding: the production leg does NOT
#: converge on M6 medium at 8 threads) in gate_results/task3_m6_triage_8t.csv.
os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,            # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
#: imported, never re-typed: the recipe and the band/experiment machinery are the single source
from run_m3_budget import (ALPHA, B_SEMI, BANDS, ETAS, M6_NEWTON_KW, M_INF,  # noqa: E402
                           band_rms, parse_experiment)
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

CSV = os.path.join(HERE, "gate_results",
                   f"task3_m6_triage_{os.environ['NUMBA_NUM_THREADS']}t.csv")
SEEDS = (0, 5, 12)
#: ★ addendum #2: MEDIUM, not coarse. The committed gate_results/m3_budget.csv says level = medium
#: on all four rows -- I had read the `main(levels=("coarse",))` DEFAULT instead of the EVIDENCE.
LEVEL = "medium"
#: ★ addendum #2: the reference is THIS pipeline's own committed ON/pressure row, not P14's anchor.
#: P14 is a different pipeline and this script has no committed coarse row at all, so the original
#: guard compared against something this leg never produced ("quote the one your pipeline ran").
REF_ROW = dict(cl_p=0.276527, cl_kj=0.281179,
               rms_LE_upper=0.235615, rms_LE_lower=0.092155)
CL_TOL_PCT, LE_TOL_PCT = 0.5, 2.5
LEG_GATE_S = 600.0


def production_solve(mc, wc, seed, ramp):
    """The (entropy=True, kutta='pressure') production leg, with the seed at the PROBE solve.

    ★ addendum #1: seeding the probe is what changes the draw; ramp level 0 keeps
    n_picard_seed = 0 verbatim because phi_init makes it inert there.
    """
    kw = dict(NEWTON_M6_RECIPE)
    for k, v in M6_NEWTON_KW.items():
        assert kw["newton_kw"][k] == v, (
            f"recipe drift: newton_kw[{k}] = {kw['newton_kw'][k]} != recorded {v}")
    kw["newton_kw"] = dict(kw["newton_kw"], entropy_correction=True)
    probe_kw = dict(M6_NEWTON_KW, n_picard_seed=seed)
    r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,
                              entropy_correction=True, **probe_kw)
    kw["newton_kw"].update(kutta_estimator="pressure", phi_init=r0["phi"],
                           gamma_init=r0["gamma"], n_picard_seed=0)
    if ramp:
        return solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw), r0
    #: arm B: the same state, no Mach ramp. GS3.3a measured this as slow and non-converged,
    #: which is an ALLOWED outcome recorded as UNDEFINED -- never as "no non-uniqueness here".
    direct = dict(kw["newton_kw"])
    return solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **direct), r0


def _f(x, nd=6):
    return "-" if x is None else format(x, f".{nd}f")


def metrics(mc, wc, r, s_ref, exp):
    phi = np.asarray(r["phi"])
    gamma = np.atleast_1d(np.asarray(r["gamma"]))
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    o = np.argsort(wc.station_z)
    clkj = float(cl_kj_3d(gamma[o], wc.station_z[o], s_ref, B_SEMI))
    #: ★ band_rms keys are per-band-per-SIDE (LE_upper / LE_lower / ...), matching the committed
    #: CSV's column names -- read from its output rather than assumed from BANDS, which is what the
    #: first run got wrong (KeyError 'LE_upper').
    curves, acc = {}, {}
    for eta in ETAS:
        try:
            curves[eta] = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        except Exception:                                        # noqa: BLE001
            continue
    for eta in list(curves):
        if eta not in exp:
            continue
        for name, (ss, nn) in band_rms(curves, exp, eta).items():
            a = acc.setdefault(name, [0.0, 0])
            a[0] += ss
            a[1] += nn
    out = dict(cl_p=float(f["cl"]), cl_kj=clkj)
    for name, (ss, nn) in sorted(acc.items()):
        out[f"rms_{name}"] = (float(np.sqrt(ss / nn)) if nn else None)
    return out


def main():
    exp = parse_experiment()
    path = os.path.join(REPO, "cases", "meshes", "onera_m6", f"{LEVEL}.msh")
    if not os.path.exists(path):
        print(f"★ mesh missing: {path} -- regenerate via generate_onera_m6.py")
        return 1
    mc, wc = cut_wake(read_mesh(path))
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    print(f"M6 triage   {LEVEL}   M{M_INF}/alpha {ALPHA}   seeds {SEEDS}   "
          f"threads {os.environ['NUMBA_NUM_THREADS']}\n")

    rows = []
    #: ★ addendum #2: arm B (direct medium) is DEFERRED -- a medium ramp is already ~515 s, close to
    #: the 600 s per-leg gate, and a direct medium solve would likely burn 80 steps for a
    #: RECORDED-only return. Whether to spend it is decided AFTER arm A says something.
    for arm, ramp in (("A_ramp", True),):
        for seed in SEEDS:
            t0 = time.perf_counter()
            try:
                r, r0 = production_solve(mc, wc, seed, ramp)
                err = ""
            except Exception as exc:                             # noqa: BLE001
                wall = time.perf_counter() - t0
                rows.append(dict(arm=arm, seed=seed, converged=False,
                                 error=type(exc).__name__, wall_s=round(wall, 1)))
                print(f"  {arm:9} seed {seed:>2}  ★ RAISED {type(exc).__name__}: "
                      f"{str(exc)[:80]}  ({wall:.0f}s)", flush=True)
                continue
            wall = time.perf_counter() - t0
            m = metrics(mc, wc, r, s_ref, exp)
            h = list(r.get("residual_history") or [])
            fr = r.get("sigma_freeze_report") or {}
            row = dict(arm=arm, seed=seed, converged=bool(r.get("converged")),
                       res_final=(h[-1] if h else None), error=err,
                       m_final=r.get("m_final"), m_last=r.get("m_last_converged"),
                       n_limited=int(r.get("n_limited") or 0),
                       n_floored=int(r.get("n_floored") or 0),
                       sigma_min=r.get("sigma_min"), m1_max=r.get("m1_max"),
                       probe_cl=None, selection_churn=fr.get("selection_churn"),
                       frozen_in_transient=fr.get("frozen_in_transient"),
                       wall_s=round(wall, 1), **m)
            rows.append(row)
            print(f"  {arm:9} seed {seed:>2}  conv={str(row['converged']):5} "
                  f"|R|={(row['res_final'] if row['res_final'] is not None else float('nan')):.2e} "
                  f"m_final={row['m_final']} cl_p={row['cl_p']:.6f} cl_KJ={row['cl_kj']:.6f} "
                  f"rms_LE(up/lo)={_f(row.get('rms_LE_upper'))}/{_f(row.get('rms_LE_lower'))}"
                  f" ({wall:.0f}s)", flush=True)
            if wall > LEG_GATE_S:
                print("    ★ cost gate exceeded -- stop"); _write(rows); return 1
    _write(rows)
    return _read(rows)


def _read(rows):
    A = [r for r in rows if r["arm"] == "A_ramp"]
    conv = [r for r in A if r.get("converged")]

    print("\n=== G-R: the seed-0 medium leg must reproduce THIS pipeline's own committed ON/pressure row ===")
    m = [r for r in A if r["seed"] == 0 and r.get("converged")]
    if not m:
        print("  ★ seed 0 did not converge -- cannot check the instrument. STOP (kill clause 1).")
        return 1
    devs = {}
    for k, want in REF_ROW.items():
        got = m[0].get(k)
        devs[k] = None if got is None else 100.0 * abs(got - want) / abs(want)
        print(f"  {k:14} {_f(got)} vs {want}  = "
              f"{'-' if devs[k] is None else format(devs[k], '.3f') + ' %'}")
    worst_cl = max(v for k, v in devs.items() if k.startswith("cl") and v is not None)
    print(f"  worst cl deviation {worst_cl:.3f} %   (tol {CL_TOL_PCT} %)")
    if worst_cl > CL_TOL_PCT:
        print("  -> ★ G-R FAIL: this is not the production path. Kill clause 1 -- stop.")
        return 1
    print("  -> G-R PASS")

    print("\n=== kill clause 4: bit-identical legs are SUSPICIOUS, not passing ===")
    _lk = sorted(k for k in conv[0] if k.startswith("rms_LE")) if conv else []
    ident = (len(conv) >= 2 and all(r["cl_p"] == conv[0]["cl_p"] for r in conv)
             and all(all(r[k] == conv[0][k] for k in _lk) for r in conv))
    print(f"  arm A converged legs: {len(conv)}/3   all bit-identical: {ident}")
    if ident:
        print("  -> ★ the seed never reached the solve. Check the instrument before reading "
              "anything (kill clause 4).")
        return 1

    print("\n=== T1-T4 on the balanced panel (arm A, converged legs only) ===")
    if len(conv) < 2:
        print(f"  -> T4  only {len(conv)} converged leg(s) -> UNDEFINED (NOT 'a small spread')")
        return 0

    def spread_pct(key):
        v = [r[key] for r in conv if r.get(key) is not None]
        if len(v) < 2:
            return None
        return 100.0 * (max(v) - min(v)) / max(abs(float(np.mean(v))), 1e-12)

    s_p, s_k = spread_pct("cl_p"), spread_pct("cl_kj")
    #: ★ LE has an upper and a lower side; the committed 4-variant CSV showed 0.3 % and 2.5 %, so
    #: BOTH are reported and the LARGER one binds -- reporting only the quiet side would be the
    #: "quote the one your pipeline actually ran" mistake.
    le_keys = sorted(k for k in rows[0] if k.startswith("rms_LE"))
    le_vals = {k: spread_pct(k) for k in le_keys}
    s_le = max((v for v in le_vals.values() if v is not None), default=None)
    for key in ("cl_p", "cl_kj"):
        v = s_p if key == "cl_p" else s_k
        print(f"  {key:14} spread {'UNDEF' if v is None else format(v, '.3f') + ' %'}")
    for k in sorted(kk for kk in rows[0] if kk.startswith("rms_")):
        v = spread_pct(k)
        print(f"  {k:14} spread {'UNDEF' if v is None else format(v, '.3f') + ' %'}"
              f"{'   <- LE binding' if v == s_le and k.startswith('rms_LE') else ''}")
    cl_bad = max(x for x in (s_p, s_k) if x is not None) > CL_TOL_PCT
    le_bad = s_le is not None and s_le > LE_TOL_PCT
    if le_bad:
        band = "T3"
        why = (f"LE band spread {s_le:.3f} % > {LE_TOL_PCT} % ⇒ the 70 % needs re-basing, and "
               f"'LE is 70 % of the budget' is meaningless at this condition")
    elif cl_bad:
        band = "★ T2"
        why = ("lift spreads but the LE band does not ⇒ the protocol's hypothesis is CONFIRMED: "
               "the lift carries the non-uniqueness and the LE deficit is separable")
    else:
        band = "★ T1"
        why = ("both within tolerance ⇒ the committed M6 coarse anchors are reproducible under "
               "the production recipe; the LE 70 % stands")
    print(f"\n  -> {band}   {why}")

    print("\n=== arm B (RECORDED -- it answers Q-B, which does NOT decide anchor credibility) ===")
    B = [r for r in rows if r["arm"] == "B_direct"]
    bc = [r for r in B if r.get("converged")]
    for r in B:
        print(f"  seed {r['seed']:>2}  conv={str(r.get('converged')):5} "
              f"{'err=' + r['error'] if r.get('error') else ''} "
              f"cl_p={'-' if 'cl_p' not in r else format(r['cl_p'], '.6f')} "
              f"({r['wall_s']:.0f}s)")
    if len(bc) < 2:
        print("  -> Q-B is NOT measurable with this recipe (fewer than two converged direct legs).")
        print("     ★ This may NOT be read as 'arm B shows no non-uniqueness'.")
    else:
        v = [r["cl_p"] for r in bc]
        sp = 100.0 * (max(v) - min(v)) / max(abs(float(np.mean(v))), 1e-12)
        print(f"  -> direct-arm cl_p spread {sp:.3f} % over {len(bc)} legs")
        if sp > CL_TOL_PCT and not cl_bad:
            print("     ★ B scatters while A is stable ⇒ THE RAMP PINS THE SELECTION.")
            print("       That is the supported sentence; 'M6 has no non-uniqueness' is NOT.")
    return 0


def _write(rows):
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")


if __name__ == "__main__":
    sys.exit(main())
