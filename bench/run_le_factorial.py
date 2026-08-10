"""Task 4: is the M6 LE Cp deficit LE-LOCAL, or the surrounding mesh's pollution floor?

Pre-registered in docs/dev_phase_two/20260805-1600-le-mechanism-prereg.md, committed
BEFORE this file was written -- read that for the legs, the criteria, the five guards
and the cost gates. The short version:

GS2.1 refined the LE only and got a non-monotone -8.3 %, and read that as "the LE
error is not resolution-controlled". But all-scales refinement moves the same band
monotonically at order 0.725. Responding to global refinement while not responding to
local refinement is the signature P11 MEASURED on the sphere: a single-variable sweep
hits the fixed-bulk-mesh pollution floor. GS2.1's own all-scales leg (L3p) was
supposed to cover that and did not converge, so the worry was never tested.

So this runs the full 2x2 GS2.1 lacked -- baseline / bulk-only / LE-only / all-scales
-- which is what makes a null result attributable. At the CHEAP end (h_wall
0.030 -> 0.020, where the all-scales effect is already measured at -30.6 % and the
solves cost 7 s and 61 s), then repeated one level finer to check the answer is not
specific to the coarse end.

Everything about the Cp extraction is imported from run_m3_budget (guard G4): same
seven stations, same five unmasked, same interpolation, same band split. The band
decomposition is asserted to reconstruct the station RMS, exactly as there.

Outputs (TRACKED): bench/gate_results/le_factorial.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import (B_SEMI, MAC, chord_at,           # noqa: E402
                                   onera_m6_wing_mesh, x_le)
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from run_le_response import le_geometry                             # noqa: E402
from run_m3_budget import (ALPHA, BANDS, ETAS, M6_NEWTON_KW,        # noqa: E402
                           M_INF, N_UNMASKED, band_rms,
                           parse_experiment, station_rms)

CSV = os.path.join(HERE, "gate_results", "le_factorial.csv")
SCRATCH = os.environ.get(
    "PYFP3D_SCRATCH",
    "/tmp/claude-1000/-home-lrz-codes-UP3D/"
    "3c5b43c4-b62c-4a09-b4da-9b9c7128d43e/scratchpad")

#: round tip + production taper, per the registration's declared deviations from
#: GS2.1's flat baseline (user ruling 2026-08-04; T0/T1 measured the taper to be what
#: makes an all-scales mesh solvable at all).
BASE = dict(r_far=15.0 * MAC, tip_cap="round", embed_wake=True)
TAPER = ("vanish_smooth", 0.05)

#: (tag, h_wall, h_le) per factorial. h_te / h_wake / h_far follow h_wall.
#: h_le is passed SEPARATELY from h_te -- GS2.1 addendum #1: `h_edge` sizes the LE
#: *and* TE, so not splitting them reads the P13 tip free-edge singularity instead.
FACTORIALS = (
    ("F1", 0.030, 0.020, -30.6),      # coarse end; all-scales effect already measured
    ("F2", 0.020, 0.015, -12.8),      # the same question one level finer
)
#: registration sec 5
LEG_GATE_S = 600.0
TOTAL_GATE_S = 1800.0
GUARD_RTOL = 0.05                     # G1 / G2: "unchanged" means within 5 %


def le_face_count(mc):
    """LE-band wall triangles, by the SAME x/c rule le_geometry uses (0 <= x/c <
    0.15). G1 needs a count as well as a spacing: a leg could keep the median
    spacing and still change how much of the LE is resolved."""
    n = 0
    for face in mc.boundary_faces["wall"]:
        c = mc.nodes[face].mean(axis=0)
        xc = (c[0] - x_le(c[2])) / chord_at(c[2])
        if 0.0 <= xc < 0.15:
            n += 1
    return n


def build(h_wall, h_le):
    return onera_m6_wing_mesh(h_wall=h_wall, h_edge=h_le, h_te=0.5 * h_wall,
                              h_wake=3.0 * h_wall, h_far=120.0 * h_wall,
                              **BASE)


def solve_taper_pressure(mc, wc):
    """The T0/T1 recipe: production taper + pressure Kutta, probe-seeded.

    The drift guard from run_m3_budget.solve is kept and runs BEFORE any intentional
    override, for the reason recorded there -- a deliberate deviation must not trip
    the check that exists to catch unintended drift.
    """
    from run_m3_budget import NEWTON_M6_RECIPE
    kw = dict(NEWTON_M6_RECIPE)
    for k, v in M6_NEWTON_KW.items():
        assert kw["newton_kw"][k] == v, (
            f"the P14 recipe's newton_kw[{k}] = {kw['newton_kw'][k]} no longer "
            f"matches the recorded {v} -- the comparison basis moved")
    taper = tip_taper_factors(wc.station_z, B_SEMI, TAPER[0], TAPER[1] * B_SEMI)
    r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,
                              tip_taper=taper, **M6_NEWTON_KW)
    kw["newton_kw"] = dict(kw["newton_kw"], tip_taper=taper,
                           kutta_estimator="pressure", phi_init=r0["phi"],
                           gamma_init=r0["gamma"])
    return solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)


def bands_of(mc, phi, exp):
    """Band RMS over the five unmasked stations, run_m3_budget's accumulation
    verbatim including its reconstruction assert (guard G4)."""
    curves = {eta: section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI,
                                    m_inf=M_INF) for eta in ETAS}
    acc, per = {}, {}
    for eta in ETAS[:N_UNMASKED]:
        per[eta] = station_rms(curves, exp, eta)[0]
        b = band_rms(curves, exp, eta)
        for kk, (ss, nn) in b.items():
            a0, n0 = acc.get(kk, (0.0, 0))
            acc[kk] = (a0 + ss, n0 + nn)
        tot_ss = sum(v[0] for v in b.values())
        tot_n = sum(v[1] for v in b.values())
        assert abs((tot_ss / max(tot_n, 1)) ** 0.5 - per[eta]) < 1e-12, (
            f"band decomposition does not reconstruct the station RMS at "
            f"eta={eta} -- it is measuring something else")
    out = {}
    for name, _, _ in BANDS:
        for side in ("upper", "lower"):
            ss, nn = acc[f"{name}_{side}"]
            out[f"{name}_{side}"] = (ss / max(nn, 1)) ** 0.5 if nn else np.nan
    all_ss = sum(v[0] for v in acc.values())
    all_n = sum(v[1] for v in acc.values())
    out["pooled"] = (all_ss / max(all_n, 1)) ** 0.5
    return out


def one_leg(tag, h_wall, h_le, exp, rows):
    t0 = time.perf_counter()
    mesh = build(h_wall, h_le)
    mc, wc = cut_wake(mesh)
    ht, hn, aniso = le_geometry(mc)
    n_le = le_face_count(mc)
    n_tet = len(mc.elements)
    t_gen = time.perf_counter() - t0
    t1 = time.perf_counter()
    r = solve_taper_pressure(mc, wc)
    wall = time.perf_counter() - t1
    hist = np.asarray(r.get("residual_history", []), dtype=float)
    res = float(hist[-1]) if len(hist) else float("nan")
    nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
    # G3: a non-converged leg's RMS is NOT a reading (the L3p lesson)
    usable = bool(r.get("converged")) and res < 1e-9 and nlim == 0 and nflr == 0
    b = bands_of(mc, np.asarray(r["phi"]), exp) if usable else {}
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    row = dict(tag=tag, h_wall=h_wall, h_le=h_le, n_tet=n_tet, n_le_faces=n_le,
               le_ht=round(ht, 8), le_hn=round(hn, 8), le_aniso=round(aniso, 4),
               converged=bool(r.get("converged")), res_final=res,
               n_limited=nlim, n_floored=nflr, usable=usable,
               m_max=round(float(np.sqrt(r["mach2_max"])), 5),
               accept_reason=r.get("accept_reason"), s_ref=round(s_ref, 6),
               t_gen_s=round(t_gen, 1), wall_s=round(wall, 1),
               **{f"band_{k}": (round(v, 6) if v == v else None)
                  for k, v in b.items()})
    rows.append(row)
    print(f"  {tag:18} tets {n_tet:>8}  LE faces {n_le:>5}  h_t {ht:.5f}  "
          f"conv={str(row['converged']):5} |R|={res:.2e} lim/flr={nlim}/{nflr}  "
          f"LE_up {row.get('band_LE_upper')}  pooled {row.get('band_pooled')}  "
          f"({t_gen:.0f}+{wall:.0f}s)", flush=True)
    if not usable:
        print(f"    ★ NOT USABLE (G3): its RMS is not a reading", flush=True)
    return row


def main():
    exp = parse_experiment()
    for eta, e in exp.items():                      # run_m3_budget's W2 guard
        if e["x"][int(np.argmax(e["cp"]))] >= 0.05:
            raise RuntimeError(f"W2 guard: eta={eta} max-Cp not at the LE")
    os.makedirs(SCRATCH, exist_ok=True)
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    rows, spent = [], 0.0
    print("LE factorial: is the LE deficit LE-local or the bulk mesh's pollution "
          "floor?\n(pre-registered 20260805-1600; round tip + production taper + "
          "pressure Kutta)\n")
    for fac, h0, h1, allscale_pct in FACTORIALS:
        print(f"--- factorial {fac}: h_wall {h0} -> {h1} "
              f"(all-scales effect already measured {allscale_pct:+.1f} %) ---")
        legs = ((f"{fac}_00_base", h0, 0.5 * h0),
                (f"{fac}_10_bulkonly", h1, 0.5 * h0),     # h_le HELD
                (f"{fac}_01_leonly", h0, 0.5 * h1),       # bulk HELD
                (f"{fac}_11_allscales", h1, 0.5 * h1))
        for tag, hw, hl in legs:
            t = time.perf_counter()
            one_leg(tag, hw, hl, exp, rows)
            dt = time.perf_counter() - t
            spent += dt
            if dt > LEG_GATE_S:
                print(f"  ★ leg gate: {dt:.0f} s > {LEG_GATE_S:.0f} s -- "
                      f"stopping factorial {fac} (registration sec 5)")
                break
        if spent > TOTAL_GATE_S:
            print(f"★ total gate: {spent:.0f} s > {TOTAL_GATE_S:.0f} s -- "
                  f"stopping after {fac}; later factorials NOT run")
            break
    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    _verdict(rows)
    return 0


def _verdict(rows):
    by = {r["tag"]: r for r in rows}
    print("\n=== verdict (criteria fixed in the registration, sec 3) ===")
    for fac, _h0, _h1, allscale_pct in FACTORIALS:
        need = {k: by.get(f"{fac}_{k}") for k in
                ("00_base", "10_bulkonly", "01_leonly", "11_allscales")}
        if any(v is None for v in need.values()):
            print(f"  {fac}: incomplete ({sum(v is None for v in need.values())} "
                  f"leg(s) missing) -- nothing read")
            continue
        base = need["00_base"]
        # --- guards first: a leg that failed a guard is not read -------------
        bad = [k for k, v in need.items() if not v["usable"]]
        if bad:
            print(f"  {fac}: G3 FAIL on {bad} -- their RMS are not readings")
        g1 = need["10_bulkonly"]
        d_ht = abs(g1["le_ht"] - base["le_ht"]) / base["le_ht"]
        d_nf = abs(g1["n_le_faces"] - base["n_le_faces"]) / base["n_le_faces"]
        g1_ok = d_ht < GUARD_RTOL and d_nf < GUARD_RTOL
        print(f"  {fac} G1 (bulk-only kept the LE fixed): h_t {100*d_ht:+.2f} %, "
              f"LE faces {100*d_nf:+.2f} %  -> {'PASS' if g1_ok else 'FAIL'}")
        g2 = need["01_leonly"]
        d_nt = abs(g2["n_tet"] - base["n_tet"]) / base["n_tet"]
        print(f"  {fac} G2 (LE-only kept the bulk fixed): tets "
              f"{100*d_nt:+.2f} %  -> {'(bulk grows with LE cells; RECORDED)'}")
        if bad or not g1_ok:
            print(f"  {fac}: guards not clean -> RECORDED, no conclusion")
            continue
        b0 = base["band_LE_upper"]
        d10 = 100.0 * (need["10_bulkonly"]["band_LE_upper"] - b0) / b0
        d01 = 100.0 * (need["01_leonly"]["band_LE_upper"] - b0) / b0
        d11 = 100.0 * (need["11_allscales"]["band_LE_upper"] - b0) / b0
        print(f"  {fac} LE-upper RMS: base {b0:.4f} | bulk-only {d10:+.1f} % | "
              f"le-only {d01:+.1f} % | all-scales {d11:+.1f} % "
              f"(prior all-scales {allscale_pct:+.1f} %)")
        thresh = 0.60 * abs(d11)
        if abs(d10) >= thresh:
            print(f"  ⇒ P: bulk-only reproduces {abs(d10)/max(abs(d11),1e-9):.0%} "
                  f"of all-scales (>= 60 %) -- the LE deficit is POLLUTION-"
                  f"DOMINATED, not LE-local. S2's target moves to the bulk mesh.")
        elif abs(d10) < 10.0:
            print(f"  ⇒ N: bulk-only only {abs(d10):.1f} % (< 10 %) -- the "
                  f"pollution hypothesis is REFUTED too; the 0.725 order comes "
                  f"from a third mechanism and needs its own round.")
        else:
            print(f"  ⇒ I: {abs(d10):.1f} % is between the thresholds "
                  f"({thresh:.1f} % and 10 %) -- RECORDED, no conclusion.")
        print(f"     additivity (recorded, not a criterion): "
              f"|d10|+|d01| = {abs(d10)+abs(d01):.1f} % vs |d11| = "
              f"{abs(d11):.1f} %")


if __name__ == "__main__":
    sys.exit(main())
