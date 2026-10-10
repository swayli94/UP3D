"""Step 0 of the LE-mechanism factorial: where is the window in which a SHOCK is
present AND refinement is usable?

Pre-registered in phases/p2/docs/dev_phase_two/20260806-0000-le-mechanism-prereg2.md sec 1b
(addendum #1, user 2026-08-06, committed before this file was written).

Why this step exists: the first-version factorial found that at M0.8395 every
refinement direction drives the peak Mach into m_cap, so no refinement arm has a
reading. My redesign escaped that by going subcritical -- and the user pointed out
that this downgrades the QUESTION rather than the condition, because with no
supersonic cell anywhere the upwinding is inert and the entropy correction is
bit-identical ON and OFF (tests/test_s1b_entropy.py::
test_scope_boundary_no_shock_means_bit_identical). At M0.50 what gets measured is the
pure P1 wall error P11 already characterised, not the mechanism behind M3's 69.6 % LE
share. And M0.70 cannot be adopted by assumption either -- NEWTON_M6_RECIPE's comment
says the ramp's first level at 0.70 is subcritical.

So the condition is MEASURED, not chosen. A Mach is in-window iff, per the
registration: both legs converge with zero clamps, both carry a supersonic zone, and at
least one station yields a locatable shock. The factorial then runs at the HIGHEST
in-window Mach. If the window is EMPTY that is the result, and the registration forbids
retreating to subsonic to manufacture a runnable experiment.

★ Only TWO meshes are needed for all six Machs -- the mesh depends on
(h_wall, h_le, h_far) and not on M_inf -- which also makes the G1 surface-identity
check a one-off rather than per-leg.

Outputs (TRACKED): bench/gate_results/le_window.csv
"""

import csv
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import (B_SEMI, MAC,                     # noqa: E402
                                   onera_m6_wing_mesh)
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from run_le_response import le_geometry                             # noqa: E402
from run_le_factorial import le_face_count                          # noqa: E402
from run_m3_budget import ALPHA, M6_NEWTON_KW, NEWTON_M6_RECIPE     # noqa: E402
from run_seed_exposure import clamp_map                             # noqa: E402

CSV = os.path.join(_GATE, "le_window.csv")
#: ★ 0.60 REMOVED after measurement: solve_newton_transonic raises when
#: m_inf < m_start, and the recipe's m_start is 0.70, so both M0.60 legs died
#: without producing a reading. No information is lost -- M0.70 measured M_max
#: 1.5358 with 214 shock cells, so 0.60 could not have been the shock-free
#: control the scan would have wanted it for.
MACHS = (0.70, 0.75, 0.78, 0.80, 0.8395)
H_WALL, H_LE = 0.020, 0.010
#: (tag, h_far). The far field is the ONLY thing that differs -- that is the point of
#: using it as the bulk arm: it does not touch the surface mesh at all.
FARS = (("base", 2.4), ("far_x2", 1.2))
TAPER = ("vanish_smooth", 0.05)
ETAS_SHOCK = (0.44, 0.65, 0.90)
#: registration sec 6 sets 400 s/leg and 1200 s for the FOUR-leg factorial. Step 0 has
#: twelve legs, so the per-leg gate is carried over unchanged and the total is stated
#: here as an explicit extension rather than silently ignored.
LEG_GATE_S = 400.0
TOTAL_GATE_S = 2400.0


def build(h_far):
    return onera_m6_wing_mesh(h_wall=H_WALL, h_edge=H_LE, h_te=0.5 * H_WALL,
                              h_wake=3.0 * H_WALL, h_far=h_far,
                              r_far=15.0 * MAC, tip_cap="round",
                              embed_wake=True)


def solve_at(mc, wc, m_inf):
    """The production transonic recipe at an arbitrary target: taper + pressure Kutta,
    probe-seeded, ramping from the recipe's own m_start. Same as run_le_factorial's
    solve_taper_pressure except the target Mach is a parameter, and the drift guard on
    M6_NEWTON_KW is kept for the reason recorded there."""
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
    return solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=ALPHA, **kw)


def main():
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    print("LE window scan (prereg 20260806-0000 sec 1b): where is a shock present AND "
          "refinement usable?\n")
    meshes, geom = {}, {}
    for tag, h_far in FARS:
        t0 = time.perf_counter()
        mc, wc = cut_wake(build(h_far))
        ht, hn, aniso = le_geometry(mc)
        geom[tag] = dict(n_tet=len(mc.elements),
                         n_wall=len(mc.boundary_faces["wall"]),
                         n_le=le_face_count(mc), le_ht=ht)
        meshes[tag] = (mc, wc)
        print(f"  mesh {tag:7} h_far {h_far:<5} tets {geom[tag]['n_tet']:>8} "
              f"wall tris {geom[tag]['n_wall']:>6} LE faces {geom[tag]['n_le']:>5} "
              f"h_t {ht:.8f}  ({time.perf_counter()-t0:.0f}s)", flush=True)
    # ---- G1, once: refining h_far must leave the SURFACE mesh untouched -------
    b, f = geom["base"], geom["far_x2"]
    d_ht = abs(f["le_ht"] - b["le_ht"]) / b["le_ht"]
    g1 = (f["n_wall"] == b["n_wall"] and f["n_le"] == b["n_le"] and d_ht < 1e-12)
    print(f"\n  G1 (h_far must not touch the surface): wall tris {b['n_wall']} vs "
          f"{f['n_wall']}, LE faces {b['n_le']} vs {f['n_le']}, "
          f"h_t rel {d_ht:.2e}  -> {'PASS' if g1 else 'FAIL'}")
    if not g1:
        print("  ★ G1 FAIL: h_far is NOT a clean bulk-only knob either. Per the "
              "registration this must be reported, not worked around -- the scan "
              "below is still run and RECORDED, but no leg may be read as "
              "'bulk-only'.", flush=True)
    print(f"  bulk cells: {b['n_tet']} -> {f['n_tet']} "
          f"({100.0 * (f['n_tet'] - b['n_tet']) / b['n_tet']:+.1f} %)\n")

    rows, spent = [], 0.0
    for m in MACHS:
        for tag, h_far in FARS:
            mc, wc = meshes[tag]
            t0 = time.perf_counter()
            try:
                r = solve_at(mc, wc, m)
            except Exception as exc:                               # noqa: BLE001
                rows.append(dict(m_inf=m, leg=tag, h_far=h_far, converged=False,
                                 note=f"{type(exc).__name__}: {exc}"))
                print(f"  M{m} {tag:7} DIED {type(exc).__name__}", flush=True)
                continue
            dt = time.perf_counter() - t0
            spent += dt
            hist = np.asarray(r.get("residual_history", []), dtype=float)
            nlim = int(r.get("n_limited") or 0)
            nflr = int(r.get("n_floored") or 0)
            m_ok = r.get("m_last_converged", r.get("m_final"))
            loc = clamp_map(mc, np.asarray(r["phi"]), m)
            shocks = {}
            for eta in ETAS_SHOCK:
                c = section_cp_curve(mc, np.asarray(r["phi"]), eta=eta,
                                     b_semi=B_SEMI, m_inf=m)
                shocks[eta] = shock_report(c, m)["upper"].get("x_shock")
            n_shock = r.get("n_shock_cells")
            mmax = float(np.sqrt(r["mach2_max"]))
            row = dict(
                m_inf=m, leg=tag, h_far=h_far, n_tet=len(mc.elements),
                converged=bool(r.get("converged")),
                res_final=(float(hist[-1]) if len(hist) else None),
                m_last_converged=m_ok, m_max=round(mmax, 5),
                n_shock_cells=n_shock, n_limited=nlim, n_floored=nflr,
                fallback_fired=bool(r.get("seed_fallback", {}).get("fired")),
                accept_reason=r.get("accept_reason"),
                **{f"x_shock_{e}": shocks[e] for e in ETAS_SHOCK},
                **loc, wall_s=round(dt, 1), note="")
            rows.append(row)
            got = [f"{e}:{'-' if shocks[e] is None else f'{shocks[e]:.3f}'}"
                   for e in ETAS_SHOCK]
            print(f"  M{m} {tag:7} conv={str(row['converged']):5} "
                  f"|R|={row['res_final']:.1e} M_max={mmax:.4f} "
                  f"shock_cells={n_shock} lim/flr={nlim}/{nflr} "
                  f"x_shock {' '.join(got)} fb={row['fallback_fired']} "
                  f"({dt:.0f}s)", flush=True)
            if nlim or nflr:
                print(f"        clamps: LE {loc['clamp_frac_LE_0_15']:.2f} / "
                      f"MID {loc['clamp_frac_MID_15_85']:.2f} / "
                      f"TE {loc['clamp_frac_TE_85_100']:.2f} / "
                      f"off {loc['clamp_frac_offbody']:.2f}", flush=True)
            if dt > LEG_GATE_S:
                print(f"  ★ leg gate {dt:.0f} s > {LEG_GATE_S:.0f} s -- stopping",
                      flush=True)
                spent = TOTAL_GATE_S + 1
                break
        if spent > TOTAL_GATE_S:
            print(f"★ total gate {spent:.0f} s > {TOTAL_GATE_S:.0f} s -- stopping; "
                  f"higher Machs NOT run", flush=True)
            break

    keys = sorted({k for r in rows for k in r})
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    # ---- the window, by the registration's own three conditions --------------
    print("\n=== the window (registration sec 1b: both legs converge unclamped, "
          "both supersonic, a shock locatable) ===")
    inwin = []
    for m in MACHS:
        #: ★ only rows that actually carry readings -- a leg that DIED has just
        #: m_inf/leg/converged/note, and the first version of this loop indexed
        #: n_shock_cells on those and crashed AFTER the CSV was written. The data
        #: survived; the verdict did not.
        legs = [r for r in rows if r.get("m_inf") == m and "n_shock_cells" in r]
        if len(legs) < 2:
            print(f"  M{m}: {len(legs)}/2 legs -- not assessable")
            continue
        c_ok = all(r["converged"] and not r["n_limited"] and not r["n_floored"]
                   for r in legs)
        sup = all((r["n_shock_cells"] or 0) > 0 or r["m_max"] > 1.0 for r in legs)
        sh = any(r[f"x_shock_{e}"] is not None for r in legs for e in ETAS_SHOCK)
        ok = c_ok and sup and sh
        print(f"  M{m}: converged-unclamped {c_ok} | supersonic {sup} | "
              f"shock locatable {sh}  -> {'IN WINDOW' if ok else 'out'}")
        if ok:
            inwin.append(m)
    if inwin:
        print(f"\n  ⇒ window = {inwin};  the factorial runs at M = {max(inwin)}")
    else:
        print("\n  ⇒ ★★ THE WINDOW IS EMPTY. Per the registration this is the result:")
        print("     on this solver 'has a shock' and 'can be refined' are mutually")
        print("     exclusive, so the LE deficit cannot be attributed by refinement")
        print("     until the m_cap / tip side is fixed. Do NOT retreat to subsonic.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
