"""GS2.1 step two: cure the tip, then read P11's all-scales lever.

Pre-registered in docs/dev_phase_two/20260801-0800-gs21-tip-prereg.md, committed before
this file was written. Read that for the legs, for G1 (the tip-cure gate that must pass
before any refinement reading is attributable) and for G2's same-recipe baseline rule.

The short version: step one measured LE resolution dead and the geometry clean, and left
all-scales refinement UNREAD because it died on the wing-tip singularity that the
baseline already carries at M_max 2.1005. B31/B32's tip_taper plus the pressure Kutta is
the existing fix; this measures whether it makes the lever readable, and then reads it.

Outputs (TRACKED): bench/gate_results/tip_allscales.csv
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
from pyfp3d.mesh.metrics import (compute_aspect_ratios,             # noqa: E402
                                 compute_min_dihedral_angles,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh, write_mesh               # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,           # noqa: E402
                                 wall_force_coefficients)
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from run_le_response import SCRATCH, build, le_geometry             # noqa: E402
from run_m3_budget import (ALPHA, BANDS, ETAS, M6_NEWTON_KW,        # noqa: E402
                           M_INF, N_UNMASKED, band_rms,
                           parse_experiment, station_rms)
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

OUT = os.path.join(HERE, "gate_results")
os.makedirs(OUT, exist_ok=True)

#: B32's adopted production taper (cases/demo/b18_wingbody_transonic CONF_TAPER)
TAPER_FORM, TAPER_RC_FRAC = "vanish_smooth", 0.05
#: (tag, h_wall, h_le, h_te). T1 reuses L3p's cached mesh by construction.
LEGS = (
    ("Tm1_coarse_taper", 0.030, 0.015, 0.015),
    ("T0_baseline_taper", 0.015, 0.0075, 0.0075),
    ("T1_allscales_taper", 0.010, 0.00375, 0.005),
    ("T2_allscales2_taper", 0.0075, 0.001875, 0.00375),
)
#: L0p's readings, for G1's tip-cure comparison
L0P = dict(m_max=2.1005, top8_in_tip=8, le_upper=0.2361, cl_p=0.263528)
COST_GATE_S = 1800.0            # pre-registered; stops the sweep
T2_GATE_S = 1200.0              # T2 only runs if T1 comes in under this
#: L3p's mesh was cached under that tag by the step-one sweep
MESH_ALIAS = {"T0_baseline_taper": "L0p_split_control",
              "T1_allscales_taper": "L3p_allscales"}


def solve_tapered(mc, wc):
    """B32's adopted production configuration: pressure Kutta + tip_taper, seeded
    from an M0.70 probe solve per the P14 recipe (the quadratic Kutta row has the
    smaller basin and the transonic driver only cold-seeds level 0)."""
    taper = tip_taper_factors(wc.station_z, B_SEMI, TAPER_FORM,
                              TAPER_RC_FRAC * B_SEMI)
    r0 = solve_newton_lifting(mc, wc, m_inf=0.70, alpha_deg=ALPHA,
                              **M6_NEWTON_KW)
    kw = dict(NEWTON_M6_RECIPE)
    for k, v in M6_NEWTON_KW.items():
        assert kw["newton_kw"][k] == v, (
            f"the P14 recipe's newton_kw[{k}] moved -- comparison basis")
    kw["newton_kw"] = dict(kw["newton_kw"], kutta_estimator="pressure",
                           tip_taper=taper, phi_init=r0["phi"],
                           gamma_init=r0["gamma"], n_picard_seed=0)
    return solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)


def load_prior():
    """Rows recorded by an earlier invocation, so a leg already measured need not be
    re-solved just to complete a sequence. Merged rows are LABELLED
    (from_prior_run=True) -- the pre-registration's order is "run the cheap legs
    first and evaluate the prediction against the RECORDED T1", not "quietly reuse
    numbers"."""
    path = os.path.join(OUT, "tip_allscales.csv")
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if not row.get("band_LE_upper"):
                continue
            d = {}
            for k, v in row.items():
                if v in ("", None):
                    d[k] = None
                elif v in ("True", "False"):
                    d[k] = v == "True"
                else:
                    try:
                        d[k] = float(v)
                    except ValueError:
                        d[k] = v
            d["leg"] = row["leg"]
            d["from_prior_run"] = True
            out[row["leg"]] = d
    return out


def main():
    exp = parse_experiment()
    only = [t for t in os.environ.get("PYFP3D_TIP_LEGS", "").split(",") if t]
    prior = load_prior()
    rows, stop = [], None
    for tag, h_wall, h_le, h_te in LEGS:
        if only and tag not in only:
            if tag in prior:
                print(f"\n=== {tag}: taken from the PRIOR run "
                      f"(LE upper {prior[tag]['band_LE_upper']:.4f}, "
                      f"labelled from_prior_run) ===")
                rows.append(prior[tag])
            else:
                print(f"\n=== {tag}: not selected and no prior row ===")
                rows.append(dict(leg=tag, note="not selected this run"))
            continue
        if stop:
            print(f"\n=== {tag}: SKIPPED -- {stop}")
            rows.append(dict(leg=tag, note=f"skipped: {stop}"))
            continue
        print(f"\n=== {tag}: h_wall={h_wall} h_le={h_le} h_te={h_te} "
              f"+ taper({TAPER_FORM}, {TAPER_RC_FRAC} b) + pressure ===",
              flush=True)
        cache = os.path.join(SCRATCH, f"{MESH_ALIAS.get(tag, tag)}.msh")
        t0 = time.perf_counter()
        try:
            if os.path.exists(cache):
                mesh = read_mesh(cache)
                print(f"  [cached mesh] {cache}", flush=True)
            else:
                mesh = build(h_wall, h_le, h_te)
                try:
                    write_mesh(mesh, cache)
                except Exception as exc:                          # noqa: BLE001
                    print(f"  (mesh not cached: {exc})")
        except Exception as exc:                                  # noqa: BLE001
            print(f"  MESH FAILED: {exc}")
            rows.append(dict(leg=tag, note=f"mesh failed: {exc}"))
            continue
        t_gen = time.perf_counter() - t0
        mc, wc = cut_wake(mesh)
        ar = compute_aspect_ratios(mc.nodes, mc.elements)
        dih = compute_min_dihedral_angles(mc.nodes, mc.elements)
        ht, hn, aniso = le_geometry(mc)
        print(f"  {len(mc.nodes)} nodes / {len(mc.elements)} tets (gen "
              f"{t_gen:.0f}s)  LE h_t {ht:.5f} h_n {hn:.5f} aniso {aniso:.3f} "
              f" max AR {ar.max():.1f} min dih {dih.min():.1f}", flush=True)

        t0 = time.perf_counter()
        try:
            r = solve_tapered(mc, wc)
        except Exception as exc:                                  # noqa: BLE001
            print(f"  SOLVE FAILED: {exc}")
            rows.append(dict(leg=tag, n_tets=len(mc.elements),
                             gen_s=round(t_gen, 1),
                             note=f"solve failed: {exc}"))
            continue
        wall = time.perf_counter() - t0
        phi = np.asarray(r["phi"])
        conv = bool(r.get("converged"))
        res = float(r.get("residual_history", [np.nan])[-1])
        print(f"  solve: conv={conv} |R|={res:.2e} n={r['n_newton']} "
              f"n_limited={r['n_limited']} n_floored={r['n_floored']} "
              f"({wall:.0f}s)", flush=True)

        # ---- G1: the tip-cure gate, read BEFORE any RMS is quoted ----------
        Bg, Vg = precompute_element_geometry(mc.nodes, mc.elements)
        gg = np.einsum("eaj,ea->ej", Bg, phi[mc.elements])
        m2 = mach_number_squared(np.einsum("ej,ej->e", gg, gg), M_INF)
        cent = mc.nodes[mc.elements].mean(axis=1)
        top8 = np.argsort(m2)[-8:][::-1]
        zc = np.clip(cent[top8, 2], 0.0, B_SEMI)
        xc8 = ((cent[top8, 0] - np.array([x_le(z) for z in zc]))
               / np.array([chord_at(z) for z in zc]))
        n_tip = int((cent[top8, 2] > 0.95 * B_SEMI).sum())
        m_max = float(np.sqrt(m2.max()))
        # ---- G1'a (20260801-1400): the check that actually matters is whether
        # a singularity sits INSIDE the measured region. M3 reads eta 0.20-0.90;
        # the tip singularity lives at eta > 0.95. So restrict M_max to eta<=0.92.
        inb = cent[:, 2] <= 0.92 * B_SEMI
        m_max_meas = float(np.sqrt(m2[inb].max())) if inb.any() else float("nan")
        print(f"  G1'a contamination: M_max inside the measured region "
              f"(eta<=0.92) = {m_max_meas:.4f}   [global {m_max:.4f}]",
              flush=True)
        # state saved so a later diagnostic never re-pays the solve (discipline #2)
        try:
            np.savez(os.path.join(SCRATCH, f"{tag}_state.npz"),
                     phi=phi, gamma=np.asarray(r["gamma"]))
        except Exception as exc:                                  # noqa: BLE001
            print(f"  (state not saved: {exc})")
        g1 = (m_max < 1.6) or (n_tip <= 2)
        print(f"  G1 tip cure: M_max {m_max:.4f} (L0p {L0P['m_max']}), "
              f"{n_tip}/8 in the tip band (L0p {L0P['top8_in_tip']}), peak at "
              f"x/c {xc8[0]:+.4f} eta {zc[0]/B_SEMI:.3f} -> "
              f"{'PASS' if g1 else 'FAIL'}", flush=True)

        curves = {e: section_cp_curve(mc, phi, eta=e, b_semi=B_SEMI,
                                      m_inf=M_INF) for e in ETAS}
        per = {e: station_rms(curves, exp, e)[0] for e in ETAS}
        acc = {}
        for e in ETAS[:N_UNMASKED]:
            for k, (ss, nn) in band_rms(curves, exp, e).items():
                a0, n0 = acc.get(k, (0.0, 0))
                acc[k] = (a0 + ss, n0 + nn)
        band = {k: (v[0] / max(v[1], 1)) ** 0.5 for k, v in acc.items()}
        pooled = float(np.mean([per[e] for e in ETAS[:N_UNMASKED]]))
        allss = sum(v[0] for v in acc.values())
        alln = sum(v[1] for v in acc.values())
        s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], phi,
                                    alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        o = np.argsort(wc.station_z)
        clkj = float(cl_kj_3d(np.atleast_1d(r["gamma"])[o], wc.station_z[o],
                              s_ref, B_SEMI))
        print(f"  ★ LE upper RMS {band['LE_upper']:.4f}   pooled {pooled:.4f}"
              f"   cl_p {f['cl']:.6f} (L0p {L0P['cl_p']}, "
              f"{100*(f['cl']/L0P['cl_p']-1):+.2f} %)   cl_KJ {clkj:.6f}",
              flush=True)
        for name, _, _ in BANDS:
            print(f"      {name:3s} upper {band[f'{name}_upper']:.4f}   "
                  f"lower {band[f'{name}_lower']:.4f}")

        rows.append(dict(
            leg=tag, h_wall=h_wall, h_le=h_le, h_te=h_te,
            n_nodes=len(mc.nodes), n_tets=len(mc.elements),
            gen_s=round(t_gen, 1), solve_s=round(wall, 1),
            le_h_t=round(ht, 6), le_h_n=round(hn, 6), le_aniso=round(aniso, 4),
            max_aspect=round(float(ar.max()), 2),
            min_dihedral=round(float(dih.min()), 3),
            converged=conv, res_final=res, n_newton=r["n_newton"],
            n_limited=r["n_limited"], n_floored=r["n_floored"],
            m_max=round(m_max, 5), top8_in_tip=n_tip,
            m_max_measured_region=round(m_max_meas, 5),
            peak_xc=round(float(xc8[0]), 5),
            peak_eta=round(float(zc[0] / B_SEMI), 5),
            g1_tip_cure="PASS" if g1 else "FAIL",
            pooled_rms_5=round(pooled, 6),
            allpoint_rms_5=round((allss / max(alln, 1)) ** 0.5, 6),
            **{f"band_{k}": round(v, 6) for k, v in band.items()},
            **{f"rms_eta{e:.2f}": round(per[e], 6) for e in ETAS},
            cl_p=round(f["cl"], 6), cl_kj=round(clkj, 6),
            sigma_min=r.get("sigma_min"), m1_max=r.get("m1_max"),
            over_budget=bool(wall > COST_GATE_S), note=""))
        if wall > COST_GATE_S:
            stop = (f"{tag} took {wall:.0f}s > the {COST_GATE_S:.0f}s gate")
            print(f"  ⚠ {stop}", flush=True)
        elif tag == "T1_allscales_taper" and wall > T2_GATE_S:
            stop = (f"T1 took {wall:.0f}s > the {T2_GATE_S:.0f}s T2 gate, so T2 "
                    f"is not run (pre-registered)")
            print(f"  ⚠ {stop}", flush=True)

    with open(os.path.join(OUT, "tip_allscales.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", os.path.join(OUT, "tip_allscales.csv"))

    got = {r["leg"]: r for r in rows if r.get("band_LE_upper") is not None
           and r.get("converged")}
    print("\n=== the registered reading (G1 -> G2) ===")
    t0r = got.get("T0_baseline_taper")
    if t0r is None:
        print("  T0 missing or non-converged: G1 cannot be read, and per the "
              "pre-registration no refinement reading is attributable. STOP.")
        return 0
    print(f"  G1: {t0r['g1_tip_cure']}  (M_max {t0r['m_max']} vs L0p "
          f"{L0P['m_max']}, tip-band {t0r['top8_in_tip']}/8 vs "
          f"{L0P['top8_in_tip']}/8)")
    if t0r["g1_tip_cure"] != "PASS":
        print("  => the taper does NOT cure this singularity on the M6 wing "
              "alone. Per the pre-registration the refinement legs are NOT "
              "readable as P11's lever; that is this round's conclusion.")
        return 0
    base = t0r["band_LE_upper"]
    print(f"  G2 baseline (SAME recipe, T0): LE upper {base:.4f}   "
          f"[step one's no-taper L0p was {L0P['le_upper']}]")
    seq = [got[t]["band_LE_upper"] for t, *_ in LEGS if t in got]
    for tag, *_ in LEGS:
        r = got.get(tag)
        if r is None:
            nc = next((x for x in rows if x["leg"] == tag), None)
            print(f"  {tag:22s} MISSING/non-converged "
                  f"({'' if nc is None else nc.get('note', '')})")
            continue
        v = r["band_LE_upper"]
        print(f"  {tag:22s} LE upper {v:.4f} ({100*(v/base-1):+.1f} % vs T0)  "
              f"{r['n_tets']} tets  {r['solve_s']:.0f}s  cl_p {r['cl_p']}")
    if len(seq) < 2:
        print("  only T0 is readable -- P11's lever is STILL unread. Recorded.")
        return 0
    best = min(seq)
    mono = all(b <= a + 1e-9 for a, b in zip(seq, seq[1:]))
    rel = abs(best - base) / base
    print(f"\n  best {best:.4f} vs T0 {base:.4f} ({100*(best/base-1):+.1f} %), "
          f"monotone: {mono}")
    # ---- P1 / P2: the pre-registered OUT-OF-SAMPLE predictions -----------
    tm1 = got.get("Tm1_coarse_taper")
    if tm1 is not None:
        p1 = tm1["band_LE_upper"] >= 0.26
        print(f"\n  P1 (coarse LE upper >= 0.26): "
              f"{tm1['band_LE_upper']:.4f} -> {'HOLDS' if p1 else 'FALSIFIED'}")
        tri = [got[t]["band_LE_upper"] for t in ("Tm1_coarse_taper",
                                                "T0_baseline_taper",
                                                "T1_allscales_taper")
               if t in got]
        p2 = len(tri) >= 2 and all(b <= a + 1e-9 for a, b in zip(tri, tri[1:]))
        print(f"  P2 (coarse -> T0 -> T1 monotone): "
              f"{' -> '.join(f'{v:.4f}' for v in tri)} -> "
              f"{'HOLDS' if p2 else 'FALSIFIED'}")
        if p1 and p2:
            print("  => the out-of-sample prediction HOLDS: all-scales "
                  "refinement behaves like a genuine convergence trend, so "
                  "step one's negative was the fixed-bulk pollution floor.")
        else:
            print("  => the out-of-sample prediction is FALSIFIED: all-scales "
                  "refinement is NOT the lever. Combined with step one, the "
                  "only remaining suspect is the intrinsic P1 capability and "
                  "D2's KILL CRITERION FIRES -- the semi-structured decision "
                  "stays the user's.")
    if best <= 0.15 and mono:
        print("  => the LE error IS resolution-controlled once the tip is "
              "cured: step one's negative was P11's fixed-bulk pollution floor. "
              "S2 has a route (all-scales refinement).")
    elif rel < 0.10:
        print("  => P11's lever fails on this wing too. The only remaining "
              "suspect is the intrinsic P1 capability at the suction peak, and "
              "D2's KILL CRITERION FIRES. S2 closes as a measured conclusion; "
              "the semi-structured decision is the user's, not this script's.")
    else:
        print("  => RECORDED (between the bands). D2 is not decided here.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
