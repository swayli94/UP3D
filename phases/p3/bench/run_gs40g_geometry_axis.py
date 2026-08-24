"""
Does k* survive a change of GEOMETRY?

Pre-registration: phases/p4/docs/dev_phase_four/20260817-2300-geometry-axis-prereg.md
Previous round:   phases/p4/docs/dev_phase_four/20260817-2100-condition-axis-verdict.md (W-DRIFT)

Wing-body -> M6 WING, at conditions MATCHED to the states already cached, so
geometry is the single variable. The criterion is not "is k* == 1" but "is k*
still equal to ITS MATCHED wing-body leg" -- which inherits last round's 3.6%
exception at M0.82 instead of pretending it away.

★ THE CONFOUND IS REGISTERED, NOT DISCOVERED: changing family changes the
production RECIPE too. A CTRL arm runs the wing at M0.88 under the WING-BODY
recipe to bound it.

★ Absolute C may NOT be compared across geometries (last round measured C
co-varying with shock strength). The registered cross-geometry quantity is the
RATIO C(P1)/C(P0).
"""

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "bench"))

import capability_matrix as cap                                 # noqa: E402
import run_gs40c_coherence_ladder as L                              # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import (build_face_adjacency,              # noqa: E402
                                 compute_edge_lengths,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

WING = REPO / "cases/meshes/onera_m6/coarse.msh"
WB = REPO / "cases/meshes/onera_m6_wingbody_conforming/coarse.msh"
NACA = REPO / "cases/meshes/naca0012_2.5d/coarse.msh"
CACHE = REPO / "bench/gate_results/gs40g_states"
WB_CACHE = REPO / "bench/gate_results/gs40f_states"
WB_ANCHOR = REPO / "bench/gate_results/gs40d_levels/coarse.npz"
OUT = REPO / "bench/gate_results/gs40g_geometry_axis.csv"
GAMMA, ALPHA = 1.4, 3.06

LADDER, ARMS = ("P0", "P1", "P2", "P3"), ("P0", "P1", "P2", "P3", "S1")
PEAK_FACTOR, FLAT_RATIO, MIN_CROSS = L.PEAK_FACTOR, L.FLAT_RATIO, L.MIN_CROSS
GATE_S = 25 * 60

#: matched pairs: wing M -> the wing-body cache holding the SAME condition
MATCH = {0.82: WB_CACHE / "M0.82.npz", 0.86: WB_CACHE / "M0.86.npz",
         0.88: WB_ANCHOR}
#: last round's published wing-body k*, per condition -- the comparison target
WB_KSTAR = {0.82: "P2", 0.86: "P1", 0.88: "P1"}
WB_C = {0.82: (0.114, 0.281, 0.291, 0.154), 0.86: (0.243, 0.457, 0.301, 0.222),
        0.88: (0.221, 0.479, 0.396, 0.291)}


def guard_recipe_m6():
    """G-RECIPE-M6: source-compare, the way run_m3_budget does it."""
    src = (REPO / "tests/test_p8_newton.py").read_text()
    need = ("dm=0.05, dm_min=0.01, freeze_tol=1e-6, intermediate_tol=1e-5",
            'precond="amg"', "freeze_refresh_max=8", "farfield_spanwise_gamma=True")
    miss = [n for n in need if n not in src]
    assert not miss, f"G-RECIPE-M6: NEWTON_M6_RECIPE drifted: {miss}"
    assert NEWTON_M6_RECIPE["freeze_tol"] == 1e-6, "G-RECIPE-M6: live value drifted"
    print(f"G-RECIPE-M6  NEWTON_M6_RECIPE source-compared ({len(need)} fragments)  PASS")


def solve_wing(tag, m_inf, wb_recipe=False):
    npz = CACHE / f"{tag}.npz"
    mc, wc = cut_wake(read_mesh(WING))
    if npz.exists():
        d = np.load(npz)
        return mc, np.asarray(d["phi"], float), bool(d["conv"]), int(d["nlim"]), \
            int(d["nflr"]), str(d["note"]), 0.0
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", 0.05 * B_SEMI)
    t0 = time.perf_counter()
    try:
        if wb_recipe:                     # ★ CTRL: the WING-BODY recipe on the wing
            seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART,
                                        alpha_deg=ALPHA, **cap.CONF_SEED_KW)
            nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                      phi_init=seed["phi"], gamma_init=seed["gamma"],
                      n_picard_seed=0, tip_taper=taper)
            r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=ALPHA,
                                       m_start=cap.WB_MSTART, dm=cap.DM,
                                       dm_min=0.01, freeze_tol=1e-5,
                                       intermediate_tol=1e-4, newton_kw=nk)
        else:                             # production M6 recipe
            kw = dict(NEWTON_M6_RECIPE)
            kw["newton_kw"] = dict(kw["newton_kw"], tip_taper=taper)
            r = solve_newton_transonic(mc, wc, m_inf=m_inf, alpha_deg=ALPHA, **kw)
    except Exception as exc:                                       # noqa: BLE001
        return mc, None, False, -1, -1, f"{type(exc).__name__}: {exc}", \
            time.perf_counter() - t0
    wall = time.perf_counter() - t0
    note = ""
    if not r["converged"]:
        h = np.asarray(r["residual_history"], float)
        mono = len(h) >= 3 and h[-1] < h[-2] < h[-3]
        note = (f"accept_reason={r.get('accept_reason')}; clamps="
                f"{int(r['n_limited'])}/{int(r['n_floored'])}; "
                f"m_last_converged={r.get('m_last_converged')}; "
                f"tail={h[-3:] if len(h) >= 3 else h}; "
                f"{'MONOTONE tail => budget_limited, NOT a ceiling' if mono else 'non-monotone tail'}")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz, phi=r["phi"], conv=bool(r["converged"]),
                        nlim=int(r["n_limited"]), nflr=int(r["n_floored"]), note=note,
                        mesh_sha=mesh_fingerprint(WING)["sha256"])
    return (mc, np.asarray(r["phi"], float), bool(r["converged"]),
            int(r["n_limited"]), int(r["n_floored"]), note, wall)


def ladder(mc, phi, m_inf):
    nodes, elements = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, elements)
    face_nb, _ = build_face_adjacency(elements)
    fa, fb = L.faces(face_nb)
    P, p1e = {}, None
    for k in (0, 1, 2, 3):
        po, pi = L.khop_patches(elements, face_nb, k)
        P[f"P{k}"] = L.patch_nodes(elements, po, pi)
        if k == 1:
            p1e = (po, pi)
    gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
    element_velocity_q2(elements, B, phi, gA, q2)
    m2A = mach_number_squared(q2, m_inf, GAMMA)
    g = {n: L.refit(nodes, elements, phi, *P[n]) for n in LADDER}
    g["S1"] = L.average_over_patch(nodes, elements, gA, *p1e)
    assert np.max(np.abs(g["P0"] - gA)) < 1e-10, "G-A"
    a = np.random.default_rng(0).normal(size=3)
    worst = max(float(np.max(np.abs(L.refit(nodes, elements, nodes @ a, *P[n]) - a)))
                for n in LADDER[1:])
    assert worst < 1e-11, f"G-CONS {worst:.2e}"
    C, res = {}, {}
    for n in ARMS:
        m2 = mach_number_squared((g[n] * g[n]).sum(axis=1), m_inf, GAMMA)
        r = L.coherence(np.sqrt(np.maximum(m2, 0.0)), fa, fb)
        C[n] = r["C"]; res[n] = r
    cent = nodes[elements].mean(axis=1)
    off, idx = P["P1"]
    smp = range(0, len(elements), max(1, len(elements) // 3000))
    rad = float(np.median([np.median(np.linalg.norm(
        nodes[idx[off[e]:off[e + 1]]] - cent[e], axis=1)) for e in smp]))
    h = float(np.median(compute_edge_lengths(nodes, elements)))
    return C, res, dict(sup_frac=float((m2A >= 1.0).mean()),
                        m_max=float(np.sqrt(m2A.max())), radius_over_h=rad / h,
                        n_tets=len(elements), worst_cons=worst)


def kstar(C):
    vals = [C[k] for k in LADDER]
    if max(vals) / max(min(vals), 1e-30) < FLAT_RATIO:
        return None, "NO-PEAK (flat)"
    if max(C["P1"], C["P2"]) < PEAK_FACTOR * max(C["P0"], C["P3"]):
        return None, "NO-PEAK (no interior max)"
    k = min(LADDER, key=lambda x: (-round(C[x], 6), LADDER.index(x)))
    return k, f"PEAK k*={k}"


def main():
    t0 = time.perf_counter()
    guard_recipe_m6()
    for p in (WING, NACA):
        man = read_manifest(p)
        if man:
            assert mesh_fingerprint(p)["sha256"] == man["sha256"], f"G-MESH {p.name}"
            print(f"G-MESH  {p.parent.name}/{p.name}  sha {man['sha256'][:12]}  PASS")
        else:
            print(f"G-MESH  {p.parent.name}/{p.name}  no manifest (tracked mesh)  noted")

    # ---- G-ANCHOR: the matched wing-body legs, from cache, zero cost ----------
    print()
    for m, npz in MATCH.items():
        d = np.load(npz)
        mcb, _ = cut_wake(read_mesh(WB))
        Cw, _r, _x = ladder(mcb, np.asarray(d["phi"], float), m)
        bad = {k: (round(Cw[k], 3), WB_C[m][i]) for i, k in enumerate(LADDER)
               if abs(Cw[k] - WB_C[m][i]) > 1e-3}
        ks, _ = kstar(Cw)
        print(f"G-ANCHOR wing-body M{m}: k*={ks} (published {WB_KSTAR[m]})  "
              f"{'PASS' if not bad and ks == WB_KSTAR[m] else '** FAIL ** ' + str(bad)}")
        assert not bad and ks == WB_KSTAR[m], "G-ANCHOR: kill criterion 2"

    rows, wing_k, ratios = [], {}, {}
    print()
    for m in (0.82, 0.86, 0.88):
        tag = f"W{m}"
        mc, phi, conv, nlim, nflr, note, wall = solve_wing(tag, m)
        if phi is None or not conv or nlim or nflr:
            print(f"[{tag}] NOT USABLE conv={conv} clamps={nlim}/{nflr}\n     {note}")
            rows.append(dict(leg=tag, geom="wing", m_inf=m, status="NOT_CONVERGED",
                             note=note, wall_s=round(wall, 1)))
            continue
        C, res, x = ladder(mc, phi, m)
        n0 = res["P0"]["n_cross"]
        counted = n0 >= MIN_CROSS
        ks, code = kstar(C)
        if counted and ks:
            wing_k[m] = ks
        ratios[tag] = (C["P1"] / max(C["P0"], 1e-30),
                       abs(C["P1"] - C["S1"]) / max(C["S1"], 1e-30))
        print(f"[{tag}] {wall:5.1f}s  {x['n_tets']} tets  sup_frac {x['sup_frac']:.4f} "
              f"M_max {x['m_max']:.3f}  n_cross {n0:>5}  radius/h {x['radius_over_h']:.3f}  "
              f"{'COUNTED' if counted else '★ EXCLUDED-SUBCRITICAL'}")
        print(f"       C = " + "  ".join(f"{n} {C[n]:.3f}" for n in ARMS)
              + f"   -> {code}   [wing-body was {WB_KSTAR[m]}]")
        for n in ARMS:
            rows.append(dict(leg=tag, geom="wing", m_inf=m, status="OK", arm=n,
                             C=round(C[n], 4), n_cross=res[n]["n_cross"],
                             max_comp=res[n]["max_comp"], n_ge50=res[n]["n_usable"],
                             sup_frac=round(x["sup_frac"], 5),
                             m_max=round(x["m_max"], 4),
                             radius_over_h=round(x["radius_over_h"], 3),
                             n_tets=x["n_tets"], wall_s=round(wall, 1)))

    # ---- CTRL: bound the recipe confound -------------------------------------
    print()
    mc, phi, conv, nlim, nflr, note, wall = solve_wing("CTRL", 0.88, wb_recipe=True)
    ctrl_k = None
    if phi is not None and conv and not nlim and not nflr:
        C, res, x = ladder(mc, phi, 0.88)
        ctrl_k, code = kstar(C)
        print(f"[CTRL] wing @ M0.88 with the WING-BODY recipe: {wall:.1f}s  "
              f"n_cross {res['P0']['n_cross']}  -> {code}")
        print("       C = " + "  ".join(f"{n} {C[n]:.3f}" for n in ARMS))
        for n in ARMS:
            rows.append(dict(leg="CTRL", geom="wing", m_inf=0.88, status="OK", arm=n,
                             C=round(C[n], 4), note="wing-body recipe"))
    else:
        print(f"[CTRL] NOT USABLE -- the recipe confound CANNOT be bounded\n     {note}")
        rows.append(dict(leg="CTRL", geom="wing", m_inf=0.88,
                         status="NOT_CONVERGED", note=note))

    # ---- verdict --------------------------------------------------------------
    print("\n" + "=" * 76)
    if len(wing_k) < 2:
        verdict = f"X-UNDEF (only {len(wing_k)} counting wing leg(s)) -- kill criterion 3"
    else:
        no_peak = [m for m in wing_k if wing_k[m] is None]
        mism = {m: (wing_k[m], WB_KSTAR[m]) for m in wing_k if wing_k[m] != WB_KSTAR[m]}
        if no_peak:
            verdict = f"X-VANISH (no interior maximum at {no_peak}) -- scope cut to wing-body"
        elif mism:
            verdict = f"X-SHIFT (k* differs from the matched wing-body leg at {mism})"
        else:
            verdict = "X-SAME (k* equals its matched wing-body leg on every counting leg)"
    for m in sorted(wing_k):
        print(f"  M{m}: wing k*={wing_k[m]}   wing-body k*={WB_KSTAR[m]}   "
              f"{'MATCH' if wing_k[m] == WB_KSTAR[m] else '★ DIFFER'}")
    print(f"  CTRL (wing-body recipe on the wing) k*={ctrl_k}  vs W0.88 "
          f"k*={wing_k.get(0.88)}  -> "
          + ("recipe confound BOUNDED" if ctrl_k and ctrl_k == wing_k.get(0.88)
             else "★ recipe confound NOT bounded"))
    print("=" * 76)
    print(f"VERDICT: {verdict}")
    print("\nRECORDED ratios (the registered cross-geometry quantity; wing-body was "
          "C(P1)/C(P0) 1.88-2.51, vs-S1 1.30-2.57):")
    for t, (r1, r2) in ratios.items():
        print(f"  {t}: C(P1)/C(P0) = {r1:.2f}   |C_P1-C_S1|/C_S1 = {r2:.2f}")
    print("\n★ Absolute C is NOT comparable across geometries; COHERENT IS NOT\n"
          "  CORRECT; and this says nothing about behaviour inside the solver.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict=verdict, wing_kstar={str(k): v for k, v in wing_k.items()},
             ctrl_kstar=ctrl_k, wall_s=round(time.perf_counter() - t0, 1)),
        indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 5"


if __name__ == "__main__":
    main()
