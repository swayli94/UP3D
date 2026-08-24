"""
Is the k*=P2 exception a STRENGTH effect or a GEOMETRY effect? -- and the 2.5-D
arm that was registered last round and never run.

Pre-registration: phases/p4/docs/dev_phase_four/20260818-0300-strength-and-2p5d-prereg.md
Previous round:   phases/p4/docs/dev_phase_four/20260818-0100-geometry-axis-verdict.md (X-SHIFT)

★ NEW PRIMARY QUANTITY, and it is a new QUESTION, not a relaxed old one: five
rounds scored the discrete k* = argmax C, and the whole dispute is one leg at a
3.6% margin. Here the primary is the CONTINUOUS margin R21 = C(P2)/C(P1), which
carries the same sign and also how much. k* is still printed for every leg and no
earlier verdict is re-judged.

★ Part A (N-2.5D) is RECORDED-ONLY -- it changes geometry, dimensionality,
condition and recipe at once, and its face-adjacency topology differs, so the
meaning of k is not fully comparable. G-2P5D prints that distribution so the
caveat carries numbers.
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

import run_capability_matrix as cap                                 # noqa: E402
import run_gs40c_coherence_ladder as L                              # noqa: E402
import run_gs40g_geometry_axis as G                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import build_face_adjacency                # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                # noqa: E402
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

WING = REPO / "cases/meshes/onera_m6/coarse.msh"
WB = REPO / "cases/meshes/onera_m6_wingbody_conforming/coarse.msh"
NACA = REPO / "cases/meshes/naca0012_2.5d"
CACHE = REPO / "bench/gate_results/gs40h_states"
OUT = REPO / "bench/gate_results/gs40h_strength_2p5d.csv"
LADDER = ("P0", "P1", "P2", "P3")

#: NEW wing legs, fixed in advance and NOT re-tuned to land on a chosen sup_frac
NEW_WING_M = (0.76, 0.78, 0.80)
#: existing legs: (geom, label, cache, m_inf) -- all already published
EXISTING = (
    ("wingbody", "WB M0.78", REPO / "bench/gate_results/gs40f_states/M0.78.npz", 0.78),
    ("wingbody", "WB M0.82", REPO / "bench/gate_results/gs40f_states/M0.82.npz", 0.82),
    ("wingbody", "WB M0.86", REPO / "bench/gate_results/gs40f_states/M0.86.npz", 0.86),
    ("wingbody", "WB M0.88", REPO / "bench/gate_results/gs40d_levels/coarse.npz", 0.88),
    ("wingbody", "WB a1.53", REPO / "bench/gate_results/gs40f_states/A1.53.npz", 0.88),
    ("wing", "W M0.82", REPO / "bench/gate_results/gs40g_states/W0.82.npz", 0.82),
    ("wing", "W M0.86", REPO / "bench/gate_results/gs40g_states/W0.86.npz", 0.86),
    ("wing", "W M0.88", REPO / "bench/gate_results/gs40g_states/W0.88.npz", 0.88),
)
#: G-ANCHOR: published C(P0..P3) for every existing leg
PUB = {"WB M0.78": (0.057, 0.143, 0.134, 0.030),
       "WB M0.82": (0.114, 0.281, 0.291, 0.154),
       "WB M0.86": (0.243, 0.457, 0.301, 0.222),
       "WB M0.88": (0.221, 0.479, 0.396, 0.291),
       "WB a1.53": (0.208, 0.419, 0.309, 0.247),
       "W M0.82": (0.102, 0.301, 0.224, 0.169),
       "W M0.86": (0.220, 0.504, 0.349, 0.257),
       "W M0.88": (0.206, 0.482, 0.369, 0.245)}
NACA_KW = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
               precond="direct", direct_refactor_every=4, n_newton_max=80,
               n_picard_seed=0)
GATE_S = 25 * 60


def guard_recipe_m1():
    src = (REPO / "bench/run_m1_gate.py").read_text().replace(" ", "")
    need = ("m_crit=0.95", "freeze_tol=1e-6", "freeze_refresh_max=8",
            'precond="direct"', "direct_refactor_every=4", "n_newton_max=80",
            "M_INF,ALPHA=0.80,1.25")
    miss = [n for n in need if n.replace(" ", "") not in src]
    assert not miss, f"G-RECIPE-M1: run_m1_gate drifted: {miss}"
    print(f"G-RECIPE-M1  {len(need)} solve parameters source-compared  PASS")


def wing_leg(m):
    tag = f"W{m}"
    npz = CACHE / f"{tag}.npz"
    mc, wc = cut_wake(read_mesh(WING))
    if npz.exists():
        d = np.load(npz)
        return mc, np.asarray(d["phi"], float), bool(d["conv"]), int(d["nlim"]), \
            int(d["nflr"]), 0.0
    taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", 0.05 * B_SEMI)
    kw = dict(NEWTON_M6_RECIPE)
    kw["newton_kw"] = dict(kw["newton_kw"], tip_taper=taper)
    from pyfp3d.solve.newton import solve_newton_transonic
    t0 = time.perf_counter()
    r = solve_newton_transonic(mc, wc, m_inf=m, alpha_deg=3.06, **kw)
    wall = time.perf_counter() - t0
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz, phi=r["phi"], conv=bool(r["converged"]),
                        nlim=int(r["n_limited"]), nflr=int(r["n_floored"]),
                        mesh_sha=mesh_fingerprint(WING)["sha256"])
    return (mc, np.asarray(r["phi"], float), bool(r["converged"]),
            int(r["n_limited"]), int(r["n_floored"]), wall)


def main():
    t0 = time.perf_counter()
    G.guard_recipe_m6(); guard_recipe_m1()
    for p in (WING, WB, NACA / "coarse.msh"):
        man = read_manifest(p)
        if man:
            assert mesh_fingerprint(p)["sha256"] == man["sha256"], f"G-MESH {p}"
    print("G-MESH  all meshes match their manifests (where stamped)  PASS\n")

    rows, legs = [], []

    # ---- existing legs: G-ANCHOR, zero cost --------------------------------
    for geom, label, npz, m in EXISTING:
        msh = WB if geom == "wingbody" else WING
        mc, _ = cut_wake(read_mesh(msh))
        d = np.load(npz)
        C, res, x = G.ladder(mc, np.asarray(d["phi"], float), m)
        bad = {k: (round(C[k], 3), PUB[label][i]) for i, k in enumerate(LADDER)
               if abs(C[k] - PUB[label][i]) > 1e-3}
        assert not bad, f"G-ANCHOR {label}: {bad} -- kill criterion 2"
        legs.append(dict(geom=geom, label=label, m_inf=m, C=C,
                         sup=x["sup_frac"], n_cross=res["P0"]["n_cross"]))
    print(f"G-ANCHOR  all {len(EXISTING)} existing legs reproduce their published "
          f"C(P0..P3) to +-0.001  PASS\n")

    # ---- Part B: new wing legs ---------------------------------------------
    for m in NEW_WING_M:
        mc, phi, conv, nlim, nflr, wall = wing_leg(m)
        if not conv or nlim or nflr:
            print(f"[W{m}] NOT USABLE conv={conv} clamps={nlim}/{nflr}")
            rows.append(dict(geom="wing", label=f"W M{m}", status="NOT_CONVERGED"))
            continue
        C, res, x = G.ladder(mc, phi, m)
        n0 = res["P0"]["n_cross"]
        if n0 < L.MIN_CROSS:
            print(f"[W{m}] {wall:5.1f}s  sup_frac {x['sup_frac']:.4f}  n_cross {n0}"
                  f"  ★ EXCLUDED-SUBCRITICAL")
            rows.append(dict(geom="wing", label=f"W M{m}", status="EXCLUDED",
                             sup_frac=round(x["sup_frac"], 5), n_cross=n0))
            continue
        legs.append(dict(geom="wing", label=f"W M{m}", m_inf=m, C=C,
                         sup=x["sup_frac"], n_cross=n0))
        print(f"[W{m}] {wall:5.1f}s  sup_frac {x['sup_frac']:.4f}  n_cross {n0:>5}  "
              f"C = " + "  ".join(f"{k} {C[k]:.3f}" for k in LADDER))

    # ---- table: R21 per leg -------------------------------------------------
    print(f"\n{'leg':12}{'geom':10}{'sup_frac':>10}{'C(P1)':>8}{'C(P2)':>8}"
          f"{'R21':>8}{'k*':>5}")
    for lg in sorted(legs, key=lambda z: z["sup"]):
        C = lg["C"]; r21 = C["P2"] / max(C["P1"], 1e-30)
        lg["R21"] = r21
        lg["kstar"] = "P2" if r21 > 1.0 else G.kstar(C)[0]
        print(f"{lg['label']:12}{lg['geom']:10}{lg['sup']:10.4f}{C['P1']:8.3f}"
              f"{C['P2']:8.3f}{r21:8.3f}{str(lg['kstar']):>5}")
        rows.append(dict(geom=lg["geom"], label=lg["label"], status="OK",
                         m_inf=lg["m_inf"], sup_frac=round(lg["sup"], 5),
                         n_cross=lg["n_cross"], R21=round(r21, 4),
                         kstar=str(lg["kstar"]),
                         **{f"C_{k}": round(lg["C"][k], 4) for k in LADDER}))

    # ---- primary criterion: strength vs geometry ---------------------------
    print("\n" + "=" * 76)
    byg = {g: [l for l in legs if l["geom"] == g] for g in ("wingbody", "wing")}
    fits, ov = {}, None
    for g, ls in byg.items():
        if len(ls) < 4:
            print(f"{g}: only {len(ls)} legs -- fit not attempted")
            continue
        x = np.log10([l["sup"] for l in ls]); y = np.array([l["R21"] for l in ls])
        p = np.polyfit(x, y, 1)
        resid = float(np.sqrt(np.mean((y - np.polyval(p, x)) ** 2)))
        fits[g] = dict(p=p, resid=resid, lo=x.min(), hi=x.max(), n=len(ls))
        print(f"{g:10} n={len(ls)}  R21 = {p[0]:+.3f}*log10(sup) {p[1]:+.3f}   "
              f"resid RMS {resid:.4f}   log10(sup) in [{x.min():.2f}, {x.max():.2f}]")
    if len(fits) < 2:
        verdict = "S-UNDEF (a geometry has fewer than 4 counting legs)"
    else:
        lo = max(f["lo"] for f in fits.values()); hi = min(f["hi"] for f in fits.values())
        if lo >= hi:
            verdict = "S-UNDEF (the strength ranges do not overlap)"
        else:
            mid = 0.5 * (lo + hi)
            off = abs(np.polyval(fits["wing"]["p"], mid)
                      - np.polyval(fits["wingbody"]["p"], mid))
            pooled = float(np.sqrt(np.mean([f["resid"] ** 2 for f in fits.values()])))
            ov = (lo, hi, mid, off, pooled)
            print(f"overlap log10(sup) [{lo:.2f}, {hi:.2f}], midpoint {mid:.2f}: "
                  f"|dR21| = {off:.4f}   pooled resid RMS = {pooled:.4f}   "
                  f"ratio {off / max(pooled, 1e-30):.2f}")
            if off <= pooled:
                verdict = "S-STRENGTH (geometry adds nothing beyond shock strength)"
            elif off >= 2 * pooled:
                verdict = "S-GEOM (a systematic geometry offset survives)"
            else:
                verdict = "S-MIX (no direction claimed)"
    print("=" * 76)
    print(f"VERDICT: {verdict}")

    p2 = [l["label"] for l in legs if l["R21"] > 1.0]
    print(f"\nD-P2: legs with R21 > 1 (k*=P2): {p2 if p2 else 'NONE'}")
    if p2 == ["WB M0.82"]:
        print("      => P2-ISOLATED: the exception is unique across "
              f"{len(legs)} legs / 2 geometries / 3 levels.")
        print("      ★ This does NOT overturn X-SHIFT; it is the context that must "
              "be quoted with it.")

    # ---- Part A: N-2.5D, RECORDED-ONLY -------------------------------------
    print("\n" + "-" * 76)
    print("Part A -- N-2.5D, RECORDED-ONLY (geometry + dimensionality + condition"
          " + recipe all change)")
    for level in ("coarse", "medium"):
        p = NACA / f"{level}.msh"
        if not p.exists():
            print(f"  {level}: mesh absent"); continue
        mc, wc = cut_wake(read_mesh(p))
        fn, _ = build_face_adjacency(mc.elements)
        nnb = (fn >= 0).sum(axis=1)
        t1 = time.perf_counter()
        r = solve_newton_lifting(mc, wc, m_inf=0.80, alpha_deg=1.25, **NACA_KW)
        w = time.perf_counter() - t1
        ok = bool(r["converged"]) and not r["n_limited"] and not r["n_floored"]
        print(f"  [{level}] {w:5.1f}s conv={r['converged']} "
              f"clamps={r['n_limited']}/{r['n_floored']}  "
              f"G-2P5D face-neighbours per tet: mean {nnb.mean():.2f} "
              f"(3-D wing mean ~3.7)")
        if not ok:
            rows.append(dict(geom="naca2.5d", label=f"NACA {level}",
                             status="NOT_CONVERGED")); continue
        C, res, x = G.ladder(mc, np.asarray(r["phi"], float), 0.80)
        n0 = res["P0"]["n_cross"]
        r21 = C["P2"] / max(C["P1"], 1e-30)
        ks = "P2" if r21 > 1 else G.kstar(C)[0]
        print(f"          sup_frac {x['sup_frac']:.4f} n_cross {n0} radius/h "
              f"{x['radius_over_h']:.3f}  C = "
              + "  ".join(f"{k} {C[k]:.3f}" for k in LADDER)
              + f"  R21 {r21:.3f}  k*={ks}"
              + ("" if n0 >= L.MIN_CROSS else "  ★ THIN"))
        rows.append(dict(geom="naca2.5d", label=f"NACA {level}", status="RECORDED",
                         m_inf=0.80, sup_frac=round(x["sup_frac"], 5), n_cross=n0,
                         R21=round(r21, 4), kstar=str(ks),
                         face_nb_mean=round(float(nnb.mean()), 3),
                         **{f"C_{k}": round(C[k], 4) for k in LADDER}))
    print("  ★ RECORDED-ONLY: M0.80/alpha1.25 is a known multi-solution condition\n"
          "    (one draw), and 2.5-D is one prism layer cut into tets, so the\n"
          "    meaning of k is not fully comparable to the 3-D legs.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w_.writeheader(); w_.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict=verdict, p2_legs=p2, overlap=None if ov is None else
             dict(lo=ov[0], hi=ov[1], offset=ov[3], pooled=ov[4]),
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 4"


if __name__ == "__main__":
    main()
