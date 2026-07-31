"""GS2.1 step one: is the M6 LE Cp error resolution-controlled at all?

Pre-registered in docs/dev_phase_two/20260801-0200-gs21-le-response-prereg.md and
committed before this file was written -- read that for the four legs and the criteria.
The short version: the expensive part of GS2.1 (true 3-D anisotropic prism layers) is
machinery to build, so first sweep the generator's EXISTING single-variable LE knob
(h_edge) and see whether the LE band's Cp error responds. L3 refines h_wall as well,
because P11 measured that a single-variable refinement can hit the bulk pollution
floor and a null result would otherwise be unattributable.

Meshes are written to a scratch directory, not into cases/meshes/ -- they are probe
meshes for one question, not a committed family.

Outputs (TRACKED): bench/gate_results/le_response.csv
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

from pyfp3d.mesh.metrics import (compute_aspect_ratios,             # noqa: E402
                                 compute_min_dihedral_angles)
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, onera_m6_wing_mesh   # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,           # noqa: E402
                                 wall_force_coefficients)
from run_m3_budget import (BANDS, ETAS, M_INF, ALPHA, N_UNMASKED,    # noqa: E402
                           band_rms, parse_experiment, solve,
                           station_rms)

OUT = os.path.join(HERE, "gate_results")
os.makedirs(OUT, exist_ok=True)

#: the committed medium's own parameters (cases/meshes/onera_m6/
#: generate_onera_m6.py::_level_params at h_wall = 0.015), so L0 reproduces it
H_WALL_MED = 0.015
BASE = dict(r_far=15.0 * MAC, tip_cap="flat", embed_wake=True)

#: (tag, h_wall, h_edge). h_wake / h_far follow _level_params from h_wall.
LEGS = (
    ("L0_baseline", H_WALL_MED, 0.5 * H_WALL_MED),
    ("L1_le_x2", H_WALL_MED, 0.25 * H_WALL_MED),
    ("L2_le_x4", H_WALL_MED, 0.125 * H_WALL_MED),
    ("L3_allscales", 0.010, 0.25 * H_WALL_MED),
)
#: the probe-leg baseline this sweep is measured against (20260801-0030)
LE_UPPER_BASE = 0.2361


def build(h_wall, h_edge):
    return onera_m6_wing_mesh(h_wall=h_wall, h_edge=h_edge,
                              h_wake=3.0 * h_wall, h_far=120.0 * h_wall,
                              **BASE)


def le_geometry(mc):
    """(h_t, h_n, anisotropy) medians in the LE band -- the 20260801-0030 read,
    repeated per leg so the sweep shows what it actually changed."""
    from pyfp3d.meshgen.wing3d import chord_at, x_le
    key2tet = {}
    for e, tet in enumerate(mc.elements):
        for f in ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)):
            key2tet[tuple(sorted(tet[list(f)]))] = e
    hn, ht = [], []
    for face in mc.boundary_faces["wall"]:
        e = key2tet.get(tuple(sorted(face)))
        if e is None:
            continue
        opp = [n for n in mc.elements[e] if n not in set(face)]
        if len(opp) != 1:
            continue
        P = mc.nodes[face]
        nv = np.cross(P[1] - P[0], P[2] - P[0])
        A2 = float(np.linalg.norm(nv))
        if A2 == 0.0:
            continue
        c = P.mean(axis=0)
        xc = (c[0] - x_le(c[2])) / chord_at(c[2])
        if not (0.0 <= xc < 0.15):
            continue
        hn.append(abs(float(np.dot(mc.nodes[opp[0]] - P[0], nv / A2))))
        ht.append(np.sqrt(A2 / 2.0))
    hn, ht = np.array(hn), np.array(ht)
    if not len(hn):
        return (np.nan,) * 3
    return (float(np.median(ht)), float(np.median(hn)),
            float(np.median(ht / hn)))


def main():
    exp = parse_experiment()
    rows = []
    for tag, h_wall, h_edge in LEGS:
        print(f"\n=== {tag}: h_wall={h_wall} h_edge={h_edge} ===", flush=True)
        t0 = time.perf_counter()
        try:
            mesh = build(h_wall, h_edge)
        except Exception as exc:                                  # noqa: BLE001
            print(f"  MESH GENERATION FAILED: {exc}")
            rows.append(dict(leg=tag, h_wall=h_wall, h_edge=h_edge,
                             note=f"mesh gen failed: {exc}"))
            continue
        t_gen = time.perf_counter() - t0
        mc, wc = cut_wake(mesh)
        ar = compute_aspect_ratios(mc.nodes, mc.elements)
        dih = compute_min_dihedral_angles(mc.nodes, mc.elements)
        ht, hn, aniso = le_geometry(mc)
        print(f"  {len(mc.nodes)} nodes / {len(mc.elements)} tets  "
              f"(gen {t_gen:.0f}s)   LE band h_t {ht:.5f} h_n {hn:.5f} "
              f"aniso {aniso:.3f}   max AR {ar.max():.1f} min dih "
              f"{dih.min():.1f} deg", flush=True)

        t0 = time.perf_counter()
        try:
            r = solve(mc, wc, True)                       # entropy at default ON
        except Exception as exc:                                  # noqa: BLE001
            print(f"  SOLVE FAILED: {exc}")
            rows.append(dict(leg=tag, h_wall=h_wall, h_edge=h_edge,
                             n_nodes=len(mc.nodes), n_tets=len(mc.elements),
                             gen_s=round(t_gen, 1),
                             note=f"solve failed: {exc}"))
            continue
        wall = time.perf_counter() - t0
        phi = np.asarray(r["phi"])
        res = float(r.get("residual_history", [np.nan])[-1])
        conv = bool(r.get("converged"))
        print(f"  solve: conv={conv} |R|={res:.2e} n_newton={r['n_newton']} "
              f"n_limited={r['n_limited']} n_floored={r['n_floored']} "
              f"({wall:.0f}s)", flush=True)

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
        d = band["LE_upper"] - LE_UPPER_BASE
        print(f"  ★ LE upper RMS {band['LE_upper']:.4f}  "
              f"(baseline {LE_UPPER_BASE}: {d:+.4f}, "
              f"{100*d/LE_UPPER_BASE:+.1f} %)   pooled {pooled:.4f}   "
              f"cl_p {f['cl']:.6f}", flush=True)
        for name, _, _ in BANDS:
            print(f"      {name:3s} upper {band[f'{name}_upper']:.4f}   "
                  f"lower {band[f'{name}_lower']:.4f}")

        rows.append(dict(
            leg=tag, h_wall=h_wall, h_edge=h_edge,
            n_nodes=len(mc.nodes), n_tets=len(mc.elements),
            gen_s=round(t_gen, 1), solve_s=round(wall, 1),
            le_h_t=round(ht, 6), le_h_n=round(hn, 6),
            le_aniso=round(aniso, 4),
            max_aspect=round(float(ar.max()), 2),
            min_dihedral=round(float(dih.min()), 3),
            converged=conv, res_final=res, n_newton=r["n_newton"],
            n_limited=r["n_limited"], n_floored=r["n_floored"],
            sigma_min=r.get("sigma_min"), m1_max=r.get("m1_max"),
            pooled_rms_5=round(pooled, 6),
            allpoint_rms_5=round((allss / max(alln, 1)) ** 0.5, 6),
            **{f"band_{k}": round(v, 6) for k, v in band.items()},
            **{f"rms_eta{e:.2f}": round(per[e], 6) for e in ETAS},
            xshock_eta0p44=shock_report(curves[0.44], M_INF)["upper"].get(
                "x_shock"),
            cl_p=round(f["cl"], 6), cl_kj=round(clkj, 6), note=""))

    with open(os.path.join(OUT, "le_response.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", os.path.join(OUT, "le_response.csv"))

    print("\n=== the registered reading ===")
    got = {r["leg"]: r for r in rows if r.get("band_LE_upper") is not None}
    for tag, _, _ in LEGS:
        r = got.get(tag)
        if r is None:
            print(f"  {tag}: MISSING (recorded, never substituted)")
            continue
        v = r["band_LE_upper"]
        print(f"  {tag:14s} LE upper {v:.4f} "
              f"({100*(v-LE_UPPER_BASE)/LE_UPPER_BASE:+.1f} %)  "
              f"h_t {r['le_h_t']:.5f} aniso {r['le_aniso']:.3f}  "
              f"{r['n_tets']} tets  {r['solve_s']:.0f}s")
    vals = [got[t]["band_LE_upper"] for t, _, _ in LEGS if t in got]
    if len(vals) < 2:
        print("  too few legs to read")
        return 0
    best = min(vals)
    rel = abs(best - LE_UPPER_BASE) / LE_UPPER_BASE
    mono = all(b <= a + 1e-9 for a, b in zip(vals, vals[1:]))
    print(f"\n  best LE upper {best:.4f} vs baseline {LE_UPPER_BASE} "
          f"({100*(best-LE_UPPER_BASE)/LE_UPPER_BASE:+.1f} %), "
          f"monotone across the legs as ordered: {mono}")
    if best <= 0.15 and mono:
        print("  => LE error IS resolution-controlled: the anisotropic-layer "
              "route is worth building")
    elif rel < 0.10:
        print("  => LE error is NOT controlled by LE resolution. MAJOR "
              "NEGATIVE: prisms will not help either; the remaining suspects "
              "are the intrinsic P1 capability (P11) and the geometry, and "
              "D2's kill criterion comes up.")
    else:
        print("  => RECORDED (between the bands). Compare L3 against L2: if L3 "
              "is clearly better, it is the bulk pollution floor and the "
              "conclusion is that ALL SCALES must be refined, not the wall "
              "alone.")
    l2, l3 = got.get("L2_le_x4"), got.get("L3_allscales")
    if l2 and l3:
        print(f"  floor check: L2 {l2['band_LE_upper']:.4f} vs L3 "
              f"{l3['band_LE_upper']:.4f} "
              f"(L3 better by {l2['band_LE_upper']-l3['band_LE_upper']:+.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
