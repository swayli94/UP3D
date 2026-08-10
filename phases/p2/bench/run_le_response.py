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
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
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

OUT = os.path.join(_GATE)
SCRATCH = os.environ.get(
    "PYFP3D_SCRATCH",
    "/tmp/claude-1000/-home-lrz-codes-UP3D/3c5b43c4-b62c-4a09-b4da-9b9c7128d43e/scratchpad")
os.makedirs(OUT, exist_ok=True)

#: the committed medium's own parameters (cases/meshes/onera_m6/
#: generate_onera_m6.py::_level_params at h_wall = 0.015), so L0 reproduces it
H_WALL_MED = 0.015
BASE = dict(r_far=15.0 * MAC, tip_cap="flat", embed_wake=True)

#: (tag, h_wall, h_le, h_te). h_wake / h_far follow _level_params from h_wall.
#: ★ REVISED by addendum #1 (20260801-0400): the first version passed only h_edge,
#: which sizes the LE *and TE* curves, so "LE x2" refined the tip trailing edge with
#: it and read the P13 free-edge singularity instead (M_max 0.9557 -> 1.3945 at
#: M_inf 0.5, 7 of the top 8 cells in the tip band). These legs hold h_te at its
#: baseline so only the LE moves. L0p passes h_te explicitly and therefore takes the
#: SPLIT-field code path at the baseline size -- it is the control for the split
#: itself, not a re-run of L0.
LEGS = (
    ("L0p_split_control", H_WALL_MED, 0.5 * H_WALL_MED, 0.5 * H_WALL_MED),
    ("L1p_le_x2", H_WALL_MED, 0.25 * H_WALL_MED, 0.5 * H_WALL_MED),
    ("L2p_le_x4", H_WALL_MED, 0.125 * H_WALL_MED, 0.5 * H_WALL_MED),
    ("L3p_allscales", 0.010, 0.25 * H_WALL_MED, 0.5 * 0.010),
)
#: addendum #1 A4: a leg past this wall time stops the SWEEP (the finer legs can only
#: be worse). Deviation from A4's letter, recorded: the driver has no time-budget hook
#: to abort a leg mid-solve, so the bound is enforced between legs and by capping the
#: Newton steps below.
COST_GATE_S = 1200.0
#: also addendum #1 A4: L1's 2613 s was 60 capped Newton steps at several levels. 30
#: still gives >2x the 14 steps the converging baseline needs, and halves a failing
#: leg's cost. Recorded as a deviation from "the P14 recipe verbatim".
N_NEWTON_MAX = 30
#: the probe-leg baseline this sweep is measured against (20260801-0030)
LE_UPPER_BASE = 0.2361


def build(h_wall, h_edge, h_te=None):
    return onera_m6_wing_mesh(h_wall=h_wall, h_edge=h_edge, h_te=h_te,
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
    stop = None
    for tag, h_wall, h_le, h_te in LEGS:
        if stop:
            print(f"\n=== {tag}: SKIPPED -- {stop}")
            rows.append(dict(leg=tag, h_wall=h_wall, h_edge=h_le,
                             h_te=h_te, note=f"skipped: {stop}"))
            continue
        print(f"\n=== {tag}: h_wall={h_wall} h_le={h_le} "
              f"h_te={h_te} ===", flush=True)
        t0 = time.perf_counter()
        # cache probe meshes: L2p alone costs 133 s to generate and a re-run
        # should not pay it again. Scratch, never cases/meshes/.
        cache = os.path.join(SCRATCH, f"{tag}.msh")
        try:
            if os.path.exists(cache):
                from pyfp3d.mesh.reader import read_mesh
                mesh = read_mesh(cache)
                print(f"  [cached mesh] {cache}", flush=True)
            else:
                mesh = build(h_wall, h_le, h_te)
                try:
                    from pyfp3d.mesh.reader import write_mesh
                    write_mesh(mesh, cache)
                except Exception as exc:                          # noqa: BLE001
                    print(f"  (mesh not cached: {exc})")
        except Exception as exc:                                  # noqa: BLE001
            print(f"  MESH GENERATION FAILED: {exc}")
            rows.append(dict(leg=tag, h_wall=h_wall, h_edge=h_le, h_te=h_te,
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
            r = solve(mc, wc, True, n_newton_max=N_NEWTON_MAX)
        except Exception as exc:                                  # noqa: BLE001
            print(f"  SOLVE FAILED: {exc}")
            rows.append(dict(leg=tag, h_wall=h_wall, h_edge=h_le, h_te=h_te,
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

        # --- addendum #1 A3: M_max AND where the peak sits, every leg ------
        from pyfp3d.mesh.metrics import precompute_element_geometry
        from pyfp3d.meshgen.wing3d import chord_at, x_le
        from pyfp3d.physics.isentropic import mach_number_squared
        Bg, Vg = precompute_element_geometry(mc.nodes, mc.elements)
        gg = np.einsum("eaj,ea->ej", Bg, phi[mc.elements])
        m2 = mach_number_squared(np.einsum("ej,ej->e", gg, gg), M_INF)
        cent = mc.nodes[mc.elements].mean(axis=1)
        top8 = np.argsort(m2)[-8:][::-1]
        zc = np.clip(cent[top8, 2], 0.0, B_SEMI)
        xc8 = (cent[top8, 0] - np.array([x_le(z) for z in zc])) / np.array(
            [chord_at(z) for z in zc])
        n_le = int((np.abs(xc8) < 0.05).sum())
        n_tip = int((cent[top8, 2] > 0.95 * B_SEMI).sum())
        print(f"  A3 peak: M_max {np.sqrt(m2.max()):.4f} at x/c "
              f"{xc8[0]:+.4f} eta {zc[0]/B_SEMI:.3f};  of the top 8: "
              f"{n_le} near the LE, {n_tip} in the tip band", flush=True)

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
            leg=tag, h_wall=h_wall, h_edge=h_le, h_te=h_te,
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
            m_max=round(float(np.sqrt(m2.max())), 5),
            peak_xc=round(float(xc8[0]), 5),
            peak_eta=round(float(zc[0] / B_SEMI), 5),
            top8_near_le=n_le, top8_in_tip=n_tip,
            cl_p=round(f["cl"], 6), cl_kj=round(clkj, 6),
            over_budget=bool(wall > COST_GATE_S), note=""))
        if wall > COST_GATE_S:
            stop = (f"{tag} took {wall:.0f}s > the {COST_GATE_S:.0f}s cost gate "
                    f"(addendum #1 A4); finer legs can only be worse")
            print(f"  ⚠ {stop}", flush=True)

    with open(os.path.join(OUT, "le_response.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", os.path.join(OUT, "le_response.csv"))

    print("\n=== the registered reading ===")
    got = {r["leg"]: r for r in rows if r.get("band_LE_upper") is not None}
    for tag, *_ in LEGS:
        r = got.get(tag)
        if r is None:
            print(f"  {tag}: MISSING (recorded, never substituted)")
            continue
        v = r["band_LE_upper"]
        print(f"  {tag:14s} LE upper {v:.4f} "
              f"({100*(v-LE_UPPER_BASE)/LE_UPPER_BASE:+.1f} %)  "
              f"h_t {r['le_h_t']:.5f} aniso {r['le_aniso']:.3f}  "
              f"{r['n_tets']} tets  {r['solve_s']:.0f}s")
    vals = [got[t]["band_LE_upper"] for t, *_ in LEGS if t in got]
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
    l2, l3 = got.get("L2p_le_x4"), got.get("L3p_allscales")
    if l2 and l3:
        print(f"  floor check: L2 {l2['band_LE_upper']:.4f} vs L3 "
              f"{l3['band_LE_upper']:.4f} "
              f"(L3 better by {l2['band_LE_upper']-l3['band_LE_upper']:+.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
