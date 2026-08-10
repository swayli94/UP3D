"""Registered item 0: a truly far-field-only arm, via the domain RADIUS `r_far`.

Pre-registered in docs/dev_phase_two/20260809-2000-farfield-arm-prereg.md, committed
before this file was written. Purpose is CLOSE-OUT -- settle whether the far field
itself controls the LE band so the capability boundary can state it -- not capability.

Why it is needed: the third LE factorial's criterion P was met, but the same-day erratum
found `h_far` is not a far-field knob. It moved the median cell size at r = 1-2 MAC by
+24.04 % and the LE band's normal spacing by 1.525 %, while the guard only looked at the
surface. So "does the FAR FIELD ITSELF control the LE band" has never been measured
cleanly, and `r_far` is the one knob that touches no surface geometry and no size field.

B5's prior (2026-07-12): the same far-field option (Dirichlet+vortex) is DOMAIN-ROBUST
over R in {15,30,60,120} c, Gamma within 0.45 %/1.09 %. So the expectation is "the LE
band barely moves" -- a prediction, which is what makes this a test. B5 measured Gamma
(an integral) in 2.5-D though, so it is a prior and not the answer.

★ The guard is what the previous erratum bought: surface AND near-body volume, because
only both together mean "far-field-only". If either fails, the registered reading is
branch N -- no clean far-field arm exists with the available knobs, a bounded negative
for the boundary document, with generator work handed to phase three.

★★ Order of operations follows from that: the guards are MESH-ONLY, and branch N needs
no Cp at all. So both meshes are built and judged FIRST, and the two solves run only if
the arm is clean. A dirty arm costs minutes instead of the leg budget, and no solve is
performed whose reading the registration would refuse to use.

Outputs (TRACKED): bench/gate_results/farfield_arm.csv
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
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, onera_m6_wing_mesh   # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402
from run_le_factorial import le_face_count                          # noqa: E402
from run_le_factorial2 import band_disp, band_err, cp_at_exp_points  # noqa: E402
from run_le_g1_volume import SHELLS, cell_sizes                     # noqa: E402
from run_le_response import le_geometry                             # noqa: E402
from run_le_window import H_LE, H_WALL, solve_at                    # noqa: E402
from run_m3_budget import parse_experiment                          # noqa: E402
from run_seed_exposure import clamp_map                             # noqa: E402

CSV = os.path.join(HERE, "gate_results", "farfield_arm.csv")
M_INF = 0.75                       # the third factorial's condition, verbatim
H_FAR = 2.4                        # ★ FIXED on both legs -- r_far is the only variable
LEGS = (("R00_base", 15.0), ("R10_domain_x2", 30.0))       # r_far, in MAC
GEN_GATE_S, LEG_GATE_S = 300.0, 600.0
NEAR_PCT = 1.0                     # G1a / G1b threshold in %, fixed in the registration
#: the unclean `h_far` arm's LE displacement, the reference the criteria are stated against
HFAR_ARM_DELTA_LE = 0.09811


def build(r_far_mac):
    """run_le_factorial2._build's baseline leg with r_far as the parameter.

    Everything else is that leg verbatim: h_wall/h_le from run_le_window, h_te = h_wall/2,
    h_wake = 3 h_wall, h_far = 2.4, round tip, wake embedded.
    """
    return onera_m6_wing_mesh(h_wall=H_WALL, h_edge=H_LE, h_te=0.5 * H_WALL,
                              h_wake=3.0 * H_WALL, h_far=H_FAR,
                              r_far=r_far_mac * MAC, tip_cap="round",
                              embed_wake=True)


def mesh_stats(mc):
    ht, hn, aniso = le_geometry(mc)
    h, r = cell_sizes(mc)
    g = dict(n_tet=len(mc.elements), n_wall=len(mc.boundary_faces["wall"]),
             n_le=le_face_count(mc), le_ht=ht, le_hn=hn, le_aniso=aniso)
    for lo, hi in SHELLS:
        m = (r >= lo) & (r < hi)
        k = shell_key(lo, hi)
        g[k] = float(np.median(h[m])) if m.any() else None
        g[k.replace("h_med", "n_cells")] = int(m.sum())
    return g


def shell_key(lo, hi):
    return f"h_med_r{lo:g}_{'inf' if hi > 1e8 else f'{hi:g}'}"


def guards(b, a):
    """G1a (surface), G1b (near body r < 2 MAC), G1c (the outer field really changed)."""
    d_ht = 100.0 * (a["le_ht"] - b["le_ht"]) / b["le_ht"]
    d_hn = 100.0 * (a["le_hn"] - b["le_hn"]) / b["le_hn"]
    print("\n=== G1a: the SURFACE must not move ===")
    print(f"  wall tris {b['n_wall']} -> {a['n_wall']}    LE faces {b['n_le']} -> {a['n_le']}")
    print(f"  LE h_t {d_ht:+.3f} %    LE h_n {d_hn:+.3f} %    (threshold {NEAR_PCT:g} %)")
    g1a = (a["n_wall"] == b["n_wall"] and a["n_le"] == b["n_le"]
           and abs(d_ht) < NEAR_PCT and abs(d_hn) < NEAR_PCT)
    print(f"  -> {'PASS' if g1a else 'FAIL'}")

    print("\n=== G1b: the NEAR BODY must not move   "
          "[the h_far arm moved r 1-2 MAC by +24.04 %] ===")
    print(f"  {'shell (MAC)':>12}{'base h_med':>13}{'r_far x2':>12}{'change':>10}"
          f"{'base cells':>12}{'x2 cells':>11}")
    worst = 0.0
    for lo, hi in SHELLS:
        k = shell_key(lo, hi)
        if b[k] is None or a[k] is None:
            continue
        d = 100.0 * (a[k] - b[k]) / b[k]
        if hi <= 2.0:
            worst = max(worst, abs(d))
        nk = k.replace("h_med", "n_cells")
        lbl = f"{lo:g}-{hi:g}" if hi < 1e8 else f"{lo:g}+"
        print(f"  {lbl:>12}{b[k]:>13.6f}{a[k]:>12.6f}{d:>9.2f}%{b[nk]:>12}{a[nk]:>11}")
    g1b = worst < NEAR_PCT
    print(f"  worst near-body (r < 2 MAC) change {worst:.2f} %  -> "
          f"{'PASS' if g1b else 'FAIL'}")

    nk = shell_key(8.0, 1e9).replace("h_med", "n_cells")
    g1c = a[nk] > 1.05 * b[nk]
    print("\n=== G1c (positive control): the OUTER field must actually change ===")
    print(f"  r > 8 MAC cells {b[nk]} -> {a[nk]}  -> "
          f"{'PASS' if g1c else 'FAIL -- this leg moved nothing anywhere'}")
    return g1a, g1b, g1c, worst, d_ht, d_hn


def main():
    exp = parse_experiment()
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    print("registered item 0: the domain-RADIUS arm (prereg 20260809-2000)")
    print("  guards are mesh-only, so they run BEFORE any solve\n")

    geom, rows = {}, []
    for tag, rf in LEGS:
        t0 = time.perf_counter()
        mc, wc = cut_wake(build(rf))
        t_gen = time.perf_counter() - t0
        geom[tag] = mesh_stats(mc)
        geom[tag].update(mesh=(mc, wc), t_gen_s=t_gen, r_far_mac=rf)
        g = geom[tag]
        print(f"  {tag:14} r_far {rf:g}*MAC  tets {g['n_tet']:>8}  wall tris "
              f"{g['n_wall']:>6}  LE faces {g['n_le']:>5}  h_t {g['le_ht']:.6f}  "
              f"h_n {g['le_hn']:.6f}  ({t_gen:.0f}s)", flush=True)
        if t_gen > GEN_GATE_S:
            print(f"  ★ mesh gate: {t_gen:.0f} s > {GEN_GATE_S:.0f} s -- stopping, "
                  f"no budget added"); return _write(rows, geom, 1)

    b, a = geom["R00_base"], geom["R10_domain_x2"]
    g1a, g1b, g1c, worst, d_ht, d_hn = guards(b, a)

    if not (g1a and g1b):
        print("\n=== reading (criteria fixed in the registration sec 5) ===")
        print("  ⇒ N: NO CLEAN FAR-FIELD ARM EXISTS with the available knobs.")
        print("     Pushing the outer boundary out also moves the surface and/or the")
        print("     near body, so 'far field vs near body' is UNTESTABLE in phase two.")
        print("     Recorded as a bounded negative in the capability boundary; decoupling")
        print("     the size fields in the generator is phase three's work. No solve is")
        print("     run -- the registration would refuse its displacement as a reading.")
        return _write(rows, geom, 0)
    if not g1c:
        print("\n  ⇒ VOID: the arm changed nothing outside either -- no reading.")
        return _write(rows, geom, 0)

    print("\n  arm is CLEAN (G1a+G1b PASS, G1c PASS) -> solving both legs\n")
    cps = {}
    for tag, rf in LEGS:
        mc, wc = geom[tag]["mesh"]
        t1 = time.perf_counter()
        r = solve_at(mc, wc, M_INF)
        wall = time.perf_counter() - t1
        hist = np.asarray(r.get("residual_history", []), dtype=float)
        rr = float(hist[-1]) if len(hist) else float("nan")
        nlim, nflr = int(r.get("n_limited") or 0), int(r.get("n_floored") or 0)
        usable = bool(r.get("converged")) and rr < 1e-9 and nlim == 0 and nflr == 0
        phi = np.asarray(r["phi"])
        cps[tag] = cp_at_exp_points(mc, phi, exp) if usable else None
        loc = clamp_map(mc, phi, M_INF)
        sh = {}
        for eta in (0.44, 0.65, 0.90):
            c = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
            sh[eta] = shock_report(c, M_INF)["upper"].get("x_shock")
        g = geom[tag]
        rows.append(dict(tag=tag, r_far_mac=rf, h_far=H_FAR,
                         **{k: (round(v, 9) if isinstance(v, float) else v)
                            for k, v in g.items() if k != "mesh"},
                         converged=bool(r.get("converged")), res_final=rr,
                         n_limited=nlim, n_floored=nflr, usable=usable,
                         m_max=round(float(np.sqrt(r["mach2_max"])), 5),
                         n_shock_cells=r.get("n_shock_cells"),
                         fallback_fired=bool(r.get("seed_fallback", {}).get("fired")),
                         **{f"x_shock_{e}": sh[e] for e in sh},
                         x_of_peak_q2=loc["x_of_peak_q2"], wall_s=round(wall, 1)))
        print(f"  {tag:14} conv={str(r.get('converged')):5} |R|={rr:.1e} "
              f"lim/flr={nlim}/{nflr} M_max={rows[-1]['m_max']:.4f} "
              f"shock={rows[-1]['n_shock_cells']}  ({wall:.0f}s)", flush=True)
        if not usable:
            print("    ★ NOT USABLE (G2): its displacement is not a reading", flush=True)
        if wall > LEG_GATE_S:
            print(f"  ★ leg gate: {wall:.0f} s > {LEG_GATE_S:.0f} s -- stopping, "
                  f"no budget added"); break

    _read(rows, cps, exp)
    return _write(rows, geom, 0)


def _read(rows, cps, exp):
    print("\n=== G2/G3 ===")
    bad = [r["tag"] for r in rows if not r["usable"]]
    print(f"  unusable legs: {bad or 'none'};  shock cells "
          f"{[r['n_shock_cells'] for r in rows]};  "
          f"x_shock(0.44) {[r['x_shock_0.44'] for r in rows]}")
    print("\n=== reading (criteria fixed in the registration sec 5) ===")
    if len(rows) < 2 or bad:
        print("  ⇒ RECORDED: a missing or unusable leg means there is no displacement")
        print("     the registration would accept as a reading."); return
    d = band_disp(cps["R10_domain_x2"], cps["R00_base"])
    z = band_disp(cps["R00_base"], cps["R00_base"])
    assert all(v == 0.0 for v in z.values()), f"zero test failed: {z}"
    eb, ea = band_err(cps["R00_base"], exp), band_err(cps["R10_domain_x2"], exp)
    for r in rows:
        src = eb if r["tag"] == "R00_base" else ea
        r.update({f"err_{k}": round(v, 6) for k, v in src.items()})
        if r["tag"] != "R00_base":
            r.update({f"disp_{k}": round(v, 6) for k, v in d.items()})
    dl = d["LE_upper"]
    print("  zero test (base minus itself): exactly 0 -> PASS")
    print(f"  Delta_LE (displacement) = {dl:.5f}    "
          f"[the unclean h_far arm gave {HFAR_ARM_DELTA_LE:.5f}]")
    print(f"  MID_up {d['MID_upper']:.5f}   TE_up {d['TE_upper']:.5f}   "
          f"pooled {d['pooled']:.5f}")
    print(f"  [RECORDED, no criterion] LE error vs EXPERIMENT "
          f"{eb['LE_upper']:.5f} -> {ea['LE_upper']:.5f} "
          f"({100 * (ea['LE_upper'] - eb['LE_upper']) / eb['LE_upper']:+.1f} %)")
    if dl < 0.02:
        print("\n  ⇒ C: the FAR FIELD ITSELF barely controls the LE band. So h_far's")
        print("     0.098 was mostly NEAR-BODY GRADING and the erratum's rewrite is")
        print("     CONFIRMED -- the boundary document says the LE deficit is dominated")
        print("     by near-body bulk grading, and B5's domain-robustness prior carries")
        print("     from a 2.5-D integral to a 3-D local quantity.")
    elif dl >= 0.05:
        print("\n  ⇒ D: the far field itself ALSO controls the LE band. The rewrite is")
        print("     INCOMPLETE -- both near-body grading and the outer boundary matter --")
        print("     and B5's prior does NOT carry to a 3-D local quantity.")
    else:
        print("\n  ⇒ I: between 0.02 and 0.05 -> RECORDED, no conclusion. Per sec 7 a")
        print("     further r_far step needs its own registration.")


def _write(rows, geom, code):
    if not rows:            # guard-only exit: still commit the mesh evidence
        rows = [dict(tag=t, **{k: (round(v, 9) if isinstance(v, float) else v)
                               for k, v in g.items() if k != "mesh"})
                for t, g in geom.items()]
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    return code


if __name__ == "__main__":
    sys.exit(main())
