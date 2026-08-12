"""Is the all-scales LE gain bought by ALLOCATION or by CELL COUNT?

Pre-registered in docs/dev_phase_three/20260814-0900-fixed-budget-allocation-prereg.md,
with addendum #1 (docs/dev_phase_three/20260814-1000-fixed-budget-addendum1.md) committed before this
file was written.

The question: proportional all-scales refinement was measured to reduce the LE-band error, but it costs
25x the cells. If the gain comes from the PROPORTION between scales, grading can reach the accuracy
inside the cost budget; if it comes from the NUMBER of cells, allocation has no solution -- and that is
the pre-registered condition under which raising the element ORDER (isoparametric P2) would get its
first measurement-backed justification. Both answers are named bands with opposite consequences.

★★★ addendum #1: every mesh here is tip_cap="round". The four legs that produced the -54.9 % headline
were built through run_le_response.build, whose BASE hard-codes tip_cap="flat", while production is
round -- and P13/G13.3 measured the flat cap DIVERGING under refinement, so the project's standing rule
is that any refinement-based claim on a flat-cap mesh has a false premise. This script therefore does
NOT reuse that builder, re-measures its own reference curve, and asserts the cap (guard G-C) -- the
guard that did not exist, which is exactly what got missed.

★ The primary reading is not "two numbers compared" but "does a re-proportioned mesh land ON this
round's own error-vs-cell-count curve or BELOW it". Cell count is an OUTPUT of the mesh generator, so
"hold it fixed" would repeat this season's fourth criterion defect (treating an output as an input);
comparing against the curve at each leg's own count absorbs the residual mismatch instead.

Outputs (TRACKED): bench/gate_results/task3_fixed_budget.csv
"""

import csv
import os
import sys
import time

#: 16, matching run_m3_budget's own setdefault, so the LE readings stay same-provenance
os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from pyfp3d.mesh.metrics import (compute_aspect_ratios,               # noqa: E402
                                 compute_min_dihedral_angles,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh, write_mesh                  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                             # noqa: E402
from pyfp3d.meshgen.wing3d import (B_SEMI, MAC, chord_at,             # noqa: E402
                                   onera_m6_wing_mesh, x_le)
from pyfp3d.physics.isentropic import mach_number_squared             # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                  # noqa: E402
from pyfp3d.post.surface import (cl_kj_3d, planform_area,             # noqa: E402
                                 wall_force_coefficients)
#: imported, never re-typed: the recipe and the band machinery are the single source
from run_m3_budget import (ALPHA, ETAS, M_INF, N_UNMASKED,            # noqa: E402
                           band_rms, parse_experiment, solve)

CSV = os.path.join(HERE, "gate_results", "task3_fixed_budget.csv")
#: ★ prereg §7: NOT cases/meshes/onera_m6/ -- that directory is glob-ingested wholesale by
#: tests/test_p2_wake_cut.py (hard rule 7), so dropping meshes there changes what the suite eats.
SCRATCH = os.environ.get(
    "PYFP3D_SCRATCH",
    "/tmp/claude-1000/-home-lrz-codes-UP3D/3c5b43c4-b62c-4a09-b4da-9b9c7128d43e/scratchpad/meshes")
#: ★★ G-C: production geometry. addendum #1 -- the whole reason this script exists standalone.
TIP_CAP = "round"
BASE = dict(r_far=15.0 * MAC, tip_cap=TIP_CAP, embed_wake=True)

#: the reference curve: two points on the production family's own ray (h_far left at its default,
#: matching what the ladder did -- min(2.5, 120 h_wall), so h_wall 0.030 is CLAMPED at 2.5 and
#: therefore off-ray, a defect _level_params has recorded since 2026-07-13 / P13-G13.3)
CURVE = (("C_coarse", 0.030, 0.015, 0.015, None),
         ("C_medium", 0.015, 0.0075, 0.0075, None))
#: allocation arms at the coarse budget. (tag, h_far, h_edge_mode, h_wall seed)
#:   h_edge_mode "prop" keeps h_edge = 0.5 h_wall (scales with the wall);
#:   "abs:<v>" pins h_edge at an absolute value so cells move INTO the LE band.
ARMS = (("A1_bulk_heavy", 1.25, "prop", 0.036),
        ("A2_wall_heavy", 5.00, "prop", 0.026),
        ("A3_le_heavy", None, "abs:0.0075", 0.034))
BUDGET_TOL = 0.10
MAX_GENS = 4
LEG_GATE_S, TOTAL_GATE_S = 900.0, 3600.0
#: G-Q is RELATIVE to the control's own mesh: coarse already sits near AR 15, so an absolute
#: gate would kill the control leg. (hex v1: the guard checked orientation only, AR was 1305.)
AR_FACTOR, DIH_FACTOR = 2.0, 0.5
SHELLS = ((0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0))


def build(h_wall, h_edge, h_te, h_far):
    """Production geometry, LE/TE split. h_far None = the generator's own default."""
    kw = dict(h_wall=h_wall, h_edge=h_edge, h_te=h_te, h_wake=3.0 * h_wall, **BASE)
    if h_far is not None:
        kw["h_far"] = h_far
    assert kw["tip_cap"] == "round", "G-C: this round is production geometry only"
    return onera_m6_wing_mesh(**kw)


def cached_mesh(tag, h_wall, h_edge, h_te, h_far):
    os.makedirs(SCRATCH, exist_ok=True)
    key = f"rnd_{h_wall:g}_{h_edge:g}_{h_te:g}_{'def' if h_far is None else f'{h_far:g}'}"
    path = os.path.join(SCRATCH, key + ".msh")
    if os.path.exists(path):
        print(f"  [cached mesh] {key}", flush=True)
        return read_mesh(path)
    t0 = time.perf_counter()
    mesh = build(h_wall, h_edge, h_te, h_far)
    print(f"  [generated] {key}  {time.perf_counter() - t0:.0f}s", flush=True)
    try:
        write_mesh(mesh, path)
    except Exception as exc:                                          # noqa: BLE001
        print(f"  (not cached: {exc})")
    return mesh


def le_geometry(mc):
    """(h_t, h_n, anisotropy) medians in the LE band -- run_le_response's read, reused verbatim
    in definition so the numbers are comparable to the phase-two guard rows."""
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
        return (float("nan"),) * 3
    return float(np.median(ht)), float(np.median(hn)), float(np.median(ht / hn))


def cell_sizes(mc):
    p = mc.nodes[mc.elements]
    v = np.abs(np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                         np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))) / 6.0
    return np.cbrt(np.maximum(v, 1e-300)), np.linalg.norm(p.mean(axis=1), axis=1)


def mesh_row(mc):
    """G-M + G-Q: what the mesh actually IS -- skin AND volume, per the phase-two lesson that a
    guard reading only the surface called an arm clean while the volume had moved 24 %."""
    ht, hn, aniso = le_geometry(mc)
    h, r = cell_sizes(mc)
    row = dict(n_tet=int(len(mc.elements)), n_node=int(len(mc.nodes)),
               n_wall=int(len(mc.boundary_faces["wall"])),
               le_ht=ht, le_hn=hn, le_aniso=aniso,
               ar_max=float(compute_aspect_ratios(mc.nodes, mc.elements).max()),
               dih_min=float(compute_min_dihedral_angles(mc.nodes, mc.elements).min()))
    for lo, hi in SHELLS:
        m = (r >= lo) & (r < hi)
        row[f"h_med_r{lo:g}_{hi:g}"] = float(np.median(h[m])) if m.any() else None
    return row


def flow_row(mc, wc, r, s_ref, exp):
    phi = np.asarray(r["phi"])
    gamma = np.atleast_1d(np.asarray(r["gamma"]))
    Bg, _Vg = precompute_element_geometry(mc.nodes, mc.elements)
    gg = np.einsum("eaj,ea->ej", Bg, phi[mc.elements])
    m2 = mach_number_squared(np.einsum("ej,ej->e", gg, gg), M_INF)
    f = wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
    o = np.argsort(wc.station_z)
    cent = mc.nodes[mc.elements].mean(axis=1)
    inb = cent[:, 2] <= 0.92 * B_SEMI
    curves, acc = {}, {}
    for eta in ETAS[:N_UNMASKED]:
        try:
            curves[eta] = section_cp_curve(mc, phi, eta=eta, b_semi=B_SEMI, m_inf=M_INF)
        except Exception:                                             # noqa: BLE001
            continue
    per_station = {}
    for eta in list(curves):
        for name, (ss, nn) in band_rms(curves, exp, eta).items():
            a = acc.setdefault(name, [0.0, 0])
            a[0] += ss
            a[1] += nn
            if name == "LE_upper":
                per_station[eta] = float(np.sqrt(ss / nn)) if nn else None
    out = dict(cl_p=float(f["cl"]),
               cl_kj=float(cl_kj_3d(gamma[o], wc.station_z[o], s_ref, B_SEMI)),
               m_max=float(np.sqrt(m2.max())),
               m_max_inboard=float(np.sqrt(m2[inb].max())) if inb.any() else None,
               m1_max=r.get("m1_max"), sigma_min=r.get("sigma_min"))
    for name, (ss, nn) in sorted(acc.items()):
        out[f"rms_{name}"] = float(np.sqrt(ss / nn)) if nn else None
    for eta, v in sorted(per_station.items()):
        out[f"le_up_eta{eta:g}"] = v
    return out


def run_leg(tag, h_wall, h_edge, h_te, h_far, exp, elapsed):
    if elapsed[0] > TOTAL_GATE_S:
        print(f"  ★ total gate {elapsed[0]:.0f}s exceeded -- {tag} NOT run (kill clause 4)")
        return None
    mc, wc = cut_wake(cached_mesh(tag, h_wall, h_edge, h_te, h_far))
    mrow = mesh_row(mc)
    s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
    t0 = time.perf_counter()
    r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=0.05)
    wall = time.perf_counter() - t0
    elapsed[0] += wall
    conv = bool(r.get("converged"))
    res = float(r.get("residual_history", [float("nan")])[-1])
    print(f"  {tag:16} tets {mrow['n_tet']:>8}  conv={conv} |R|={res:.2e} "
          f"lim={r.get('n_limited')} flr={r.get('n_floored')}  ({wall:.0f}s)", flush=True)
    row = dict(leg=tag, tip_cap=TIP_CAP, h_wall=h_wall, h_edge=h_edge, h_te=h_te,
               h_far=("default" if h_far is None else h_far), converged=conv,
               res_final=res, n_limited=r.get("n_limited"), n_floored=r.get("n_floored"),
               solve_s=round(wall, 1), over_leg_gate=bool(wall > LEG_GATE_S), **mrow)
    if conv:
        row.update(flow_row(mc, wc, r, s_ref, exp))
    return row


def search_budget(tag, h_far, mode, seed, target, h_te_ref):
    """Cell count is an OUTPUT: search h_wall, never claim it was set. Up to MAX_GENS
    generations (kill clause 3 forbids buying more)."""
    h = seed
    hist = []
    for k in range(MAX_GENS):
        h_edge = 0.5 * h if mode == "prop" else float(mode.split(":")[1])
        mesh = cached_mesh(f"{tag}_g{k}", h, h_edge, h_te_ref, h_far)
        n = len(mesh.elements)
        hist.append((h, h_edge, n))
        rel = n / target - 1.0
        print(f"    search {tag} gen{k}: h_wall {h:.5f} h_edge {h_edge:.5f} "
              f"-> {n} tets ({100 * rel:+.1f} % of budget)", flush=True)
        if abs(rel) <= BUDGET_TOL:
            return h, h_edge, n, True, hist
        #: n ~ h^-3 locally; damped so a bad local exponent cannot fling the search
        h *= (n / target) ** (1.0 / 3.0)
    h, h_edge, n = min(hist, key=lambda t: abs(t[2] / target - 1.0))
    return h, h_edge, n, False, hist


def main():
    print("resolved threads: " + ", ".join(
        f"{k}={os.environ.get(k)}" for k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS",
                                             "OPENBLAS_NUM_THREADS")))
    print(f"load average: {os.getloadavg()}")
    print(f"★ G-C: tip_cap = {TIP_CAP} on every leg (addendum #1)\n")
    exp = parse_experiment()
    rows, elapsed = [], [0.0]

    print("=== reference curve (HEAD, round cap, production leg) ===")
    for tag, hw, he, hte, hf in CURVE:
        row = run_leg(tag, hw, he, hte, hf, exp, elapsed)
        if row is not None:
            rows.append(row)
    _write(rows)

    ctrl = next((x for x in rows if x["leg"] == "C_coarse"), None)
    tgt = next((x for x in rows if x["leg"] == "C_medium"), None)
    if ctrl is None or not ctrl["converged"]:
        print("\n★ control leg C_coarse did not converge -- STOP (kill clause 2).")
        return 1
    if tgt is None or not tgt["converged"]:
        print("\n★ C_medium did not converge -- the curve has one point. STOP (kill clause 1).")
        return 1
    G = tgt and (ctrl["rms_LE_upper"] - tgt["rms_LE_upper"])
    print(f"\n=== the gain G on the production family = {ctrl['rms_LE_upper']:.6f} - "
          f"{tgt['rms_LE_upper']:.6f} = {G:+.6f} "
          f"({100 * G / ctrl['rms_LE_upper']:+.1f} %) ===")
    if G <= 0:
        print("  -> ★ G <= 0: on ROUND caps, coarse is not worse than medium. That REFUTES the")
        print("     flat-cap ladder's reading rather than confirming it. STOP (kill clause 1),")
        print("     and report it as the finding -- do not switch quantity and retry.")
        _write(rows)
        return 1
    budget = ctrl["n_tet"]
    print(f"  budget = C_coarse's own {budget} tets (+-{100 * BUDGET_TOL:.0f} %); "
          f"target error = C_medium's {tgt['rms_LE_upper']:.6f}")

    print("\n=== allocation arms at the coarse budget ===")
    for tag, h_far, mode, seed in ARMS:
        h, h_edge, n, ok, _hist = search_budget(tag, h_far, mode, seed, budget, 0.015)
        if not ok:
            print(f"  ★ {tag}: budget NOT met after {MAX_GENS} generations "
                  f"({n} vs {budget}) -- recorded, read on the curve only (kill clause 3)")
        row = run_leg(tag, h, h_edge, 0.015, h_far, exp, elapsed)
        if row is not None:
            row["budget_met"] = ok
            rows.append(row)
            _write(rows)
    _write(rows)
    return report(rows, ctrl, tgt, G, budget)


def report(rows, ctrl, tgt, G, budget):
    arms = [x for x in rows if x["leg"].startswith("A")]
    print("\n=== G-C / G-B / G-Q ===")
    caps = {x["tip_cap"] for x in rows}
    print(f"  G-C tip_cap set = {caps}  -> {'PASS' if caps == {'round'} else '★ FAIL'}")
    good = []
    for x in arms:
        inb = abs(x["n_tet"] / budget - 1.0) <= BUDGET_TOL
        q = (x["ar_max"] <= AR_FACTOR * ctrl["ar_max"]
             and x["dih_min"] >= DIH_FACTOR * ctrl["dih_min"])
        print(f"  {x['leg']:16} tets {x['n_tet']:>7} ({100 * (x['n_tet'] / budget - 1):+5.1f} %) "
              f"budget {'OK ' if inb else '★NO'}  AR {x['ar_max']:7.2f} "
              f"(ctrl {ctrl['ar_max']:.2f}) dih {x['dih_min']:6.3f} "
              f"(ctrl {ctrl['dih_min']:.3f}) quality {'OK' if q else '★EXCLUDED'}  "
              f"conv={x['converged']}")
        if q and x["converged"]:
            good.append(x)
    print("\n=== G-M: what each leg actually changed (skin AND volume) ===")
    print(f"  {'leg':16}{'LE h_t':>10}{'LE h_n':>10}{'aniso':>8}"
          + "".join(f"{f'h@r{lo:g}-{hi:g}':>12}" for lo, hi in SHELLS))
    for x in [ctrl, tgt] + arms:
        print(f"  {x['leg']:16}{x['le_ht']:>10.6f}{x['le_hn']:>10.6f}{x['le_aniso']:>8.3f}"
              + "".join(f"{(x.get(f'h_med_r{lo:g}_{hi:g}') or float('nan')):>12.5f}"
                        for lo, hi in SHELLS))
    print("\n=== G-S: the solution moves too (RECORDED, not noise) ===")
    for x in [ctrl, tgt] + [a for a in arms if a["converged"]]:
        print(f"  {x['leg']:16} M_max {x['m_max']:.4f} (inboard {x['m_max_inboard']:.4f})  "
              f"m1_max {x.get('m1_max')}  sigma_min {x.get('sigma_min')}  "
              f"cl_p {x['cl_p']:.6f}")

    print("\n=== B1 / B2 / B3 (binding = LE upper RMS) ===")
    if len(good) < 2:
        print(f"  -> ★ only {len(good)} converged+quality arm(s): UNDEFINED, no direction claimed "
              f"(standing question 2 -- a spread over fewer than two legs is not 'small').")
        return 0
    if len({round(x["rms_LE_upper"], 12) for x in good} | {round(ctrl["rms_LE_upper"], 12)}) == 1:
        print("  -> ★ every leg bit-identical: NOT B2, SUSPICION -- the knob never reached the")
        print("     solve. Check the instrument first (kill clause 5).")
        return 1
    best = min(good, key=lambda x: x["rms_LE_upper"])
    r = (ctrl["rms_LE_upper"] - best["rms_LE_upper"]) / G
    print(f"  control C_coarse {ctrl['rms_LE_upper']:.6f}  ->  best {best['leg']} "
          f"{best['rms_LE_upper']:.6f}")
    for x in sorted(good, key=lambda y: y["rms_LE_upper"]):
        print(f"    {x['leg']:16} {x['rms_LE_upper']:.6f}  "
              f"({100 * (x['rms_LE_upper'] / ctrl['rms_LE_upper'] - 1):+6.2f} % vs control, "
              f"r = {(ctrl['rms_LE_upper'] - x['rms_LE_upper']) / G:+.3f})")
    #: ★ two-sidedness of the DOSE: bulk-heavy and wall-heavy must straddle the control,
    #: else "proportion is the axis" is not established even if one arm looks good.
    a1 = next((x for x in good if x["leg"].startswith("A1")), None)
    a2 = next((x for x in good if x["leg"].startswith("A2")), None)
    straddle = (a1 is not None and a2 is not None
                and (a1["rms_LE_upper"] - ctrl["rms_LE_upper"])
                * (a2["rms_LE_upper"] - ctrl["rms_LE_upper"]) < 0)
    print(f"\n  dose two-sidedness (A1 and A2 on opposite sides of control): "
          f"{'YES' if straddle else '★ NO'}"
          + ("" if a1 is not None and a2 is not None else "   (one arm missing)"))
    if r >= 1.0 / 3.0 and straddle:
        print(f"  -> ★ B1  r = {r:.3f} >= 1/3 AND the dose is two-sided ⇒ ALLOCATION IS A LEVER:")
        print("     the same cell count buys materially better accuracy ⇒ next spend is a")
        print("     graded/adaptive recipe, and raising the element ORDER is DEFERRED.")
    elif all(abs(x["rms_LE_upper"] / ctrl["rms_LE_upper"] - 1) < 0.10 for x in good):
        print(f"  -> ★★ B2  every arm within +-10 % of control (best r = {r:.3f}) ⇒ the gain is")
        print("     CELL COUNT, not allocation ⇒ the pre-registered condition for isoparametric")
        print("     P2 is MET for the first time (it becomes a measurement-backed candidate).")
        print("     ★ Still not an authorisation: penetration 4/4 and the prohibition stand.")
    else:
        print(f"  -> B3  r = {r:.3f}, straddle={straddle} ⇒ RECORDED, no direction claimed.")
        print("     ★ a one-sided dose is itself informative: bulk fine/coarse is not the axis.")
    return 0


def _write(rows):
    if not rows:
        return
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
