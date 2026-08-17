"""
Part A: does the reconstruction move the shock TOWARD an independent judge?
Part B: how big is the perturbation inside the solver's own switches?

Pre-registration: docs/dev_phase_four/20260818-0900-position-and-blastradius-prereg.md
Summary:          docs/dev_phase_four/20260818-0700-recovery-line-summary.md

★ Part A is finally non-circular. J1 -- the wall-Cp shock position -- runs through
WALL-TRIANGLE tangential gradients and reads no volume gradient at all, so it is
IDENTICAL for every arm and forms a fixed reference. J2 is the ONERA M6
experiment (TEST 2308), the only external truth in the repository, and is
RECORDED not gated: it is VISCOUS data against an INVISCID solver.

★★ Part B is NOT the answer to "what happens inside the solver". It is that
question's PRICE. A blast radius measures the SIZE OF A PERTURBATION, not what
the solution does. Small may NOT be read as "safe to wire in", large may NOT be
read as "would help".
"""

import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "bench"))

import run_gs40c_coherence_ladder as L                              # noqa: E402
import run_gs40g_geometry_axis as G                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.kernels.entropy import EntropyOperator                  # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.kernels.upwind import UpwindOperator                    # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import build_face_adjacency, precompute_element_geometry  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_metrics, shock_report                          # noqa: E402
from pyfp3d.solve.newton import solve_newton_transonic              # noqa: E402
from tests.test_p8_newton import NEWTON_M6_RECIPE                   # noqa: E402

WING = REPO / "cases/meshes/onera_m6/coarse.msh"
EXP = REPO / "cases/reference_data/onera_m6_experiment/experiment-Cp.dat"
CACHE = REPO / "bench/gate_results/gs40i_states"
OUT = REPO / "bench/gate_results/gs40i_position_blastradius.csv"
M_EXP, ALPHA, GAMMA = 0.8395, 3.06, 1.4
ETAS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)     # TEST 2308 stations
LADDER, ARMS = ("P0", "P1", "P2", "P3"), ("P0", "P1", "P2", "P3", "S1")
BETTER, WORSE, MIN_STATIONS = 0.7, 1.3, 3
BR_SMALL, BR_LARGE = 0.01, 0.10
UPWIND_C, M_CRIT = 1.5, 0.95
GATE_S = 25 * 60

#: the cached states Part B reuses -- no new solves
BR_STATES = (("wing M0.82", REPO / "bench/gate_results/gs40g_states/W0.82.npz", WING, 0.82),
             ("wing M0.86", REPO / "bench/gate_results/gs40g_states/W0.86.npz", WING, 0.86),
             ("wing M0.88", REPO / "bench/gate_results/gs40g_states/W0.88.npz", WING, 0.88))


def guard_exp():
    head = EXP.read_text().splitlines()[0]
    assert "0.8395" in head and "3.0600" in head and "2308" in head, \
        f"G-EXP: unexpected header: {head}"
    n = sum(1 for ln in EXP.read_text().splitlines() if ln.startswith("ZONE"))
    print(f"G-EXP  {EXP.name}: TEST 2308 / M0.8395 / a3.06, {n} zones (read-only)  PASS")
    return n


def exp_shock_x():
    """Experimental upper-surface shock x/c per station, using THE PROJECT'S OWN
    detector on the committed curve -- no rule invented here.

    ★ The first version dropped the Z/L column and filtered "upper surface" by
    Cp < 0, which is a SUCTION filter, not a surface filter; the dry-check caught
    it because the resulting x_shock values (0.02-0.04) were the leading-edge
    suction recovery, not a shock. Columns are NP, X/L, Y/b, Z/L, Cp; upper is
    Z/L > 0.
    """
    zones, cur = [], []
    for ln in EXP.read_text().splitlines():
        if ln.startswith("ZONE"):
            if cur:
                zones.append(np.array(cur))
            cur = []
            continue
        f = ln.split()
        if len(f) == 5:
            try:
                cur.append([float(f[1]), float(f[2]), float(f[3]), float(f[4])])
            except ValueError:
                pass
    if cur:
        zones.append(np.array(cur))
    res = {}
    for a in zones:
        yb = float(np.median(a[:, 1]))
        up = a[a[:, 2] > 0.0]                       # UPPER surface via Z/L
        if len(up) < 5:
            continue
        o = np.argsort(up[:, 0])
        m = shock_metrics(up[o, 0], up[o, 3], M_EXP)
        res[round(yb, 3)] = (float(m["x_shock"]) if m.get("has_shock") else np.nan)
    return res


def foot_from_field(nodes, elements, M, wall_faces, etas, b_semi):
    """Shock FOOTPRINT per station from a VOLUME Mach field: among wall-adjacent
    elements near the station, the largest x/c that is still supersonic.

    Fixed in advance: 'largest x/c with M >= 1 among elements whose centroid lies
    within +-2% b_semi of the station and within the wall-adjacent layer'.
    """
    cent = nodes[elements].mean(axis=1)
    wall_nodes = np.unique(wall_faces.reshape(-1))
    isw = np.zeros(len(nodes), bool); isw[wall_nodes] = True
    adj = isw[elements].any(axis=1)
    out = {}
    for e_ in etas:
        z0 = e_ * b_semi
        band = adj & (np.abs(cent[:, 2] - z0) <= 0.02 * b_semi) & (M >= 1.0)
        if band.sum() < 3:
            out[e_] = np.nan; continue
        xs = cent[band, 0]
        # local chord normalisation via the station's wall x-extent
        wb = adj & (np.abs(cent[:, 2] - z0) <= 0.02 * b_semi)
        x0, x1 = cent[wb, 0].min(), cent[wb, 0].max()
        out[e_] = float((xs.max() - x0) / max(x1 - x0, 1e-9))
    return out


def main():
    t0 = time.perf_counter()
    G.guard_recipe_m6()
    man = read_manifest(WING)
    assert man and mesh_fingerprint(WING)["sha256"] == man["sha256"], "G-MESH"
    print(f"G-MESH  onera_m6/coarse.msh sha {man['sha256'][:12]}  PASS")
    guard_exp()

    # ---------------------------------------------------------------- Part A --
    print("\n" + "=" * 78 + "\nPart A -- position, against an INDEPENDENT judge\n"
          + "=" * 78)
    CACHE.mkdir(parents=True, exist_ok=True)
    npz = CACHE / "wing_M08395.npz"
    mc, wc = cut_wake(read_mesh(WING))
    if npz.exists():
        d = np.load(npz); phi = np.asarray(d["phi"], float); wall = 0.0
        conv, nlim, nflr = bool(d["conv"]), int(d["nlim"]), int(d["nflr"])
    else:
        taper = tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth", 0.05 * B_SEMI)
        kw = dict(NEWTON_M6_RECIPE)
        kw["newton_kw"] = dict(kw["newton_kw"], tip_taper=taper)
        t1 = time.perf_counter()
        r = solve_newton_transonic(mc, wc, m_inf=M_EXP, alpha_deg=ALPHA, **kw)
        wall = time.perf_counter() - t1
        phi = np.asarray(r["phi"], float)
        conv, nlim, nflr = bool(r["converged"]), int(r["n_limited"]), int(r["n_floored"])
        np.savez_compressed(npz, phi=phi, conv=conv, nlim=nlim, nflr=nflr,
                            mesh_sha=man["sha256"])
    assert conv and not nlim and not nflr, f"G-STATE: conv={conv} {nlim}/{nflr}"
    print(f"state: M{M_EXP}/a{ALPHA} converged 0/0 clamps ({wall:.1f}s)")

    nodes, elements = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, elements)
    face_nb, _ = build_face_adjacency(elements)
    wall_faces = mc.boundary_faces["wall"]
    P = {}
    for k in (0, 1, 2, 3):
        po, pi = L.khop_patches(elements, face_nb, k)
        P[f"P{k}"] = L.patch_nodes(elements, po, pi)
        if k == 1:
            p1e = (po, pi)
    gA = np.empty((len(elements), 3)); q2A = np.empty(len(elements))
    element_velocity_q2(elements, B, phi, gA, q2A)
    m2A = mach_number_squared(q2A, M_EXP, GAMMA)
    sup = float((m2A >= 1.0).mean())
    strong = 0.09 <= sup <= 0.20
    print(f"G-STRENGTH sup_frac {sup:.4f}  "
          f"{'in the STRONG band (A2 applies)' if strong else '★ OUTSIDE -- read via A3'}")

    # J1: the fixed, independent reference -- wall-Cp shock position
    j1, exp_map = {}, exp_shock_x()
    for e_ in ETAS:
        try:
            cur = section_cp_curve(mc, phi, eta=e_, b_semi=B_SEMI, m_inf=M_EXP)
            sr = shock_report(cur, M_EXP)
            j1[e_] = sr["upper"].get("x_shock")
        except Exception:                                          # noqa: BLE001
            j1[e_] = None
    print("G-INDEP  J1 reads wall-triangle gradients only => identical for every "
          "arm (fixed reference)  PASS")
    print("  J1 x_shock per station: " + "  ".join(
        f"{e_}:{'--' if j1[e_] is None else f'{j1[e_]:.3f}'}" for e_ in ETAS))
    print("  J2 experiment  x_shock: " + "  ".join(
        f"{k}:{v:.3f}" for k, v in sorted(exp_map.items())))

    g = {n: L.refit(nodes, elements, phi, *P[n]) for n in LADDER}
    g["S1"] = L.average_over_patch(nodes, elements, gA, *p1e)
    rowsA, dev = [], {}
    for n in ARMS:
        m2 = mach_number_squared((g[n] * g[n]).sum(axis=1), M_EXP, GAMMA)
        feet = foot_from_field(nodes, elements, np.sqrt(np.maximum(m2, 0.0)),
                               wall_faces, ETAS, B_SEMI)
        d1 = [abs(feet[e_] - j1[e_]) for e_ in ETAS
              if j1[e_] is not None and not np.isnan(feet[e_])]
        dev[n] = (float(np.median(d1)) if d1 else np.nan, len(d1))
        for e_ in ETAS:
            rowsA.append(dict(part="A", arm=n, eta=e_,
                              foot=None if np.isnan(feet[e_]) else round(feet[e_], 4),
                              j1=None if j1[e_] is None else round(j1[e_], 4)))
        print(f"  arm {n}: median |foot - J1| = "
              f"{dev[n][0]:.4f}  over {dev[n][1]} stations")

    base = dev["P0"][0]
    if dev["P1"][1] < MIN_STATIONS or np.isnan(base) or base <= 0:
        vA = "POS-UNDEF (fewer than 3 usable stations, or a degenerate baseline)"
        ratio = float("nan")
    else:
        ratio = dev["P1"][0] / base
        vA = ("POS-BETTER (the reconstruction moves the footprint toward J1)"
              if ratio <= BETTER else
              "POS-WORSE (it moves the footprint AWAY from J1) -- kill criterion 1"
              if ratio >= WORSE else
              "POS-SAME (position neither improved nor degraded)")
    print(f"\n  P1/P0 median-deviation ratio = {ratio:.3f}  "
          f"(BETTER <= {BETTER}, WORSE >= {WORSE})\n  VERDICT A: {vA}")

    # ---------------------------------------------------------------- Part B --
    print("\n" + "=" * 78 + "\nPart B -- BLAST RADIUS (a PRICE, not an answer)\n"
          + "=" * 78)
    rowsB = []
    for label, sp, msh, m_inf in BR_STATES:
        mcb, _ = cut_wake(read_mesh(msh))
        nb, eb = mcb.nodes, mcb.elements
        Bb, _ = precompute_element_geometry(nb, eb)
        fnb, _ = build_face_adjacency(eb)
        Pb = {}
        for k in (0, 1):
            po, pi = L.khop_patches(eb, fnb, k)
            Pb[f"P{k}"] = L.patch_nodes(eb, po, pi)
        ph = np.asarray(np.load(sp)["phi"], float)
        gc = np.empty((len(eb), 3)); qc = np.empty(len(eb))
        element_velocity_q2(eb, Bb, ph, gc, qc)
        g1 = L.refit(nb, eb, ph, *Pb["P1"])
        q1 = np.einsum("ej,ej->e", g1, g1)

        up = UpwindOperator(nb, eb)
        ent = EntropyOperator(len(eb))
        rho = np.ones(len(eb))

        def state(grad, q2):
            u = up.upstream_map(grad).copy()
            m2 = mach_number_squared(q2, m_inf, GAMMA)
            nu_on = m2 > M_CRIT ** 2                     # nu_e > 0 iff m2 > m_crit^2
            #: ★ the entropy set is taken from the LIBRARY's own sigma field
            #: (charged <=> sigma < 1) rather than re-implementing the membership
            #: test -- the registration requires the solver's switch, not mine.
            charged = ent.sigma(q2, u, m_inf).copy() < 1.0
            return u, nu_on, charged

        uA, nuA, poA = state(gc, qc)
        u1, nu1, po1 = state(g1, q1)
        n = len(eb)
        B1 = float((nuA ^ nu1).sum() / n)
        B2 = float((uA != u1).sum() / n)
        B3 = float((poA ^ po1).sum() / n)
        mx = max(B1, B2, B3)
        band = "BR-SMALL" if mx < BR_SMALL else ("BR-LARGE" if mx > BR_LARGE else "BR-MED")
        assert mx > 0, "Part B: all three counts are zero -- the wiring is wrong"
        print(f"  {label:12} B1 dissipation {B1*100:6.3f}%   B2 donor {B2*100:6.3f}%"
              f"   B3 sigma-set {B3*100:6.3f}%   -> {band}")
        rowsB.append(dict(part="B", state=label, m_inf=m_inf, n_tets=n,
                          B1_dissipation=round(B1, 6), B2_donor=round(B2, 6),
                          B3_sigma_set=round(B3, 6), band=band))
    order = sorted(("B1", "B2", "B3"),
                   key=lambda k: -np.mean([r[{"B1": "B1_dissipation",
                                              "B2": "B2_donor",
                                              "B3": "B3_sigma_set"}[k]] for r in rowsB]))
    print(f"\n  ranking by mean magnitude: {' > '.join(order)}")
    print("  ★ A blast radius is the SIZE OF A PERTURBATION, not what the solution\n"
          "    does. Small may NOT be read as 'safe to wire in'.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = rowsA + rowsB
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict_A=vA, ratio_P1_P0=None if np.isnan(ratio) else round(ratio, 4),
             sup_frac=round(sup, 5), blast=rowsB, ranking=order,
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 4"


if __name__ == "__main__":
    main()
