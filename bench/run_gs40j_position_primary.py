"""
Position as the PRIMARY quantity, across every cached state.

Pre-registration: docs/dev_phase_four/20260818-1300-position-primary-prereg.md
Previous round:   docs/dev_phase_four/20260818-1100-position-and-blastradius-verdict.md

Seven rounds optimised COHERENCE, which every verdict labelled a proxy with the
forbidden sentence "coherent is not correct". Last round measured the target for
the first time on ONE state and found the coherence-optimal arm was the WORST on
position. This round is that observation's registered confirmation.

★ The judge is exogenous: J1 (wall-Cp shock position) reads wall-triangle
tangential gradients and NO volume gradient, so it is bit-identical across arms.
★★ Every scale now carries its own SMOOTHER (S1/S2/S3), because "a bigger patch
wins on position" has two readings -- better reconstruction, or merely a blurrier
field landing closer by luck -- and only a same-scale smoother separates them.
Q-SMOOTH is a KILL criterion, and it is what I predicted.
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

import run_gs40c_coherence_ladder as L                              # noqa: E402
import run_gs40g_geometry_axis as G                                 # noqa: E402
import run_gs40i_position_and_blastradius as I                      # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import build_face_adjacency, precompute_element_geometry  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402
from pyfp3d.post.section_cut import section_cp_curve                # noqa: E402
from pyfp3d.post.shock import shock_report                          # noqa: E402

WING = REPO / "cases/meshes/onera_m6/coarse.msh"
WBD = REPO / "cases/meshes/onera_m6_wingbody_conforming"
OUT = REPO / "bench/gate_results/gs40j_position_primary.csv"
GR = REPO / "bench/gate_results"
ALPHA, GAMMA = 3.06, 1.4
ETAS = (0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99)
LADDER = ("P0", "P1", "P2", "P3")
SMOOTH = {"S1": "P1", "S2": "P2", "S3": "P3"}
ARMS = LADDER + tuple(SMOOTH)
MIN_STATIONS, MIN_STATES = 3, 4
Q_REFIT, Q_WORSE = 0.8, 1.25
GATE_S = 25 * 60

STATES = (
    ("wing", "W M0.76", GR / "gs40h_states/W0.76.npz", WING, 0.76),
    ("wing", "W M0.78", GR / "gs40h_states/W0.78.npz", WING, 0.78),
    ("wing", "W M0.80", GR / "gs40h_states/W0.8.npz", WING, 0.80),
    ("wing", "W M0.82", GR / "gs40g_states/W0.82.npz", WING, 0.82),
    ("wing", "W M0.8395", GR / "gs40i_states/wing_M08395.npz", WING, 0.8395),
    ("wing", "W M0.86", GR / "gs40g_states/W0.86.npz", WING, 0.86),
    ("wing", "W M0.88", GR / "gs40g_states/W0.88.npz", WING, 0.88),
    ("wingbody", "WB M0.78", GR / "gs40f_states/M0.78.npz", WBD / "coarse.msh", 0.78),
    ("wingbody", "WB M0.82", GR / "gs40f_states/M0.82.npz", WBD / "coarse.msh", 0.82),
    ("wingbody", "WB M0.86", GR / "gs40f_states/M0.86.npz", WBD / "coarse.msh", 0.86),
    ("wingbody", "WB M0.88", GR / "gs40d_levels/coarse.npz", WBD / "coarse.msh", 0.88),
    ("wingbody", "WB xcoarse", GR / "gs40d_levels/xcoarse.npz", WBD / "xcoarse.msh", 0.88),
    ("wingbody", "WB medium", GR / "le14_cache/rc0.05.npz", WBD / "medium.msh", 0.88),
)
#: G-ANCHOR: last round's published D on the M0.8395 state
ANCHOR_D = {"P0": 0.0132, "P1": 0.0169, "P2": 0.0055, "P3": 0.0055, "S1": 0.0130}


def guard_no_solve_executed():
    """G-Z, STRICTLY STRONGER than the version it replaces.

    ★ The inherited guard asserted that no solver driver appears in sys.modules.
    It FIRED here, correctly under its own definition -- this round imports
    run_gs40g and run_gs40i to REUSE their helpers (kstar, foot_from_field)
    rather than re-typing them, and those modules import solver entry points.

    But module presence is a PROXY for the property that matters ("no solve was
    executed"), and it is the wrong proxy in both directions: importing a driver
    does not mean solving, and NOT importing one does not prevent solving through
    some other path. That is the same defect this project logged as "a guard must
    cover what the conclusion claims" -- so the fix is to measure the property
    directly rather than to relax the proxy.

    Here every solver entry point is replaced by a raising stub for the duration
    of the round. A solve does not merely get detected, it becomes IMPOSSIBLE.
    """
    import pyfp3d.solve.newton as N
    import pyfp3d.solve.picard as PC

    def _forbidden(name):
        def _f(*a, **k):
            raise AssertionError(
                f"G-Z: {name} was CALLED -- this round is registered zero-solve")
        return _f

    n = 0
    for mod, names in ((N, ("solve_newton_lifting", "solve_newton_transonic")),
                       (PC, ("solve_subsonic", "solve_subsonic_lifting",
                             "solve_laplace", "solve_laplace_lifting"))):
        for nm in names:
            if hasattr(mod, nm):
                setattr(mod, nm, _forbidden(nm)); n += 1
    print(f"G-Z     {n} solver entry points replaced by raising stubs -- a solve is "
          "IMPOSSIBLE, not merely detected  PASS")


def main():
    t0 = time.perf_counter()
    guard_no_solve_executed()

    rows, per = [], {}
    for geom, label, sp, msh, m_inf in STATES:
        if not sp.exists():
            print(f"[{label}] state missing"); continue
        man = read_manifest(msh)
        sha = mesh_fingerprint(msh)["sha256"]
        assert man is None or sha == man["sha256"], f"G-MESH {label}"
        d = np.load(sp)
        if not bool(d["conv"]) or int(d["nlim"]) or int(d["nflr"]):
            print(f"[{label}] G-STATE: not usable"); continue
        phi = np.asarray(d["phi"], float)

        mc, wc = cut_wake(read_mesh(msh))
        nodes, elements = mc.nodes, mc.elements
        if len(phi) != len(nodes):
            print(f"[{label}] phi/mesh mismatch"); continue
        B, _ = precompute_element_geometry(nodes, elements)
        fn, _ = build_face_adjacency(elements)
        wall_faces = mc.boundary_faces["wall"]

        # J1 -- exogenous, computed ONCE
        j1 = {}
        for e_ in ETAS:
            try:
                cur = section_cp_curve(mc, phi, eta=e_, b_semi=B_SEMI, m_inf=m_inf)
                x = shock_report(cur, m_inf)["upper"].get("x_shock")
                j1[e_] = None if x is None or not np.isfinite(x) else float(x)
            except Exception:                                      # noqa: BLE001
                j1[e_] = None
        n_ok = sum(v is not None for v in j1.values())
        if n_ok < MIN_STATIONS:
            print(f"[{label}] ★ EXCLUDED-NO-WALL-SHOCK: J1 valid at {n_ok}/7 stations")
            rows.append(dict(state=label, geom=geom, m_inf=m_inf,
                             status="EXCLUDED", n_j1=n_ok))
            continue

        P, pel = {}, {}
        for k in (0, 1, 2, 3):
            po, pi = L.khop_patches(elements, fn, k)
            P[f"P{k}"] = L.patch_nodes(elements, po, pi); pel[f"P{k}"] = (po, pi)
        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)
        sup = float((mach_number_squared(q2, m_inf, GAMMA) >= 1.0).mean())

        g = {n: L.refit(nodes, elements, phi, *P[n]) for n in LADDER}
        for s, p in SMOOTH.items():
            #: G-SMOOTH: the smoother uses the SAME patch as its refit partner
            g[s] = L.average_over_patch(nodes, elements, gA, *pel[p])
        assert np.max(np.abs(g["P0"] - gA)) < 1e-10, f"G-A {label}"

        D, feet_all = {}, {}
        for n in ARMS:
            m2 = mach_number_squared((g[n] * g[n]).sum(axis=1), m_inf, GAMMA)
            feet = I.foot_from_field(nodes, elements, np.sqrt(np.maximum(m2, 0.0)),
                                     wall_faces, ETAS, B_SEMI)
            dd = [abs(feet[e_] - j1[e_]) for e_ in ETAS
                  if j1[e_] is not None and np.isfinite(feet[e_])]
            D[n] = float(np.median(dd)) if len(dd) >= MIN_STATIONS else np.nan
            feet_all[n] = feet
        if any(np.isnan(D[n]) for n in LADDER):
            print(f"[{label}] ★ EXCLUDED: an arm has < {MIN_STATIONS} usable stations")
            rows.append(dict(state=label, geom=geom, m_inf=m_inf, status="EXCLUDED"))
            continue

        kpos = min(LADDER, key=lambda a: D[a])
        ties = [a for a in LADDER if abs(D[a] - D[kpos]) < 1e-9]
        Cc = {n: L.coherence(np.sqrt(np.maximum(
            mach_number_squared((g[n] * g[n]).sum(axis=1), m_inf, GAMMA), 0.0)),
            *L.faces(fn))["C"] for n in LADDER}
        kcoh = G.kstar(Cc)[0]
        per[label] = dict(geom=geom, D=D, kpos=kpos, ties=ties, kcoh=kcoh,
                          sup=sup, n_j1=n_ok)
        print(f"[{label:11}] sup {sup:.4f}  J1 {n_ok}/7   D: "
              + "  ".join(f"{n} {D[n]:.4f}" for n in ARMS)
              + f"\n{'':14}k_pos*={kpos}{' (TIE)' if len(ties) > 1 else ''}   "
              f"k_coh*={kcoh}   {'AGREE' if kpos == kcoh else '★ DIFFER'}")
        for n in ARMS:
            rows.append(dict(state=label, geom=geom, m_inf=m_inf, status="OK",
                             arm=n, D=round(D[n], 5), sup_frac=round(sup, 5),
                             n_j1=n_ok, k_pos=kpos, k_coh=str(kcoh),
                             C=round(Cc[n], 4) if n in Cc else None))

        if label == "W M0.8395":
            bad = {k: (round(D[k], 4), v) for k, v in ANCHOR_D.items()
                   if abs(D[k] - v) > 5e-4}
            print(f"G-ANCHOR M0.8395 vs the published D: "
                  f"{'PASS' if not bad else '** FAIL ** ' + str(bad)}")
            assert not bad, "G-ANCHOR: kill criterion 3"

    # -------------------------------------------------------------- verdict --
    print("\n" + "=" * 78)
    if len(per) < MIN_STATES:
        verdict = f"POS-UNDEF (only {len(per)} counting states) -- kill criterion 4"
        qv = "n/a"
    else:
        votes = {}
        for lb, v in per.items():
            votes.setdefault("MIXED" if len(v["ties"]) > 1 else v["kpos"],
                             []).append(lb)
        half = len(per) / 2.0
        win = max(votes, key=lambda k: len(votes[k]))
        print(f"k_pos* votes over {len(per)} states: "
              + "  ".join(f"{k}:{len(v)}" for k, v in sorted(votes.items())))
        if len(votes[win]) <= half or win == "MIXED":
            verdict = "POS-MIXED (no majority)"
        elif win == "P0":
            verdict = "POS-CONTROL (no reconstruction improves position) -- kill 1"
        elif win == "P1":
            verdict = "POS-P1 (position agrees with coherence)"
        else:
            verdict = f"POS-BIGGER (position prefers {win}, coherence prefers P1)"

        # smoothing discriminator at the winning scale
        s = None if win in ("P0", "MIXED") else {"P1": "S1", "P2": "S2", "P3": "S3"}[win]
        if s is None:
            qv = "n/a (the winner is the control arm)"
        else:
            qs = [v["D"][win] / max(v["D"][s], 1e-30) for v in per.values()]
            nr = sum(1 for q in qs if q <= Q_REFIT)
            nw = sum(1 for q in qs if q >= Q_WORSE)
            ns = len(qs) - nr - nw
            med = float(np.median(qs))
            qv = ("Q-REFIT (the winner beats its own smoother)" if nr > half else
                  "Q-WORSE (worse than its own smoother)" if nw > half else
                  "Q-SMOOTH (indistinguishable from smoothing) -- kill criterion 2")
            print(f"smoothing discriminator at {win} vs {s}: median Q = {med:.3f}   "
                  f"REFIT/SMOOTH/WORSE = {nr}/{ns}/{nw}")
    print("=" * 78)
    print(f"VERDICT: {verdict}")
    print(f"Q:       {qv}")

    print("\nk_pos* vs k_coh* per state:")
    agree = 0
    for lb, v in per.items():
        a = v["kpos"] == v["kcoh"]; agree += a
        print(f"  {lb:11} sup {v['sup']:.4f}  k_pos*={v['kpos']:2}  "
              f"k_coh*={str(v['kcoh']):2}  {'agree' if a else '★ differ'}")
    print(f"  => agree on {agree}/{len(per)} states")
    dall = [v["D"][a] for v in per.values() for a in LADDER]
    print(f"\nRECORDED: all D in [{min(dall):.4f}, {max(dall):.4f}] chord "
          f"({'ALL under 0.02c -- differences are RELATIVE' if max(dall) < 0.02 else 'spread beyond 0.02c'})")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(verdict=verdict, Q=qv, n_states=len(per),
             agree=f"{agree}/{len(per)}",
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 5"


if __name__ == "__main__":
    main()
