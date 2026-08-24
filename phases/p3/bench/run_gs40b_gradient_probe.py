"""
Registered item 1 -- the ZERO-SOLVE kill test for interior gradient reconstruction.

Pre-registration: phases/p4/docs/dev_phase_four/20260817-0300-gradient-probe-prereg.md
Addendum #1:      phases/p4/docs/dev_phase_four/20260817-0400-gradient-probe-addendum1.md

THE QUESTION. Hold a CONVERGED phi fixed and swap ONLY the recovery operator.
Does the heavy tail of the single-face Mach jump collapse -- and does it collapse
on the SLIVER elements without taking the SHOCK with it?

WHY IT CAN BE THIS CHEAP. phi is fixed, so the only independent variable is the
recovery operator. That is the clean single-variable experiment this project
could never construct on the mesh side (every mesh knob was measured acting
outside its name); here it is constructive.

★ WHAT THIS CANNOT SHOW (the forbidden sentence, from the registration): it
answers "can a patch reconstruction produce a cleaner Mach FIELD from the same
phi", NOT "does using it inside the solver produce a better phi". The first is a
NECESSARY condition for the second, never a sufficient one. No result here may
be written as "interior gradient reconstruction works".

ARMS (four scored + one recorded):
    A   per-element constant gradient          -- the control, = today
    B1  linear patch LSQ on the node 1-ring    -- addendum #1
    B2  quadratic patch LSQ on the node 1-ring -- addendum #1
    S   volume-weighted node-averaged gradient -- ★ a DELIBERATE fake
    B0  linear LSQ on the face-neighbour patch -- RECORDED, the registered form

★ Arm S is the design's core. The blocking conclusion (GS1b.10 §9.4) says a
sliver's LOCAL Mach signature IS a shock's, so ANY operator that flattens both
shrinks the tail. S is the operator that only flattens. If B1/B2 cannot be told
apart from S, the route is smoothing, not reconstruction, and it dies.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.coloring import node_to_element_csr                # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import (build_face_adjacency,              # noqa: E402
                                 compute_aspect_ratios,
                                 compute_tet_volumes,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402

MESH = REPO / "cases/meshes/onera_m6_wingbody_conforming/medium.msh"
STATE = REPO / "bench/gate_results/le14_cache/rc0.05.npz"
OUT = REPO / "bench/gate_results/gs40b_gradient_probe.csv"
M_INF, GAMMA = 0.88, 1.4

#: registered thresholds -- NOT to be moved after seeing a number
R_SLIVER_PASS, R_SHOCK_KEEP, R_RATIO, S_INDIST = 0.5, 0.8, 0.6, 0.20
AR_PCT, SHOCK_MIN_FACES, MIN_N_FACES = 99.0, 50, 100
QUAD_MIN_PTS = 10          # 3-D quadratic has 10 coefficients
G_CONS_BAR = 1e-12         # the REGISTERED bar (prereg 1.2.3); B2 exempt per addendum #2
GATE_S = 20 * 60           # kill criterion 4


# ---------------------------------------------------------------- guards ----
def guard_zero_solve():
    """G-Z: this script may not import a solver driver.

    ★ The forbidden list is ASSEMBLED AT RUNTIME rather than written as literals,
    because a guard that scans source for its own forbidden words matches itself
    -- that has now happened three times in this project (`pgrep -f`, the G-Z
    forbidden-word list, the load-average guard). Here we inspect sys.modules,
    not source text, so there is nothing to match.
    """
    stem = "solve_newton" if True else ""
    banned = {stem + "_lifting", stem + "_transonic",
              "solve" + "_subsonic", "solve" + "_transonic_lifting"}
    hit = [m for m in sys.modules
           if any(b in dir(sys.modules[m]) for b in banned
                  if sys.modules[m] is not None)]
    assert not hit, f"G-Z: a solver driver was imported: {hit}"
    return sorted(banned)


def guard_mesh_identity():
    """G-MESH: pin the mesh by its GS4.0 manifest, not by its filename."""
    man = read_manifest(MESH)
    assert man is not None, f"G-MESH: no manifest beside {MESH} -- run GS4.0's stamping"
    live = mesh_fingerprint(MESH)
    assert live["sha256"] == man["sha256"], (
        f"G-MESH: {MESH} no longer matches its committed manifest "
        f"({live['sha256'][:12]} vs {man['sha256'][:12]}) -- a mesh file is part of "
        "the provenance, exactly like the thread count")
    return man


def guard_state(d):
    """G-Z (state half): the probe is only meaningful on a genuine solution."""
    assert bool(d["conv"]), "G-Z: the cached state is NOT converged"
    assert int(d["nlim"]) == 0 and int(d["nflr"]) == 0, (
        f"G-Z: clamped state ({int(d['nlim'])}/{int(d['nflr'])}) -- GS1.4 says a "
        "clamped state is not a converged flow")


# ------------------------------------------------------------ recovery ------
def arm_A(nodes, elements, B, phi):
    """The control: today's per-element constant gradient."""
    g = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
    element_velocity_q2(elements, B, phi, g, q2)
    return g


def _patch_lsq(nodes, elements, phi, patch_off, patch_idx, order):
    """Least-squares fit of phi over each NODE's patch; return the node gradients.

    order=1 fits a + g.d ; order=2 fits a + g.d + 1/2 d'Hd and returns (g, H).
    Coordinates are scaled by the patch radius so the normal equations stay
    conditioned; the gradient is unscaled on the way out.

    ★ Registered in advance (addendum #1 §2.2): a node whose patch has fewer than
    QUAD_MIN_PTS points falls back to LINEAR, and the fallback fraction is
    reported per partition. Not a 2-ring: a 2-ring grows the patch most exactly
    where the slivers are, which would manufacture the effect being measured.
    """
    n_node = len(nodes)
    g = np.zeros((n_node, 3)); H = np.zeros((n_node, 3, 3))
    fell_back = np.zeros(n_node, dtype=bool)
    for n in range(n_node):
        idx = patch_idx[patch_off[n]:patch_off[n + 1]]
        if len(idx) < 4:
            continue
        d = nodes[idx] - nodes[n]
        h = np.max(np.abs(d))
        if h <= 0.0:
            continue
        ds = d / h
        f = phi[idx] - phi[n]
        use_quad = order == 2 and len(idx) >= QUAD_MIN_PTS
        if use_quad:
            X = np.column_stack([ds, 0.5 * ds[:, 0] ** 2, 0.5 * ds[:, 1] ** 2,
                                 0.5 * ds[:, 2] ** 2, ds[:, 0] * ds[:, 1],
                                 ds[:, 0] * ds[:, 2], ds[:, 1] * ds[:, 2]])
            c, _res, rank, _sv = np.linalg.lstsq(X, f, rcond=None)
            #: ★ addendum #2: point COUNT is not the same test as RANK. A patch can
            #: hold 10+ points and still be rank-deficient for the quadratic (nearly
            #: coplanar 1-rings are common next to slivers -- exactly where this probe
            #: looks), and lstsq then returns the MINIMUM-NORM solution, which spreads
            #: a genuine linear gradient into H. Caught by G-CONS on a synthetic linear
            #: field BEFORE any measurement: 6.5e-02 against a 1e-8 bar.
            #: Addendum #1 already registered "fall back to linear when the quadratic
            #: is underdetermined" -- rank deficiency IS underdetermination, so this
            #: fixes the implementation to match the registered intent rather than
            #: changing the rule. Strictly conservative: it can only make B2 *less*
            #: different from B1.
            if rank < X.shape[1]:
                use_quad = False
        if order == 2 and not use_quad:
            fell_back[n] = True
        if not use_quad:
            c, *_ = np.linalg.lstsq(ds, f, rcond=None)
        g[n] = c[:3] / h
        if use_quad:
            H[n] = np.array([[c[3], c[6], c[7]],
                             [c[6], c[4], c[8]],
                             [c[7], c[8], c[5]]]) / (h * h)
    return g, H, fell_back


def _node_patch(elements, n_node):
    """Node 1-ring: for each node, the distinct nodes of its incident elements."""
    off, ne = node_to_element_csr(elements)
    starts = np.zeros(n_node + 1, dtype=np.int64)
    buf = []
    for n in range(n_node):
        els = ne[off[n]:off[n + 1]]
        u = np.unique(elements[els].reshape(-1)) if len(els) else np.empty(0, np.int64)
        buf.append(u)
        starts[n + 1] = starts[n] + len(u)
    return starts, (np.concatenate(buf) if buf else np.empty(0, np.int64))


def arm_B(nodes, elements, phi, patch, order):
    """Node-patch LSQ (linear or quadratic), evaluated at the element centroid."""
    g_n, H_n, fb = _patch_lsq(nodes, elements, phi, patch[0], patch[1], order)
    cent = nodes[elements].mean(axis=1)
    g_e = np.zeros((len(elements), 3))
    for k in range(4):
        nk = elements[:, k]
        d = cent - nodes[nk]
        g_e += g_n[nk] + (np.einsum("eij,ej->ei", H_n[nk], d) if order == 2 else 0.0)
    return g_e / 4.0, fb


def arm_S(nodes, elements, g_A):
    """★ THE DELIBERATE FAKE: volume-weighted node average of the CONSTANT element
    gradients, pushed back to elements. It never refits phi -- it only flattens.
    """
    vol = compute_tet_volumes(nodes, elements)
    n_node = len(nodes)
    num = np.zeros((n_node, 3)); den = np.zeros(n_node)
    for k in range(4):
        np.add.at(num, elements[:, k], g_A * vol[:, None])
        np.add.at(den, elements[:, k], vol)
    g_n = num / np.maximum(den, 1e-300)[:, None]
    return g_n[elements].mean(axis=1)


def arm_B0(nodes, elements, phi, face_nb):
    """RECORDED: the registered face-neighbour linear patch (addendum #1 §2.3)."""
    cent = nodes[elements].mean(axis=1)
    g_e = np.zeros((len(elements), 3))
    for e in range(len(elements)):
        nb = face_nb[e][face_nb[e] >= 0]
        idx = np.unique(np.concatenate([elements[e], elements[nb].reshape(-1)])
                        if len(nb) else elements[e])
        d = nodes[idx] - cent[e]
        h = np.max(np.abs(d))
        if h <= 0 or len(idx) < 4:
            continue
        #: ★ addendum #2: an EXPLICIT intercept. The first version centred phi on
        #: `phi[idx].mean()` while centring d on the element CENTROID -- two different
        #: origins, so the constant term was not eliminated and a linear field came
        #: back with a 4.0e-01 gradient error. G-CONS caught it on synthetic data
        #: before any measurement. Fitting c0 explicitly makes the origin irrelevant.
        X = np.column_stack([np.ones(len(idx)), d / h])
        c, *_ = np.linalg.lstsq(X, phi[idx], rcond=None)
        g_e[e] = c[1:] / h
    return g_e


# ------------------------------------------------------------ scoring -------
def face_list(face_nb):
    """Internal faces as unique (e, e') pairs -- the sample set, fixed once."""
    e = np.repeat(np.arange(len(face_nb)), 4)
    o = face_nb.reshape(-1)
    m = (o >= 0) & (o > e)
    return e[m], o[m]


def connected_components(fa, fb, sel, face_nb):
    """Sizes of the connected components of the SELECTED faces.

    Faces are connected when they share an element. A shock is a SURFACE; a
    sliver is a POINT -- that is the only discriminator available here that does
    not lean on a local statistic, which the blocking conclusion says cannot work.
    """
    from collections import defaultdict, deque
    idx = np.flatnonzero(sel)
    by_elem = defaultdict(list)
    for i in idx:
        by_elem[fa[i]].append(i)
        by_elem[fb[i]].append(i)
    seen, sizes = set(), {}
    for i in idx:
        if i in seen:
            continue
        q, comp = deque([i]), []
        seen.add(i)
        while q:
            j = q.popleft(); comp.append(j)
            for e in (fa[j], fb[j]):
                for k in by_elem[e]:
                    if k not in seen:
                        seen.add(k); q.append(k)
        for j in comp:
            sizes[j] = len(comp)
    out = np.zeros(len(fa), dtype=np.int64)
    for j, s in sizes.items():
        out[j] = s
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="coarse mesh + its own cached state is NOT available; "
                         "this only shrinks the printout")
    args = ap.parse_args()
    t0 = time.perf_counter()

    banned = guard_zero_solve()
    print(f"G-Z  no solver driver imported (checked {len(banned)} names)")
    man = guard_mesh_identity()
    print(f"G-MESH  {MESH.name}  sha {man['sha256'][:12]}  "
          f"{man['n_nodes']} nodes / {man['n_tets']} tets  PASS")

    d = np.load(STATE)
    guard_state(d)
    phi = np.asarray(d["phi"], dtype=np.float64)
    print(f"G-Z  state {STATE.name}: conv={bool(d['conv'])} "
          f"lim/flr={int(d['nlim'])}/{int(d['nflr'])} |R|={float(d['res']):.3e}  PASS")

    mc, _wc = cut_wake(read_mesh(MESH))
    nodes, elements = mc.nodes, mc.elements
    assert len(phi) == len(nodes), (
        f"phi has {len(phi)} entries but the cut mesh has {len(nodes)} nodes")
    print(f"cut mesh: {len(nodes)} nodes / {len(elements)} tets  "
          f"({time.perf_counter() - t0:.1f} s)")

    B, _ = precompute_element_geometry(nodes, elements)
    face_nb, _ = build_face_adjacency(elements)
    ar = compute_aspect_ratios(nodes, elements)
    fa, fb = face_list(face_nb)
    print(f"internal faces: {len(fa)}")

    # ---- arms -------------------------------------------------------------
    g = {}
    g["A"] = arm_A(nodes, elements, B, phi)
    # G-A: the control must BIT-reproduce the library kernel
    gg = np.empty_like(g["A"]); qq = np.empty(len(elements))
    element_velocity_q2(elements, B, phi, gg, qq)
    assert np.array_equal(g["A"], gg), "G-A: control arm does not reproduce the kernel"
    print("G-A  control bit-reproduces element_velocity_q2  PASS")

    patch = _node_patch(elements, len(nodes))
    psz = np.diff(patch[0])
    print(f"node 1-ring: min {psz.min()} / p05 {np.percentile(psz, 5):.0f} / "
          f"median {np.median(psz):.0f} / max {psz.max()}; "
          f"below {QUAD_MIN_PTS}: {100 * (psz < QUAD_MIN_PTS).mean():.2f}%")

    g["B1"], _ = arm_B(nodes, elements, phi, patch, 1)
    g["B2"], fb_quad = arm_B(nodes, elements, phi, patch, 2)
    g["S"] = arm_S(nodes, elements, g["A"])
    g["B0"] = arm_B0(nodes, elements, phi, face_nb)
    print(f"arms built ({time.perf_counter() - t0:.1f} s)")

    # ---- G-CONS: linear reproduction --------------------------------------
    #: the REGISTERED bar is 1e-12 (the code briefly carried a laxer 1e-8; corrected
    #: to the registration). ★ B2 is exempt from the ASSERT and is not scored --
    #: addendum #2: it fails at the registered bar because that bar is below the
    #: floating-point floor of a quadratic LSQ (max cond 7.2e8 => cond*eps = 1.6e-7,
    #: measured worst 2.4e-7, median error 8.9e-16). I do NOT relax the bar to a
    #: cond-aware one, because that is moving a threshold after seeing the number.
    #: B2 therefore stays in the printout and the CSV with scored=False.
    rng = np.random.default_rng(0)
    a = rng.normal(size=3)
    gcons = {}
    for name, order in (("B1", 1), ("B2", 2), ("B0", 0)):
        lin = nodes @ a
        gl = (arm_B0(nodes, elements, lin, face_nb) if name == "B0"
              else arm_B(nodes, elements, lin, patch, order)[0])
        err = float(np.max(np.abs(gl - a)))
        gcons[name] = err
        scored = name != "B2"
        tag = ("PASS" if err < G_CONS_BAR else "** FAIL **") + (
            "" if scored else "  [NOT SCORED -- addendum #2]")
        print(f"G-CONS {name}: linear reproduction max err {err:.3e}  {tag}")
        if scored:
            assert err < G_CONS_BAR, (
                f"G-CONS: {name} does not reproduce a linear field "
                f"({err:.3e} >= {G_CONS_BAR}) -- kill criterion 2")

    # ---- Mach fields and face jumps ---------------------------------------
    M = {k: np.sqrt(mach_number_squared(np.einsum("ej,ej->ej", v, v).sum(axis=1),
                                        M_INF, GAMMA))
         for k, v in g.items()}
    J = {k: np.abs(m[fa] - m[fb]) for k, m in M.items()}

    # ---- partitions: computed ONCE on the control (G-P) --------------------
    ar_thr = np.percentile(ar, AR_PCT)
    sliver = (ar[fa] >= ar_thr) | (ar[fb] >= ar_thr)
    MA = M["A"]
    crosses = ((MA[fa] >= 1.0) & (MA[fb] < 1.0)) | ((MA[fb] >= 1.0) & (MA[fa] < 1.0))
    comp = connected_components(fa, fb, crosses, face_nb)
    shock = crosses & (comp >= SHOCK_MIN_FACES)
    bulk = ~sliver & ~shock
    parts = {"SLIVER": sliver, "SHOCK": shock, "BULK": bulk}
    print("\npartitions (on arm A, shared by every arm -- G-P):")
    for p, m in parts.items():
        print(f"  {p:7} n_faces {int(m.sum()):>8}  ({100 * m.mean():5.2f}%)")
    print(f"  AR p{AR_PCT:.0f} = {ar_thr:.2f}; sonic-crossing faces {int(crosses.sum())}, "
          f"of which in components >= {SHOCK_MIN_FACES}: {int(shock.sum())}")

    # 8th question: does any partition make a verdict unreachable?
    undefined = [p for p, m in parts.items() if int(m.sum()) < MIN_N_FACES]
    if undefined:
        print(f"\n★ UNDEFINED partitions (n < {MIN_N_FACES}): {undefined} "
              "-- the registration says report, do not relax")

    # ---- ratios -----------------------------------------------------------
    rows = []
    print(f"\np99 of |dM| per partition, and R = p99(arm)/p99(A):")
    hdr = f"{'partition':9}" + "".join(f"{k:>12}" for k in ("A", "B1", "B2", "S", "B0"))
    print(hdr)
    Rs = {}
    for p, m in parts.items():
        if not m.any():
            continue
        base = np.percentile(J["A"][m], 99)
        line = f"{p:9}{base:12.5f}"
        Rs[p] = {}
        for k in ("B1", "B2", "S", "B0"):
            v = np.percentile(J[k][m], 99)
            r = v / base if base > 0 else float("nan")
            Rs[p][k] = r
            line += f"{r:11.3f}x"
            rows.append(dict(partition=p, arm=k, n_faces=int(m.sum()),
                             p99_A=round(float(base), 6),
                             p99_arm=round(float(v), 6), R=round(float(r), 4)))
        print(line)

    # quadratic fallback per partition (addendum #1 §2.2)
    print("\nB2 quadratic fallback-to-linear, per partition (addendum #1):")
    for p, m in parts.items():
        if not m.any():
            continue
        nds = np.unique(elements[np.unique(np.concatenate([fa[m], fb[m]]))].reshape(-1))
        print(f"  {p:7} {100 * fb_quad[nds].mean():5.2f}% of its nodes")

    # ---- verdict ----------------------------------------------------------
    print("\n" + "=" * 72)
    verdicts = {}
    #: ★ The first version short-circuited to UNDEFINED whenever ANY partition was
    #: missing. That is NOT what the registration says, and the difference matters:
    #:   G-REAL  needs SLIVER *and* SHOCK           -> unreachable without SHOCK
    #:   G-NULL  needs SLIVER only
    #:   G-SMEAR is a DISJUNCTION -- "R(SHOCK) < 0.8  OR  |R_B - R_S|/R_S < 20%" --
    #:           and its second branch needs only SLIVER and the S arm.
    #: The registration's own rationale for that branch never mentions the shock:
    #: "if B cannot be told apart from S, reconstruction bought nothing, it is just
    #: smoothing." So evaluating it with SHOCK undefined follows the registered text.
    #: This is code being corrected to match the registration, NOT a threshold move:
    #: every number in the conditions below is the registered one.
    for k in ("B1",):          # ★ B2 not scored -- addendum #2
        if "SLIVER" not in Rs:
            verdicts[k] = "UNDEFINED (SLIVER partition empty)"
            print(f"{k}: UNDEFINED -- SLIVER empty")
            continue
        rsl = Rs["SLIVER"][k]
        rsh = Rs["SHOCK"][k] if "SHOCK" in Rs else None
        s_gap = abs(rsl - Rs["SLIVER"]["S"]) / max(Rs["SLIVER"]["S"], 1e-30)
        shock_bad = (rsh is not None) and rsh < R_SHOCK_KEEP
        if rsl > R_SHOCK_KEEP:
            v = "G-NULL  (the heavy tail did not collapse)"
        elif (rsl <= R_SLIVER_PASS and rsh is not None
              and rsh >= R_SHOCK_KEEP and rsl <= R_RATIO * rsh):
            v = "G-REAL  (tail collapsed on slivers, shock kept)"
        elif rsl <= R_SLIVER_PASS and (shock_bad or s_gap < S_INDIST):
            v = "G-SMEAR (it is smoothing, not reconstruction)"
        elif rsh is None:
            v = "UNDEFINED (SHOCK empty; G-REAL unreachable, G-SMEAR not triggered)"
        else:
            v = "G-MIX   (no direction claimed)"
        verdicts[k] = v
        print(f"{k}: R(SLIVER)={rsl:.3f}  "
              f"R(SHOCK)={'UNDEFINED' if rsh is None else f'{rsh:.3f}'}  "
              f"|R_B-R_S|/R_S={s_gap:.3f}  (S_INDIST={S_INDIST})\n    =>  {v}")
    print("=" * 72)
    if "SHOCK" not in Rs:
        print("★ SHOCK was UNDEFINED, so the G-REAL branch was never testable: this "
              "round CANNOT say 'the shock is preserved'.\n"
              "  The registered 50-face threshold is NOT relaxed.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    import csv
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    meta = dict(mesh_sha=man["sha256"], state=STATE.name, m_inf=M_INF,
                ar_threshold=float(ar_thr), verdicts=verdicts,
                wall_s=round(time.perf_counter() - t0, 1))
    OUT.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"\nwrote {OUT.name} + .json;  total {time.perf_counter() - t0:.1f} s "
          f"(gate {GATE_S} s)")
    assert time.perf_counter() - t0 < GATE_S, "kill criterion 4: over the time gate"


if __name__ == "__main__":
    main()
