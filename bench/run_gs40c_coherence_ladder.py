"""
Registered item 1, re-opening measurement: sonic-set COHERENCE across a patch-size ladder.

Pre-registration: docs/dev_phase_four/20260817-0900-coherence-ladder-prereg.md
Previous round:   docs/dev_phase_four/20260817-0700-gradient-probe-verdict.md (G-SMEAR)

THE TARGET IS LAST ROUND'S FAILURE MODE. There, the SHOCK partition came back
EMPTY -- 18746 sonic-crossing faces in 6273 components, median 2, largest 47,
none >= 50 -- so the G-REAL branch was never testable. The deliverable here is to
make that test possible again.

THE AXIS. Face-hop distance k, with the METHOD held fixed at a linear refit,
because last round measured that changing the PATCH moves the answer 3.5x more
than changing the METHOD. Patch scale is therefore the single independent
variable.

THE CRITERION IS A SHAPE, NOT A LEVEL. C(k) with an interior maximum means small
patches are a genuine optimum; C monotone decreasing means last round's B0
reading was just "closer to the control" and this line closes.

★ THREE FORBIDDEN SENTENCES (carried + one new):
  1. a post-processing probe is a NECESSARY condition, never sufficient;
  2. zero solves => this says nothing about behaviour INSIDE the solver;
  3. ★ COHERENT IS NOT CORRECT -- there is no external truth for the sonic
     surface in cases/reference_data/, so coherence is a structural necessary
     condition and may never be read as "the shock is correctly located".
"""

import csv
import json
import sys
import time
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from pyfp3d.mesh.coloring import node_to_element_csr                # noqa: E402
from pyfp3d.mesh.manifest import mesh_fingerprint, read_manifest    # noqa: E402
from pyfp3d.mesh.metrics import (build_face_adjacency,              # noqa: E402
                                 compute_tet_volumes,
                                 precompute_element_geometry)
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2             # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared           # noqa: E402

MESH = REPO / "cases/meshes/onera_m6_wingbody_conforming/medium.msh"
CACHE = REPO / "bench/gate_results/le14_cache"
OUT = REPO / "bench/gate_results/gs40c_coherence_ladder.csv"
M_INF, GAMMA = 0.88, 1.4

DISCOVERY = "rc0.05"                    # the state the hypothesis came from
BINDING = ("rc0.025", "rc0.045")        # out-of-sample (another taper radius only)

#: registered -- not to be moved after seeing a number
COMP_MIN = 10          # C's definition (from last round's published table)
USABLE_MIN = 50        # D-USABLE; last round's registered value, NOT the measured 47
PEAK_FACTOR = 1.30     # C-PEAK: interior max >= this x the endpoints
FLAT_RATIO = 1.15      # below this spread => C-FLAT => UNDEFINED
MIN_CROSS = 100        # below this many crossing faces => report only
GATE_S = 25 * 60       # kill criterion 5

#: G-REPRO: P1 must reproduce last round's arm B0 on the discovery state
REPRO = dict(n_cross=13295, max_comp=118, C=0.428)


# ------------------------------------------------------------------ guards --
def guard_zero_solve():
    """G-Z: no solver driver. Inspects sys.modules, never its own source, so it
    cannot match itself -- the self-reference trap that has fired three times."""
    banned = {"solve" + "_newton_lifting", "solve" + "_newton_transonic",
              "solve" + "_subsonic", "solve" + "_transonic_lifting"}
    hit = [m for m, mod in sys.modules.items()
           if mod is not None and any(b in dir(mod) for b in banned)]
    assert not hit, f"G-Z: a solver driver was imported: {hit}"


def guard_mesh():
    man = read_manifest(MESH)
    assert man is not None, "G-MESH: no manifest -- run GS4.0's stamping"
    assert mesh_fingerprint(MESH)["sha256"] == man["sha256"], "G-MESH: mesh moved"
    return man


def load_state(tag):
    """G-Z (state half). ★ The two le10 states are UNUSABLE and the reason is
    registered rather than the states being dropped silently: their npz carry no
    nlim/nflr, so the clamp assertion cannot be satisfied at all."""
    d = np.load(CACHE / f"{tag}.npz")
    assert bool(d["conv"]), f"{tag}: not converged"
    assert int(d["nlim"]) == 0 and int(d["nflr"]) == 0, f"{tag}: clamped"
    return np.asarray(d["phi"], dtype=np.float64)


# ------------------------------------------------------------- patches ------
def khop_patches(elements, face_nb, k):
    """Elements within k face-hops, as CSR. Mesh-only => built once, reused for
    every state. k=0 is the element itself."""
    n = len(elements)
    if k == 0:
        return np.arange(n + 1), np.arange(n)
    rows = np.repeat(np.arange(n), 4)
    cols = face_nb.reshape(-1)
    m = cols >= 0
    A = sp.coo_matrix((np.ones(m.sum(), np.int8), (rows[m], cols[m])),
                      shape=(n, n)).tocsr()
    A = A + sp.eye(n, dtype=np.int8, format="csr")
    R = A.copy()
    for _ in range(k - 1):
        R = (R @ A)
        R.data[:] = 1
    R = R.tocsr()
    return R.indptr, R.indices


def patch_nodes(elements, poff, pidx):
    """Distinct nodes of each element's patch (CSR)."""
    n = len(elements)
    starts = np.zeros(n + 1, dtype=np.int64)
    buf = []
    for e in range(n):
        u = np.unique(elements[pidx[poff[e]:poff[e + 1]]].reshape(-1))
        buf.append(u)
        starts[e + 1] = starts[e] + len(u)
    return starts, np.concatenate(buf)


def node_ring_patch(elements, n_node):
    """N1: the node 1-ring, per ELEMENT (union over the element's 4 nodes)."""
    off, ne = node_to_element_csr(elements)
    n = len(elements)
    starts = np.zeros(n + 1, dtype=np.int64)
    buf = []
    for e in range(n):
        els = np.concatenate([ne[off[v]:off[v + 1]] for v in elements[e]])
        u = np.unique(elements[np.unique(els)].reshape(-1))
        buf.append(u); starts[e + 1] = starts[e] + len(u)
    return starts, np.concatenate(buf)


# ------------------------------------------------------------- recovery -----
def refit(nodes, elements, phi, noff, nidx, cond_out=None):
    """Linear LSQ of phi over each element's patch, gradient at the centroid.

    ★ EXPLICIT INTERCEPT (addendum #2): the first version of this operator centred
    phi and d on two different origins and returned a 0.4 gradient error on a
    linear field. Fitting c0 makes the origin irrelevant.
    """
    cent = nodes[elements].mean(axis=1)
    g = np.zeros((len(elements), 3))
    kmax = 0.0
    for e in range(len(elements)):
        idx = nidx[noff[e]:noff[e + 1]]
        if len(idx) < 4:
            continue
        d = nodes[idx] - cent[e]
        h = np.max(np.abs(d))
        if h <= 0:
            continue
        X = np.column_stack([np.ones(len(idx)), d / h])
        if cond_out is not None:
            kmax = max(kmax, np.linalg.cond(X))
        c, *_ = np.linalg.lstsq(X, phi[idx], rcond=None)
        g[e] = c[1:] / h
    if cond_out is not None:
        cond_out.append(kmax)
    return g


def average_over_patch(nodes, elements, g_const, poff, pidx):
    """S1: volume-weighted AVERAGE of the constant element gradients over the SAME
    patch as P1 -- the pure smoother, kept scale-for-scale so 'refit vs average'
    is testable at every k."""
    vol = compute_tet_volumes(nodes, elements)
    g = np.zeros((len(elements), 3))
    for e in range(len(elements)):
        els = pidx[poff[e]:poff[e + 1]]
        w = vol[els]
        g[e] = (g_const[els] * w[:, None]).sum(axis=0) / w.sum()
    return g


# ------------------------------------------------------------- scoring ------
def faces(face_nb):
    e = np.repeat(np.arange(len(face_nb)), 4)
    o = face_nb.reshape(-1)
    m = (o >= 0) & (o > e)
    return e[m], o[m]


def coherence(M, fa, fb):
    """C = share of sonic-crossing faces sitting in components of >= COMP_MIN.

    A SHARE, not a count: an arm that simply produces more crossing faces must
    not score higher for that reason. Counts are printed alongside anyway (8th
    question: a share can hide a nearly empty numerator)."""
    cr = ((M[fa] >= 1) & (M[fb] < 1)) | ((M[fb] >= 1) & (M[fa] < 1))
    idx = np.flatnonzero(cr)
    if len(idx) == 0:
        return dict(n_cross=0, n_comp=0, max_comp=0, n_usable=0, C=float("nan"),
                    biggest=np.empty(0, np.int64))
    by = defaultdict(list)
    for i in idx:
        by[fa[i]].append(i); by[fb[i]].append(i)
    seen, comps = set(), []
    for i in idx:
        if i in seen:
            continue
        q, cur = deque([i]), []
        seen.add(i)
        while q:
            j = q.popleft(); cur.append(j)
            for el in (fa[j], fb[j]):
                for kk in by[el]:
                    if kk not in seen:
                        seen.add(kk); q.append(kk)
        comps.append(cur)
    sizes = np.array([len(c) for c in comps])
    big = comps[int(np.argmax(sizes))]
    return dict(n_cross=len(idx), n_comp=len(sizes), max_comp=int(sizes.max()),
                n_usable=int((sizes >= USABLE_MIN).sum()),
                C=float(sizes[sizes >= COMP_MIN].sum() / len(idx)),
                biggest=np.array(big))


def pca_shape(nodes, elements, fa, fb, big):
    """RECORDED only: is the largest component sheet-like, filament-like or a blob?"""
    if len(big) < 10:
        return None, None
    c = nodes[elements].mean(axis=1)
    p = 0.5 * (c[fa[big]] + c[fb[big]])
    p = p - p.mean(axis=0)
    ev = np.sort(np.linalg.eigvalsh(np.cov(p.T)))[::-1]
    return float(ev[1] / ev[0]), float(ev[2] / ev[0])


def main():
    t0 = time.perf_counter()
    guard_zero_solve()
    man = guard_mesh()
    print(f"G-Z   no solver driver imported            PASS")
    print(f"G-MESH {MESH.name} sha {man['sha256'][:12]}  PASS")

    mc, _ = cut_wake(read_mesh(MESH))
    nodes, elements = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, elements)
    face_nb, _ = build_face_adjacency(elements)
    fa, fb = faces(face_nb)
    print(f"mesh: {len(nodes)} nodes / {len(elements)} tets / {len(fa)} internal "
          f"faces  ({time.perf_counter() - t0:.1f} s)")

    # patches: mesh-only, built ONCE
    P = {}
    for k in (0, 1, 2, 3):
        po, pi = khop_patches(elements, face_nb, k)
        P[f"P{k}"] = patch_nodes(elements, po, pi)
        if k == 1:
            p1_elems = (po, pi)
        print(f"  patch P{k}: mean {np.diff(P[f'P{k}'][0]).mean():.1f} nodes "
              f"({time.perf_counter() - t0:.1f} s)")
    P["N1"] = node_ring_patch(elements, len(nodes))
    print(f"  patch N1: mean {np.diff(P['N1'][0]).mean():.1f} nodes "
          f"({time.perf_counter() - t0:.1f} s)")

    rows, results = [], {}
    for tag in (DISCOVERY,) + BINDING:
        phi = load_state(tag)
        assert len(phi) == len(nodes)
        gA = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
        element_velocity_q2(elements, B, phi, gA, q2)

        g = {}
        conds = []
        for name in ("P0", "P1", "P2", "P3", "N1"):
            g[name] = refit(nodes, elements, phi, *P[name],
                            cond_out=conds if name == "P1" else None)
        g["S1"] = average_over_patch(nodes, elements, gA, *p1_elems)

        if tag == DISCOVERY:
            # G-A: the k=0 refit is the tet's own linear function == the kernel
            err = np.max(np.abs(g["P0"] - gA))
            print(f"\nG-A   P0 (k=0 refit) vs element_velocity_q2: {err:.3e}"
                  f"  {'PASS' if err < 1e-10 else '** FAIL **'}")
            assert err < 1e-10, "G-A failed"
            # G-CONS with a CONDITION-DERIVED bar (defect #10's fix)
            a = np.random.default_rng(0).normal(size=3)
            kappa = conds[0]
            bar = max(1e-12, 10 * kappa * np.finfo(float).eps)
            for name in ("P1", "P2", "P3", "N1"):
                gl = refit(nodes, elements, nodes @ a, *P[name])
                e_lin = float(np.max(np.abs(gl - a)))
                ok = e_lin < bar
                print(f"G-CONS {name}: {e_lin:.3e}  bar {bar:.3e} "
                      f"(kappa {kappa:.2e}, measured not guessed)  "
                      f"{'PASS' if ok else '** FAIL **'}")
                assert ok, f"G-CONS {name}: kill criterion 4"

        res = {}
        for name in ("P0", "P1", "P2", "P3", "N1", "S1"):
            m2 = mach_number_squared((g[name] * g[name]).sum(axis=1), M_INF, GAMMA)
            r = coherence(np.sqrt(np.maximum(m2, 0.0)), fa, fb)
            l2l1, l3l1 = pca_shape(nodes, elements, fa, fb, r["biggest"])
            res[name] = r
            rows.append(dict(state=tag, arm=name,
                             mean_patch_nodes=round(float(np.diff(P[name][0]).mean()), 1)
                             if name in P else None,
                             n_cross=r["n_cross"], n_comp=r["n_comp"],
                             max_comp=r["max_comp"], n_ge50=r["n_usable"],
                             C=round(r["C"], 4),
                             pca_l2l1=None if l2l1 is None else round(l2l1, 3),
                             pca_l3l1=None if l3l1 is None else round(l3l1, 3)))
        results[tag] = res

        if tag == DISCOVERY:
            p1 = res["P1"]
            ok = (p1["n_cross"] == REPRO["n_cross"]
                  and p1["max_comp"] == REPRO["max_comp"]
                  and abs(p1["C"] - REPRO["C"]) < 5e-4)
            print(f"\nG-REPRO P1 vs last round's B0: n_cross {p1['n_cross']} "
                  f"(want {REPRO['n_cross']}), max {p1['max_comp']} "
                  f"(want {REPRO['max_comp']}), C {p1['C']:.4f} "
                  f"(want {REPRO['C']})  {'PASS' if ok else '** FAIL **'}")
            assert ok, "G-REPRO: kill criterion 3 -- instrument, not a finding"

        print(f"\n[{tag}]  {'arm':4}{'nodes':>8}{'cross':>8}{'comps':>7}"
              f"{'max':>6}{'>=50':>6}{'C':>8}{'l2/l1':>7}{'l3/l1':>7}")
        for name in ("P0", "P1", "P2", "P3", "N1", "S1"):
            r = res[name]
            l2, l3 = pca_shape(nodes, elements, fa, fb, r["biggest"])
            mp = np.diff(P[name][0]).mean() if name in P else float("nan")
            print(f"       {name:4}{mp:8.1f}{r['n_cross']:>8}{r['n_comp']:>7}"
                  f"{r['max_comp']:>6}{r['n_usable']:>6}{r['C']:8.3f}"
                  f"{(l2 if l2 else float('nan')):7.2f}{(l3 if l3 else float('nan')):7.2f}")

    # ------------------------------------------------------------ verdict ---
    print("\n" + "=" * 74)
    verdict = {}
    for tag in BINDING:
        C = {k: results[tag][k]["C"] for k in ("P0", "P1", "P2", "P3")}
        few = [k for k in C if results[tag][k]["n_cross"] < MIN_CROSS]
        if few or any(np.isnan(v) for v in C.values()):
            verdict[tag] = f"UNDEFINED (thin/absent sonic set: {few})"
        elif max(C.values()) / max(min(C.values()), 1e-30) < FLAT_RATIO:
            verdict[tag] = "C-FLAT => UNDEFINED (the ladder separated nothing)"
        elif max(C["P1"], C["P2"]) >= PEAK_FACTOR * max(C["P0"], C["P3"]):
            verdict[tag] = "C-PEAK  (small patch is a genuine interior optimum)"
        elif C["P0"] >= C["P1"] >= C["P2"] >= C["P3"]:
            verdict[tag] = "C-MONO-DOWN (B0 was 'closer to control' => line closes)"
        elif C["P0"] <= C["P1"] <= C["P2"] <= C["P3"]:
            verdict[tag] = "C-MONO-UP (opposite of last round; recorded, not claimed)"
        else:
            verdict[tag] = "C-MIX (no direction claimed)"
        print(f"{tag}: C(P0..P3) = " + " ".join(f"{C[k]:.3f}" for k in
              ("P0", "P1", "P2", "P3")) + f"   => {verdict[tag]}")
    # D-USABLE and the smoother check, on the binding legs
    for tag in BINDING:
        r = results[tag]
        best = max(("P0", "P1", "P2", "P3", "N1"), key=lambda k: r[k]["C"])
        gap = abs(r["P1"]["C"] - r["S1"]["C"]) / max(r["S1"]["C"], 1e-30)
        print(f"{tag}: D-USABLE best={best} components>={USABLE_MIN}: "
              f"{r[best]['n_usable']}  ({'MET' if r[best]['n_usable'] >= 1 else 'not met'})"
              f"   |C_P1-C_S1|/C_S1 = {gap:.3f}")
    print("=" * 74)
    print("★ COHERENT IS NOT CORRECT: there is no external truth for the sonic\n"
          "  surface, so this is a structural necessary condition only.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    OUT.with_suffix(".json").write_text(json.dumps(
        dict(mesh_sha=man["sha256"], verdict=verdict,
             wall_s=round(time.perf_counter() - t0, 1)), indent=2) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nwrote {OUT.name} + .json;  {dt:.1f} s (gate {GATE_S} s)")
    assert dt < GATE_S, "kill criterion 5: over the time gate"


if __name__ == "__main__":
    main()
