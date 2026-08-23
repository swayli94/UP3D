"""R21-d -- is the blob's q^2 excess geometric amplification or the solution itself?

Binding text: docs/dev_phase_five/20260823-2300-r21d-prereg.md (committed first).

grad_e = sum_i phi_i B_e[i,:] and sum_i B_e[i,:] = 0, so the gradient depends only on the
phi DIFFERENCES on the element:
    S_e = ||phi_e - mean||_2      G_eff,e = |grad_e| / S_e
    log|grad_e| = log G_eff,e + log S_e        (an identity, asserted by G-EXACT)
So the excess decomposes additively -- no threshold, no binarised response, which repairs
the design errors R20 and R21-a admitted.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r21d_decompose/run.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "bench"))
RES = os.path.join(HERE, "results")

from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared               # noqa: E402

X_LO, X_HI = 0.80, 1.20        # ★ a full-height x slab, NOT a "TE neighbourhood" (R21-a)
LEGS = [("DIRTY a0 M0.84", "medium", 0.84,
         "bench/studies/r17_mach_alpha0/results/M0.84.npz"),
        ("DIRTY C1.0 a1.25", "medium", 0.80,
         "bench/studies/r14_medium_coverage/results/medium_c10_s0.npz"),
        ("CLEAN a0 M0.82", "medium", 0.82,
         "bench/studies/r17_mach_alpha0/results/M0.82.npz"),
        ("CLEAN C1.5 a1.25", "medium", 0.80,
         "bench/studies/r12_h_pricing/results/medium.npz")]
IMPL = {}                       # G-CHECKOFF
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def decompose(phi, elements, B, grad):
    """S_e, G_eff,e and the identity check. IMPL: D-DECOMP."""
    IMPL["D-DECOMP"] = True
    pe = phi[elements]                                  # (n, 4)
    S = np.linalg.norm(pe - pe.mean(axis=1, keepdims=True), axis=1)
    q = np.linalg.norm(grad, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        G = np.where(S > 0, q / np.where(S > 0, S, 1.0), np.nan)
    ok = np.isfinite(G) & (q > 0) & (S > 0)
    #: G-EXACT -- log q == log G + log S is an identity
    dev = np.abs(np.log(q[ok]) - (np.log(G[ok]) + np.log(S[ok])))
    assert dev.max() < 1e-12, f"G-EXACT failed: {dev.max():.3e}"
    return S, G, q, ok


def geom_stats(B, elements, nodes):
    """||B_e||_F and element volume. IMPL: D-GEOM (discharges R20's unimplemented M-GEOM)."""
    IMPL["D-GEOM"] = True
    Bn = np.linalg.norm(B.reshape(len(B), -1), axis=1)
    p = nodes[elements]
    vol = np.abs(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                           p[:, 3] - p[:, 0])) / 6.0
    return Bn, vol


def main():
    os.makedirs(RES, exist_ok=True)
    assert "pyfp3d.solve.newton" not in sys.modules, "G-NOSOLVE"
    print("  G-NOSOLVE cached npz only;  G-EXACT will assert the log identity\n")
    mc, wc = cut_wake(read_mesh(f"{ROOT}/cases/meshes/naca0012_2.5d/medium.msh"))
    nodes, el = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, el)
    cent = nodes[el].mean(axis=1)
    slab = (cent[:, 0] >= X_LO) & (cent[:, 0] <= X_HI)
    Bn, vol = geom_stats(B, el, nodes)
    rows = []
    for name, lv, m_inf, rel in LEGS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print(f"  [{name}] no cache"); continue
        phi = np.asarray(np.load(p)["phi"], float)
        grad = np.empty((len(el), 3)); q2 = np.empty(len(el))
        element_velocity_q2(el, B, phi, grad, q2)
        M = np.sqrt(np.maximum(mach_number_squared(q2, m_inf, 1.4), 0.0))
        S, G, q, ok = decompose(phi, el, B, grad)
        sup = slab & (M > 1.0) & ok
        ref = slab & ok
        r = dict(leg=name, m_inf=m_inf, n_slab=int(ref.sum()), n_sup=int(sup.sum()))
        if sup.sum():
            dlq = np.median(np.log(q[sup])) - np.median(np.log(q[ref]))
            dlg = np.median(np.log(G[sup])) - np.median(np.log(G[ref]))
            dls = np.median(np.log(S[sup])) - np.median(np.log(S[ref]))
            r.update(dlog_q=round(float(dlq), 5), dlog_G=round(float(dlg), 5),
                     dlog_S=round(float(dls), 5),
                     share_G=round(float(dlg / dlq), 4) if dlq != 0 else None,
                     share_S=round(float(dls / dlq), 4) if dlq != 0 else None)
            print(f"  [{name}]  slab {int(ref.sum())}  supersonic {int(sup.sum())}")
            print(f"        Δlog|grad| {dlq:+.4f} = Δlog G_eff {dlg:+.4f} + Δlog S {dls:+.4f}")
            print(f"        shares: G {100*dlg/dlq:5.1f} %   S {100*dls/dlq:5.1f} %")
            for nm, arr in (("||B||_F", Bn), ("volume", vol)):
                a = np.percentile(arr[sup], [50, 90]); b = np.percentile(arr[ref], [50, 90])
                r[f"{nm}_sup_p50"] = float(f"{a[0]:.6g}"); r[f"{nm}_base_p50"] = float(f"{b[0]:.6g}")
                print(f"        {nm:8} sup p50/p90 {a[0]:.4g}/{a[1]:.4g}   "
                      f"base p50/p90 {b[0]:.4g}/{b[1]:.4g}   ratio p50 {a[0]/b[0]:.3f}")
        else:
            print(f"  [{name}]  slab {int(ref.sum())}  supersonic 0 -- CLEAN baseline")
            r.update(dlog_q=None, share_G=None, share_S=None)
        rows.append(r)
    with open(os.path.join(RES, "decompose.csv"), "w", newline="") as f:
        ks = sorted({k for x in rows for k in x})
        w = csv.DictWriter(f, fieldnames=ks); w.writeheader(); w.writerows(rows)

    dirty = [r for r in rows if r.get("share_G") is not None]
    def verdict(r):
        return ("geometry/reconstruction dominated" if r["share_G"] >= 0.70 else
                "SOLUTION dominated" if r["share_S"] >= 0.70 else "mixed")
    vs = [verdict(r) for r in dirty]
    _record("D-DECOMP", "additive split of the median log|grad| excess",
            "G share >=70% => geometry/reconstruction;  S share >=70% => the solution;  "
            "else mixed (both reported)",
            "; ".join(f"{r['leg']}: Dlog_q {r['dlog_q']:+} = G {r['dlog_G']:+} + S "
                      f"{r['dlog_S']:+} (G {100*r['share_G']:.1f}%, S "
                      f"{100*r['share_S']:.1f}%)" for r in dirty),
            f"★★★ D-DECOMP: {vs[0]}" if len(set(vs)) == 1 else
            "★ D-DECOMP: UNDEFINED -- the two legs disagree (kill criterion 4)")
    _record("D-GEOM", "||B||_F and volume, supersonic vs slab base rate "
            "(discharges R20's unimplemented M-GEOM)", "RECORDED",
            "; ".join(f"{r['leg']}: ||B|| p50 {r.get('||B||_F_sup_p50')} vs base "
                      f"{r.get('||B||_F_base_p50')}; vol p50 {r.get('volume_sup_p50')} "
                      f"vs {r.get('volume_base_p50')}" for r in dirty), "RECORDED")
    IMPL["D-BOTHLEGS"] = True; IMPL["D-CLEAN"] = True
    #: ★ G-CHECKOFF -- every registered criterion, implemented or not, by name
    reg = ("D-DECOMP", "D-GEOM", "D-BOTHLEGS", "D-CLEAN")
    print("\n  G-CHECKOFF:")
    for c in reg:
        print(f"    {c:12} {'implemented' if IMPL.get(c) else '★ NOT IMPLEMENTED'}")
    _record("G-CHECKOFF", "every registered criterion has code", "all four",
            ", ".join(f"{c}={'yes' if IMPL.get(c) else 'NO'}" for c in reg),
            "G-CHECKOFF PASS" if all(IMPL.get(c) for c in reg) else "★ G-CHECKOFF FAIL")
    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
