"""R21-a -- is the TE supersonic set one blob or salt-and-pepper? With a null.

Binding text: phases/p5/docs/dev_phase_five/20260823-2100-r21a-prereg.md (committed first).

2146 of 6652 band elements is 32.3% occupancy and the 3-D site-percolation threshold is
about 0.31, so a random set of the same size would very likely connect too -- "it is one
blob" says nothing without a same-size null. The primary statistic is therefore the
boundary-to-volume ratio, which is not sitting on the threshold; component counting is
secondary and declared in advance to be weak.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r21a_coherence/run.py
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
from pyfp3d.mesh.metrics import build_face_adjacency, precompute_element_geometry  # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared               # noqa: E402

TE_LO, TE_HI = 0.80, 1.20
SEED, N_NULL = 20260823, 200           # G-SEED
LEGS = [("DIRTY a0 M0.84", "medium", 0.84,
         "bench/studies/r17_mach_alpha0/results/M0.84.npz"),
        ("DIRTY C1.0 a1.25", "medium", 0.80,
         "bench/studies/r14_medium_coverage/results/medium_c10_s0.npz")]
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def pv_ratio(sel, fn, pool):
    """boundary faces per selected element. Faces to OUTSIDE the pool count as boundary
    too, so the statistic is comparable between the measured set and the null (both are
    subsets of the same pool)."""
    idx = np.where(sel)[0]
    if idx.size == 0:
        return np.nan
    nb = fn[idx]                                  # (k, 4) face neighbours, -1 = none
    inside = np.zeros(len(fn), bool); inside[idx] = True
    ext = (nb < 0) | ~inside[np.clip(nb, 0, None)]
    return float(ext.sum() / idx.size)


def components(sel, fn):
    idx = np.where(sel)[0]
    if idx.size == 0:
        return 0, 0.0
    pos = -np.ones(len(fn), np.int64); pos[idx] = np.arange(idx.size)
    parent = np.arange(idx.size)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for k, e in enumerate(idx):
        for nb in fn[e]:
            if nb >= 0 and pos[nb] >= 0:
                ra, rb = find(k), find(pos[nb])
                if ra != rb:
                    parent[ra] = rb
    roots = np.array([find(i) for i in range(idx.size)])
    _, cnt = np.unique(roots, return_counts=True)
    return int(cnt.size), float(cnt.max() / idx.size)


def main():
    os.makedirs(RES, exist_ok=True)
    assert "pyfp3d.solve.newton" not in sys.modules, "G-NOSOLVE"
    print(f"  G-NOSOLVE cached npz only;  G-SEED seed={SEED}, {N_NULL} null draws\n")
    meshes, rows = {}, []
    for name, lv, m_inf, rel in LEGS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print(f"  [{name}] no cache"); continue
        if lv not in meshes:
            mc, wc = cut_wake(read_mesh(f"{ROOT}/cases/meshes/naca0012_2.5d/{lv}.msh"))
            B, _ = precompute_element_geometry(mc.nodes, mc.elements)
            fn, _ = build_face_adjacency(mc.elements)
            meshes[lv] = (mc, wc, B, fn, mc.nodes[mc.elements].mean(axis=1))
        mc, wc, B, fn, cent = meshes[lv]
        d = np.load(p)
        g = np.empty((len(mc.elements), 3)); q2 = np.empty(len(mc.elements))
        element_velocity_q2(mc.elements, B, np.asarray(d["phi"], float), g, q2)
        M = np.sqrt(np.maximum(mach_number_squared(q2, m_inf, 1.4), 0.0))
        band = (cent[:, 0] >= TE_LO) & (cent[:, 0] <= TE_HI)
        sup = band & (M > 1.0)
        pool = np.where(band)[0]; k = int(sup.sum())
        if not k:
            continue
        #: G-NULL -- same pool, same size
        rng = np.random.default_rng(SEED)
        pv_null, comp_null, frac_null = [], [], []
        for _ in range(N_NULL):
            r = np.zeros(len(fn), bool)
            r[rng.choice(pool, size=k, replace=False)] = True
            pv_null.append(pv_ratio(r, fn, pool))
            c, f_ = components(r, fn)
            comp_null.append(c); frac_null.append(f_)
        pv = pv_ratio(sup, fn, pool)
        nc, fmax = components(sup, fn)
        q = np.percentile(pv_null, [5, 50, 95])
        qc = np.percentile(comp_null, [5, 50, 95])
        qf = np.percentile(frac_null, [5, 50, 95])
        xs = cent[sup, 0]; ys = np.abs(cent[sup, 1])
        rows.append(dict(leg=name, level=lv, m_inf=m_inf, n_pool=len(pool), n_sup=k,
                         occupancy=round(k / len(pool), 4),
                         pv=round(pv, 4), pv_p5=round(q[0], 4), pv_p50=round(q[1], 4),
                         pv_p95=round(q[2], 4), n_comp=nc,
                         comp_p5=qc[0], comp_p50=qc[1], comp_p95=qc[2],
                         frac_max=round(fmax, 4), frac_p5=round(qf[0], 4),
                         frac_p50=round(qf[1], 4), frac_p95=round(qf[2], 4),
                         x_lo=round(float(xs.min()), 4), x_hi=round(float(xs.max()), 4),
                         y_max=round(float(ys.max()), 4)))
        r_ = rows[-1]
        print(f"  [{name}]  pool {len(pool)}  sup {k}  occupancy {r_['occupancy']}")
        print(f"        P/V      measured {pv:.4f}   null p5/p50/p95 "
              f"{q[0]:.4f}/{q[1]:.4f}/{q[2]:.4f}")
        print(f"        n_comp   measured {nc:6d}   null {qc[0]:.0f}/{qc[1]:.0f}/{qc[2]:.0f}")
        print(f"        max frac measured {fmax:.4f}   null {qf[0]:.4f}/{qf[1]:.4f}/{qf[2]:.4f}",
              flush=True)

    with open(os.path.join(RES, "coherence.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    verdicts = []
    for r in rows:
        v = ("field structure" if r["pv"] < r["pv_p5"] else
             "MORE scattered than chance" if r["pv"] > r["pv_p95"] else
             "indistinguishable from a same-size random set")
        verdicts.append(v)
    _record("C-PV", "boundary/volume of the supersonic set vs a same-size null",
            "below null p5 => field structure;  inside [p5,p95] => 'one blob' unsupported;"
            "  above p95 => more scattered than chance",
            "; ".join(f"{r['leg']}: P/V {r['pv']} vs null [{r['pv_p5']}, {r['pv_p95']}]"
                      f" -> {v}" for r, v in zip(rows, verdicts)),
            f"★ C-PV: {verdicts[0]}" if len(set(verdicts)) == 1 else
            "★ C-PV: UNDEFINED -- the two legs disagree (kill criterion 4)")
    _record("C-COMP", "components and largest-component share vs the null",
            "declared in advance to be WEAK -- 32% occupancy sits on the 3-D percolation "
            "threshold, so this is secondary",
            "; ".join(f"{r['leg']}: {r['n_comp']} comps (null {r['comp_p5']:.0f}-"
                      f"{r['comp_p95']:.0f}), max frac {r['frac_max']} (null "
                      f"{r['frac_p5']}-{r['frac_p95']})" for r in rows), "RECORDED")
    _record("C-EXT", "extent of the set", "RECORDED",
            "; ".join(f"{r['leg']}: x/c {r['x_lo']}-{r['x_hi']}, max |y| {r['y_max']}"
                      for r in rows), "RECORDED")
    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
