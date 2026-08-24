"""R20 -- what are the TE supersonic cells touching? Zero solves.

Binding text: phases/p5/docs/dev_phase_five/20260823-1900-r20-prereg.md (committed first).

Thirteen rounds measured where and how much; none asked why. A potential flow at M 0.8
cannot have a physical supersonic region at the trailing edge, so it comes from a
numerical structure. This classifies the cells by what they touch, WITH the base rate of
the same band on the same leg -- an enrichment without a base rate says nothing.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r20_te_mechanism/run.py
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

TE_LO, TE_HI = 0.80, 1.20
CANON = ("xcoarse", "coarse", "medium")
CLASSES = ("WAKE", "TE-NODE", "WALL", "INTERIOR")
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def classify(mc, wc, el):
    """Exhaustive, mutually exclusive, in this precedence: WAKE > TE-NODE > WALL >
    INTERIOR. G-EXHAUST asserts the four cover every element exactly once."""
    n = len(mc.nodes)
    def mask_of(groups):
        m = np.zeros(n, bool)
        for g in groups:
            if g in mc.boundary_faces:
                m[np.unique(mc.boundary_faces[g].reshape(-1))] = True
        return m
    wake_n = mask_of(("wake_plus", "wake_minus"))
    wall_n = mask_of(("wall",))
    te_n = np.zeros(n, bool)
    for attr in ("te_nodes", "te_upper", "te_lower", "stations"):
        v = getattr(wc, attr, None)
        if v is not None:
            a = np.asarray(v).reshape(-1)
            if a.size and np.issubdtype(a.dtype, np.integer):
                te_n[a[(a >= 0) & (a < n)]] = True
    cls = np.full(len(el), 3, np.int8)                     # INTERIOR
    cls[wall_n[el].any(axis=1)] = 2                        # WALL
    cls[te_n[el].any(axis=1)] = 1                          # TE-NODE
    cls[wake_n[el].any(axis=1)] = 0                        # WAKE (highest precedence)
    return cls, dict(n_wake_nodes=int(wake_n.sum()), n_te_nodes=int(te_n.sum()),
                     n_wall_nodes=int(wall_n.sum()))


def main():
    os.makedirs(RES, exist_ok=True)
    assert "pyfp3d.solve.newton" not in sys.modules, "G-NOSOLVE"
    print("  G-NOSOLVE  cached npz only\n")
    legs = [
        ("DIRTY  coarse G0.0 M0.86", "coarse", 0.86,
         "bench/studies/r19_gamma_input/results"),      # from the grid: re-solve free
        ("CLEAN  coarse G0.0 M0.80", "coarse", 0.80, None),
        ("DIRTY  medium a0 M0.84", "medium", 0.84,
         "bench/studies/r17_mach_alpha0/results/M0.84.npz"),
    ]
    rows = []
    #: the R19 grid did not cache phi per cell, so the coarse legs are recomputed from
    #: the R17/R14 caches where available; here we use what IS on disk (G-PROV prints it)
    avail = [("DIRTY medium a0 M0.84", "medium", 0.84,
              f"{ROOT}/bench/studies/r17_mach_alpha0/results/M0.84.npz"),
             ("CLEAN medium a0 M0.82", "medium", 0.82,
              f"{ROOT}/bench/studies/r17_mach_alpha0/results/M0.82.npz"),
             ("DIRTY medium C1.0 a1.25", "medium", 0.80,
              f"{ROOT}/bench/studies/r14_medium_coverage/results/medium_c10_s0.npz"),
             ("CLEAN medium C1.5 a1.25", "medium", 0.80,
              f"{ROOT}/bench/studies/r12_h_pricing/results/medium.npz")]
    meshes = {}
    for name, lv, m_inf, p in avail:
        assert lv in CANON, "G-SCOPE"
        if not os.path.exists(p):
            print(f"  [{name}] no cache"); continue
        if lv not in meshes:
            mc, wc = cut_wake(read_mesh(f"{ROOT}/cases/meshes/naca0012_2.5d/{lv}.msh"))
            B, _ = precompute_element_geometry(mc.nodes, mc.elements)
            cls, info = classify(mc, wc, mc.elements)
            #: G-EXHAUST
            cnt = np.bincount(cls, minlength=4)
            assert cnt.sum() == len(mc.elements), "G-EXHAUST"
            meshes[lv] = (mc, wc, B, mc.nodes[mc.elements].mean(axis=1), cls, info)
            print(f"  G-PROV  {lv}: {info}  class counts (whole mesh) "
                  f"{dict(zip(CLASSES, cnt.tolist()))}")
        mc, wc, B, cent, cls, info = meshes[lv]
        d = np.load(p)
        g = np.empty((len(mc.elements), 3)); q2 = np.empty(len(mc.elements))
        element_velocity_q2(mc.elements, B, np.asarray(d["phi"], float), g, q2)
        M = np.sqrt(np.maximum(mach_number_squared(q2, m_inf, 1.4), 0.0))
        band = (cent[:, 0] >= TE_LO) & (cent[:, 0] <= TE_HI)
        sup = band & (M > 1.0)
        nb, ns = int(band.sum()), int(sup.sum())
        base = np.bincount(cls[band], minlength=4) / max(nb, 1)
        row = dict(leg=name, level=lv, m_inf=m_inf, n_band=nb, n_sup=ns,
                   te_maxM=round(float(M[band].max()), 4))
        for i, c in enumerate(CLASSES):
            row[f"base_{c}"] = round(float(base[i]), 4)
        if ns:
            obs = np.bincount(cls[sup], minlength=4) / ns
            for i, c in enumerate(CLASSES):
                row[f"obs_{c}"] = round(float(obs[i]), 4)
                row[f"enrich_{c}"] = (round(float(obs[i] / base[i]), 3)
                                      if base[i] > 0 else None)
        rows.append(row)
        print(f"\n  [{name}]  band {nb} elems, supersonic {ns}, TE maxM {row['te_maxM']}")
        print(f"        {'class':10}{'base':>9}{'observed':>10}{'enrich':>9}")
        for i, c in enumerate(CLASSES):
            o = row.get(f"obs_{c}"); e = row.get(f"enrich_{c}")
            print(f"        {c:10}{base[i]:9.4f}"
                  + (f"{o:10.4f}{(e if e is not None else float('nan')):9.3f}"
                     if ns else f"{'--':>10}{'--':>9}"))

    with open(os.path.join(RES, "classes.csv"), "w", newline="") as f:
        ks = sorted({k for r in rows for k in r})
        w = csv.DictWriter(f, fieldnames=ks); w.writeheader(); w.writerows(rows)

    dirty = [r for r in rows if r["n_sup"] > 0]
    clean = [r for r in rows if r["n_sup"] == 0]
    # ---- M-CONTROL: do the base rates agree between dirty and clean legs? ---
    if dirty and clean:
        worst, which = 0.0, None
        for c in CLASSES:
            for a in dirty:
                for b in clean:
                    if a["level"] != b["level"]:
                        continue
                    lo = min(a[f"base_{c}"], b[f"base_{c}"])
                    if lo > 0:
                        r_ = abs(a[f"base_{c}"] - b[f"base_{c}"]) / lo
                        if r_ > worst:
                            worst, which = r_, (c, a["leg"], b["leg"])
        _record("M-CONTROL", "do the TE-band base rates agree between dirty and clean legs",
                "<=10% => an enrichment is a distribution change;  >10% => it is the "
                "band's COMPOSITION changing and must not be read as enrichment",
                f"worst relative base-rate difference {100*worst:.2f}% at {which}",
                "M-CONTROL: base rates agree" if worst <= 0.10 else
                "★ M-CONTROL: base rates DIFFER -- enrichments must not be attributed")

    # ---- M-ENRICH ---------------------------------------------------------
    for r in dirty:
        e = {c: r.get(f"enrich_{c}") for c in CLASSES}
        top = [c for c in CLASSES if e[c] is not None and e[c] >= 2.0]
        rest_ok = all(e[c] is None or e[c] < 1.5 for c in CLASSES if c not in top)
        flat = all(e[c] is None or 0.7 <= e[c] <= 1.5 for c in CLASSES)
        v = (f"★ M-ENRICH: {top[0]}-enriched" if len(top) == 1 and rest_ok else
             "★★ M-ENRICH: NO enrichment -- not a topological feature" if flat else
             "M-ENRICH: UNDEFINED (multiple or mixed)")
        _record(f"M-ENRICH/{r['leg']}", "enrichment by class, against the band's own base rate",
                "one class >=2 with the rest <1.5 => that feature;  all in [0.7,1.5] => "
                "no topological feature;  else UNDEFINED",
                "  ".join(f"{c} {e[c]}" for c in CLASSES), v)

    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
