"""What does the level-set LE alpha-dependence actually LOOK like, against conforming?

The alpha discriminator found the LE-band R degrading with alpha on BOTH level-set 3-D
geometries (m6wing 0.825 -> 0.988 = +18 %, wingbody 0.934 -> 1.184 = +24 %) while every
conforming case stayed flat within the instrument's ~10 %. R is a ratio of two norms, so the
thing to draw is not the raw Cp but the TWO REFINEMENT DIFFERENCES whose norms form it:

    d1(x) = Cp(coarse) - Cp(xcoarse)      the first refinement's movement
    d2(x) = Cp(medium) - Cp(coarse)       the second refinement's movement

Converging means d2 is visibly smaller than d1. R ~ 1 means the second refinement moved the
answer as much as the first -- the picture of a solution that is not settling.

Two figures per geometry:
  A  d1 and d2 over the LE band at one station, 2x2 = path x alpha. The contrast is the point.
  B  per-STATION LE-band d1/d2/R, so the alpha-dependence is located spanwise rather than
     only known as a pooled number. If it concentrates at the root station the wake-junction
     classes are implicated; if it is uniform, the sheet's global influence is.

Also writes the section-Cp CSVs per (geom, path, alpha, level) -- the next-phase reference
material, and the raw data behind both figures.

Outputs (TRACKED): bench/gate_results/capability/le_alpha_<geom>_curves.png
                   bench/gate_results/capability/le_alpha_<geom>_spanwise.png
                   bench/gate_results/capability/le_alpha_<geom>_cp.csv
                   bench/gate_results/le_alpha_stations.csv
"""

import csv
import math
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.post.unified import section_cp                          # noqa: E402

OUT = os.path.join(HERE, "gate_results")
ART = os.path.join(OUT, "capability")
os.makedirs(ART, exist_ok=True)

M_INF = 0.50
ALPHAS = (0.0, 3.06)
ETAS = (0.20, 0.44, 0.65, 0.80, 0.90)
LE_LO, LE_HI, NQ = 0.0, 0.15, 400
XQ = np.linspace(LE_LO, LE_HI, NQ)

GEOMS = [
    ("m6wing", [("conforming", "onera_m6", cap.conf_wing),
                ("level-set", "onera_m6_wakefree", cap.ls_wing)],
     ("xcoarse_ss", "coarse_ss", "medium")),
    ("wingbody", [("conforming", "onera_m6_wingbody_conforming", cap.conf_wingbody),
                  ("level-set", "onera_m6_wingbody", cap.ls_wingbody)],
     ("xcoarse", "coarse", "medium")),
]


def resample(sec, side):
    x, cp = np.asarray(sec[f"x_{side}"]), np.asarray(sec[f"cp_{side}"])
    o = np.argsort(x)
    return np.interp(XQ, x[o], cp[o])


def main():
    station_rows = []
    for geom, paths, levels in GEOMS:
        #: cp[(path, alpha, level)][eta][side] -> Cp resampled on XQ
        cp = {}
        for path, mdir, fn in paths:
            for a in ALPHAS:
                for lv in levels:
                    mp = os.path.join(REPO, "cases", "meshes", mdir, f"{lv}.msh")
                    if not os.path.exists(mp):
                        print(f"  missing {mp}"); continue
                    try:
                        mesh, op, r, phi, mvop = fn(mp, M_INF, a)
                    except Exception as exc:                       # noqa: BLE001
                        print(f"  {geom}/{path} a{a} {lv}: "
                              f"{type(exc).__name__}: {exc}", flush=True)
                        continue
                    kw = (dict(phi=phi) if mvop is None
                          else dict(mvop=mvop, phi_ext=phi))
                    m_eff = M_INF if mvop is None else r.get("m_final", M_INF)
                    d = {}
                    for eta in ETAS:
                        sec = section_cp(mesh, eta=eta, b_semi=cap.B_SEMI,
                                         m_inf=m_eff, **kw)
                        d[eta] = {s: resample(sec, s) for s in ("upper", "lower")}
                    cp[(path, a, lv)] = d
                    print(f"  {geom}/{path} a{a} {lv} done", flush=True)

        # ---- CSV of every resampled LE curve (the raw data behind both figures) ----
        p = os.path.join(ART, f"le_alpha_{geom}_cp.csv")
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["path", "alpha", "level", "eta", "side", "x_over_c", "cp"])
            for (path, a, lv), d in cp.items():
                for eta, sides in d.items():
                    for s, v in sides.items():
                        for xv, cv in zip(XQ, v):
                            w.writerow([path, a, lv, f"{eta:.2f}", s,
                                        f"{xv:.5f}", f"{cv:.6f}"])
        print(f"wrote {p}")

        lx, lc, lm = levels

        def diffs(path, a, eta, side):
            k = [(path, a, l) for l in levels]
            if not all(kk in cp for kk in k):
                return None, None
            return (cp[k[1]][eta][side] - cp[k[0]][eta][side],
                    cp[k[2]][eta][side] - cp[k[1]][eta][side])

        # ---------------- figure A: the two difference curves ------------------
        fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0), sharex=True)
        ETA_SHOW = 0.44
        for i, (path, _m, _f) in enumerate(paths):
            for j, a in enumerate(ALPHAS):
                ax = axes[i][j]
                d1, d2 = diffs(path, a, ETA_SHOW, "upper")
                if d1 is None:
                    ax.text(0.5, 0.5, "missing", ha="center"); continue
                n1, n2 = (math.sqrt(float(np.trapz(d1 ** 2, XQ))),
                          math.sqrt(float(np.trapz(d2 ** 2, XQ))))
                ax.plot(XQ, d1, lw=1.8, color="tab:blue",
                        label=f"d1 = coarse − xcoarse   ‖·‖={n1:.4f}")
                ax.plot(XQ, d2, lw=1.8, color="tab:red",
                        label=f"d2 = medium − coarse   ‖·‖={n2:.4f}")
                ax.axhline(0, color="k", lw=0.6)
                ax.set_title(f"{path}   α = {a}   →  R = {n2 / n1:.3f}"
                             if n1 else f"{path} α={a}", fontsize=10.5)
                ax.grid(alpha=0.3)
                ax.legend(fontsize=8, loc="lower right")
                if i == 1:
                    ax.set_xlabel("x/c  (LE band)")
                if j == 0:
                    ax.set_ylabel(r"$\Delta C_p$  (upper surface)")
        fig.suptitle(f"{geom}: what the LE alpha-dependence IS — the two refinement "
                     f"movements at η={ETA_SHOW}, M{M_INF}\n"
                     "converging = red visibly below blue; R→1 = the second refinement "
                     "moved the answer as much as the first", fontsize=11)
        fig.tight_layout()
        fa = os.path.join(ART, f"le_alpha_{geom}_curves.png")
        fig.savefig(fa, dpi=130); plt.close(fig)
        print(f"wrote {fa}")

        # ---------------- figure B: spanwise localisation ----------------------
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
        for j, a in enumerate(ALPHAS):
            ax = axes[j]
            for path, _m, _f in paths:
                Rs = []
                for eta in ETAS:
                    tot1 = tot2 = 0.0
                    for side in ("upper", "lower"):
                        d1, d2 = diffs(path, a, eta, side)
                        if d1 is None:
                            continue
                        tot1 += float(np.trapz(d1 ** 2, XQ))
                        tot2 += float(np.trapz(d2 ** 2, XQ))
                    n1, n2 = math.sqrt(tot1), math.sqrt(tot2)
                    R = n2 / n1 if n1 else float("nan")
                    Rs.append(R)
                    station_rows.append(dict(geom=geom, path=path, alpha=a, eta=eta,
                                             d1=round(n1, 8), d2=round(n2, 8),
                                             R=round(R, 5)))
                ax.plot(ETAS, Rs, "-o", lw=1.8, ms=6, label=path)
            r_first = 0.5 if geom == "m6wing" else 1.071
            ax.axhline(r_first, color="g", ls="--", lw=1.2,
                       label=f"first order (R={r_first})")
            ax.axhline(1.0 if geom == "m6wing" else 1.810, color="r", ls=":", lw=1.2,
                       label="falsification ceiling")
            ax.set_title(f"α = {a}", fontsize=11)
            ax.set_xlabel(r"$\eta = z/b_{semi}$"); ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        axes[0].set_ylabel("LE-band R  per station")
        fig.suptitle(f"{geom}: is the LE alpha-dependence localised spanwise? "
                     f"(root η=0.20 → tip η=0.90)", fontsize=11)
        fig.tight_layout()
        fb = os.path.join(ART, f"le_alpha_{geom}_spanwise.png")
        fig.savefig(fb, dpi=130); plt.close(fig)
        print(f"wrote {fb}")

    if station_rows:
        p = os.path.join(OUT, "le_alpha_stations.csv")
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(station_rows[0]))
            w.writeheader(); w.writerows(station_rows)
        print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
