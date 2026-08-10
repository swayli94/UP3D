"""LE-2: is the M6 wing-alone's LE deficit localised at the ROOT?

LE-1 found the LE-band order p = 0.27-0.37 on the M6 wing-alone (both wake paths) while the
geometrically MORE complex wing-body is healthy there (p = 1.19-1.42) and the 2.5-D airfoil
sits at 0.75-0.85. Both 3-D cases share the tip cap, the 30-degree sweep, the taper and a 3-D
wake. The one thing only the wing-alone has is a root section TERMINATING ON THE SYMMETRY
PLANE -- the wing-body fairs that root into the fuselage.

So: refine LE-1's banding from x/c alone to x/c x eta and re-read. If the deficit concentrates
at the root stations, the root/symmetry-plane termination is the site. If it is spanwise
uniform, that hypothesis is dead and the remaining candidates are sweep / taper / LE-line
curvature.

Same instrument as LE-1 (band-integrated lift, no sampling, no section cutting), and the same
guard extended one level down: the eta bins of a chordwise band must sum to that band's own
total, which is itself bit-identical to the library's cl. A binning slip therefore cannot pass.

★ Read the eta profile, not a single number. The previous round's failure was reading a pooled
scalar whose station scatter was the size of the effect; here every bin is an INTEGRAL over its
own triangles, so there is no sampling noise -- but bins with few triangles or near-cancelling
d1 still give unreadable R, so the triangle count and d1 magnitude are reported alongside every
R and small-d1 bins are flagged rather than quoted.

Outputs (TRACKED): bench/gate_results/le2_spanwise.csv
                   bench/gate_results/capability/le2_spanwise.png
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
#: ★ archive-move fix (2026-08-10): `bench/gate_results/` STAYED at the repo's bench/
#: -- the 7 kept scripts write there and the capability boundary cites those CSVs by
#: path -- so an archived script must reach ACROSS to it, not look below itself.
_GATE = str(__import__('pathlib').Path(__file__).resolve().parents[3]
            / 'bench' / 'gate_results')
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le            # noqa: E402
from pyfp3d.post.surface import _pressure_force, planform_area      # noqa: E402
from pyfp3d.post.unified import (_cp_from_q2,                       # noqa: E402
                                 _conforming_wall_state,
                                 _d11_wall_state, wall_forces)

OUT = os.path.join(_GATE)
ART = os.path.join(OUT, "capability")
CSV = os.path.join(OUT, "le2_spanwise.csv")
M_INF, ALPHA = 0.50, 3.06
LE_LO, LE_HI = 0.0, 0.15
ETA_EDGES = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0001)
#: a bin holding too few triangles cannot be read at the coarsest level -- flagged, not quoted
MIN_TRI = 8

CASES = [
    ("m6wing", "conforming", "onera_m6", cap.conf_wing,
     ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015)),
    ("m6wing", "level-set", "onera_m6_wakefree", cap.ls_wing,
     ("xcoarse_ss", "coarse_ss", "medium"), (0.060, 0.030, 0.015)),
    ("wingbody", "conforming", "onera_m6_wingbody_conforming", cap.conf_wingbody,
     ("xcoarse", "coarse", "medium"), (0.044, 0.030, 0.015)),
    ("wingbody", "level-set", "onera_m6_wingbody", cap.ls_wingbody,
     ("xcoarse", "coarse", "medium"), (0.044, 0.030, 0.015)),
]


def le_bins(mesh, *, phi=None, mvop=None, phi_ext=None, alpha_deg=0.0,
            m_inf=0.0, s_ref=1.0, wall_tag="wall", u_inf=1.0):
    """Per-eta-bin lift contribution of the LE chordwise band, plus that band's total."""
    wall = np.asarray(mesh.boundary_faces[wall_tag], dtype=np.int64)
    if phi is not None:
        q2, _, area, n_out = _conforming_wall_state(mesh, phi, wall, u_inf, 0)
    else:
        q2, _, area, n_out = _d11_wall_state(mesh, mvop, phi_ext, wall, u_inf, 0)
    cp = _cp_from_q2(q2, m_inf, 1.4)
    cen = mesh.nodes[wall].mean(axis=1)
    z = cen[:, 2]
    xc = (cen[:, 0] - np.array([x_le(zi) for zi in z])) / \
         np.array([chord_at(zi) for zi in z])
    eta = z / B_SEMI
    in_le = (xc >= LE_LO) & (xc < LE_HI)
    out = {}
    for i in range(len(ETA_EDGES) - 1):
        lo, hi = ETA_EDGES[i], ETA_EDGES[i + 1]
        m = in_le & (eta >= lo) & (eta < hi)
        n = int(m.sum())
        if n == 0:
            out[(lo, hi)] = (0.0, 0)
            continue
        _cf, cl, _cd = _pressure_force(cp[m], area[m], n_out[m], s_ref, alpha_deg)
        out[(lo, hi)] = (cl, n)
    #: bins outside the eta range (z < 0 or eta >= 1.0001) so the guard stays exact
    m_extra = in_le & ~((eta >= ETA_EDGES[0]) & (eta < ETA_EDGES[-1]))
    extra = 0.0
    if np.any(m_extra):
        _cf, extra, _cd = _pressure_force(cp[m_extra], area[m_extra],
                                          n_out[m_extra], s_ref, alpha_deg)
    _cf, le_tot, _cd = _pressure_force(cp[in_le], area[in_le], n_out[in_le],
                                       s_ref, alpha_deg)
    return out, extra, le_tot, int(in_le.sum())


def main():
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    centres = [0.5 * (ETA_EDGES[i] + ETA_EDGES[i + 1])
               for i in range(len(ETA_EDGES) - 1)]
    for geom, path, mdir, fn, levels, hs in CASES:
        hx, hc, hm = hs
        r_first = (hm - hc) / (hc - hx)
        r_max = math.log(hm / hc) / math.log(hc / hx)
        print(f"\n=== {geom} / {path}  (p=1 -> R {r_first:.3f}, falsified >= "
              f"{r_max:.3f}) ===", flush=True)
        vals, ntri = {}, {}
        ok = True
        for lv in levels:
            mp = os.path.join(REPO, "cases", "meshes", mdir, f"{lv}.msh")
            if not os.path.exists(mp):
                print(f"  {lv}: mesh missing"); ok = False; break
            mesh, op, r, phi, mvop = fn(mp, M_INF, ALPHA)
            kw = (dict(phi=phi) if mvop is None else dict(mvop=mvop, phi_ext=phi))
            m_eff = M_INF if mvop is None else r.get("m_final", M_INF)
            sref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
            b, extra, le_tot, n_le = le_bins(mesh, alpha_deg=ALPHA, m_inf=m_eff,
                                             s_ref=sref, **kw)
            #: the guard, one level down from LE-1's
            s = sum(v for v, _n in b.values()) + extra
            err = abs(s - le_tot)
            if err > 1e-12:
                print(f"  ★ {lv}: eta-bin GUARD FAILED, sum-vs-band {err:.2e}",
                      flush=True)
            print(f"  {lv:11s} n_LE_tri={n_le:6d}  LE_total {le_tot:+.7f}  "
                  f"guard {err:.1e}", flush=True)
            print("             bins " + "  ".join(
                f"[{lo:.1f},{hi:.1f}):{v:+.6f}(n{n})" for (lo, hi), (v, n)
                in b.items()), flush=True)
            vals[lv], ntri[lv] = b, {k: n for k, (v, n) in b.items()}
        if not ok or len(vals) < 3:
            continue
        lx, lc, lm = levels
        Rs, flags = [], []
        for k in vals[lx]:
            d1 = vals[lc][k][0] - vals[lx][k][0]
            d2 = vals[lm][k][0] - vals[lc][k][0]
            R = d2 / d1 if d1 else float("nan")
            nmin = min(ntri[lv].get(k, 0) for lv in levels)
            weak = nmin < MIN_TRI or abs(d1) < 1e-5
            Rs.append(R); flags.append(weak)
            print(f"    eta [{k[0]:.1f},{k[1]:.1f})  d1 {d1:+.7f}  d2 {d2:+.7f}  "
                  f"R {R:+.4f}  n_min {nmin:4d}"
                  f"{'   ← too weak to read' if weak else ''}", flush=True)
            rows.append(dict(geom=geom, path=path, eta_lo=k[0], eta_hi=k[1],
                             d1=round(d1, 10), d2=round(d2, 10),
                             R=(None if R != R else round(R, 5)),
                             n_tri_min=nmin, readable=not weak,
                             R_first_order=round(r_first, 4),
                             R_falsify_ceiling=round(r_max, 4)))
        ax = axes[0 if geom == "m6wing" else 1]
        good = [(c, R) for c, R, w in zip(centres, Rs, flags) if not w]
        if good:
            ax.plot([g[0] for g in good], [g[1] for g in good], "-o", lw=1.8,
                    ms=6, label=path)
        bad = [(c, R) for c, R, w in zip(centres, Rs, flags) if w and R == R]
        if bad:
            ax.plot([b[0] for b in bad], [b[1] for b in bad], "x", ms=9, mew=2,
                    color="gray", label=f"{path} (too weak)")
        ax.axhline(r_first, color="g", ls="--", lw=1.2)
        ax.axhline(r_max, color="r", ls=":", lw=1.2)
        ax.set_title(f"{geom}   (— first order {r_first:.2f}, ⋯ ceiling "
                     f"{r_max:.2f})", fontsize=10.5)
        ax.set_xlabel(r"$\eta = z/b_{semi}$   (root → tip)")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    axes[0].set_ylabel("LE-band R per spanwise bin")
    fig.suptitle("LE-2: is the wing-alone LE deficit at the ROOT? "
                 "(M0.50, α 3.06; lower R = better convergence)", fontsize=11)
    fig.tight_layout()
    p = os.path.join(ART, "le2_spanwise.png")
    fig.savefig(p, dpi=130); plt.close(fig)
    print(f"\nwrote {p}")
    if rows:
        with open(CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        print(f"wrote {CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
