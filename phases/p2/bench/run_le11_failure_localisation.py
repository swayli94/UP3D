"""LE-11: stop hypothesising -- LOCALISE the failure and the entropy correction's effect.

Three hypotheses have now been refuted by my own measurements in a row: extra artificial
dissipation (LE-9: nu_max +55 % changed nothing and erased the shock), the entropy correction
weakening the shock (LE-9: M2max 6.426 ON against 5.607 OFF, i.e. STRONGER), and supersonic-region
connectivity (LE-10: largest-component fraction 0.92-0.95 across both groups, no separation).
Continuing to guess a fourth mechanism has falling expected value, so this round diagnoses instead
of testing a hypothesis -- it asks WHERE, which needs no prior guess about WHY.

Two questions, both answerable from the cached phi at zero solve cost:

  (a) WHERE does the stalled residual live? The three failing legs stall at |R| 7.3e-06 /
      7.0e-05 / 3.2e-06, so there is a nonzero residual to localise. The converging legs sit at
      ~1e-14 and act as the control (nothing to find). Localised against geometry: spanwise eta,
      chordwise x/c, distance to the tip, and whether the element is supersonic.
  (b) WHERE do the converging and failing SOLUTIONS differ? The entropy-ON and entropy-OFF fields
      differ by 7 % in cl -- mapping that difference spatially says what the entropy correction
      actually does, rather than inferring it from a scalar.

★ Caching gap recorded: LE-10 stored only phi, not gamma, so the full Newton residual (which
includes the Kutta rows and the FROZEN upwind selection) cannot be reconstructed. What is computed
here is the mass-conservation residual with a LIVE selection -- the right quantity for "where is
it big", but not bit-identical to the solver's own frozen residual. Next cache writes gamma too.

Outputs (TRACKED): bench/gate_results/le11_localisation.csv
                   bench/gate_results/capability/le11_residual_map.png
"""

import csv
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

from pyfp3d.kernels.upwind import UpwindOperator                    # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry          # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI, chord_at, x_le            # noqa: E402
from pyfp3d.physics.isentropic import (density_field,               # noqa: E402
                                       mach_number_squared)

OUT = os.path.join(_GATE)
ART = os.path.join(OUT, "capability")
CACHE = os.path.join(OUT, "le10_cache")
CSV = os.path.join(OUT, "le11_localisation.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_INF = 0.84
LEGS = ["ref_entropyON_c1.5", "ref_entropyOFF_c1.5", "entropyOFF_c2.0",
        "entropyOFF_c3.0", "taperON_entropyOFF"]
KEYS = ["leg", "converged", "res_l2", "res_max", "argmax_eta", "argmax_xc",
        "argmax_dist_tip", "argmax_supersonic", "top1pct_med_eta",
        "top1pct_med_xc", "top1pct_frac_tip", "top1pct_frac_supersonic", "note"]


def main():
    mc, wc = cut_wake(read_mesh(MP))
    B, V = precompute_element_geometry(mc.nodes, mc.elements)
    upw = UpwindOperator(mc.nodes, mc.elements, weighted=False)
    cen = mc.nodes[mc.elements].mean(axis=1)
    eta = cen[:, 2] / B_SEMI
    #: chordwise position needs the local chord; outside the wing span it is meaningless,
    #: so it is left as NaN there rather than silently wrapping (the fuselage is in this mesh)
    xc = np.full(len(cen), np.nan)
    inwing = (cen[:, 2] > 0.0) & (cen[:, 2] < B_SEMI)
    z = cen[inwing, 2]
    xc[inwing] = (cen[inwing, 0] - np.array([x_le(zi) for zi in z])) / \
        np.array([chord_at(zi) for zi in z])
    dist_tip = np.abs(cen[:, 2] - B_SEMI)

    rows, fields = [], {}
    for leg in LEGS:
        npz = os.path.join(CACHE, f"{leg}.npz")
        if not os.path.exists(npz):
            print(f"  {leg}: no cache"); continue
        d = np.load(npz)
        phi, conv = d["phi"], bool(d["conv"])
        g = np.einsum("eaj,ea->ej", B, phi[mc.elements])
        q2 = np.einsum("ej,ej->e", g, g)
        m2 = mach_number_squared(q2, M_INF)
        rho = density_field(q2, M_INF)
        #: LIVE selection (see the caching-gap note): the frozen one the solver used is
        #: not recoverable without gamma, so this is a localisation quantity, not the
        #: solver's own residual.
        #: ★ signature and constants READ, not recalled -- this round I wrote
        #: `density_ratio` (does not exist; it is density_field) and passed NINE
        #: positional args to rho_tilde, which takes
        #: (grad, q2, rho, m_inf, C, m_crit, gamma, rho_floor). Production constants come
        #: from tests/test_p8_newton.py: upwind_c 1.5, m_crit 0.95, rho_floor 0.05.
        rho_t = upw.rho_tilde(g, q2, rho, M_INF, 1.5, 0.95, 1.4, 0.05)
        #: nodal mass-conservation residual  R_i = sum_e rho_t V B_{e,i}.grad(phi)_e
        contrib = (rho_t * V)[:, None] * np.einsum("eaj,ej->ea", B, g)
        R = np.zeros(len(mc.nodes))
        np.add.at(R, mc.elements, contrib)
        free = np.setdiff1d(np.arange(len(mc.nodes)),
                            np.unique(mc.boundary_faces["farfield"]))
        Rf = np.abs(R[free])
        #: ★ NORMALISE by element volume. The first version localised the RAW nodal
        #: residual and every leg reported its argmax at the same far-field element
        #: (eta = 8.23, d_tip = 8.65) -- because the raw residual is volume-weighted and
        #: h_far = 200*h_wall, so it only said "the biggest raw residual is in the biggest
        #: element". Five legs giving the identical location was the tell. Dividing by V
        #: gives an error DENSITY, and the wing near-field mask keeps the far field from
        #: dominating by sheer element size.
        Re = np.abs(contrib).sum(axis=1) / np.maximum(V, 1e-30)
        near = (np.abs(cen[:, 2]) < 1.3 * B_SEMI) & (cen[:, 0] > -1.0) & \
               (cen[:, 0] < 3.0) & (np.abs(cen[:, 1]) < 1.0)
        Re = np.where(near, Re, 0.0)
        i = int(np.argmax(Re))
        cut = np.quantile(Re, 0.99)
        top = Re >= cut
        fields[leg] = (Re, conv)
        row = dict(leg=leg, converged=conv,
                   res_l2=float(np.linalg.norm(Rf)), res_max=float(Rf.max()),
                   argmax_eta=round(float(eta[i]), 4),
                   argmax_xc=(None if xc[i] != xc[i] else round(float(xc[i]), 4)),
                   argmax_dist_tip=round(float(dist_tip[i]), 5),
                   argmax_supersonic=bool(m2[i] > 1.0),
                   top1pct_med_eta=round(float(np.median(eta[top])), 4),
                   top1pct_med_xc=round(float(np.nanmedian(xc[top])), 4),
                   top1pct_frac_tip=round(float(np.mean(dist_tip[top] < 0.1 * B_SEMI)), 4),
                   top1pct_frac_supersonic=round(float(np.mean(m2[top] > 1.0)), 4),
                   note="")
        print(f"  {leg:22s} conv={conv!s:5s} |R|2={row['res_l2']:.3e} "
              f"max={row['res_max']:.3e}  argmax: eta={row['argmax_eta']} "
              f"x/c={row['argmax_xc']} d_tip={row['argmax_dist_tip']} "
              f"ss={row['argmax_supersonic']}", flush=True)
        print(f"  {'':22s} top-1%: med eta={row['top1pct_med_eta']} "
              f"med x/c={row['top1pct_med_xc']} "
              f"frac within 10% of tip={row['top1pct_frac_tip']} "
              f"frac supersonic={row['top1pct_frac_supersonic']}", flush=True)
        rows.append(row)

    # ---- (b) where do the ON and OFF solutions differ? --------------------
    a = os.path.join(CACHE, "ref_entropyON_c1.5.npz")
    b = os.path.join(CACHE, "ref_entropyOFF_c1.5.npz")
    if os.path.exists(a) and os.path.exists(b):
        pa, pb = np.load(a)["phi"], np.load(b)["phi"]
        ga = np.einsum("eaj,ea->ej", B, pa[mc.elements])
        gb = np.einsum("eaj,ea->ej", B, pb[mc.elements])
        dq = np.abs(np.einsum("ej,ej->e", ga, ga)
                    - np.einsum("ej,ej->e", gb, gb))
        j = int(np.argmax(dq))
        cut = np.quantile(dq, 0.99)
        top = dq >= cut
        print(f"\n  (b) ON vs OFF |dq2| max at eta={eta[j]:.4f} "
              f"x/c={xc[j]:.4f} d_tip={dist_tip[j]:.5f}")
        print(f"      top-1% of |dq2|: med eta={np.median(eta[top]):.4f} "
              f"med x/c={np.nanmedian(xc[top]):.4f} "
              f"frac within 10% of tip={np.mean(dist_tip[top] < 0.1*B_SEMI):.4f}")
        rows.append(dict(leg="ON_minus_OFF_dq2", converged=None,
                         argmax_eta=round(float(eta[j]), 4),
                         argmax_xc=(None if xc[j] != xc[j] else round(float(xc[j]), 4)),
                         argmax_dist_tip=round(float(dist_tip[j]), 5),
                         top1pct_med_eta=round(float(np.median(eta[top])), 4),
                         top1pct_med_xc=round(float(np.nanmedian(xc[top])), 4),
                         top1pct_frac_tip=round(float(np.mean(dist_tip[top] < 0.1*B_SEMI)), 4),
                         note="solution difference, not a residual"))

    # ---- figure: residual vs eta and vs x/c for each leg -----------------
    if fields:
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8))
        for leg, (Re, conv) in fields.items():
            cutq = np.quantile(Re, 0.99)
            m = Re >= cutq
            st = "-" if conv else "--"
            axes[0].hist(eta[m], bins=25, histtype="step", ls=st, lw=1.6,
                         label=f"{leg} {'(conv)' if conv else '(FAIL)'}")
            ok = m & ~np.isnan(xc)
            axes[1].hist(xc[ok], bins=25, range=(0, 1), histtype="step", ls=st,
                         lw=1.6, label=leg)
        axes[0].set_xlabel(r"$\eta = z/b_{semi}$"); axes[0].set_ylabel("count in top-1% |R|")
        axes[1].set_xlabel("x/c")
        for ax in axes:
            ax.grid(alpha=0.3); ax.legend(fontsize=7)
        fig.suptitle("LE-11: where the top-1% mass residual sits "
                     "(solid = converged, dashed = failed)", fontsize=11)
        fig.tight_layout()
        p = os.path.join(ART, "le11_residual_map.png")
        fig.savefig(p, dpi=130); plt.close(fig)
        print(f"\nwrote {p}")
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
