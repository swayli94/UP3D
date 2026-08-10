"""LE-12: do the taper and the entropy correction act on the SAME degree of freedom?

LE-11 localised the failure: the residual-density peak sits at the tip (argmax at eta ~ 1.001,
d_tip 0.001, supersonic; 83 % of the top-1 % within 10 % of the tip) in ALL five legs, and what
separates converging from failing is not the LOCATION but the MAGNITUDE there -- 7.5e-04/8.2e-04
converging against 3.9e-03/3.9e-03/1.1e-02 failing, a factor 5-15 at the same place. So the taper
and the entropy correction both suppress one peak rather than moving it.

The ambiguity LE-11 could not resolve: "the effect is largest at the tip" is equally consistent
with acting THERE and with the tip amplifying ANY perturbation. This round separates them, at zero
solve cost, using a paired construction the cache happens to allow -- both differences share the
SAME baseline state:

    D_entropy = phi(entropy ON,  no taper) - phi(entropy OFF, no taper)
    D_taper   = phi(entropy OFF, taper)    - phi(entropy OFF, no taper)

Same baseline, one knob each. Then:

  the two difference FIELDS are spatially aligned  => they act on the same degree of freedom,
      and the shared tip localisation is mechanism, not amplification
  aligned only near the tip, uncorrelated elsewhere => the tip is amplifying two unrelated
      perturbations, and the "same mechanism" reading is wrong

Measured: Pearson correlation of the signed dq2 fields (globally, in the tip region, and outside
it), cosine similarity, and the best-fit proportionality D_entropy ~ k * D_taper with its residual
-- if one is a scalar multiple of the other they are the same direction in state space.

Outputs (TRACKED): bench/gate_results/le12_same_dof.csv
                   bench/gate_results/capability/le12_same_dof.png
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

from pyfp3d.mesh.metrics import precompute_element_geometry          # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402

OUT = os.path.join(_GATE)
ART = os.path.join(OUT, "capability")
CACHE = os.path.join(OUT, "le10_cache")
CSV = os.path.join(OUT, "le12_same_dof.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
BASE = "ref_entropyOFF_c1.5"
ENT = "ref_entropyON_c1.5"
TAP = "taperON_entropyOFF"


def main():
    for k in (BASE, ENT, TAP):
        if not os.path.exists(os.path.join(CACHE, f"{k}.npz")):
            raise SystemExit(f"missing cache: {k}")
    mc, _wc = cut_wake(read_mesh(MP))
    B, V = precompute_element_geometry(mc.nodes, mc.elements)
    cen = mc.nodes[mc.elements].mean(axis=1)
    eta = cen[:, 2] / B_SEMI

    def q2_of(name):
        phi = np.load(os.path.join(CACHE, f"{name}.npz"))["phi"]
        g = np.einsum("eaj,ea->ej", B, phi[mc.elements])
        return np.einsum("ej,ej->e", g, g)

    q_base, q_ent, q_tap = q2_of(BASE), q2_of(ENT), q2_of(TAP)
    d_ent = q_ent - q_base
    d_tap = q_tap - q_base

    #: near-field only, the same mask LE-11 needed -- far-field elements are enormous and
    #: carry no information about the tip mechanism
    near = ((np.abs(cen[:, 2]) < 1.3 * B_SEMI) & (cen[:, 0] > -1.0)
            & (cen[:, 0] < 3.0) & (np.abs(cen[:, 1]) < 1.0))
    tip = near & (np.abs(cen[:, 2] - B_SEMI) < 0.10 * B_SEMI)
    out = near & ~tip

    rows = []
    for label, m in (("near_field", near), ("tip_region", tip),
                     ("outside_tip", out)):
        a, b = d_ent[m], d_tap[m]
        if len(a) < 10 or a.std() == 0 or b.std() == 0:
            continue
        r = float(np.corrcoef(a, b)[0, 1])
        cos = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
        #: best-fit scalar multiple: if one difference is k times the other, they are the
        #: same direction in state space and the knobs move the same degree of freedom
        k = float(a @ b / (b @ b))
        resid = float(np.linalg.norm(a - k * b) / np.linalg.norm(a))
        print(f"  {label:13s} n={m.sum():7d}  corr={r:+.4f}  cos={cos:+.4f}  "
              f"k={k:+.4f}  ||D_ent - k D_tap||/||D_ent|| = {resid:.4f}",
              flush=True)
        rows.append(dict(region=label, n=int(m.sum()), corr=round(r, 6),
                         cosine=round(cos, 6), k_fit=round(k, 6),
                         rel_residual=round(resid, 6)))

    #: magnitude profile along the span, both differences, to see them by eye
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.6))
    bins = np.linspace(0.0, 1.15, 40)
    ctr = 0.5 * (bins[1:] + bins[:-1])
    for lab, d, st in (("entropy ON − OFF", d_ent, "-o"),
                       ("taper − none", d_tap, "--s")):
        prof = [np.sqrt(np.mean(d[near & (eta >= lo) & (eta < hi)] ** 2))
                if np.any(near & (eta >= lo) & (eta < hi)) else np.nan
                for lo, hi in zip(bins[:-1], bins[1:])]
        axes[0].semilogy(ctr, prof, st, ms=4, lw=1.6, label=lab)
    axes[0].axvline(1.0, color="k", ls=":", lw=1.0)
    axes[0].set_xlabel(r"$\eta = z/b_{semi}$")
    axes[0].set_ylabel(r"RMS $\Delta q^2$ in the strip")
    axes[0].grid(alpha=0.3); axes[0].legend(fontsize=9)
    axes[0].set_title("spanwise profile of each knob's effect", fontsize=10.5)
    sub = np.flatnonzero(near)
    if len(sub) > 40000:
        sub = sub[:: len(sub) // 40000]
    axes[1].scatter(d_tap[sub], d_ent[sub], s=2, alpha=0.25,
                    c=np.where(tip[sub], "tab:red", "tab:blue"))
    axes[1].set_xlabel(r"$\Delta q^2$  (taper − none)")
    axes[1].set_ylabel(r"$\Delta q^2$  (entropy ON − OFF)")
    axes[1].grid(alpha=0.3)
    axes[1].set_title("red = tip region, blue = rest of the near field",
                      fontsize=10.5)
    fig.suptitle("LE-12: do the two knobs move the same degree of freedom?",
                 fontsize=11)
    fig.tight_layout()
    p = os.path.join(ART, "le12_same_dof.png")
    fig.savefig(p, dpi=130); plt.close(fig)
    print(f"\nwrote {p}")
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"wrote {CSV}")

    print("\n=== reading ===")
    d = {r["region"]: r for r in rows}
    t, o = d.get("tip_region"), d.get("outside_tip")
    if t and o:
        print(f"  tip corr {t['corr']:+.3f}   outside-tip corr {o['corr']:+.3f}")
        if t["corr"] > 0.7 and o["corr"] > 0.7:
            print("  => SAME DEGREE OF FREEDOM everywhere: the two knobs move the")
            print("  field in the same direction, so the shared tip localisation is")
            print("  mechanism, not amplification.")
        elif t["corr"] > 0.7 >= o["corr"]:
            print("  => aligned ONLY at the tip: the tip amplifies two otherwise")
            print("  unrelated perturbations, and 'same mechanism' is the wrong")
            print("  reading -- what they share is the site, not the direction.")
        else:
            print("  => NOT aligned even at the tip: the two knobs act differently")
            print("  and their substitutability needs another explanation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
