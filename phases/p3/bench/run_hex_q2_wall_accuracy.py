"""Q2 — is the WALL SURFACE resolved better on a structured body-fitted grid?

The user's second concern, pre-registered in
phases/p3/docs/dev_phase_three/20260811-0300-hex-mesh-prereg.md §4 Q2. The cylinder is the first leg
because it has an EXACT potential solution and the unstructured family is already measured,
so this is a same-geometry A/B against external truth rather than against our own history.

Bands (fixed before running):
    W1  hex max |Cp error| <= 1/2 the unstructured value AT EQUAL WALL DOF   -> PASS
    W2  hex convergence order >= 1.5 (the unstructured family is ~1.0)       -> PASS
    W3  facet-normal error                                       RECORDED, no criterion

★ EQUAL WALL DOF is achievable exactly here, and that is a direct payoff of G0: `n_theta` is
a single-variable knob, so the hex ladder is built with n_theta = 126 / 252 / 504 to match the
unstructured wall counts node for node. Without that property the comparison would have had to
argue about "roughly comparable resolution", which is the kind of argument phase two kept
losing.

★★ SAME EXTRACTOR: both families go through `tests/mesh_utils.run_cylinder_case` -- the same
Laplace solve, the same quadratic wall-tangential gradient recovery, the same Cp definition and
the same exact reference. A wall-accuracy A/B computed by two different extractors would not be
a measurement of the mesh.

⚠ Registered risk 2 stands: an O-grid wall is still FLAT FACETS, so W3 is not expected to
improve. It is reported so it cannot later be counted as a gain.

Outputs (TRACKED): bench/gate_results/hex_q2_wall_accuracy.csv
"""

import csv
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "tests"))

from mesh_utils import run_cylinder_case                        # noqa: E402
from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.meshgen.extrude import extrude_single_layer         # noqa: E402
from pyfp3d.meshgen.structured import cylinder_o_grid_2d        # noqa: E402

CSV = os.path.join(HERE, "gate_results", "hex_q2_wall_accuracy.csv")
UNSTRUCTURED = os.path.join(REPO, "cases", "meshes", "cylinder_2.5d")
#: (level, unstructured wall-node count, hex h_wall_normal, dz).
#: ★ n_theta is set to HALF the target wall-node count, because `extrude_single_layer` puts the
#: 2-D wall nodes on BOTH planes: extruded wall nodes = 2 x n_theta. The first run asserted
#: 252 vs 126 and caught the factor -- which is what that assert exists for, and why the match
#: is asserted rather than assumed.
LADDER = (("coarse", 126, 0.04, 0.20),
          ("medium", 252, 0.02, 0.10),
          ("fine", 504, 0.01, 0.05))
R_FAR, GROWTH = 20.0, 1.15
W1_FACTOR = 0.5          # hex must be at most half the unstructured max error
W2_ORDER = 1.5           # hex order must reach this (unstructured measures ~1.0)


def facet_normal_error_deg(mesh):
    wall = np.asarray(mesh.boundary_faces["wall"], dtype=np.int64)
    p = mesh.nodes[wall]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    c = p.mean(axis=1); c[:, 2] = 0.0
    c /= np.linalg.norm(c, axis=1, keepdims=True)
    ang = np.degrees(np.arccos(np.clip(np.abs(np.einsum("ij,ij->i", n, c)), -1.0, 1.0)))
    return float(ang.max())


def order_of(errs, refine=2.0):
    """Observed order between successive levels of a factor-2 ladder."""
    return [float(np.log(errs[i] / errs[i + 1]) / np.log(refine))
            for i in range(len(errs) - 1)]


def main():
    rows, unst, hexs = [], [], []
    print("Q2: wall accuracy against the EXACT cylinder solution, same extractor both "
          "families\n")
    for level, n_wall, h_n, dz in LADDER:
        # ---- unstructured (the existing family) -----------------------------
        mu = read_mesh(os.path.join(UNSTRUCTURED, f"{level}.msh"))
        ru = run_cylinder_case(None, mesh=mu)
        eu_max, eu_mean = float(ru["error"].max()), float(ru["error"].mean())
        nu = len(np.unique(np.asarray(mu.boundary_faces["wall"])))

        # ---- structured, n_theta SET to match the wall DOF exactly ----------
        assert n_wall % 2 == 0, f"{level}: odd wall-node target {n_wall}"
        pts, tris, edges, info = cylinder_o_grid_2d(
            radius=1.0, r_far=R_FAR, n_theta=n_wall // 2, h_wall_normal=h_n, growth=GROWTH)
        mh = extrude_single_layer(pts, tris, edges, dz=dz, name=f"cyl_hex_{level}")
        rh = run_cylinder_case(None, mesh=mh)
        eh_max, eh_mean = float(rh["error"].max()), float(rh["error"].mean())
        nh = len(np.unique(np.asarray(mh.boundary_faces["wall"])))

        assert nh == nu, f"wall DOF not matched at {level}: {nh} vs {nu}"
        unst.append(eu_max); hexs.append(eh_max)
        for fam, e_max, e_mean, m, nw, res in (
                ("unstructured", eu_max, eu_mean, mu, nu, ru),
                ("structured_hex", eh_max, eh_mean, mh, nh, rh)):
            rows.append(dict(level=level, family=fam, wall_nodes=nw,
                             n_nodes=len(m.nodes), n_tets=len(m.elements),
                             cp_err_max=round(e_max, 8), cp_err_mean=round(e_mean, 8),
                             normal_err_max_deg=round(facet_normal_error_deg(m), 6),
                             residual_norm=float(res["residual_norm"]),
                             n_cg=int(res["n_cg_iterations"])))
        print(f"  {level:7} wall DOF {nu:>4} (matched)   "
              f"unstructured max |dCp| {eu_max:.5f}   hex {eh_max:.5f}   "
              f"ratio {eh_max / eu_max:.3f}")

    ou, oh = order_of(unst), order_of(hexs)
    print(f"\n  observed order (max |dCp|):  unstructured {['%.3f' % o for o in ou]}   "
          f"hex {['%.3f' % o for o in oh]}")
    print(f"  normal-error max (deg):      "
          f"unstructured {[r['normal_err_max_deg'] for r in rows if r['family']=='unstructured']}"
          f"   hex {[r['normal_err_max_deg'] for r in rows if r['family']=='structured_hex']}")

    ratios = [h / u for h, u in zip(hexs, unst)]
    w1 = all(r <= W1_FACTOR for r in ratios)
    w2 = min(oh) >= W2_ORDER
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== reading (bands fixed in the pre-registration §4 Q2) ===")
    print(f"  W1  hex/unstructured max-error ratios {['%.3f' % r for r in ratios]} "
          f"vs <= {W1_FACTOR}  -> {'PASS' if w1 else 'FAIL'}")
    print(f"  W2  hex order min {min(oh):.3f} vs >= {W2_ORDER}  -> "
          f"{'PASS' if w2 else 'FAIL'}   (unstructured min {min(ou):.3f})")
    print(f"  W3  RECORDED, no criterion: the facet-normal error is O(h) on BOTH families -- "
          f"an O-grid wall is still flat facets (registered risk 2).")
    if w1 or w2:
        print("\n  ⇒ Q2 has at least one band PASS: route (A) buys wall accuracy, so per the")
        print("     kill criterion it stays alive regardless of Q1.")
    else:
        print("\n  ⇒ Q2 FAILS BOTH bands. Per the kill criterion this alone does NOT kill (A) --")
        print("     Q1 (donor determinacy) is the other half, and (A) dies only if BOTH fail.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
