"""W2-conf -- the conforming oracle's fuselage-lift decomposition, measured
with the SAME extractor as the F1 legs (fuselage_cl_parts, bw0.06).

F1 showed the fragment clip RAISES cl_fus_out (medium a=3.06: 0.0214 ->
0.0504), overshooting the pre-registered |d| <= 20% carryover guard. But
the A-side total (0.0357) sits exactly ON the B9 conforming oracle total
(0.0356) -- so the question "is C's out-band increase physical (the root
vortex now convects along the fuselage side, as the conforming sheet has
always done) or an artifact redistribution?" needs the oracle's band/out
SPLIT, not its total. No committed artifact carries it, so measure it here:
B9's committed conf_medium.npz (pressure leg, read-only) on the conforming
medium mesh through cut_wake, with a single-valued side_potentials shim --
the same D11/Cp/pressure-force core as wall_forces/fuselage_cl_parts.

Artifacts: results/w2_conf.csv (appended to f1 evidence)
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B23 = REPO_ROOT / "cases/analysis/b23_junction_discriminator"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(B23))

from run_w2_decomp import fuselage_cl_parts  # noqa: E402
from pyfp3d.mesh.reader import read_mesh  # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake  # noqa: E402
from pyfp3d.post.surface import planform_area  # noqa: E402

CONF_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody_conforming"
B9_NPZ = REPO_ROOT / "cases/demo/b9_wingbody/results/conf_medium.npz"
ALPHA = 3.06
OUT = HERE / "results"


class _SingleValued:
    """side_potentials shim: the conforming field is single-valued, so the
    D11 upper/lower selection is a no-op."""

    def side_potentials(self, phi):
        return phi, phi


def main():
    d = np.load(B9_NPZ, allow_pickle=True)
    phi = np.asarray(d["pressure_phi"], dtype=np.float64)
    mesh, wc = cut_wake(read_mesh(str(CONF_DIR / "medium.msh")))
    assert len(phi) == len(mesh.nodes), (len(phi), len(mesh.nodes))
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    rows = fuselage_cl_parts(mesh, _SingleValued(), phi, ALPHA, s_ref)
    out = OUT / "w2_conf.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "level", "alpha", "all", "poles",
                    "band_bw0.06", "out_bw0.06"])
        w.writerow(["conforming_pressure", "medium", ALPHA, rows["all"],
                    rows["poles"], rows["band_bw0.06"], rows["out_bw0.06"]])
    print(f"conforming medium a={ALPHA}: all={rows['all']:+.5f} "
          f"poles={rows['poles']:+.5f} band={rows['band_bw0.06']:+.5f} "
          f"out={rows['out_bw0.06']:+.5f}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
