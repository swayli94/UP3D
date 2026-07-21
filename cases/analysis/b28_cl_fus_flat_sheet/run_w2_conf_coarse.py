"""B28 coarse oracle: the conforming COARSE fuselage-lift decomposition --
the level-endpoint anchor the b25 w2_conf.csv (medium-only) does not carry.

Same extractor as run_w2_conf.py (fuselage_cl_parts, bw0.06), read-only on
B9's committed conf_coarse.npz (pressure leg). No solve.

Artifacts: results/w2_conf_coarse.csv
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
B9_NPZ = REPO_ROOT / "cases/demo/b9_wingbody/results/conf_coarse.npz"
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
    mesh, wc = cut_wake(read_mesh(str(CONF_DIR / "coarse.msh")))
    assert len(phi) == len(mesh.nodes), (len(phi), len(mesh.nodes))
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    rows = fuselage_cl_parts(mesh, _SingleValued(), phi, ALPHA, s_ref)
    out = OUT / "w2_conf_coarse.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "level", "alpha", "all", "poles",
                    "band_bw0.06", "out_bw0.06"])
        w.writerow(["conforming_pressure", "coarse", ALPHA, rows["all"],
                    rows["poles"], rows["band_bw0.06"], rows["out_bw0.06"]])
    print(f"conforming coarse a={ALPHA}: all={rows['all']:+.5f} "
          f"poles={rows['poles']:+.5f} band={rows['band_bw0.06']:+.5f} "
          f"out={rows['out_bw0.06']:+.5f}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
