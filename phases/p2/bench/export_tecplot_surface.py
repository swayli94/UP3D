"""Export a mesh's boundary surfaces as Tecplot ASCII (FEPOINT / ET=TRIANGLE), one zone each.

Requested so the round-tip ONERA M6 wall can be opened and inspected directly. FEPOINT with
ET=TRIANGLE is Tecplot's standard finite-element surface format, so an unstructured triangulated
surface loads without any structured-grid assumption.

Each boundary group becomes its own ZONE, and each zone is re-indexed to only the nodes it uses,
so the wall zone does not drag the whole 80k-node volume mesh along -- Tecplot then shows the wing
surface and the wake sheet as separate, independently toggleable objects, which is the point when
what you want to look at is how the sheet terminates at the tip.

`farfield` is skipped by default: it is a coarse outer sphere that would swamp the view scale.

Usage:
    python bench/export_tecplot_surface.py <mesh.msh> [out.dat] [--groups wall,wake]

Outputs (TRACKED): whatever path is given, default alongside the mesh as <stem>_surface.dat
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)

from pyfp3d.mesh.reader import read_mesh                            # noqa: E402

SKIP_DEFAULT = ("farfield",)


def write_tecplot(mesh, out_path, groups):
    with open(out_path, "w") as fh:
        fh.write('TITLE = "%s boundary surfaces"\n'
                 % os.path.basename(out_path).replace("_surface.dat", ""))
        fh.write('VARIABLES = "X", "Y", "Z"\n')
        for name in groups:
            tris = np.asarray(mesh.boundary_faces[name], dtype=np.int64)
            if tris.size == 0:
                continue
            #: re-index to this zone's own nodes -- otherwise every zone would carry all
            #: 80k volume nodes and Tecplot would show a cloud of unused points
            used = np.unique(tris.reshape(-1))
            remap = -np.ones(len(mesh.nodes), dtype=np.int64)
            remap[used] = np.arange(len(used))
            local = remap[tris]
            assert local.min() >= 0, "re-index produced an unused reference"
            pts = mesh.nodes[used]
            fh.write(f'ZONE T="{name}", N={len(pts)}, E={len(local)}, '
                     f'F=FEPOINT, ET=TRIANGLE\n')
            for x, y, z in pts:
                fh.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
            for a, b, c in local + 1:            # Tecplot is 1-based
                fh.write(f"{a} {b} {c}\n")
            print(f"  zone {name:10s} {len(pts):7d} nodes  {len(local):7d} tris")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__); return 2
    mesh_path = args[0]
    out = args[1] if len(args) > 1 else \
        os.path.splitext(mesh_path)[0] + "_surface.dat"
    sel = None
    for a in sys.argv[1:]:
        if a.startswith("--groups"):
            sel = a.split("=", 1)[1].split(",") if "=" in a else None
    mesh = read_mesh(mesh_path)
    groups = sel or [g for g in mesh.boundary_faces if g not in SKIP_DEFAULT]
    print(f"{mesh_path}: {len(mesh.nodes)} nodes / {len(mesh.elements)} tets")
    write_tecplot(mesh, out, groups)
    #: validate what was written rather than trusting the writer
    n_zone = sum(1 for line in open(out) if line.startswith("ZONE"))
    size = os.path.getsize(out) / 1e6
    print(f"wrote {out}  ({n_zone} zones, {size:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
