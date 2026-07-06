"""
Gmsh generator for the spherical-shell validation case (gate G1.2: an
incompressible sphere, Cp vs 1 - 9/4 sin^2(theta)).

Two concentric spheres (r=1 "wall", r=r_out "farfield") are cut with the
OpenCASCADE kernel to leave the shell as the fluid volume; mesh size is
graded fine near the wall and coarse near the far field via a Distance +
Threshold field, per docs/roadmap.md Sec 0.1 tooling policy (gmsh CLI/API,
-format msh4).

Usage:
    python generate_sphere_shell.py --level medium
    python generate_sphere_shell.py --level coarse --level medium --level fine
"""

import argparse
from pathlib import Path

import gmsh

LEVELS = {
    # h_min, h_max, r_out, dist_min, dist_max
    "coarse": dict(h_min=0.20, h_max=3.0, r_out=20.0, dist_min=0.3, dist_max=8.0),
    "medium": dict(h_min=0.08, h_max=3.0, r_out=20.0, dist_min=0.3, dist_max=10.0),
    "fine": dict(h_min=0.03, h_max=2.0, r_out=15.0, dist_min=0.3, dist_max=8.0),
}


def generate_sphere_shell(
    out_path: Path,
    r_in: float = 1.0,
    r_out: float = 20.0,
    h_min: float = 0.08,
    h_max: float = 3.0,
    dist_min: float = 0.3,
    dist_max: float = 10.0,
) -> Path:
    """Build and mesh the shell between two concentric spheres; write a msh4 file."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("sphere_shell")

        inner = gmsh.model.occ.addSphere(0, 0, 0, r_in)
        outer = gmsh.model.occ.addSphere(0, 0, 0, r_out)
        gmsh.model.occ.cut([(3, outer)], [(3, inner)])
        gmsh.model.occ.synchronize()

        wall_tag = farfield_tag = None
        for dim, tag in gmsh.model.getEntities(2):
            bbox = gmsh.model.getBoundingBox(dim, tag)
            size = max(bbox[3] - bbox[0], bbox[4] - bbox[1], bbox[5] - bbox[2])
            if abs(size / 2 - r_in) < 0.1:
                wall_tag = tag
            elif abs(size / 2 - r_out) < 0.5:
                farfield_tag = tag
        if wall_tag is None or farfield_tag is None:
            raise RuntimeError("Could not identify wall/farfield surfaces after boolean cut")

        volumes = gmsh.model.getEntities(3)
        gmsh.model.addPhysicalGroup(2, [wall_tag], name="wall")
        gmsh.model.addPhysicalGroup(2, [farfield_tag], name="farfield")
        gmsh.model.addPhysicalGroup(3, [volumes[0][1]], name="fluid")

        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [wall_tag])
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", h_min)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", h_max)
        gmsh.model.mesh.field.setNumber(2, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(2, "DistMax", dist_max)
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(3)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(out_path))
    finally:
        gmsh.finalize()

    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--level", action="append", choices=sorted(LEVELS), default=None,
        help="Resolution level(s) to generate (default: medium)",
    )
    args = parser.parse_args()
    levels = args.level or ["medium"]

    out_dir = Path(__file__).parent
    for level in levels:
        out_path = out_dir / f"{level}.msh"
        generate_sphere_shell(out_path, **LEVELS[level])
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
