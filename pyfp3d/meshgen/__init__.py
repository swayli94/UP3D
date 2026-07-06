"""
meshgen: scripted quasi-2D mesh generation (docs/roadmap.md Track M).

M0 deliverable: single-layer extruded quasi-2D tet meshes built from a
vanilla-Gmsh 2D triangulation (planar.py) and a globally consistent
prism -> 3-tet split (extrude.py). Node duplication for the wake stays
solver-side (mesh/wake_cut.py, P2) per agent-rules hard rule 8.
"""

from pyfp3d.meshgen.extrude import extrude_single_layer, assert_quad_split_consistency

__all__ = ["extrude_single_layer", "assert_quad_split_consistency"]
