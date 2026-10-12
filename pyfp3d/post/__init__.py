"""Post-processing: output writers and visualization.

  - vtk_out.py: ParaView-compatible .vtu writer

Reference: design.md §7 (Architecture)
"""

from . import vtk_out

__all__ = ["vtk_out"]
