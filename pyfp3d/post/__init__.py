"""Post-processing: output writers and visualization.

  - vtk_out.py: ParaView-compatible .vtu writer
  - artifacts.py: Headless PNG/CSV artifact generation for gate checks

Reference: design.md §7 (Architecture)
"""

from . import vtk_out

__all__ = ["vtk_out"]
