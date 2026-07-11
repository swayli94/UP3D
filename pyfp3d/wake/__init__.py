"""
Track B level-set embedded wake (docs/design_track_b.md; roadmap Track B).

Parallel to the conforming path (mesh/wake_cut.py + constraints/wake.py) --
nothing here is imported by the shipped solver paths.
"""

from pyfp3d.wake.levelset import WakeLevelSet
from pyfp3d.wake.cut_elements import CutElementMap
from pyfp3d.wake.multivalued import MultivaluedOperator

__all__ = ["WakeLevelSet", "CutElementMap", "MultivaluedOperator"]
