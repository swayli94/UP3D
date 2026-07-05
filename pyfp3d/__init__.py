"""
pyFP3D: 3D Unstructured-mesh Full-Potential Transonic Flow Solver

A Python + Numba implementation of the full-potential (FP) method for steady 
external flows over wings at transonic Mach numbers (0.3–0.87).

Core modules:
  - mesh: I/O, topology, metrics, coloring
  - physics: Isentropic constitutive relations
  - kernels: Residual assembly, artificial compressibility
  - solve: Linear and nonlinear solvers
  - post: VTK output, visualization

Reference: docs/design.md
Roadmap: docs/roadmap.md
"""

__version__ = "0.0.1"
__author__ = "swayli94"

import os

# Development mode flag: PYFP3D_NOJIT=1 disables Numba JIT for debugging
NOJIT = os.environ.get("PYFP3D_NOJIT", "0") == "1"
