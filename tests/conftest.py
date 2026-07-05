"""
Pytest configuration and shared fixtures for pyFP3D test suite.

Gate-level tests should use these fixtures to ensure consistent reference data,
mesh sets, and artifact storage.
"""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def artifacts_dir():
    """Temporary directory for test artifacts (PNG, CSV, VTU files).
    
    In CI, these should be captured and stored for gate validation.
    """
    with tempfile.TemporaryDirectory(prefix="pyfp3d_artifacts_") as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def gate_artifacts_dir(artifacts_dir, request):
    """
    Create a subdirectory for the current test's gate artifacts.
    
    Naming convention: artifacts_dir / test_gate_id /
    """
    gate_id = request.node.name.replace("test_", "").upper()
    gate_dir = artifacts_dir / gate_id
    gate_dir.mkdir(parents=True, exist_ok=True)
    return gate_dir


@pytest.fixture
def reference_mesh_dir():
    """Return path to cases/reference_data (immutable reference meshes and data)."""
    repo_root = Path(__file__).parent.parent
    ref_dir = repo_root / "cases" / "reference_data"
    ref_dir.mkdir(parents=True, exist_ok=True)
    return ref_dir


@pytest.fixture
def mesh_dir():
    """Return path to cases/meshes (where coarse/medium/fine families are stored)."""
    repo_root = Path(__file__).parent.parent
    mesh_dir = repo_root / "cases" / "meshes"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    return mesh_dir


@pytest.fixture
def disable_numba_jit(monkeypatch):
    """
    Fixture to disable Numba JIT during tests for debugging.
    
    Usage:
        def test_something(disable_numba_jit):
            # Numba @njit will run in object mode, allowing prints/pdb
            ...
    """
    monkeypatch.setenv("PYFP3D_NOJIT", "1")
