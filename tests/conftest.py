"""
Pytest configuration and shared fixtures for pyFP3D test suite.

Gate-level tests should use these fixtures to ensure consistent reference data,
mesh sets, and artifact storage.
"""

import pytest
import os
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent


@pytest.fixture
def artifacts_dir():
    """Persistent root for gate artifacts (PNG, CSV, VTU files).

    Defaults to <repo>/artifacts (gitignored) so the headless evidence each
    visual gate produces survives the test run and can be inspected -- the
    CLAUDE.md/roadmap workflow requires every visual gate to leave PNG+CSV
    artifacts behind. (This used to be a tempfile.TemporaryDirectory, which
    deleted every artifact at teardown and left artifacts/ permanently
    empty.) Set PYFP3D_ARTIFACTS_DIR to redirect, e.g. to a CI upload dir.
    """
    env = os.environ.get("PYFP3D_ARTIFACTS_DIR")
    base = Path(env) if env else REPO_ROOT / "artifacts"
    base.mkdir(parents=True, exist_ok=True)
    return base


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
