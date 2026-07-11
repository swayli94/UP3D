"""
Track B / B4.5 gate: far-field A/B -- Dirichlet+vortex (option a) vs the
Lopez-style Neumann outlet (option b). docs/roadmap.md Track B B4.5;
docs/design_track_b.md section 5.4.

The heavy Lopez-style domain-size re-calibration (Gamma vs far-field radius R,
both mesh families) lives in the demo cases/demo/b4p5_farfield/. These are the
cheap regression locks on the committed 15-chord coarse meshes:

  - option a (vortex) is accurate at 15c: Gamma within 2% of the conforming
    solver (the vortex correction is why pyFP3D's compact domain works);
  - option b (neumann) truncates the O(Gamma/r) far-field tail, so at 15c its
    Gamma sits a few percent BELOW option a -- the signature that it needs a
    larger domain (the demo shows it converging up to a as R grows);
  - the far-field split into inflow (Dirichlet) / outflow (Neumann) is sane and
    the freestream net flux through the closed boundary is ~0.

VERDICT: option a stays the default (solve_multivalued_lifting(farfield=
"vortex")); option b is validated but domain-hungry.
"""

import csv
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.picard import solve_subsonic_lifting
from pyfp3d.solve.picard_ls import (
    _farfield_split, _neumann_outlet_rhs, solve_multivalued_lifting,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO = Path(__file__).parent.parent
M0_DIR = REPO / "cases" / "meshes" / "naca0012_2.5d"
M3_DIR = REPO / "cases" / "meshes" / "naca0012_wakefree_2.5d"
ALPHA, M_INF = 2.0, 0.5

_CACHE = {}


def _mesh(directory):
    p = directory / "coarse.msh"
    if not p.exists():
        pytest.skip(f"{p} not generated (gitignored)")
    return read_mesh(p)


def _solve_b(mesh, mode):
    key = (id(mesh), mode)
    if key not in _CACHE:
        z = mesh.nodes[:, 2]
        wls = WakeLevelSet(
            np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
            direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
        _CACHE[key] = solve_multivalued_lifting(
            mvop, mesh, M_INF, alpha_deg=ALPHA, farfield=mode)
    return _CACHE[key]


def _conforming_gamma(mesh):
    key = (id(mesh), "conf")
    if key not in _CACHE:
        mc, wc = cut_wake(mesh)
        ref = solve_subsonic_lifting(mc, wc, M_INF, alpha_deg=ALPHA)
        _CACHE[key] = float(np.mean(ref["gamma"]))
    return _CACHE[key]


class TestFarfieldSplit:
    """Unit checks on the inflow/outflow split and the Neumann RHS -- pure
    geometry/assembly, no solve."""

    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    def test_split_is_sane(self, directory):
        mesh = _mesh(directory)
        faces, area, u_dot_n, inflow = _farfield_split(mesh, ALPHA, 1.0)
        assert np.all(area > 0.0), "degenerate far-field face area"
        assert (u_dot_n < 0).any() and (u_dot_n > 0).any(), \
            "far field must have both inflow and outflow faces"
        assert len(inflow) > 0

    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    def test_freestream_net_flux_through_closed_boundary_is_zero(self, directory):
        # oint u.n dS = 0 over a closed surface for the divergence-free
        # freestream: outward normals + areas are consistent iff this holds.
        mesh = _mesh(directory)
        faces, area, u_dot_n, _ = _farfield_split(mesh, ALPHA, 1.0)
        net = float(np.sum(u_dot_n * area))
        scale = float(np.sum(np.abs(u_dot_n) * area))
        assert abs(net) / scale < 1e-3, f"net freestream flux {net:.3e} (rel {abs(net)/scale:.1e})"

    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    def test_neumann_rhs_only_on_outflow(self, directory):
        mesh = _mesh(directory)
        faces, area, u_dot_n, _ = _farfield_split(mesh, ALPHA, 1.0)
        n_total = len(mesh.nodes)   # >= main dofs; extra aux beyond it here
        b = _neumann_outlet_rhs(mesh, ALPHA, 1.0, n_total)
        assert np.sum(b) > 0.0, "outflow flux RHS should be net positive"
        # the RHS support is a subset of outflow-face nodes
        support = np.flatnonzero(b != 0.0)
        outflow_nodes = np.unique(faces[u_dot_n > 0].ravel())
        assert set(support.tolist()).issubset(set(outflow_nodes.tolist()))


class TestOptionAIsDefaultAndAccurate:
    def test_vortex_matches_conforming_M0(self):
        mesh = _mesh(M0_DIR)
        ga = _solve_b(mesh, "vortex")
        assert ga["converged"]
        gc = _conforming_gamma(mesh)
        e = abs(ga["gamma"] - gc) / gc
        assert e < 0.02, f"option a Gamma {ga['gamma']:.4f} vs conf {gc:.4f} ({e*100:.2f}%)"

    def test_vortex_is_the_default(self):
        # regression pin: the shipped default is option a
        import inspect
        sig = inspect.signature(solve_multivalued_lifting)
        assert sig.parameters["farfield"].default == "vortex"


class TestOptionBTruncatesAtCompactDomain:
    """At 15c option b (Neumann, no vortex) truncates O(Gamma/r): its Gamma
    sits a few percent BELOW option a. This is the signature that it needs
    the domain grown (the demo shows the convergence up to a as R grows)."""

    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    def test_neumann_below_vortex_by_a_few_percent(self, directory):
        mesh = _mesh(directory)
        ga = _solve_b(mesh, "vortex")
        gb = _solve_b(mesh, "neumann")
        assert gb["converged"]
        gap = (ga["gamma"] - gb["gamma"]) / ga["gamma"]
        assert 0.02 < gap < 0.08, (
            f"{directory.name}: option b Gamma {gb['gamma']:.4f} vs a "
            f"{ga['gamma']:.4f}, truncation gap {gap*100:.2f}% "
            f"(expected ~4% at 15c)")


if __name__ == "__main__":
    pytest.main([__file__, "-x", "-q"])
