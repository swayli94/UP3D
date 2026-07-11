"""
Track B / B3 gate: lifting solve with IMPLICIT Kutta on the level-set path
(docs/roadmap.md Track B B3; docs/design_track_b.md §9).

No Γ secant and no master–slave Γ constraint: the TE jump is carried by the
multivalued aux DOFs, the g₁+g₂ wake LS convects it downstream, and its VALUE
is set by the B4 nonlinear TE pressure-equality (Bernoulli) Kutta. Γ is a
result, not an unknown of an outer loop.

The reference bracket is the project's committed corrected-panel data
(`cases/reference_data/naca0012_m05/cl_reference.csv`, the same file the
conforming G3.2 gate reads): the exact subcritical full-potential cl is
bracketed by the Prandtl–Glauert and Kármán–Tsien corrections. Gate values are
read from the CSV, never embedded here (roadmap reference-data discipline).

Dual-mesh rule (roadmap Track B working rules): every gate runs on BOTH the
wake-EMBEDDED M0 family (which also permits a strict same-mesh A/B against the
conforming solver) and the wake-FREE M3 family (generic cuts, no `wake` tag at
all — the actual Track B workflow form, where no conforming counterpart exists).
"""

import csv
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.picard import solve_laplace_lifting, solve_subsonic_lifting
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
MESHES = REPO_ROOT / "cases" / "meshes"
M0_DIR = MESHES / "naca0012_2.5d"                 # wake-embedded
M3_DIR = MESHES / "naca0012_wakefree_2.5d"        # wake-free (workflow form)
ALPHA = 2.0
M_INF = 0.5


def _band():
    path = (REPO_ROOT / "cases" / "reference_data" / "naca0012_m05"
            / "cl_reference.csv")
    with open(path) as f:
        for row in csv.DictReader(f):
            if abs(float(row["alpha_deg"]) - ALPHA) < 1e-9:
                return float(row["cl_pg"]), float(row["cl_kt"])
    raise LookupError(f"alpha={ALPHA} not in {path}")


def _load(directory, level):
    path = directory / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _solve_b(mesh, m_inf):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),      # chord plane, design.md §4
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return solve_multivalued_lifting(mvop, mesh, m_inf, alpha_deg=ALPHA)


class TestV3LiftBracket:
    """V3: M0.5 α=2° cl inside the [PG, KT] bracket — on BOTH mesh families."""

    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    def test_cl_inside_pg_kt(self, directory):
        cl_pg, cl_kt = _band()
        mesh = _load(directory, "medium")
        r = _solve_b(mesh, M_INF)
        assert r["converged"], f"{directory.name}: not converged"
        cl = r["cl_kj"]
        assert cl_pg <= cl <= cl_kt, (
            f"{directory.name}: cl_KJ {cl:.4f} outside [PG {cl_pg:.4f}, "
            f"KT {cl_kt:.4f}]"
        )


class TestSameMeshAB:
    """Strict same-mesh A/B against the conforming solver (only possible on the
    wake-EMBEDDED family). B4's TE Kutta brings this inside 1%; the B3 TE row
    (lower-side mass conservation) was 42% out.

    Run on COARSE only: it pins the A/B cheaply, while the medium gate is
    carried by the EXTERNAL [PG, KT] bracket above (which needs no conforming
    solve) and by B4's medium Γ-vs-conforming gate. Keeps the suite inside the
    G8.3 CI budget."""

    @pytest.mark.parametrize("level", ["coarse"])
    @pytest.mark.parametrize("m_inf", [0.0, M_INF])
    def test_gamma_matches_conforming(self, level, m_inf):
        mesh = _load(M0_DIR, level)
        mesh_cut, wc = cut_wake(mesh)
        ref = (solve_laplace_lifting(mesh_cut, wc, alpha_deg=ALPHA) if m_inf == 0
               else solve_subsonic_lifting(mesh_cut, wc, m_inf, alpha_deg=ALPHA))
        g_ref = float(np.mean(ref["gamma"]))
        g_b = _solve_b(mesh, m_inf)["gamma"]
        err = abs(g_b - g_ref) / abs(g_ref)
        assert err < 0.02, (
            f"{level} M{m_inf}: Γ_B {g_b:.4f} vs conforming {g_ref:.4f} "
            f"({err*100:.1f}%)"
        )


class TestWakeFreeAgreesWithEmbedded:
    """The workflow payoff: a mesh whose topology knows NOTHING about the wake
    (no `wake` tag, generic cuts) reproduces the embedded-mesh circulation."""

    def test_wakefree_matches_embedded(self):
        m0 = _load(M0_DIR, "medium")
        m3 = _load(M3_DIR, "medium")
        assert "wake" in m0.boundary_faces
        assert "wake" not in m3.boundary_faces
        g0 = _solve_b(m0, M_INF)["gamma"]
        g3 = _solve_b(m3, M_INF)["gamma"]
        assert abs(g3 - g0) / g0 < 0.02, f"embedded {g0:.4f} vs free {g3:.4f}"


class TestGammaEmerges:
    """Γ is a RESULT, not an outer-loop unknown: no secant, no master–slave."""

    def test_no_gamma_constraint_and_jump_is_convected(self):
        mesh = _load(M0_DIR, "coarse")
        z = mesh.nodes[:, 2]
        wls = WakeLevelSet(
            np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
            direction=(1.0, 0.0, 0.0),
        )
        cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                           wall_nodes=np.unique(mesh.boundary_faces["wall"]))
        mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
        r = solve_multivalued_lifting(mvop, mesh, 0.0, alpha_deg=ALPHA)

        # the g2 condition convects the jump: it stays ~constant along the wake
        cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
        jump = mvop.node_jump(r["phi_ext"], cut_nodes)
        d = cm.d[cut_nodes]
        near = jump[(d > 0.05) & (d < 1.0)]
        far = jump[d > 3.0]
        assert len(near) and len(far)
        assert abs(far.mean() - near.mean()) / abs(near.mean()) < 0.10, (
            "the wake jump is not being convected (it decays downstream) — "
            "check the far-field aux DOFs are FREE, not pinned to the vortex"
        )
        assert abs(near.mean() - r["gamma"]) / abs(r["gamma"]) < 0.10


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
