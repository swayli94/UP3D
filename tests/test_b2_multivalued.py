"""
Track B / B2 gate: multivalued (CutFEM-style) FE assembly on the level-set
cut mesh (docs/roadmap.md Track B B2; docs/design_track_b.md sections
2.1/2.5, D6).

B2 delivers the extended-DOF assembly (kernels/cut_assembly.py,
wake/multivalued.py) and proves it is CONSISTENT: a cut element is the same
P1 element matrix assembled twice with two DOF numberings, and with the B2
continuity ("weld") closure the extended system reduces EXACTLY to the
single-valued one. So the three non-lifting gates must reproduce the
single-valued answer to machine precision on BOTH mesh families (dual-mesh
rule):

  V0  freestream (phi = U.x, full Dirichlet) captured to < 1e-12,
  V1  MMS L2 convergence slope >= 1.9 on a cube cut in generic position,
  a=0 Laplace: TE jump ~ 0 (weld forbids a jump) => cl ~ 0, and the main
      potential matches the single-valued solve_laplace oracle.

The physical wake jump ([phi] != 0, implicit Kutta) is B3 -- not exercised
here. Dual-mesh: the wake-embedded M0/M1 meshes (every wake-plane node on
the level set, the eps side-shift at scale) and the wake-free M3/M4 meshes
(generic cuts, the workflow form). The M6 (3D) families are gitignored, so
those parametrizations skip unless the meshes were generated locally.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.dirichlet import freestream_phi
from pyfp3d.mesh.metrics import compute_tet_volumes
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.wing3d import B_SEMI, x_te
from pyfp3d.solve.picard import solve_laplace
from pyfp3d.solve.picard_ls import solve_multivalued_laplace
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

from .mesh_utils import cube_boundary_mask, generate_structured_cube_mesh
from .test_laplace_mms import (
    _consistent_load_vector,
    _lumped_nodal_volumes,
    phi_exact_fn,
)

REPO_ROOT = Path(__file__).parent.parent
M0_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_2.5d"
M3_DIR = REPO_ROOT / "cases" / "meshes" / "naca0012_wakefree_2.5d"
M1_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"
M4_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wakefree"


def _wall_nodes(mesh):
    return np.unique(mesh.boundary_faces["wall"])


def _all_boundary_nodes(mesh):
    return np.unique(
        np.concatenate([f.ravel() for f in mesh.boundary_faces.values()])
    )


def _naca_levelset(mesh, alpha_deg=0.0):
    z = mesh.nodes[:, 2]
    a = np.radians(alpha_deg)
    te = np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]])
    return WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))


def _m6_levelset(alpha_deg=0.0):
    a = np.radians(alpha_deg)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    return WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))


def _require(path):
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored; see B1 test header)")
    return read_mesh(path)


def _mvop(mesh, wls):
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=_wall_nodes(mesh))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm)


# ---------------------------------------------------------------------------
# Synthetic consistency pin (no mesh files): the welded double-assembly IS
# the single-valued stiffness matrix.
# ---------------------------------------------------------------------------

class TestWeldReducesToSingleValued:
    def _two_tet_strip(self):
        # wake plane y = 0, TE at (1,0,0), direction +x -- both tets cut.
        nodes = np.array([
            [1.5, -0.5, 0.0], [2.5, -0.5, 0.0], [1.5, 0.5, 0.0],
            [1.5, 0.0, 1.0], [2.5, 0.5, 1.0],
        ])
        elements = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
        return nodes, elements

    def test_extended_matrix_folds_to_stiffness(self):
        from pyfp3d.kernels.jacobian import PicardOperator

        nodes, elements = self._two_tet_strip()
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(nodes, elements, wls)
        mvop = MultivaluedOperator(nodes, elements, cm)
        A_ext = mvop.assemble_matrix().toarray()

        # Fold every aux column/row back onto its main node (aux_k = main_j)
        # and drop the aux rows: must equal the single-valued stiffness.
        n_main = cm.n_main
        node_of_aux = np.full(cm.n_total_dofs, -1, dtype=np.int64)
        cut_nodes = np.flatnonzero(cm.ext_dof_of_node >= 0)
        node_of_aux[cm.ext_dof_of_node[cut_nodes]] = cut_nodes

        A_std = PicardOperator(nodes, elements).assemble_matrix().toarray()
        folded = A_ext[:n_main, :n_main].copy()
        for k in range(n_main, cm.n_total_dofs):
            folded[:n_main, node_of_aux[k]] += A_ext[:n_main, k]
        assert np.allclose(folded, A_std, atol=1e-13)

    def test_aux_rows_are_the_weld(self):
        nodes, elements = self._two_tet_strip()
        wls = WakeLevelSet([1.0, 0.0, 0.0], direction=(1.0, 0.0, 0.0))
        cm = CutElementMap(nodes, elements, wls)
        mvop = MultivaluedOperator(nodes, elements, cm)
        A = mvop.assemble_matrix().tocsr()
        for j in np.flatnonzero(cm.ext_dof_of_node >= 0):
            k = cm.ext_dof_of_node[j]
            row = A[k].toarray().ravel()
            assert row[k] == pytest.approx(1.0)
            assert row[j] == pytest.approx(-1.0)
            assert np.count_nonzero(row) == 2  # weld: aux_k - main_j = 0


# ---------------------------------------------------------------------------
# V0 -- freestream preservation on the cut mesh (both families, both alpha).
# ---------------------------------------------------------------------------

class TestV0Freestream:
    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    @pytest.mark.parametrize("level", ["coarse", "medium"])
    @pytest.mark.parametrize("alpha", [0.0, 4.0])
    def test_naca_2p5d(self, directory, level, alpha):
        mesh = _require(directory / f"{level}.msh")
        mvop = _mvop(mesh, _naca_levelset(mesh, alpha))
        b = _all_boundary_nodes(mesh)
        r = solve_multivalued_laplace(mvop, b, freestream_phi(mesh.nodes[b], alpha))
        err = np.max(np.abs(r["phi"] - freestream_phi(mesh.nodes, alpha)))
        assert err < 1e-12, f"{directory.name} a={alpha}: V0 err {err:.2e}"
        # the weld carries no jump on a single-valued field
        assert np.max(np.abs(r["te_jump"])) < 1e-12

    @pytest.mark.parametrize("directory", [M1_DIR, M4_DIR])
    def test_onera_m6_3d(self, directory):
        mesh = _require(directory / "coarse.msh")
        mvop = _mvop(mesh, _m6_levelset())
        b = _all_boundary_nodes(mesh)
        r = solve_multivalued_laplace(mvop, b, freestream_phi(mesh.nodes[b], 0.0))
        err = np.max(np.abs(r["phi"] - freestream_phi(mesh.nodes, 0.0)))
        assert err < 1e-12, f"{directory.name}: 3D V0 err {err:.2e}"
        assert np.max(np.abs(r["te_jump"])) < 1e-12


# ---------------------------------------------------------------------------
# V1 -- MMS convergence on a cube cut in generic position.
# ---------------------------------------------------------------------------

class TestV1MMS:
    def _run(self, n):
        nodes, elements = generate_structured_cube_mesh(n=n, L=1.0)
        dn = np.where(cube_boundary_mask(nodes, L=1.0))[0]
        pe = phi_exact_fn(nodes[:, 0], nodes[:, 1], nodes[:, 2])
        vol = compute_tet_volumes(nodes, elements)
        load = _consistent_load_vector(nodes, elements, vol)
        # tilted wake half-plane -> cuts pass through elements in general
        # position (not along the structured grid planes).
        te = np.array([[0.3, 0.5, 0.0], [0.3, 0.5, 1.0]])
        a = np.radians(8.0)
        wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
        cm = CutElementMap(nodes, elements, wls)
        mvop = MultivaluedOperator(nodes, elements, cm)
        r = solve_multivalued_laplace(mvop, dn, pe[dn], body_source_rhs=load)
        lump = _lumped_nodal_volumes(elements, vol, len(nodes))
        l2 = float(np.sqrt(np.sum(lump * (r["phi"] - pe) ** 2) / np.sum(lump)))
        return 1.0 / n, l2, len(cm.cut_elems), r["residual_norm"]

    def test_mms_slope(self):
        levels = [self._run(n) for n in (4, 8, 16)]
        for h, l2, n_cut, res in levels:
            assert n_cut > 0, "cut set empty -- MMS would not test the cut path"
            assert res < 1e-8, f"h={h}: linear solve not converged ({res:.2e})"
        h = np.array([lvl[0] for lvl in levels])
        err = np.array([lvl[1] for lvl in levels])
        slope = np.polyfit(np.log(h), np.log(err), 1)[0]
        assert slope >= 1.9, f"MMS slope {slope:.3f} < 1.9 (errors {err})"


# ---------------------------------------------------------------------------
# alpha = 0 -- non-lifting Laplace gives no jump (cl ~ 0) and matches the
# single-valued oracle.
# ---------------------------------------------------------------------------

class TestLaplaceAlpha0:
    @pytest.mark.parametrize("directory", [M0_DIR, M3_DIR])
    @pytest.mark.parametrize("level", ["coarse", "medium"])
    def test_no_lift_matches_single_valued(self, directory, level):
        mesh = _require(directory / f"{level}.msh")
        mvop = _mvop(mesh, _naca_levelset(mesh, 0.0))
        ff = np.unique(mesh.boundary_faces["farfield"])
        phi_ff = freestream_phi(mesh.nodes[ff], 0.0)
        r = solve_multivalued_laplace(mvop, ff, phi_ff)

        gamma = r["te_jump"]
        cl_kj = 2.0 * float(np.mean(gamma))  # chord = 1
        assert np.max(np.abs(gamma)) < 1e-10, "a=0 must carry no circulation"
        assert abs(cl_kj) < 1e-10

        oracle = solve_laplace(mesh.nodes, mesh.elements, ff, phi_ff,
                               rtol=1e-12, maxiter=5000)["phi"]
        assert np.max(np.abs(r["phi"] - oracle)) < 1e-9


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
