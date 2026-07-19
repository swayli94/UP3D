"""
P11 curved wall-adjacent elements (solve/curved_wall.py; design.md Sec 5.1.3,
roadmap track_p.md Sec P11) -- structural locks + the measured NEGATIVE.

The route's outcome, measured 2026-07-19 (demo cases/demo/p11_curved_walls/):
mapped-P1 curved elements move the medium-sphere Cp error only 11.56% ->
11.33% (= the G1.4 Option A oracle ceiling), because the pre-registered
superparametric-consistency risk fires -- the mapped-P1 basis on quadratic
geometry carries an O(h) gradient-consistency error of the same order as the
O(h) facet-normal error it removes. The deeper finding (E5/E6/E8 in the
demo): the G1.6 error was never dominated by the wall's geometric crime in
the first place -- the h_min-sweep order collapse is the fixed-bulk-mesh
pollution floor, and the medium mesh's 11.6% is the intrinsic P1-field
capability at h = 0.08 (structured-shell control: ~11%).

These tests lock:
  1. the null test: a planar wall + planar projection produces a delta that
     is EXACTLY zero (both stiffness paths share one code path);
  2. quadrature correctness: the flat-geometry quadrature stiffness matches
     the independent P1 reference assembly to machine precision;
  3. quadrature-rule insensitivity (deg2 vs deg3);
  4. the superparametric premise: linear fields are NOT reproduced on curved
     elements, with an O(h)-sized gradient deviation (the mechanism of the
     negative -- if this ever reads near machine zero, the basis changed and
     the P11 verdict must be re-examined);
  5. the measured anchors: coarse/medium curved-path Cp error, and the fact
     that the curved route does NOT close the < 2% gate (the G1.6 xfail in
     test_laplace_sphere.py stays authoritative for the gate itself);
  6. the structured-shell control: the icosphere family's wall phi error
     converges at ~second order with FLAT facets -- the measurement that
     re-attributes G1.6 away from the boundary variational crime.
"""

import numpy as np
import pytest

from pyfp3d.kernels.residual import assemble_stiffness_matrix_reference
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.surface import wall_tangential_gradient_quadratic
from pyfp3d.solve.curved_wall import (
    _L,
    _RULES,
    _tet10_dshape,
    assemble_curved_stiffness_delta,
    curved_volumes,
    curved_wall_geometry,
    element_stiffness_tet10,
    plane_closest_point_normal,
)
from pyfp3d.solve.picard import solve_laplace
from pyfp3d.solve.wall_correction import sphere_closest_point_normal

from .mesh_utils import generate_sphere_shell_mesh, generate_structured_cube_mesh

A = 1.0


def _boundary_tris_on_plane(nodes, elements, axis, value, tol=1e-12):
    from collections import defaultdict

    count = defaultdict(int)
    face_defs = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    for tet in np.asarray(elements):
        for fd in face_defs:
            tri = tuple(sorted((int(tet[fd[0]]), int(tet[fd[1]]), int(tet[fd[2]]))))
            count[tri] += 1
    tris = np.array([t for t, c in count.items() if c == 1], dtype=np.int64)
    on = np.all(np.abs(nodes[tris][:, :, axis] - value) < tol, axis=1)
    return tris[on]


def _run_sphere(mesh_path, curved, rule="deg2"):
    mesh = read_mesh(mesh_path)
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])
    r = np.linalg.norm(nodes, axis=1)
    phi_exact = nodes[:, 0] * (1.0 + 0.5 * A**3 / r**3)

    delta = None
    if curved:
        geo = curved_wall_geometry(
            nodes, elements, wall_faces,
            lambda p: sphere_closest_point_normal(p, radius=A),
        )
        delta = assemble_curved_stiffness_delta(len(nodes), elements, geo, rule=rule)
    result = solve_laplace(
        nodes, elements, farfield_nodes, phi_exact[farfield_nodes],
        stiffness_delta=delta, rtol=1e-11, maxiter=3000,
    )
    phi = result["phi"]
    grad_wall = wall_tangential_gradient_quadratic(nodes, wall_faces, phi)
    q2 = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
    cp = 1.0 - q2
    cos_t = nodes[wall_nodes, 0] / r[wall_nodes]
    cp_exact = 1.0 - 2.25 * (1.0 - cos_t**2)
    return float(np.abs(cp - cp_exact).max()), result


class TestStructuralLocks:
    def test_null_delta_on_planar_wall(self):
        """Planar wall + planar projection: midpoints project onto themselves
        bitwise, so dA must be EXACTLY zero (curved and flat stiffness share
        one code path). Also covers mixed tet orientations: the structured
        cube mesh carries no orientation guarantee."""
        nodes, elements = generate_structured_cube_mesh(6)
        wall = _boundary_tris_on_plane(nodes, elements, axis=2, value=0.0)
        cb = plane_closest_point_normal((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
        geo = curved_wall_geometry(nodes, elements, wall, cb)
        assert geo["max_offset"] == 0.0
        dA = assemble_curved_stiffness_delta(len(nodes), elements, geo)
        assert dA.nnz == 0 or np.abs(dA.data).max() == 0.0

    def test_flat_quadrature_matches_p1_reference(self, mesh_dir):
        """The tet10 quadrature stiffness on STRAIGHT geometry equals the
        independent serial P1 triplet assembly (constant integrand: every
        rule is exact)."""
        import scipy.sparse as sp

        mesh = read_mesh(mesh_dir / "sphere_shell" / "coarse.msh")
        geo = curved_wall_geometry(
            mesh.nodes, mesh.elements, mesh.boundary_faces["wall"],
            lambda p: sphere_closest_point_normal(p, radius=A),
        )
        K_flat = element_stiffness_tet10(geo["geom10_flat"])
        sub = np.asarray(mesh.elements)[geo["curved_tets"]]
        uniq, sub_l = np.unique(sub, return_inverse=True)
        sub_l = sub_l.reshape(sub.shape).astype(np.int32)
        A_ref = assemble_stiffness_matrix_reference(mesh.nodes[uniq], sub_l)
        rows = np.repeat(sub_l.astype(np.int64), 4, axis=1).ravel()
        cols = np.tile(sub_l.astype(np.int64), (1, 4)).ravel()
        A_q = sp.coo_matrix((K_flat.ravel(), (rows, cols)), shape=A_ref.shape).tocsr()
        rel = np.abs(A_q - A_ref).max() / np.abs(A_ref).max()
        assert rel < 1e-12

    def test_quadrature_rule_insensitivity(self, mesh_dir):
        """deg2 vs deg3 element stiffness on the curved coarse-sphere layer
        differ only at quadrature-truncation level (the integrand is smooth
        and the map perturbation O(h^2))."""
        mesh = read_mesh(mesh_dir / "sphere_shell" / "coarse.msh")
        geo = curved_wall_geometry(
            mesh.nodes, mesh.elements, mesh.boundary_faces["wall"],
            lambda p: sphere_closest_point_normal(p, radius=A),
        )
        K2 = element_stiffness_tet10(geo["geom10"], rule="deg2")
        K3 = element_stiffness_tet10(geo["geom10"], rule="deg3")
        scale = np.abs(K2).max()
        assert np.abs(K2 - K3).max() / scale < 1e-4

    def test_curved_layer_removes_sliver_volume(self, mesh_dir):
        """The curved layer's volume DECREASES vs flat (the fluid-side
        slivers between the chordal facets and the sphere are removed), by
        an amount consistent with the O(h^2) gap: coarse ~0.060, well under
        1% of the layer volume at medium."""
        for level, lo, hi in (("coarse", 0.04, 0.09), ("medium", 0.005, 0.02)):
            mesh = read_mesh(mesh_dir / "sphere_shell" / f"{level}.msh")
            geo = curved_wall_geometry(
                mesh.nodes, mesh.elements, mesh.boundary_faces["wall"],
                lambda p: sphere_closest_point_normal(p, radius=A),
            )
            removed = float(
                curved_volumes(geo["geom10_flat"]).sum()
                - curved_volumes(geo["geom10"]).sum()
            )
            assert lo < removed < hi, f"{level}: removed sliver volume {removed}"


class TestSuperparametricPremise:
    def test_linear_field_not_reproduced_on_curved_elements(self, mesh_dir):
        """THE mechanism of the P11 negative, locked as a premise: the
        mapped-P1 basis on curved elements does NOT reproduce u = x -- the
        gradient deviation is O(h)-sized (measured 0.138 max on coarse), the
        same order as the facet-normal error the curving removes. If this
        assert ever fails LOW, the field basis changed (e.g. an isoparametric
        P2 layer landed) and the P11 curved-route verdict must be
        re-examined."""
        mesh = read_mesh(mesh_dir / "sphere_shell" / "coarse.msh")
        geo = curved_wall_geometry(
            mesh.nodes, mesh.elements, mesh.boundary_faces["wall"],
            lambda p: sphere_closest_point_normal(p, radius=A),
        )
        bary, _ = _RULES["deg2"]
        dN = _tet10_dshape(bary)
        J = np.einsum("cai,aqj->cqij", geo["geom10"], dN)
        Jinv = np.linalg.inv(J)
        G = np.einsum("cqji,kj->cqki", Jinv, _L)
        u = geo["geom10"][:, :4, 0]  # u = x at the vertices
        grad_u = np.einsum("cqki,ck->cqi", G, u)
        dev = np.linalg.norm(grad_u - np.array([1.0, 0.0, 0.0]), axis=-1)
        assert 0.05 < dev.max() < 0.5
        assert 0.005 < dev.mean() < 0.1


class TestMeasuredAnchors:
    def test_sphere_coarse_curved_anchor(self, mesh_dir):
        """Curved-path coarse Cp error: measured 0.2035 (flat: 0.2305) --
        an improvement, but far from the 2% gate."""
        err, result = _run_sphere(mesh_dir / "sphere_shell" / "coarse.msh", curved=True)
        assert result["residual_norm"] < 1e-6
        assert abs(err - 0.2035) < 0.005

    def test_sphere_medium_curved_negative(self, mesh_dir):
        """The P11 negative, locked: curved wall elements land at 11.33% on
        medium -- the same value as the G1.4 Option A oracle ceiling (0.1133)
        and nowhere near the 2% gate. The G1.6 strict xfail in
        test_laplace_sphere.py remains the gate's lock."""
        err, _ = _run_sphere(mesh_dir / "sphere_shell" / "medium.msh", curved=True)
        assert abs(err - 0.1133) < 0.003
        assert err > 0.02  # the gate the route was built for stays open


class TestStructuredShellControl:
    def test_icosphere_flat_facets_converge_second_order(self):
        """The re-attribution measurement: an icosphere-extruded shell with
        the SAME flat-facet wall 'crime' has wall phi error converging at
        ~2nd order (measured 1.86 between s3 and s4) -- the boundary
        variational crime is NOT what limits the gmsh-family sphere case.
        (Full 3-level family + curved A/B in the P11 demo.)"""
        errs = []
        for subdiv, n_layers in ((3, 24), (4, 48)):
            nodes, elements, wall_nodes, farfield_nodes = generate_sphere_shell_mesh(
                subdivisions=subdiv, n_layers=n_layers,
                r_inner=1.0, r_outer=25.0, grading=1.5,
            )
            r = np.linalg.norm(nodes, axis=1)
            phi_exact = nodes[:, 0] * (1.0 + 0.5 * A**3 / r**3)
            result = solve_laplace(
                nodes, elements, farfield_nodes, phi_exact[farfield_nodes],
                rtol=1e-11, maxiter=3000,
            )
            errs.append(np.abs(result["phi"] - phi_exact)[wall_nodes].max())
        order = np.log(errs[0] / errs[1]) / np.log(2.0)
        assert order > 1.6, f"structured-shell phi order {order:.2f}"


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
