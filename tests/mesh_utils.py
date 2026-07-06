"""
Shared synthetic mesh generators and case helpers for solver validation
gates (P1+).

The generators are dependency-free (no Gmsh) so P1 gates don't block on
mesh generation tooling: a structured cube (Kuhn triangulation, for MMS
convergence studies) and a spherical shell (icosphere extruded radially,
for the incompressible-sphere Cp gate). The cylinder_2.5d case helpers
(analytic solution + end-to-end solve/Cp recovery) are shared between the
M0 pipeline validation and the G1.3 cylinder oracle pre-study.
"""

import numpy as np


def generate_structured_cube_mesh(n: int, L: float = 1.0) -> tuple:
    """
    Structured n x n x n cube mesh, each sub-cube split into 6 tets via the
    Kuhn (Freudenthal) triangulation.

    Using the same global axis order (x, then y, then z) in every sub-cube
    makes the triangulation automatically face-conforming across the whole
    grid -- no cracks or T-junctions.

    Args:
        n: number of subdivisions per side
        L: cube edge length

    Returns:
        (nodes, elements): nodes (( n+1)**3, 3), elements (6*n**3, 4)
    """
    lin = np.linspace(0.0, L, n + 1)
    X, Y, Z = np.meshgrid(lin, lin, lin, indexing="ij")
    nodes = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1).astype(np.float64)

    def idx(i, j, k):
        return i * (n + 1) * (n + 1) + j * (n + 1) + k

    unit_offsets = {0: np.array([1, 0, 0]), 1: np.array([0, 1, 0]), 2: np.array([0, 0, 1])}
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]

    elements = np.empty((6 * n * n * n, 4), dtype=np.int32)
    e = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                base = np.array([i, j, k])
                for perm in perms:
                    v0 = base
                    v1 = v0 + unit_offsets[perm[0]]
                    v2 = v1 + unit_offsets[perm[1]]
                    v3 = v2 + unit_offsets[perm[2]]
                    elements[e] = [idx(*v0), idx(*v1), idx(*v2), idx(*v3)]
                    e += 1

    return nodes, elements


def cube_boundary_mask(nodes: np.ndarray, L: float = 1.0, tol: float = 1e-9) -> np.ndarray:
    """Boolean mask of nodes lying on the boundary of a [0, L]^3 cube."""
    lo = np.any(nodes < tol, axis=1)
    hi = np.any(nodes > (L - tol), axis=1)
    return lo | hi


def _icosahedron() -> tuple:
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array(
        [
            [-1, phi, 0], [1, phi, 0], [-1, -phi, 0], [1, -phi, 0],
            [0, -1, phi], [0, 1, phi], [0, -1, -phi], [0, 1, -phi],
            [phi, 0, -1], [phi, 0, 1], [-phi, 0, -1], [-phi, 0, 1],
        ],
        dtype=np.float64,
    )
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)
    faces = np.array(
        [
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
        ],
        dtype=np.int64,
    )
    return verts, faces


def icosphere(subdivisions: int) -> tuple:
    """Unit-sphere triangulation: icosahedron + recursive midpoint subdivision."""
    verts, faces = _icosahedron()
    verts = list(verts)

    for _ in range(subdivisions):
        midpoint_cache = {}

        def get_midpoint(i1, i2):
            key = (min(i1, i2), max(i1, i2))
            if key in midpoint_cache:
                return midpoint_cache[key]
            v = (verts[i1] + verts[i2]) / 2.0
            v /= np.linalg.norm(v)
            verts.append(v)
            new_idx = len(verts) - 1
            midpoint_cache[key] = new_idx
            return new_idx

        new_faces = []
        for (a, b, c) in faces:
            ab = get_midpoint(a, b)
            bc = get_midpoint(b, c)
            ca = get_midpoint(c, a)
            new_faces.append((a, ab, ca))
            new_faces.append((b, bc, ab))
            new_faces.append((c, ca, bc))
            new_faces.append((ab, bc, ca))
        faces = np.array(new_faces, dtype=np.int64)

    return np.array(verts, dtype=np.float64), faces


def generate_sphere_shell_mesh(
    subdivisions: int = 2,
    n_layers: int = 16,
    r_inner: float = 1.0,
    r_outer: float = 25.0,
    grading: float = 1.5,
) -> tuple:
    """
    Tetrahedral mesh of the shell between two concentric spheres.

    An icosphere surface pattern is extruded radially in graded layers; each
    prism between consecutive layers is split into 3 tets by sorting the
    triangle's surface-vertex indices (p < q < r) and cutting
    (p,q,r,r'), (p,q,r',q'), (p,q',r',p'). Because the diagonal used on any
    shared vertical quad face depends only on the (unordered) pair of surface
    indices at its corners -- not on which of the two adjacent prisms
    generated it -- this is watertight by construction.

    Returns:
        (nodes, elements, wall_nodes, farfield_nodes)
    """
    verts, faces = icosphere(subdivisions)
    n_v = len(verts)

    k = np.arange(n_layers + 1)
    radii = r_inner + (r_outer - r_inner) * (k / n_layers) ** grading

    nodes = np.vstack([verts * r for r in radii]).astype(np.float64)

    elements = np.empty((len(faces) * n_layers * 3, 4), dtype=np.int32)
    e = 0
    for layer in range(n_layers):
        off0 = layer * n_v
        off1 = (layer + 1) * n_v
        for (a, b, c) in faces:
            p, q, r = sorted((int(a), int(b), int(c)))
            P, Q, R = p + off0, q + off0, r + off0
            Pp, Qp, Rp = p + off1, q + off1, r + off1
            elements[e] = [P, Q, R, Rp]
            elements[e + 1] = [P, Q, Rp, Qp]
            elements[e + 2] = [P, Qp, Rp, Pp]
            e += 3

    wall_nodes = np.arange(n_v, dtype=np.int64)
    farfield_nodes = np.arange(n_layers * n_v, (n_layers + 1) * n_v, dtype=np.int64)

    return nodes, elements, wall_nodes, farfield_nodes


# ---------------------------------------------------------------------------
# Shared cylinder_2.5d case helpers (used by test_m0_cylinder.py and the
# G1.3 cylinder oracle pre-study, test_wall_correction_cylinder.py)
# ---------------------------------------------------------------------------

CYLINDER_RADIUS = 1.0


def element_gradients_all(nodes, elements, phi):
    """Constant P1 gradient of phi in every tet, vectorized."""
    p, el = nodes, elements
    e = np.stack(
        [p[el[:, 1]] - p[el[:, 0]],
         p[el[:, 2]] - p[el[:, 0]],
         p[el[:, 3]] - p[el[:, 0]]], axis=1
    )
    d = np.stack(
        [phi[el[:, 1]] - phi[el[:, 0]],
         phi[el[:, 2]] - phi[el[:, 0]],
         phi[el[:, 3]] - phi[el[:, 0]]], axis=1
    )
    return np.linalg.solve(e, d[:, :, None])[:, :, 0]


def cylinder_phi_exact(nodes, a=CYLINDER_RADIUS):
    """phi = U x (1 + a^2/r^2) for unit freestream past a unit cylinder."""
    r2 = nodes[:, 0] ** 2 + nodes[:, 1] ** 2
    return nodes[:, 0] * (1.0 + a**2 / r2)


def cylinder_grad_exact(points, a=CYLINDER_RADIUS):
    """Analytic gradient of cylinder_phi_exact at arbitrary points, (m, 3)."""
    x, y = points[:, 0], points[:, 1]
    r2 = x**2 + y**2
    grad = np.zeros_like(points, dtype=np.float64)
    grad[:, 0] = 1.0 + a**2 / r2 - 2.0 * a**2 * x**2 / r2**2
    grad[:, 1] = -2.0 * a**2 * x * y / r2**2
    return grad


def run_cylinder_case(mesh_path, body_source_rhs=None, mesh=None):
    """Solve incompressible flow past the quasi-2D cylinder and recover wall
    Cp(theta); pass body_source_rhs to solve with a wall-flux correction."""
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.post.surface import wall_tangential_gradient_quadratic
    from pyfp3d.solve.picard import solve_laplace

    if mesh is None:
        mesh = read_mesh(mesh_path)
    nodes, elements = mesh.nodes, mesh.elements
    wall_faces = mesh.boundary_faces["wall"]
    wall_nodes = np.unique(wall_faces)
    farfield_nodes = np.unique(mesh.boundary_faces["farfield"])

    phi_exact = cylinder_phi_exact(nodes)

    result = solve_laplace(
        nodes, elements, farfield_nodes, phi_exact[farfield_nodes],
        body_source_rhs=body_source_rhs, rtol=1e-11, maxiter=3000,
    )
    phi = result["phi"]

    grad_wall = wall_tangential_gradient_quadratic(nodes, wall_faces, phi)
    q_squared = np.sum(grad_wall[wall_nodes] ** 2, axis=1)
    cp_numeric = 1.0 - q_squared

    r2 = nodes[:, 0] ** 2 + nodes[:, 1] ** 2
    sin2 = nodes[wall_nodes, 1] ** 2 / r2[wall_nodes]
    cp_exact = 1.0 - 4.0 * sin2

    return {
        "mesh": mesh,
        "phi": phi,
        "wall_nodes": wall_nodes,
        "sin2": sin2,
        "cp_numeric": cp_numeric,
        "cp_exact": cp_exact,
        "error": np.abs(cp_numeric - cp_exact),
        "n_cg_iterations": result["n_cg_iterations"],
        "residual_norm": result["residual_norm"],
    }
