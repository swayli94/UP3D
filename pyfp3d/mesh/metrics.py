"""
Mesh geometry metrics: volumes, gradients, face adjacency.

Numba-jitted kernels for:
  - Tetrahedral volume and Jacobian
  - Linear basis function gradients
  - Edge and face adjacency lists
  - Element-wise metrics (dihedral angle, aspect ratio)

Reference: design.md §7 (Architecture), standard FEM geometry formulas

All functions use SoA arrays and are @njit-compatible.
"""

import numba
import numpy as np
from typing import Tuple, Dict, List


@numba.njit(cache=True)
def tet_volume_and_jacobian(nodes: np.ndarray, tet: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Compute volume and reference Jacobian for a tetrahedral element.
    
    Given 4 nodes of a tet, compute:
      V = |det(J)| / 6
      J = [x1-x0, x2-x0, x3-x0]  (3x3 matrix, shape (3,3))
    
    Args:
        nodes: (n_nodes, 3) nodal coordinates
        tet: (4,) tetrahedral node indices
        
    Returns:
        (volume, J): volume (float) and Jacobian matrix (3, 3)
    """
    # Extract node coordinates
    x0 = nodes[tet[0]]
    x1 = nodes[tet[1]]
    x2 = nodes[tet[2]]
    x3 = nodes[tet[3]]
    
    # Edge vectors
    v1 = x1 - x0  # shape (3,)
    v2 = x2 - x0
    v3 = x3 - x0
    
    # Build Jacobian matrix (column-major storage for standard orientation)
    J = np.zeros((3, 3), dtype=np.float64)
    J[:, 0] = v1
    J[:, 1] = v2
    J[:, 2] = v3
    
    # Compute determinant (scalar triple product)
    # det = v1 · (v2 × v3)
    cross_product = np.zeros(3, dtype=np.float64)
    cross_product[0] = v2[1] * v3[2] - v2[2] * v3[1]
    cross_product[1] = v2[2] * v3[0] - v2[0] * v3[2]
    cross_product[2] = v2[0] * v3[1] - v2[1] * v3[0]
    
    det = v1[0] * cross_product[0] + v1[1] * cross_product[1] + v1[2] * cross_product[2]
    
    volume = abs(det) / 6.0
    
    return volume, J


@numba.njit(cache=True)
def compute_tet_volumes(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """
    Compute volumes for all tetrahedral elements.
    
    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) element connectivity
        
    Returns:
        volumes: (n_tets,) volume of each tet
    """
    n_tets = len(elements)
    volumes = np.zeros(n_tets, dtype=np.float64)
    
    for e in range(n_tets):
        vol, _ = tet_volume_and_jacobian(nodes, elements[e])
        volumes[e] = vol
    
    return volumes


@numba.njit(cache=True)
def linear_basis_gradient(J_inv: np.ndarray, basis_node: int) -> np.ndarray:
    """
    Compute gradient of linear basis function on a tet.
    
    For tet node i in reference element [0, 1, 2, 3]:
      ∇φ_i = (J^{-1})^T · ∇_ref φ_i
    
    where ∇_ref is the reference gradient:
      ∇_ref φ_0 = [-1, -1, -1]
      ∇_ref φ_i = e_i for i=1,2,3 (standard basis vectors)
    
    Args:
        J_inv: (3, 3) inverse Jacobian matrix
        basis_node: Local node number (0, 1, 2, or 3)
        
    Returns:
        grad_phi: (3,) gradient of basis function
    """
    if basis_node == 0:
        grad_ref = np.array([-1.0, -1.0, -1.0], dtype=np.float64)
    elif basis_node == 1:
        grad_ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    elif basis_node == 2:
        grad_ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    else:  # basis_node == 3
        grad_ref = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    
    # grad_phi = (J^{-T}) · grad_ref
    grad_phi = J_inv.T @ grad_ref
    
    return grad_phi


@numba.njit(cache=True)
def element_gradients(nodes: np.ndarray, elements: np.ndarray, tet_index: int) -> np.ndarray:
    """
    Compute basis function gradients for all nodes of a single tet.
    
    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) element connectivity
        tet_index: Which element to compute
        
    Returns:
        grads: (4, 3) gradient for each of 4 basis functions
    """
    tet = elements[tet_index]
    vol, J = tet_volume_and_jacobian(nodes, tet)
    
    # Compute inverse Jacobian
    det = np.linalg.det(J)
    if abs(det) < 1e-20:
        # Degenerate element; return zero gradients
        return np.zeros((4, 3), dtype=np.float64)
    
    J_inv = np.linalg.inv(J)
    
    # Compute gradients for all 4 basis functions
    grads = np.zeros((4, 3), dtype=np.float64)
    for i in range(4):
        grads[i] = linear_basis_gradient(J_inv, i)
    
    return grads


@numba.njit(cache=True)
def build_face_adjacency(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build element-to-element adjacency via faces.
    
    Each tet has 4 faces. Two tets sharing a face are neighbors.
    
    Returns:
        face_neighbors: (n_tets, 4) array where face_neighbors[e, f] is the
                       neighbor tet across face f of element e (or -1 if boundary)
        face_orientations: (n_tets, 4) orientation flag for each shared face
        
    Note: This is a simplified version; full implementation would handle
          periodic boundary conditions and mesh boundaries.
    """
    n_tets = len(elements)
    
    # For each tet, define its 4 faces as sorted node triples
    # Face 0: nodes (1, 2, 3) opposite node 0
    # Face 1: nodes (0, 2, 3) opposite node 1
    # Face 2: nodes (0, 1, 3) opposite node 2
    # Face 3: nodes (0, 1, 2) opposite node 3
    
    face_defs = np.array(
        [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]], dtype=np.int32
    )
    
    # Build a map from sorted face to (element, local_face_index)
    face_to_tet: Dict = {}  # Key: (n0, n1, n2), Value: (tet_idx, face_idx, global_n0, ...)
    
    for e in range(n_tets):
        tet = elements[e]
        for f in range(4):
            # Get global node indices of this face
            local_nodes = face_defs[f]
            global_nodes = np.array([tet[local_nodes[0]], tet[local_nodes[1]], tet[local_nodes[2]]])
            
            # Sort for consistent key
            sorted_nodes = np.sort(global_nodes)
            key = (sorted_nodes[0], sorted_nodes[1], sorted_nodes[2])
            
            if key not in face_to_tet:
                face_to_tet[key] = []
            face_to_tet[key].append((e, f))
    
    # Build neighbor arrays
    face_neighbors = np.full((n_tets, 4), -1, dtype=np.int32)
    face_orientations = np.zeros((n_tets, 4), dtype=np.int32)
    
    for key, tet_face_list in face_to_tet.items():
        if len(tet_face_list) == 2:
            (e1, f1), (e2, f2) = tet_face_list
            face_neighbors[e1, f1] = e2
            face_neighbors[e2, f2] = e1
            # Orientations can be computed from relative face orderings
            face_orientations[e1, f1] = 1
            face_orientations[e2, f2] = -1
    
    return face_neighbors, face_orientations


@numba.njit(cache=True)
def compute_edge_lengths(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """
    Compute all edge lengths in the mesh.
    
    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) element connectivity
        
    Returns:
        edge_lengths: All edge lengths (flattened)
    """
    n_tets = len(elements)
    edge_lengths_list = []
    
    # Define the 6 edges of a tet
    edges = np.array(
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32
    )
    
    for e in range(n_tets):
        tet = elements[e]
        for edge_idx in range(6):
            n0 = tet[edges[edge_idx, 0]]
            n1 = tet[edges[edge_idx, 1]]
            
            v = nodes[n1] - nodes[n0]
            length = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
            edge_lengths_list.append(length)
    
    return np.array(edge_lengths_list, dtype=np.float64)


@numba.njit(cache=True)
def compute_aspect_ratios(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """
    Compute aspect ratio (longest / shortest edge) for each tet.
    
    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) element connectivity
        
    Returns:
        aspect_ratios: (n_tets,) aspect ratio for each element
    """
    n_tets = len(elements)
    aspect_ratios = np.zeros(n_tets, dtype=np.float64)
    
    edges = np.array(
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int32
    )
    
    for e in range(n_tets):
        tet = elements[e]
        min_length = 1e10
        max_length = 1e-10
        
        for edge_idx in range(6):
            n0 = tet[edges[edge_idx, 0]]
            n1 = tet[edges[edge_idx, 1]]
            
            v = nodes[n1] - nodes[n0]
            length = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
            
            min_length = min(min_length, length)
            max_length = max(max_length, length)
        
        if min_length > 1e-20:
            aspect_ratios[e] = max_length / min_length
        else:
            aspect_ratios[e] = 1e10  # Degenerate element
    
    return aspect_ratios


if __name__ == "__main__":
    # Self-test: unit cube volume
    print("=== Mesh Metrics Self-Test ===\n")
    
    # Define unit cube: 8 vertices, 5 tets
    nodes = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ], dtype=np.float64)
    
    # 5 tets decomposing the unit cube
    elements = np.array([
        [0, 1, 3, 4],
        [1, 2, 3, 6],
        [1, 3, 4, 6],
        [3, 4, 6, 7],
        [1, 4, 5, 6],
    ], dtype=np.int32)
    
    volumes = compute_tet_volumes(nodes, elements)
    total_volume = np.sum(volumes)
    
    print(f"Unit cube volume test:")
    print(f"  Total volume: {total_volume:.6f} (expected 1.0)")
    print(f"  Per-element: {volumes}")
    print(f"  Error: {abs(total_volume - 1.0):.2e}")
    
    if abs(total_volume - 1.0) < 1e-12:
        print("✓ Volume test PASSED")
    else:
        print("✗ Volume test FAILED")
    
    # Test gradient recovery
    print(f"\nLinear field gradient test:")
    grads = element_gradients(nodes, elements, 0)
    print(f"  Element 0 basis gradients:\n{grads}")
    print("  (For linear field, should sum to zero)")
    grad_sum = np.sum(grads, axis=0)
    print(f"  Sum: {grad_sum} (expected [0, 0, 0])")
    
    if np.allclose(grad_sum, 0):
        print("✓ Gradient test PASSED")
    else:
        print("✗ Gradient test FAILED")
