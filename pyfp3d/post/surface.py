"""
Nodal post-processing derived from the element-constant P1 gradient.

Reference: design.md Sec 9 ("Nodal q^2, M, Cp via volume-weighted element
averages").
"""

import numpy as np

from pyfp3d.mesh.metrics import compute_tet_volumes, element_gradients


def nodal_gradient_recovery(nodes: np.ndarray, elements: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Recover a nodal gradient field from the element-constant P1 gradient by
    a volume-weighted average over each node's incident elements.

    Args:
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity
        phi: (n_nodes,) nodal potential

    Returns:
        grad_nodal: (n_nodes, 3) recovered nodal gradient
    """
    n_nodes = len(nodes)
    n_tets = len(elements)
    volumes = compute_tet_volumes(nodes, elements)

    grad_elem = np.empty((n_tets, 3), dtype=np.float64)
    for e in range(n_tets):
        grads = element_gradients(nodes, elements, e)
        grad_elem[e] = phi[elements[e]] @ grads

    grad_sum = np.zeros((n_nodes, 3), dtype=np.float64)
    weight_sum = np.zeros(n_nodes, dtype=np.float64)
    for i in range(4):
        node_ids = elements[:, i]
        np.add.at(grad_sum, node_ids, volumes[:, None] * grad_elem)
        np.add.at(weight_sum, node_ids, volumes)

    return grad_sum / weight_sum[:, None]
