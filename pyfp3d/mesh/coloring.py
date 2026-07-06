"""
Element graph coloring for parallel assembly via @prange.

Greedy graph coloring algorithm: ensure no two same-color elements share a node,
allowing safe parallel updates in Numba @prange loops.

Reference: design.md §7 (Numba kernel architecture, rule 2: colored assembly),
standard graph coloring for FEM assembly.

Note: this module is pure-Python preprocessing (dict/set/list based) that runs
once at startup; it produces the color partition that the @njit assembly
kernels consume. It is intentionally NOT numba-jitted itself.
"""

import numba
import numpy as np
from typing import Tuple, List, Dict


def build_element_connectivity(elements: np.ndarray) -> Dict:
    """
    Build node-to-elements adjacency list.
    
    Args:
        elements: (n_tets, 4) element connectivity
        
    Returns:
        node_to_elements: Dict[node_id] -> array of element indices sharing that node
    """
    n_tets = len(elements)
    n_nodes = np.max(elements) + 1
    
    # Build adjacency lists
    node_to_elements: Dict = {}
    for node in range(n_nodes):
        node_to_elements[node] = []
    
    for e in range(n_tets):
        for local_node in range(4):
            global_node = elements[e, local_node]
            node_to_elements[global_node].append(e)
    
    return node_to_elements


def node_to_element_csr(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """CSR-style node -> incident elements map (offsets, element ids).

    Vectorized numpy build; feeds the numba coloring kernel and any other
    node-patch consumer that needs flat arrays instead of dict/list
    adjacency (design.md Sec 7 rule 1: SoA only inside kernels).
    """
    n_tets = len(elements)
    n_nodes = int(np.max(elements)) + 1
    flat_nodes = np.asarray(elements, dtype=np.int64).reshape(-1)
    flat_elems = np.repeat(np.arange(n_tets, dtype=np.int64), 4)

    order = np.argsort(flat_nodes, kind="stable")
    counts = np.bincount(flat_nodes, minlength=n_nodes)
    offsets = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    return offsets, flat_elems[order]


@numba.njit(cache=True)
def _greedy_coloring_kernel(
    elements: np.ndarray, node_offsets: np.ndarray, node_elems: np.ndarray
) -> np.ndarray:
    """Same greedy algorithm as the original pure-Python loop (same element
    visit order, smallest available color), so the assignment is identical
    -- but O(1000x) faster on real meshes, which lets one-shot assembly
    calls afford a fresh coloring. `mark[c] == e` flags color c as forbidden
    for element e without clearing the workspace between elements."""
    n_tets = len(elements)
    colors = np.full(n_tets, -1, dtype=np.int32)
    max_colors = 256
    mark = np.full(max_colors, -1, dtype=np.int64)

    for e in range(n_tets):
        for i in range(4):
            node = elements[e, i]
            for k in range(node_offsets[node], node_offsets[node + 1]):
                c = colors[node_elems[k]]
                if c >= 0:
                    mark[c] = e
        c = 0
        while c < max_colors and mark[c] == e:
            c += 1
        if c >= max_colors:
            raise ValueError(
                "greedy coloring exceeded 256 colors -- pathological mesh "
                "connectivity (a node shared by >255 elements)"
            )
        colors[e] = c
    return colors


def greedy_coloring(elements: np.ndarray) -> Tuple[np.ndarray, int]:
    """
    Greedy graph coloring of element connectivity.

    Ensures no two elements of the same color share a node.

    Algorithm:
      1. For each element in order
      2. Find forbidden colors (colors of adjacent elements)
      3. Assign smallest available color

    Args:
        elements: (n_tets, 4) element connectivity

    Returns:
        (colors, n_colors): element colors and number of colors used

    Note:
        Upper bound on chromatic number for tetrahedral mesh: typically 5-6 in 3D.
        This greedy algorithm is not optimal but is deterministic. Since P3
        the greedy loop runs in a numba kernel with the same visit order
        (identical assignment to the original pure-Python implementation).
    """
    node_offsets, node_elems = node_to_element_csr(elements)
    colors = _greedy_coloring_kernel(
        np.asarray(elements), node_offsets, node_elems
    )
    n_colors = int(np.max(colors) + 1)
    return colors, n_colors


def color_partition_csr(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """Colored element partition in flat CSR form for the assembly kernels.

    Returns:
        (color_offsets, color_elems, n_colors): elements of color c are
        color_elems[color_offsets[c]:color_offsets[c+1]], in ascending
        element order within each color (deterministic: with a valid
        coloring no two same-color elements share a node, so per-node and
        per-nnz scatter-add order is fixed by the color sequence alone,
        independent of prange thread scheduling).
    """
    colors, n_colors = greedy_coloring(elements)
    order = np.argsort(colors, kind="stable")
    counts = np.bincount(colors, minlength=n_colors)
    color_offsets = np.zeros(n_colors + 1, dtype=np.int64)
    np.cumsum(counts, out=color_offsets[1:])
    return color_offsets, order.astype(np.int64), n_colors


def validate_coloring(elements: np.ndarray, colors: np.ndarray) -> bool:
    """
    Verify that the coloring is valid (no same-color elements share a node).
    
    Args:
        elements: (n_tets, 4) element connectivity
        colors: (n_tets,) color of each element
        
    Returns:
        True if coloring is valid
        
    Raises:
        AssertionError if any two adjacent elements have the same color
    """
    n_tets = len(elements)
    
    # Build node-to-elements adjacency
    node_to_elements = build_element_connectivity(elements)
    
    # Check: no two adjacent elements have the same color
    for node, adjacent_elements in node_to_elements.items():
        for i in range(len(adjacent_elements)):
            for j in range(i + 1, len(adjacent_elements)):
                e1 = adjacent_elements[i]
                e2 = adjacent_elements[j]
                
                assert colors[e1] != colors[e2], \
                    f"Same-color adjacency: elements {e1} and {e2} both color {colors[e1]}"
    
    return True


def partition_by_color(colors: np.ndarray) -> List[np.ndarray]:
    """
    Partition elements by color.
    
    Useful for explicit parallel loops where each color is processed independently.
    
    Args:
        colors: (n_tets,) color of each element
        
    Returns:
        color_lists: List of arrays, each containing elements of one color
    """
    n_colors = int(np.max(colors) + 1)
    n_tets = len(colors)
    
    # Count elements per color
    color_counts = np.zeros(n_colors, dtype=np.int32)
    for e in range(n_tets):
        color_counts[colors[e]] += 1
    
    # Build arrays for each color
    color_lists = []
    offsets = np.zeros(n_colors, dtype=np.int32)
    
    for color in range(n_colors):
        color_lists.append(np.zeros(color_counts[color], dtype=np.int32))
    
    # Fill arrays
    for e in range(n_tets):
        color = colors[e]
        idx = offsets[color]
        color_lists[color][idx] = e
        offsets[color] += 1
    
    return color_lists


if __name__ == "__main__":
    # Self-test: unit cube coloring
    print("=== Element Graph Coloring Self-Test ===\n")
    
    # Unit cube decomposed into 5 tets
    elements = np.array([
        [0, 1, 3, 4],
        [1, 2, 3, 6],
        [1, 3, 4, 6],
        [3, 4, 6, 7],
        [1, 4, 5, 6],
    ], dtype=np.int32)
    
    colors, n_colors = greedy_coloring(elements)
    
    print(f"Unit cube coloring:")
    print(f"  Elements: {len(elements)}")
    print(f"  Colors used: {n_colors}")
    print(f"  Assignment: {colors}")
    
    # Validate
    is_valid = validate_coloring(elements, colors)
    print(f"  Valid coloring: {is_valid}")
    
    if is_valid:
        print("✓ Coloring test PASSED")
    else:
        print("✗ Coloring test FAILED")
    
    # Partition by color
    color_lists = partition_by_color(colors)
    print(f"\nPartitioning by color:")
    for color, elements_of_color in enumerate(color_lists):
        print(f"  Color {color}: elements {elements_of_color}")
