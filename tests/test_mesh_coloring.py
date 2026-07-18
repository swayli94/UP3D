"""
Gate G0.3: Element graph coloring validity test.

Verify that no two same-color elements share a node.

This enables safe @prange parallelization in assembly loops.
"""

import pytest
import numpy as np
from pyfp3d.mesh.coloring import greedy_coloring, validate_coloring, partition_by_color


def create_unit_cube_elements():
    """5 tets decomposing unit cube."""
    return np.array([
        [0, 1, 3, 4],
        [1, 2, 3, 6],
        [1, 3, 4, 6],
        [3, 4, 6, 7],
        [1, 4, 5, 6],
    ], dtype=np.int32)


def create_simple_tet_chain():
    """Linear chain of 3 tets sharing faces."""
    return np.array([
        [0, 1, 2, 3],
        [1, 2, 3, 4],
        [2, 3, 4, 5],
    ], dtype=np.int32)


class TestElementColoring:
    """Gate G0.3: Element graph coloring validity."""
    
    def test_coloring_validity(self):
        """Verify no same-color elements share a node."""
        elements = create_unit_cube_elements()
        colors, n_colors = greedy_coloring(elements)
        
        # Manual validation
        assert validate_coloring(elements, colors), "Coloring is invalid"
    
    def test_coloring_count(self):
        """Coloring should use reasonable number of colors (≤ 6 for 3D tets)."""
        elements = create_unit_cube_elements()
        colors, n_colors = greedy_coloring(elements)
        
        # Upper bound for tetrahedral mesh in 3D
        assert n_colors <= 6, f"Too many colors: {n_colors}"
        assert n_colors >= 1, f"Invalid color count: {n_colors}"
    
    def test_simple_chain_coloring(self):
        """Test coloring on simple chain of tets."""
        elements = create_simple_tet_chain()
        colors, n_colors = greedy_coloring(elements)
        
        assert validate_coloring(elements, colors), "Coloring is invalid"

        # Chain of 3 tets sharing nodes should need at least 2 colors
        assert n_colors >= 2, f"Chain coloring used {n_colors} color(s), expected ≥2"
    
    def test_coloring_deterministic(self):
        """Coloring should be deterministic."""
        elements = create_unit_cube_elements()
        
        colors1, n_colors1 = greedy_coloring(elements)
        colors2, n_colors2 = greedy_coloring(elements)
        
        assert np.array_equal(colors1, colors2), "Coloring is non-deterministic"
        assert n_colors1 == n_colors2, "Color count changed"
    
    def test_coloring_partitioning(self):
        """Partition elements by color and verify grouping."""
        elements = create_unit_cube_elements()
        colors, n_colors = greedy_coloring(elements)
        
        color_lists = partition_by_color(colors)
        
        assert len(color_lists) == n_colors, f"Partition count mismatch"
        
        # Verify all elements are covered exactly once
        total_covered = sum(len(cl) for cl in color_lists)
        assert total_covered == len(elements), "Elements not covered exactly once"
        
        # Verify each partition contains elements of one color
        for color, color_elements in enumerate(color_lists):
            for elem_idx in color_elements:
                assert colors[elem_idx] == color, f"Element {elem_idx} in wrong color partition"


class TestColoringArtifacts:
    """Generate visual artifacts for G0.3."""
    
    def test_export_coloring_visualization(self, gate_artifacts_dir):
        """Export 3D coloring visualization."""
        from pyfp3d.post.vtk_out import write_vtu
        
        # Unit cube mesh
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
        
        elements = create_unit_cube_elements()
        colors, n_colors = greedy_coloring(elements)
        
        # Convert element colors to nodal field (average)
        node_colors = np.zeros(len(nodes), dtype=np.float64)
        node_count = np.zeros(len(nodes), dtype=np.int32)
        
        for e in range(len(elements)):
            for local_node in range(4):
                node_id = elements[e, local_node]
                node_colors[node_id] += colors[e]
                node_count[node_id] += 1
        
        node_colors /= np.maximum(node_count, 1)
        
        # Write VTU
        output_file = gate_artifacts_dir / "element_colors_3d.vtu"
        write_vtu(
            output_file,
            nodes,
            elements,
            point_data={"color_index": node_colors},
            verbose=True,
        )
        
        assert output_file.exists(), "VTU file not created"
        
        # Write CSV summary
        csv_file = gate_artifacts_dir / "summary.csv"
        with open(csv_file, 'w') as f:
            f.write("metric,value\n")
            f.write(f"n_colors,{n_colors}\n")
            f.write(f"n_elements,{len(elements)}\n")
            for color in range(n_colors):
                count = np.sum(colors == color)
                f.write(f"color_{color}_count,{count}\n")
        
        assert csv_file.exists(), "CSV file not created"


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
