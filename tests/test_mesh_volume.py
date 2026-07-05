"""
Gate G0.1: Volume conservation test.

Verify that ΣV_e = exact domain volume to machine precision.

This is the first mesh infrastructure gate.
"""

import pytest
import numpy as np
from pathlib import Path
from pyfp3d.mesh.metrics import compute_tet_volumes


def create_unit_cube_mesh():
    """Create unit cube decomposed into 5 tetrahedra."""
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
    
    elements = np.array([
        [0, 1, 3, 4],
        [1, 2, 3, 6],
        [1, 3, 4, 6],
        [3, 4, 6, 7],
        [1, 4, 5, 6],
    ], dtype=np.int32)
    
    return nodes, elements


def create_unit_octahedron_mesh():
    """Create a regular octahedron (vertex-to-center distance 1), fan-tessellated
    into 8 tets from the origin. Analytic volume is 4/3 * r^3 = 4/3."""
    nodes = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ], dtype=np.float64)

    # 8 tets: fan from origin covering both the +z and -z apex pyramids
    elements = np.array([
        [0, 1, 2, 3],
        [0, 2, 4, 3],
        [0, 4, 5, 3],
        [0, 5, 1, 3],
        [0, 1, 2, 6],
        [0, 2, 4, 6],
        [0, 4, 5, 6],
        [0, 5, 1, 6],
    ], dtype=np.int32)

    return nodes, elements


class TestVolumeConservation:
    """Gate G0.1: Mesh volume conservation."""
    
    def test_unit_cube_volume(self):
        """Unit cube should have volume = 1.0."""
        nodes, elements = create_unit_cube_mesh()
        volumes = compute_tet_volumes(nodes, elements)
        
        total_volume = np.sum(volumes)
        error = abs(total_volume - 1.0)
        
        assert error < 1e-12, f"Cube volume error: {error:.2e}, got {total_volume}"
    
    def test_unit_octahedron_volume(self):
        """Regular octahedron (vertex-to-center distance 1) has volume 4/3."""
        nodes, elements = create_unit_octahedron_mesh()
        volumes = compute_tet_volumes(nodes, elements)

        total_volume = np.sum(volumes)
        expected_volume = 4.0 / 3.0
        error = abs(total_volume - expected_volume)

        assert error < 1e-12, \
            f"Octahedron volume error: {error:.2e}, expected {expected_volume}"
    
    def test_positive_volumes(self):
        """All element volumes must be positive."""
        nodes, elements = create_unit_cube_mesh()
        volumes = compute_tet_volumes(nodes, elements)
        
        assert np.all(volumes > 0), "Negative or zero volume detected"
        assert np.all(volumes < 1.0), "Individual tet volume > 1 (impossible for unit cube)"
    
    def test_scaled_mesh_volume(self):
        """Scaling nodes by factor α should scale volume by α³."""
        nodes, elements = create_unit_cube_mesh()
        volumes_original = compute_tet_volumes(nodes, elements)
        total_original = np.sum(volumes_original)
        
        # Scale by factor 2
        scale = 2.0
        nodes_scaled = scale * nodes
        volumes_scaled = compute_tet_volumes(nodes_scaled, elements)
        total_scaled = np.sum(volumes_scaled)
        
        expected_ratio = scale**3
        actual_ratio = total_scaled / total_original
        
        assert abs(actual_ratio - expected_ratio) < 1e-12, \
            f"Scaling: expected ratio {expected_ratio}, got {actual_ratio}"


class TestVolumeArtifacts:
    """Generate visual artifacts for G0.1."""
    
    def test_export_volume_error_heatmap(self, gate_artifacts_dir):
        """Export per-element volume comparison."""
        from pyfp3d.post.vtk_out import write_vtu
        
        nodes, elements = create_unit_cube_mesh()
        volumes = compute_tet_volumes(nodes, elements)
        
        # Compute per-element error (e.g., vs uniform distribution)
        average_volume = np.mean(volumes)
        error = np.zeros(len(nodes), dtype=np.float64)
        
        # Assign element error to nodes (average)
        volume_errors = volumes - average_volume
        for e in range(len(elements)):
            for local_node in range(4):
                node_id = elements[e, local_node]
                error[node_id] += abs(volume_errors[e]) / 4.0
        
        # Write VTU
        output_file = gate_artifacts_dir / "volume_error_heatmap.vtu"
        write_vtu(
            output_file,
            nodes,
            elements,
            point_data={"volume_error": error},
            verbose=True,
        )
        
        assert output_file.exists(), "VTU file not created"
        
        # Write CSV summary
        csv_file = gate_artifacts_dir / "summary.csv"
        with open(csv_file, 'w') as f:
            f.write("metric,value\n")
            f.write(f"total_volume,{np.sum(volumes):.6e}\n")
            f.write(f"expected_volume,1.0\n")
            f.write(f"absolute_error,{abs(np.sum(volumes) - 1.0):.6e}\n")
            f.write(f"relative_error,{abs(np.sum(volumes) - 1.0):.6e}\n")
        
        assert csv_file.exists(), "CSV file not created"


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
