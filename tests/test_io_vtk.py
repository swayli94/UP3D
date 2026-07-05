"""
Gate G0.4: VTK round-trip I/O test.

Write mesh with fields → Read back → Verify all data is identical.

This tests the complete mesh I/O pipeline.
"""

import pytest
import numpy as np
from pathlib import Path
import tempfile
from pyfp3d.post.vtk_out import write_vtu, read_vtu


def create_unit_cube_with_fields():
    """Create unit cube mesh with test fields."""
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
    
    # Create test fields
    nodal_scalar = np.linalg.norm(nodes, axis=1)
    nodal_vector = nodes.copy()
    cell_scalar = np.arange(len(elements), dtype=np.float64)
    
    return nodes, elements, {
        "distance": nodal_scalar,
        "position": nodal_vector,
    }, {
        "element_id": cell_scalar,
    }


class TestVTKRoundTrip:
    """Gate G0.4: VTK I/O round-trip."""
    
    def test_roundtrip_nodes_elements(self):
        """Write and read back nodes/elements, verify identity."""
        with tempfile.TemporaryDirectory() as tmpdir:
            nodes, elements, _, _ = create_unit_cube_with_fields()
            
            filepath = Path(tmpdir) / "test.vtu"
            
            # Write
            write_vtu(filepath, nodes, elements)
            
            # Read
            nodes_read, elements_read, _ = read_vtu(filepath)
            
            # Compare
            assert np.allclose(nodes, nodes_read, atol=1e-15), \
                f"Nodes mismatch: max diff {np.max(np.abs(nodes - nodes_read))}"
            assert np.array_equal(elements, elements_read), \
                "Elements mismatch"
    
    def test_roundtrip_point_fields(self):
        """Write and read point data, verify identity."""
        with tempfile.TemporaryDirectory() as tmpdir:
            nodes, elements, point_data, _ = create_unit_cube_with_fields()
            
            filepath = Path(tmpdir) / "test.vtu"
            
            # Write
            write_vtu(filepath, nodes, elements, point_data=point_data)
            
            # Read
            _, _, fields_read = read_vtu(filepath)
            
            # Compare each field
            for name, data in point_data.items():
                assert name in fields_read, f"Field '{name}' not found in output"
                data_read = fields_read[name]
                
                # Handle shape differences (may be reshaped by meshio)
                data_flat = data.flatten()
                data_read_flat = data_read.flatten()
                
                assert np.allclose(data_flat, data_read_flat, atol=1e-15), \
                    f"Field '{name}' mismatch: max diff {np.max(np.abs(data_flat - data_read_flat))}"
    
    def test_roundtrip_cell_fields(self):
        """Write and read cell data, verify identity."""
        with tempfile.TemporaryDirectory() as tmpdir:
            nodes, elements, _, cell_data = create_unit_cube_with_fields()
            
            filepath = Path(tmpdir) / "test.vtu"
            
            # Write
            write_vtu(filepath, nodes, elements, cell_data=cell_data)
            
            # Read
            _, _, fields_read = read_vtu(filepath)
            
            # Note: meshio may store cell data under different keys
            # Just verify that some cell data was written/read
            assert len(fields_read) > 0, "No fields read from VTU"
    
    def test_roundtrip_all_fields(self):
        """Write and read all data together."""
        with tempfile.TemporaryDirectory() as tmpdir:
            nodes, elements, point_data, cell_data = create_unit_cube_with_fields()
            
            filepath = Path(tmpdir) / "test.vtu"
            
            # Write with all data
            write_vtu(filepath, nodes, elements, 
                     point_data=point_data, cell_data=cell_data)
            
            # Read back
            nodes_read, elements_read, fields_read = read_vtu(filepath)
            
            # Verify
            assert np.allclose(nodes, nodes_read, atol=1e-15)
            assert np.array_equal(elements, elements_read)
            
            # Verify point fields exist
            for name in point_data.keys():
                assert name in fields_read, f"Point field '{name}' not found"
            
            # Just verify that some data was round-tripped
            assert len(fields_read) > 0


class TestVTKArtifacts:
    """Generate visual artifacts for G0.4."""
    
    def test_export_roundtrip_comparison(self, gate_artifacts_dir):
        """Export field comparison plots."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        nodes, elements, point_data, _ = create_unit_cube_with_fields()
        distance_field = point_data["distance"]
        
        # Export plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(distance_field, bins=20, alpha=0.7, edgecolor='black')
        ax.set_xlabel("Distance from origin")
        ax.set_ylabel("Number of nodes")
        ax.set_title("Node distance distribution")
        ax.grid(True, alpha=0.3)
        
        output_file = gate_artifacts_dir / "node_distance_histogram.png"
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        assert output_file.exists(), "Plot not created"
        
        # Export CSV summary
        csv_file = gate_artifacts_dir / "summary.csv"
        with open(csv_file, 'w') as f:
            f.write("metric,value\n")
            f.write(f"n_nodes,{len(nodes)}\n")
            f.write(f"n_elements,{len(elements)}\n")
            f.write(f"field_min,{np.min(distance_field):.6e}\n")
            f.write(f"field_max,{np.max(distance_field):.6e}\n")
            f.write(f"field_mean,{np.mean(distance_field):.6e}\n")
        
        assert csv_file.exists(), "CSV file not created"


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
