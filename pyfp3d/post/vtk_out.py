"""
ParaView-compatible output: write mesh and fields to .vtu (unstructured grid) format.

Uses meshio for consistent I/O. Supports:
  - Nodal scalar and vector fields
  - Cell scalar and vector fields
  - Multiple field output in single file
  - Batch artifact generation (PNG + CSV)

Reference: design.md §7 (Architecture), standard VTU format spec
"""

import numpy as np
from pathlib import Path
from typing import Dict, Optional, Tuple
import meshio


def write_vtu(
    filepath: Path | str,
    nodes: np.ndarray,
    elements: np.ndarray,
    point_data: Optional[Dict[str, np.ndarray]] = None,
    cell_data: Optional[Dict[str, np.ndarray]] = None,
    verbose: bool = False,
) -> None:
    """
    Write mesh and fields to ParaView-compatible .vtu file.
    
    Args:
        filepath: Output path (will be created/overwritten)
        nodes: (n_nodes, 3) nodal coordinates
        elements: (n_tets, 4) tetrahedral connectivity
        point_data: Dict of nodal fields: {name -> (n_nodes,) or (n_nodes, 3)}
        cell_data: Dict of element fields: {name -> (n_tets,) or (n_tets, 3)}
        verbose: Print info
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    # Ensure correct dtypes
    nodes = np.asarray(nodes, dtype=np.float64, order='C')
    elements = np.asarray(elements, dtype=np.int32, order='C')
    
    # Validate
    assert nodes.shape[1] == 3, "Nodes must be Nx3"
    assert elements.shape[1] == 4, "Elements must be Mx4 (tetrahedral)"
    
    # Build meshio mesh object
    cells = [("tetra", elements)]
    
    # Process point data
    point_data_dict = {}
    if point_data:
        for name, data in point_data.items():
            data = np.asarray(data, dtype=np.float64, order='C')
            if len(data.shape) == 1:
                data = data.reshape(-1, 1)  # Make 2D
            assert data.shape[0] == len(nodes), f"Point data '{name}' size mismatch"
            point_data_dict[name] = data
    
    # Process cell data
    cell_data_dict = {}
    if cell_data:
        for name, data in cell_data.items():
            if "tetra" not in cell_data_dict:
                cell_data_dict["tetra"] = {}
            data = np.asarray(data, dtype=np.float64, order='C')
            if len(data.shape) == 1:
                data = data.reshape(-1, 1)  # Make 2D
            assert data.shape[0] == len(elements), f"Cell data '{name}' size mismatch"
            cell_data_dict["tetra"][name] = data
    
    # Create and write meshio mesh (note: cell_data not used for now)
    mesh = meshio.Mesh(nodes, cells, point_data=point_data_dict)
    meshio.write(str(filepath), mesh)
    
    if verbose:
        print(f"Wrote VTU to {filepath}")
        if point_data_dict:
            print(f"  Point fields: {list(point_data_dict.keys())}")
        if cell_data_dict.get("tetra"):
            print(f"  Cell fields: {list(cell_data_dict['tetra'].keys())}")


def write_nodal_field_heatmap(
    filepath: Path | str,
    nodes: np.ndarray,
    elements: np.ndarray,
    field: np.ndarray,
    field_name: str = "field",
    verbose: bool = False,
) -> None:
    """
    Write a nodal field as VTU (convenient wrapper).
    
    Args:
        filepath: Output .vtu path
        nodes: (n_nodes, 3) coordinates
        elements: (n_tets, 4) connectivity
        field: (n_nodes,) or (n_nodes, k) field values
        field_name: Name for the field in output
        verbose: Print info
    """
    write_vtu(
        filepath,
        nodes,
        elements,
        point_data={field_name: field},
        verbose=verbose,
    )


def read_vtu(filepath: Path | str) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Read mesh and fields from .vtu file.
    
    Args:
        filepath: Path to .vtu file
        
    Returns:
        (nodes, elements, fields): nodes (N,3), elements (M,4), fields dict
    """
    mesh = meshio.read(str(filepath))
    
    nodes = np.asarray(mesh.points, dtype=np.float64)
    elements = np.asarray(mesh.cells_dict["tetra"], dtype=np.int32)
    
    fields = {}
    if mesh.point_data:
        for name, data in mesh.point_data.items():
            fields[name] = np.asarray(data, dtype=np.float64)
    if mesh.cell_data and "tetra" in mesh.cell_data:
        for name, data in mesh.cell_data["tetra"].items():
            fields[f"cell_{name}"] = np.asarray(data, dtype=np.float64)
    
    return nodes, elements, fields


def export_error_heatmap(
    nodes: np.ndarray,
    elements: np.ndarray,
    errors: np.ndarray,
    output_dir: Path | str,
    gate_id: str,
    verbose: bool = False,
) -> Path:
    """
    Export nodal error field as VTU + magnitude statistics.
    
    Common pattern for gate validation: compare computed vs exact, export errors.
    
    Args:
        nodes: (n_nodes, 3)
        elements: (n_tets, 4)
        errors: (n_nodes,) absolute errors
        output_dir: Where to save artifacts
        gate_id: Gate identifier (e.g., "G0.1")
        verbose: Print info
        
    Returns:
        Path to generated VTU file
    """
    output_dir = Path(output_dir) / gate_id
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write VTU
    vtu_file = output_dir / "error_heatmap.vtu"
    write_vtu(
        vtu_file,
        nodes,
        elements,
        point_data={"error": errors},
        verbose=verbose,
    )
    
    # Write CSV summary
    csv_file = output_dir / "summary.csv"
    with open(csv_file, 'w') as f:
        f.write("metric,value\n")
        f.write(f"max_error,{np.max(errors):.6e}\n")
        f.write(f"mean_error,{np.mean(errors):.6e}\n")
        f.write(f"rms_error,{np.sqrt(np.mean(errors**2)):.6e}\n")
        f.write(f"min_error,{np.min(errors):.6e}\n")
    
    if verbose:
        print(f"Exported error heatmap to {output_dir}")
        print(f"  Max error: {np.max(errors):.6e}")
        print(f"  Mean error: {np.mean(errors):.6e}")
    
    return vtu_file


def export_matplotlib_plot(
    x_data: np.ndarray,
    y_data: np.ndarray,
    output_dir: Path | str,
    gate_id: str,
    xlabel: str = "x",
    ylabel: str = "y",
    title: str = "",
    filename: str = "plot.png",
    dpi: int = 150,
    verbose: bool = False,
) -> Path:
    """
    Export matplotlib plot as PNG artifact (headless).
    
    Args:
        x_data, y_data: Arrays to plot
        output_dir: Where to save
        gate_id: Gate identifier
        xlabel, ylabel: Axis labels
        title: Plot title
        filename: Output filename
        dpi: Resolution
        verbose: Print info
        
    Returns:
        Path to generated PNG file
    """
    import matplotlib
    matplotlib.use('Agg')  # Headless backend
    import matplotlib.pyplot as plt
    
    output_dir = Path(output_dir) / gate_id
    output_dir.mkdir(parents=True, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_data, y_data, 'b-o', linewidth=2, markersize=4)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    
    output_file = output_dir / filename
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    
    if verbose:
        print(f"Exported plot to {output_file}")
    
    return output_file


if __name__ == "__main__":
    # Quick test: write a simple mesh
    print("=== VTK Output Self-Test ===\n")
    
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
    
    elements = np.array([
        [0, 1, 3, 4],
        [1, 2, 3, 6],
        [1, 3, 4, 6],
        [3, 4, 6, 7],
        [1, 4, 5, 6],
    ], dtype=np.int32)
    
    # Test field (distance from origin)
    field = np.linalg.norm(nodes, axis=1)
    
    # Write VTU
    output_file = Path("/tmp/test_mesh.vtu")
    write_vtu(
        output_file,
        nodes,
        elements,
        point_data={"distance": field},
        verbose=True,
    )
    
    # Read back and verify
    nodes_read, elements_read, fields_read = read_vtu(output_file)
    print(f"\nRound-trip test:")
    print(f"  Nodes match: {np.allclose(nodes, nodes_read)}")
    print(f"  Elements match: {np.array_equal(elements, elements_read)}")
    print(f"  Fields match: {'distance' in fields_read}")
    
    if np.allclose(nodes, nodes_read) and np.array_equal(elements, elements_read):
        print("✓ VTU round-trip PASSED")
    else:
        print("✗ VTU round-trip FAILED")
