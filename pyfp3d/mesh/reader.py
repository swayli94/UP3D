"""
Mesh I/O: read .msh files via meshio and convert to Structure-of-Arrays format.

This module:
  - Reads .msh (Gmsh format) into meshio mesh object
  - Validates tetrahedral element type
  - Converts to SoA (Structure-of-Arrays) for efficient Numba kernels
  - Extracts and organizes boundary tags
  - Returns a Mesh object with all necessary topology information

Reference: design.md §4 (Mesh Requirements), §7 (Architecture)
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import meshio


class Mesh:
    """
    Structure-of-Arrays mesh representation optimized for Numba kernels.
    
    All arrays are C-contiguous (row-major) for Numba @njit compatibility.
    """
    
    def __init__(self):
        """Initialize empty mesh."""
        # Nodal data
        self.nodes: np.ndarray = np.empty((0, 3), dtype=np.float64)  # shape (n_nodes, 3)
        
        # Element connectivity (tetrahedral)
        self.elements: np.ndarray = np.empty((0, 4), dtype=np.int32)  # shape (n_tets, 4)
        
        # Boundary faces (triangular): indices into nodes
        self.boundary_faces: Dict[str, np.ndarray] = {}  # {"wall": (n_faces, 3), ...}
        
        # Element-to-boundary-tag mapping (for assembly)
        self.element_tags: np.ndarray = np.empty(0, dtype=np.int32)  # shape (n_tets,)
        self.tag_names: List[str] = []  # e.g. ["bulk", "farfield", ...]
        
        # Metadata
        self.name: str = "unnamed"
        self.unit: str = "m"
    
    def __repr__(self) -> str:
        return (
            f"Mesh(n_nodes={len(self.nodes)}, n_tets={len(self.elements)}, "
            f"boundary_tags={list(self.boundary_faces.keys())})"
        )
    
    def validate(self) -> bool:
        """
        Validate mesh topology consistency.
        
        Returns:
            True if all checks pass
            
        Raises:
            AssertionError if any invariant is violated
        """
        n_nodes = len(self.nodes)
        n_tets = len(self.elements)
        
        # Check node coordinates are real and finite
        assert np.all(np.isfinite(self.nodes)), "Non-finite node coordinates"
        
        # Check element connectivity
        assert np.all(self.elements >= 0), "Negative node index in element"
        assert np.all(self.elements < n_nodes), "Node index out of bounds"
        
        # Check boundary faces reference valid nodes
        for tag, faces in self.boundary_faces.items():
            assert np.all(faces >= 0), f"Negative node in {tag} boundary"
            assert np.all(faces < n_nodes), f"Out-of-bounds node in {tag} boundary"
        
        # Check element tags
        if len(self.element_tags) > 0:
            assert len(self.element_tags) == n_tets, "Element tag count mismatch"
            assert np.all(self.element_tags >= 0), "Negative element tag"
            assert np.all(self.element_tags < len(self.tag_names)), "Tag out of range"
        
        return True


def read_mesh(filepath: Path | str, verbose: bool = False) -> Mesh:
    """
    Read a .msh (Gmsh) file and convert to SoA mesh.
    
    Args:
        filepath: Path to .msh file
        verbose: Print info about mesh
        
    Returns:
        Mesh object with nodes, elements, boundary faces, and tags
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If mesh contains non-tetrahedral elements or is empty
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Mesh file not found: {filepath}")
    
    # Read with meshio
    mesh_obj = meshio.read(str(filepath))
    
    if verbose:
        print(f"Read mesh from {filepath}")
        print(f"  Points: {mesh_obj.points.shape}")
        print(f"  Cell types: {list(mesh_obj.cells_dict.keys())}")
    
    # Extract nodes
    nodes = np.asarray(mesh_obj.points, dtype=np.float64, order='C')
    
    # Extract tetrahedral elements
    if "tetra" not in mesh_obj.cells_dict:
        raise ValueError("Mesh must contain tetrahedral elements (cell type 'tetra')")
    
    tets = np.asarray(mesh_obj.cells_dict["tetra"], dtype=np.int32, order='C')
    
    if len(tets) == 0:
        raise ValueError("No tetrahedral elements found in mesh")
    
    # Extract boundary faces (triangles on domain boundary), split by physical
    # surface tag name (e.g. "wall", "farfield") when available so BCs can be
    # applied per named boundary; falls back to a single "all_triangles" group
    # for untagged meshes.
    boundary_faces_dict = {}
    if "triangle" in mesh_obj.cells_dict:
        triangles = np.asarray(mesh_obj.cells_dict["triangle"], dtype=np.int32, order='C')

        surface_tags = None
        if "gmsh:physical" in mesh_obj.cell_data_dict:
            physical_tags = mesh_obj.cell_data_dict["gmsh:physical"]
            if "triangle" in physical_tags:
                surface_tags = np.asarray(physical_tags["triangle"], dtype=np.int32)

        surface_name_by_tag = {}
        if hasattr(mesh_obj, "field_data"):
            surface_name_by_tag = {
                tag_id: name for name, (tag_id, dim) in mesh_obj.field_data.items() if dim == 2
            }

        if surface_tags is not None and surface_name_by_tag:
            for tag_id, name in surface_name_by_tag.items():
                mask = surface_tags == tag_id
                if np.any(mask):
                    boundary_faces_dict[name] = triangles[mask]
        else:
            boundary_faces_dict["all_triangles"] = triangles
    
    # Organize boundary tags (if cell_data exists)
    element_tags = np.zeros(len(tets), dtype=np.int32)
    tag_names = ["bulk"]  # Default tag for all elements
    
    if "gmsh:physical" in mesh_obj.cell_data_dict:
        physical_tags = mesh_obj.cell_data_dict["gmsh:physical"]
        if "tetra" in physical_tags:
            element_tags = np.asarray(physical_tags["tetra"], dtype=np.int32)
    
    # Build tag name map if available
    if hasattr(mesh_obj, 'field_data'):
        tag_map = mesh_obj.field_data  # Dict mapping tag name -> (tag_id, dimension)
        tag_names = [""] * (max(element_tags) + 1)
        for name, (tag_id, dim) in tag_map.items():
            if dim == 3:  # Volume region
                if tag_id < len(tag_names):
                    tag_names[tag_id] = name
    
    # Create output mesh
    mesh = Mesh()
    mesh.nodes = nodes
    mesh.elements = tets
    mesh.boundary_faces = boundary_faces_dict
    mesh.element_tags = element_tags
    mesh.tag_names = tag_names  # keep aligned with element_tags (indices are tag ids)
    mesh.name = filepath.stem
    
    # Validate
    mesh.validate()
    
    if verbose:
        print(f"  → {mesh}")
    
    return mesh


def write_mesh(mesh: Mesh, filepath: Path | str, verbose: bool = False) -> None:
    """
    Write mesh to .msh (Gmsh) format via meshio.
    
    Args:
        mesh: Mesh object
        filepath: Output path
        verbose: Print info
    """
    filepath = Path(filepath)
    
    # Build meshio mesh object
    cells = [("tetra", mesh.elements)]
    if "all_triangles" in mesh.boundary_faces:
        cells.append(("triangle", mesh.boundary_faces["all_triangles"]))
    
    mesh_obj = meshio.Mesh(mesh.nodes, cells)
    
    meshio.write(str(filepath), mesh_obj)
    
    if verbose:
        print(f"Wrote mesh to {filepath}")


def mesh_stats(mesh: Mesh) -> Dict:
    """
    Compute basic mesh statistics.
    
    Args:
        mesh: Mesh object
        
    Returns:
        Dict with: n_nodes, n_elements, volume, aspect_ratios, etc.
    """
    from pyfp3d.mesh.metrics import compute_tet_volumes, compute_edge_lengths
    
    n_nodes = len(mesh.nodes)
    n_tets = len(mesh.elements)
    
    # Compute volumes (detailed in metrics.py)
    volumes = compute_tet_volumes(mesh.nodes, mesh.elements)
    total_volume = np.sum(volumes)
    
    # Edge statistics
    edge_lengths = compute_edge_lengths(mesh.nodes, mesh.elements)
    
    stats = {
        "n_nodes": n_nodes,
        "n_elements": n_tets,
        "total_volume": total_volume,
        "min_volume": np.min(volumes),
        "max_volume": np.max(volumes),
        "mean_volume": np.mean(volumes),
        "min_edge_length": np.min(edge_lengths),
        "max_edge_length": np.max(edge_lengths),
        "mean_edge_length": np.mean(edge_lengths),
    }
    
    return stats


if __name__ == "__main__":
    # Quick test: read and validate a mesh if available
    import sys
    
    if len(sys.argv) > 1:
        mesh_file = Path(sys.argv[1])
        mesh = read_mesh(mesh_file, verbose=True)
        print(f"✓ Mesh loaded: {mesh}")
    else:
        print("Usage: python -m pyfp3d.mesh.reader <mesh.msh>")
