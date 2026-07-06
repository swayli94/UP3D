"""
Primary regression test for freestream preservation.

This test MUST pass before advancing to any new phase.
It validates that ∇²φ = 0 with φ = x (linear freestream) yields R ≈ 0.

Run with: pytest tests/test_v0_freestream.py -xvs
"""

import pytest
import numpy as np

from .mesh_utils import generate_structured_cube_mesh, cube_boundary_mask


def test_import_pyfp3d():
    """Smoke test: can we import the main package?"""
    import pyfp3d
    assert pyfp3d.NOJIT is False, "PYFP3D_NOJIT not set (OK for CI)"


def test_import_physics():
    """Smoke test: can we import physics module?"""
    from pyfp3d.physics import isentropic
    assert hasattr(isentropic, "density_isentropic")
    assert hasattr(isentropic, "pressure_coefficient")


def test_isentropic_stagnation():
    """
    Test isentropic density at stagnation (q² = 0).
    
    At stagnation, total enthalpy is conserved, so density rises:
    ρ₀ = ρ∞ * [1 + (γ-1)/2 · M∞²]^{1/(γ-1)}
    
    For M∞ = 0.5, this gives ρ₀ ≈ 1.1297.
    """
    from pyfp3d.physics.isentropic import density_isentropic
    import math
    
    M_inf = 0.5
    gamma = 1.4
    q_squared = 0.0  # Stagnation
    
    rho = density_isentropic(q_squared, M_inf)
    
    # Expected: [1 + 0.2 * 0.25]^2.5 = 1.05^2.5
    rho_expected = (1.0 + 0.5 * (gamma - 1.0) * M_inf**2) ** (1.0 / (gamma - 1.0))
    
    assert abs(rho - rho_expected) < 1e-14, f"Stagnation density off: computed {rho}, expected {rho_expected}"


def test_isentropic_freestream():
    """
    Test isentropic density at freestream (q² = 1).
    
    At q² = 1 (freestream speed), density should be reference = 1.0.
    """
    from pyfp3d.physics.isentropic import density_isentropic, mach_number_squared
    
    M_inf = 0.5
    q_squared = 1.0  # Freestream
    
    rho = density_isentropic(q_squared, M_inf)
    M_sq_local = mach_number_squared(q_squared, M_inf)
    
    # At freestream (q = U∞), local Mach should equal M∞
    assert abs(rho - 1.0) < 1e-14
    assert abs(M_sq_local - M_inf**2) < 1e-14, f"Freestream Mach off: {M_sq_local}"


def test_pressure_coefficient_bounds():
    """
    Test that computed Cp stays within physical bounds for a range of Mach and speeds.
    
    This is a sanity check before gate closures (G1–G3).
    """
    from pyfp3d.physics.isentropic import (
        pressure_coefficient, 
        validate_physics_bounds,
        critical_speed_squared
    )
    import numpy as np
    
    M_inf = 0.7
    q_sq_crit = critical_speed_squared(M_inf)
    
    # Sample speeds from subsonic to near-sonic
    q_sq_array = np.linspace(0, q_sq_crit * 1.1, 20)
    
    rho_array = np.ones_like(q_sq_array)
    q_array = np.sqrt(q_sq_array)
    M_array = np.sqrt(q_sq_array) * M_inf  # Approximate
    Cp_array = np.array([pressure_coefficient(q_sq, M_inf) for q_sq in q_sq_array])
    
    # Should not raise an error
    validate_physics_bounds(rho_array, q_array, M_array, Cp_array, M_inf)


def test_residual_freestream_preservation():
    """
    Mesh-level freestream preservation (design.md Sec 3): phi = x on any
    tet mesh must give a machine-zero assembled residual. This is the
    check agent-rules.md hard rule #1 refers to -- run this test first
    after any kernel/assembly change.

    Only interior nodes are checked: a boundary node's residual is the
    divergence-theorem flux integral of grad(phi) through its own boundary
    support, which is nonzero whenever the boundary isn't a solid wall (zero
    flux) or isn't force-free by symmetry -- that row gets overwritten by
    the BC anyway, so it isn't part of the "residual should vanish" claim.
    """
    from pyfp3d.kernels.residual import assemble_residual

    nodes, elements = generate_structured_cube_mesh(n=4, L=1.0)
    phi = nodes[:, 0].copy()  # uniform freestream aligned with x

    R = assemble_residual(nodes, elements, phi)
    interior = ~cube_boundary_mask(nodes, L=1.0)

    assert np.max(np.abs(R[interior])) < 1e-12, (
        f"Freestream residual not machine-zero: {np.max(np.abs(R[interior])):.3e}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
