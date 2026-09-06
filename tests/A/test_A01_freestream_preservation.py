"""
Primary regression test for freestream preservation.

This test MUST pass before advancing to any new phase.
It validates that ∇²φ = 0 with φ = x (linear freestream) yields R ≈ 0.

Run with: pytest tests/A/test_A01_freestream_preservation.py -xvs
"""

import pytest
import numpy as np

from tests.mesh_utils import generate_structured_cube_mesh, cube_boundary_mask


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


#: ★★★ W1.1 / H6 (2026-09-06). 直到这一轮，本文件**只测 `phi = x` 一个方向**
#: （上面那条，:116）。业界惯例（SU2 / CFL3D 的回归套件）是**任意方向均匀流**，
#: 因为它抓的是「度量 / 装配只在坐标轴对齐时才正确」那一族 bug —— 非结构网格
#: 项目的常见病，而单方向检验**在原理上看不见它**。
#:
#: 方向取三个轴 + 三个**非轴对齐**方向；后者才是新增的判别力。
#:
#: ★★ **G-TEETH 实测（2026-09-06）**：给 `precompute_element_geometry` 的 B 的
#: **y 列**加一个逐单元的 1e-3 扰动（`B[e,:,1] *= 1 + 1e-3*((e%7)-3)`）——
#: **`+x` 与 `+z` 保持绿**（a_y = 0，被扰动的那一列根本不进来），
#: **`+y` 与三条斜方向全红**。⇒ **历史上那条单方向测试（φ = x）看不见这个缺陷**，
#: 而这正是本组存在的理由。
#: ★ 那次扰动只红了 `test_compressible_assembly_...` 那一组：
#: `PicardOperator` 走 `precompute_element_geometry`（jacobian.py:250），而
#: `kernels/residual.assemble_residual` 自己调 `element_gradients`
#: （residual.py:173）—— **两条装配路径的几何来源不同**，所以下面两组不是重复，
#: 各自覆盖一条路。
_DIRECTIONS = (
    ("+x", (1.0, 0.0, 0.0)),          # 历史基线，保持
    ("+y", (0.0, 1.0, 0.0)),
    ("+z", (0.0, 0.0, 1.0)),
    ("oblique_123", (1.0, 2.0, 3.0)),
    ("oblique_1m10", (1.0, -1.0, 0.0)),
    ("oblique_irrational", (0.31, -0.72, 0.55)),
)


def _unit(v):
    a = np.asarray(v, dtype=np.float64)
    return a / np.linalg.norm(a)


@pytest.mark.parametrize("name,d", _DIRECTIONS, ids=[n for n, _ in _DIRECTIONS])
def test_freestream_preservation_is_direction_independent(name, d):
    """φ = a·x for ANY unit a must give a machine-zero interior residual.

    ★ 判据与 x 方向**同一个容差**（1e-12）—— 一个只在轴对齐时成立的装配
    在这里必须红。
    ★ G-DOMAIN：相反结果 = 斜方向也机器零 ⇒ 装配与坐标轴无关，静默通过。
    """
    from pyfp3d.kernels.residual import assemble_residual

    nodes, elements = generate_structured_cube_mesh(n=4, L=1.0)
    phi = nodes @ _unit(d)
    R = assemble_residual(nodes, elements, phi)
    interior = ~cube_boundary_mask(nodes, L=1.0)
    worst = float(np.max(np.abs(R[interior])))
    assert worst < 1e-12, (
        f"direction {name} = {_unit(d)}: interior residual {worst:.3e} "
        f"-- freestream preservation is direction-DEPENDENT")


@pytest.mark.parametrize("name,d", _DIRECTIONS, ids=[n for n, _ in _DIRECTIONS])
def test_compressible_assembly_is_direction_independent(name, d):
    """同样的检验，但走**生产的可压装配** `PicardOperator.assemble_residual`。

    ★★ 上一条测的是 `kernels/residual.py` 的 Laplace 装配；生产的跨声速路径
    走的是 `PicardOperator` + 逐元素 `rho_tilde`。均匀流 ⇒ ρ 逐元素为常数，
    所以残差仍必须机器零 —— 这一条把方向无关性从 kernel 扩到**上线的那条装配**。
    """
    from pyfp3d.kernels.jacobian import PicardOperator

    nodes, elements = generate_structured_cube_mesh(n=4, L=1.0)
    op = PicardOperator(nodes, elements)
    phi = nodes @ _unit(d)
    _, q2 = op.velocities(phi)
    assert np.allclose(q2, 1.0, rtol=1e-12), (
        f"{name}: |grad phi|^2 = {q2.min():.6f}..{q2.max():.6f}, expected 1")
    rho_t = np.full(len(elements), 0.7, dtype=np.float64)   # 任意常数 rho
    R = op.assemble_residual(phi, rho_t)
    interior = ~cube_boundary_mask(nodes, L=1.0)
    worst = float(np.max(np.abs(R[interior])))
    assert worst < 1e-12, (
        f"direction {name}: production assembly interior residual "
        f"{worst:.3e} -- direction-DEPENDENT")


@pytest.mark.parametrize("name,d", _DIRECTIONS, ids=[n for n, _ in _DIRECTIONS])
def test_direction_independence_on_an_unstructured_mesh(name, d):
    """第三条腿：**真三维非结构**网格（已提交的 `sphere_shell/coarse.msh`）。

    ★ 结构化立方体的每个 tet 都是同一族的仿射像；非结构网格才让**每个单元的
    形函数梯度都不同**。方向无关性若只在结构化网格上成立，那是运气。

    ★★★ **不要换成 2.5-D 族**（`naca0012_2.5d` / `cylinder_2.5d`）。实测：
    它们是**单层拉伸**网格，两个 z 面都进 `symmetry` 组 ⇒ **全部节点都在边界上，
    内点数恒为 0**（naca0012_2.5d coarse: 5610/5610；cylinder_2.5d coarse:
    2412/2412）。⇒ **2.5-D 网格在原理上承载不了「内点残差」这类判据。**
    我第一版就是拿 naca0012_2.5d 写的，是下面那条**前提断言**把它拦下来的
    （前提断言先于物理断言，正是为此）。`sphere_shell/coarse.msh` 有 2297 个内点。
    """
    from pathlib import Path

    from pyfp3d.mesh.reader import read_mesh
    from tests.conftest import REPO_ROOT

    p = Path(REPO_ROOT) / "cases" / "meshes" / "sphere_shell" / "coarse.msh"
    if not p.exists():                       # 与 W0.1 同一条约定
        pytest.skip(f"{p.name} not generated")
    mesh = read_mesh(p)
    from pyfp3d.kernels.residual import assemble_residual
    phi = mesh.nodes @ _unit(d)
    R = assemble_residual(mesh.nodes, mesh.elements, phi)
    bnd = np.unique(np.concatenate(
        [f.ravel() for f in mesh.boundary_faces.values()]))
    interior = np.ones(len(mesh.nodes), dtype=bool)
    interior[bnd] = False
    assert interior.sum() > 100, (
        f"premise: the mesh must have interior nodes, got "
        f"{interior.sum()} of {len(mesh.nodes)} -- a single-layer 2.5-D "
        f"mesh has NONE (see the docstring)")
    worst = float(np.max(np.abs(R[interior])))
    #: ★ 容差比结构化那条松：非结构网格的单元体积跨 ~3 个数量级，累加的
    #: 舍入随之放大。**实测**最坏 ~1e-14（见轮次记录），带留 100x。
    assert worst < 1e-12, (
        f"direction {name}: unstructured interior residual {worst:.3e}")


if __name__ == "__main__":
    pytest.main([__file__, "-xvs"])
