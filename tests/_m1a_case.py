"""M1a 算例的共享装配 —— 被 C 类（网格收敛性质）与 E 类（逐级已提交值）两半共用。

★ 2026-08-24 重编号时抽出。沿用本仓库既有惯例（`mesh_utils.py`、`v5_state.py`、`_tol.py`）：
共享代码放 `tests/_*.py`，两半 import 它，**不复制** —— 复制会分叉。
★ 已知：`_solve` 无缓存，两半各自求解 ⇒ 与拆分前**行为相同**（拆分前也跑两次）。
"""
import numpy as np
import pytest
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting


M_INF, ALPHA = 0.72, 1.25    #: 包线内的 M1a 工况（自 test_s1_m1a_envelope 抽出）

LOCK = {
    "xcoarse": dict(cl=0.222742, m_max=1.0972),
    "coarse": dict(cl=0.242984, m_max=1.1356),
    "medium": dict(cl=0.253237, m_max=1.1537),
}
LADDER = ("xcoarse", "coarse", "medium")
CL_RTOL = 2.0e-3          # 0.2 %: well inside run-to-run scatter, far below
RATIO_MAX = 0.7           # measured 0.5065 in envelope, -0.6003 out
D2_REL_MAX = 0.05         # measured 4.22 % in envelope, 16.33 % out
def _solve(mesh_dir, level):
    path = mesh_dir / "naca0012_2.5d" / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path.name} not generated "
                    "(cases/meshes/naca0012_2.5d/generate_naca0012.py)")
    mc, wc = cut_wake(read_mesh(path))
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
                             freeze_refresh_max=8, precond="direct",
                             direct_refactor_every=4, n_newton_max=400)
    dz = float(np.ptp(mc.nodes[:, 2]))
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF)
    return r, float(f["cl"])
@pytest.fixture(scope="module")
def mesh_dir():
    from tests.conftest import REPO_ROOT
    return REPO_ROOT / "cases" / "meshes"
