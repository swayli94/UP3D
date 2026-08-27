"""
Pytest configuration and shared fixtures for pyFP3D test suite.

Gate-level tests should use these fixtures to ensure consistent reference data,
mesh sets, and artifact storage.
"""

import pytest
import os
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent


@pytest.fixture
def gate_evidence_dir(request):
    """C / D 类门的**被跟踪**证据目录：`cases/gates/<门号>/`。

    ★★★ 2026-08-24：取代已删除的 `artifacts_dir` / `gate_artifacts_dir`。
    那两个 fixture 指向 gitignored 的 `<repo>/artifacts/`，于是
    **workflow 规则要求「每个视觉门留下产物」，而 .gitignore 保证它永不进 HEAD**
    —— 按纪律 3 它从来就不是证据，却有 11 处文档把它当证据引用。
    ★ 而 `mkdir(exist_ok=True)` 让它**每次 pytest 一跑就重新出现**，所以删目录不是修法，
    必须删 fixture 本身。

    **只有 C / D 类用它** —— 它们要给人肉眼检查（压力分布 / 力系数 / 网格收敛）。
    A / B / E / F 不写图：它们的判据**就是断言本身**，写出来没人看的 PNG 是纯 churn
    （实测那 5 个 A/B 的 savefig 里 4 个没有任何文档引用）。需要临时看图用 pytest 的 `tmp_path`。

    **写与不写**（避免 churn，同时保证图与断言同源）：
      - 平时跑：**不写**，断言对着已提交的 `summary.csv` ⇒ 代码一改答案就红，逼出一次有意的刷新；
      - `PYFP3D_GATE_FIGURES=1`：同一次计算重写 `summary.csv` + `<门号>/*.png`。
    ⇒ 图与断言来自**同一次计算，构造保证**；刷新走本项目已有的再基线勘误纪律。
    """
    #: 门号 = 测试文件名的 <CLASS><nn>_<stem> 段，于是「门 <-> 证据目录」可机械互查
    stem = Path(str(request.node.fspath)).stem          # test_C03_laplace_sphere
    gate = stem[len("test_"):] if stem.startswith("test_") else stem
    d = REPO_ROOT / "cases" / "gates" / gate
    if os.environ.get("PYFP3D_GATE_FIGURES"):
        #: ★ 建门目录本身，不建 `figures/` 子目录 —— 2026-08-26 实测：fixture 建了它，
        #: 而**八个门一个都没往里写**（PNG 全写在门目录根下），git 又不跟踪空目录，
        #: 于是它只是本地噪声。**规则与实现不一致，改规则跟实现走。**
        d.mkdir(parents=True, exist_ok=True)
    return d


def gate_figures_enabled():
    """True 时才写图/CSV。★ 平时为 False，所以一次普通 pytest 不会脏工作树。"""
    return bool(os.environ.get("PYFP3D_GATE_FIGURES"))


@pytest.fixture
def reference_mesh_dir():
    """Return path to cases/reference_data (immutable reference meshes and data)."""
    repo_root = Path(__file__).parent.parent
    ref_dir = repo_root / "cases" / "reference_data"
    ref_dir.mkdir(parents=True, exist_ok=True)
    return ref_dir


@pytest.fixture
def mesh_dir():
    """Return path to cases/meshes (where coarse/medium/fine families are stored)."""
    repo_root = Path(__file__).parent.parent
    mesh_dir = repo_root / "cases" / "meshes"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    return mesh_dir


@pytest.fixture
def disable_numba_jit(monkeypatch):
    """
    Fixture to disable Numba JIT during tests for debugging.
    
    Usage:
        def test_something(disable_numba_jit):
            # Numba @njit will run in object mode, allowing prints/pdb
            ...
    """
    monkeypatch.setenv("PYFP3D_NOJIT", "1")
