"""A56 -- 翼身生成器的生成-比对冒烟门（A 类；预注册
phases/p6/docs/dev_phase_six/20261011-1600-code-surface-prereg.md §2.3）。

★ 由来：`pyfp3d/meshgen/wingbody.py` 的 9 个定义（含 443 行的
`onera_m6_wingbody_mesh`）此前只被 `cases/meshes/onera_m6_wingbody*` 的两个
生成脚本调用、**零测试执行**。本条**不是覆盖**，是「生成器改动会响」的
最小保险：以 coarse 档再生成一次翼身 conforming 网格（实测 9.6 s，
铁律上限 2 min），与**已提交的 manifest** 比对。

★★ 生成是**位确定**的（实测 2026-10-11：重生成写盘的 sha256 与 committed
manifest 逐字一致），所以指纹断言是**相等**而不是带 —— manifest 自己的
note 就写着 "regeneration must reproduce this sha256"。.msh 本身被
gitignore（M1/M5 政策），**已提交的指纹就是 manifest**，门只读它。

★ G-TEETH：把 committed manifest 的 sha256 改坏一位 ⇒ 门红；恢复 ⇒ 绿。
"""
import hashlib
import json
import time

import pytest

from pyfp3d.mesh.reader import write_mesh
from pyfp3d.meshgen.fuselage import FuselageParams
from pyfp3d.meshgen.wingbody import onera_m6_wingbody_mesh
from tests.conftest import REPO_ROOT

#: 与 cases/meshes/onera_m6_wingbody_conforming/generate_onera_m6_wingbody_conforming.py
#: 的 `_level_params(LEVELS["coarse"])` 逐字一致（h_wall 一档定全家）。
COARSE_PARAMS = dict(h_wall=0.030, h_far=6.0, h_wake=0.09, h_edge=0.015,
                     h_tip=0.0075, h_junction=0.015, h_body=0.06,
                     h_body_tip=0.0075)
MANIFEST = (REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody_conforming"
            / "coarse.msh.manifest.json")


def _fingerprint_errors(mesh, fresh_sha, mani):
    """生成网格 vs manifest 的全部指纹差（空 = 一致）。抽成函数是为了
    能用**扰动副本**喂它（test 2）—— 比对逻辑内联在断言里的门，被扰动的
    参考传不进去，等于没有读者。"""
    errs = []
    if len(mesh.nodes) != mani["n_nodes"]:
        errs.append(f"n_nodes {len(mesh.nodes)} != {mani['n_nodes']}")
    if len(mesh.elements) != mani["n_tets"]:
        errs.append(f"n_tets {len(mesh.elements)} != {mani['n_tets']}")
    got_tris = {k: len(v) for k, v in mesh.boundary_faces.items()}
    if got_tris != mani["n_tris"]:
        errs.append(f"boundary tris {got_tris} != {mani['n_tris']}")
    if fresh_sha != mani["sha256"]:
        errs.append(f"sha256 {fresh_sha} != {mani['sha256']}")
    return errs


@pytest.fixture(scope="module")
def fresh(tmp_path_factory):
    t0 = time.perf_counter()
    mesh = onera_m6_wingbody_mesh(fuselage=FuselageParams(), tip_cap="round",
                                  embed_wake=True, **COARSE_PARAMS)
    dt = time.perf_counter() - t0
    out = tmp_path_factory.mktemp("wb") / "fresh_coarse.msh"
    write_mesh(mesh, out)
    return mesh, hashlib.sha256(out.read_bytes()).hexdigest(), dt


def test_coarse_regeneration_matches_the_committed_fingerprint(fresh):
    """(a) 节点/单元数与边界组拓扑指纹 == committed manifest（相等，不是带）；
    (b) 重生成写盘的 .msh sha256 == manifest 的 sha256。"""
    mesh, sha, dt = fresh
    assert dt < 120, f"coarse 生成 {dt:.0f} s 超出本门的成本铁律（120 s）"
    errs = _fingerprint_errors(mesh, sha, json.loads(MANIFEST.read_text()))
    assert not errs, (
        "the wingbody generator no longer reproduces the committed coarse "
        f"mesh: {errs}（counts 看不见排序/符号翻转,所以指纹含写盘 sha256）")


def test_the_comparison_actually_reads_the_manifest(fresh):
    """把**扰动过**的 manifest 副本喂给比对函数，差必须被报出来 ——
    否则「比对」是按构造恒真的（F06 的坑：参考从没被读过）。"""
    mesh, sha, _ = fresh
    mani = json.loads(MANIFEST.read_text())
    assert not _fingerprint_errors(mesh, sha, mani), "baseline 必须先干净"
    for key, bad in (("n_nodes", mani["n_nodes"] + 1),
                     ("n_tris", dict(mani["n_tris"], wall=mani["n_tris"]["wall"] + 1)),
                     ("sha256", ("1" if mani["sha256"][0] != "1" else "0")
                                + mani["sha256"][1:])):
        pert = dict(mani)
        pert[key] = bad
        errs = _fingerprint_errors(mesh, sha, pert)
        assert errs, f"perturbed manifest[{key}] produced NO error -- the gate is not reading it"
