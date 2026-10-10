"""B05 — 壁面 Cp 锯齿：**载荷不平滑、截面平滑一遍**（B 类）。

★★★ **这道门的定义被推翻过两次，两次都是被测量推翻的。**记在这里，因为判据的形状
就是那两次推翻的产物：

1. 我最初提议「**补一条锯齿度量**」—— **实测为假**：`cp_oscillation_metric` 早就在
   `post/section_cut.py`，`test_B03` 用 7 条锁着，而且它**比我打算写的更细**
   （只在**斜率反号**的内点取二阶差分，专门把 2 格锯齿与激波真实陡升分开；
   我打算写的「二阶差分范数」正是它 docstring 点名否掉的做法）。
   我会说错，是因为我只看了 P6 那两个测试的**名字**没读实现 —— 「**名字 ≠ 内容**」。
2. 改提议「锁生产的恢复压住了锯齿」—— **也为假**：活代码 12 处 `section_cp_curve`
   调用**没有一处开平滑**，而 G6.3（2026-07-08 实测）明写**故意不开**。
   ★ 若照建，那会是一道**与生产配置相反**的门：要么永远红，要么有人为了让它绿
   而把生产改成一个**已被测量否定**的配置 —— 后者比没有这道门糟得多。

⇒ **使用者 2026-08-25 的裁决把界划对了，而且是按「被积分的是什么」划的**：

| 量 | 平滑 | 理由 |
|---|---|---|
| **全表面总升力** | **不平滑** | 对整个曲面积分，±锯齿**天然抵消**；平滑只引入额外误差 |
| **截面 Cp 提取** | **平滑** | 截面切割是 1-D 切线，**没有配对**；展向不均匀 ⇒ **系统性偏采一号** |
| **截面升力系数** | **平滑** | ★ 它**由截面 Cp 积分而来** ⇒ 与截面 Cp 同侧（我原先漏了这一条） |

## 判据的每个数字都来自实测（M6 coarse，M0.8395 / α3.06，2026-08-25）

| 量 | sp=0 | sp=1 | sp=2 |
|---|---|---|---|
| 锯齿度量中位（7 站上表面） | **0.1742** | **0.0469** | 0.0607 |
| 7 站 pooled RMS（M3a 自己的 `station_rms`） | 0.22029 | **0.21125**（−4.10 %） | 0.22266（+1.08 %） |
| **LE 带** | 0.32452 | **0.30688**（−5.44 %） | 0.33220（**+8.3 %**） |
| MID 带 | 0.11182 | 0.11082（−0.90 %） | 0.11061 |
| TE 带 | 0.26224 | 0.25879（−1.32 %） | 0.25498 |

★★ **一遍有净收益、两遍削过头，而且同一个机理解释了两者**：收益与代价都落在 **LE 带**，
正是 G6.3 说的「平滑**抹掉前缘吸力峰**」。

## ★★★ RECORDED（不设门）：锯齿**不是** M3a 地板的成分

M3a 的地板在 **MID(79 %)/TE**，而平滑改善的是 **LE（−5.44 %）**，MID 只动 −0.90 %。
⇒ 「那个无粘模型地板 0.0516–0.0707 里掺着后处理产物」这个**推测不成立**。
★ 这条是**否定结论**，写在这里是为了不让它被重新提出来。

★ 另：本轮的 pooled 0.220 是 **coarse**，M3a 记录的 0.1007 是 **medium 全点** ——
**不是同一个数**，不得相互引用。
★★ 而更早我自己写过一版 RMS（上下两面都插到**全部**实验点上）得到 0.40，
我当时把它记成「另一种口径」放过了 —— **差 4 倍的量级不符应该当成「其中一个是错的」
去定位，而不是当成「两种约定」去解释**。那一版把真实改善从 4.10 % 稀释成 0.46 %。
"""
import numpy as np
import pytest

from bench.recipes import NEWTON_M6_RECIPE
from bench.run_m3_budget import BANDS, band_rms, parse_experiment, station_rms
from pyfp3d.constraints.wake import tip_taper_factors
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import B_SEMI
from pyfp3d.physics.isentropic import GAMMA
from pyfp3d.post.section_cut import cp_oscillation_metric, section_cp_curve
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_transonic
from tests.conftest import REPO_ROOT

M_INF, ALPHA = 0.8395, 3.06
CP_STAR = 2 / (GAMMA * M_INF ** 2) * (
    ((1 + 0.5 * (GAMMA - 1) * M_INF ** 2) / (1 + 0.5 * (GAMMA - 1)))
    ** (GAMMA / (GAMMA - 1)) - 1)

#: 实测 0.1742 —— 带取 0.05，余量 3.5x。★ 这一条是 G-DOMAIN 守卫：
#: 若未平滑时锯齿本来就 ~0，下面两条会**空过**（vacuous pass），整道门失去意义。
SAW_PRESENT_MIN = 0.05
#: 实测 0.1742 -> 0.0469 = 3.7x。带取 2.0x，余量 1.85x。
SAW_SUPPRESS_MIN = 2.0
#: 实测 -4.10 %。带取 -1.0 %，余量 4.1x。
RMS_IMPROVE_MIN = 0.01


@pytest.fixture(scope="module")
def m6_curves():
    """M6 coarse 的生产状态 + 三档平滑的 7 站截面。

    ★ 求解构造**逐字复制** `tests/E/test_E01::_m6_case`（taper + NEWTON_M6_RECIPE +
    transonic）—— 那是生产路径本身，而 `bench/run_m3_budget.py` 有一条守卫盯着它不许变。
    实测 ~8 s / conv=True / 0 钳制。
    """
    p = REPO_ROOT / "cases" / "meshes" / "onera_m6" / "coarse.msh"
    if not p.exists():
        pytest.skip("onera_m6/coarse.msh 未生成（cases/meshes/onera_m6/generate_onera_m6.py）")
    mc, wc = cut_wake(read_mesh(str(p)))
    kw = dict(NEWTON_M6_RECIPE)
    kw["newton_kw"] = dict(kw["newton_kw"], tip_taper=tip_taper_factors(
        wc.station_z, B_SEMI, "vanish_smooth", 0.05 * B_SEMI))
    r = solve_newton_transonic(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)
    assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0, (
        "M6 coarse 生产配方没有干净收敛（conv=%s lim=%s flr=%s）—— "
        "下面每一条读数都建立在这个状态上，所以这是前置条件不是附带检查"
        % (r["converged"], r["n_limited"], r["n_floored"]))
    exp = parse_experiment()
    etas = sorted(exp)
    curves = {sp: {e: section_cp_curve(mc, r["phi"], eta=e, b_semi=B_SEMI,
                                       m_inf=M_INF, smooth_passes=sp)
                   for e in etas} for sp in (0, 1, 2)}
    return exp, etas, curves, mc, r


def _pooled(curves, exp, etas):
    ss, n = 0.0, 0
    for e in etas:
        v, k = station_rms(curves, exp, e)
        ss += v * v * k
        n += k
    return (ss / max(n, 1)) ** 0.5


def _median_saw(curves, etas):
    m = [cp_oscillation_metric(np.asarray(curves[e]["x_upper"]),
                               np.asarray(curves[e]["cp_upper"]),
                               CP_STAR).get("metric", np.nan) for e in etas]
    return float(np.nanmedian(m))


def test_the_sawtooth_is_actually_there_unsmoothed(m6_curves):
    """★ G-DOMAIN 守卫：**先证明有锯齿可压**，否则下面两条是空过。

    「相反的结果会落在哪里」：若未平滑的度量本来就 < 0.05，那说明这张网格上没有
    可测的锯齿 —— 那时后两条即使全绿也**什么都没证明**，而这道门会给出一个虚假的安心。
    """
    exp, etas, curves, _, _ = m6_curves
    med = _median_saw(curves[0], etas)
    assert med > SAW_PRESENT_MIN, (
        "未平滑的锯齿度量中位只有 %.5f（要求 > %.2f）—— 没有可压的锯齿，"
        "本文件后两条会空过。检查是不是换了网格/工况，而不是放宽这个数"
        % (med, SAW_PRESENT_MIN))


def test_one_smoothing_pass_suppresses_the_sawtooth(m6_curves):
    """一遍平滑把截面锯齿压下去（实测 0.1742 -> 0.0469 = 3.7x）。"""
    exp, etas, curves, _, _ = m6_curves
    m0, m1 = _median_saw(curves[0], etas), _median_saw(curves[1], etas)
    assert m1 > 0.0, "平滑后的度量为 0 或 nan —— 度量本身失效，不是锯齿没了"
    assert m0 / m1 >= SAW_SUPPRESS_MIN, (
        "一遍平滑只把锯齿中位从 %.5f 压到 %.5f（%.2fx，要求 >= %.1fx）"
        % (m0, m1, m0 / m1, SAW_SUPPRESS_MIN))


def test_one_pass_improves_the_seven_station_comparison(m6_curves):
    """★★ 平滑**改善了对实验的对比**，不只是「变好看」（实测 -4.10 %）。

    用 M3a 自己的 `station_rms`（逐点按实验的 upper 侧标志映射），
    ★ **不是**我自己写的版本 —— 那一版把上下两面都插到全部实验点上，
    把真实的 4.10 % 稀释成 0.46 %，并让我得出「平滑没有可测收益」的错误结论。
    """
    exp, etas, curves, _, _ = m6_curves
    p0, p1 = _pooled(curves[0], exp, etas), _pooled(curves[1], exp, etas)
    gain = (p0 - p1) / p0
    assert gain >= RMS_IMPROVE_MIN, (
        "一遍平滑对 7 站 pooled RMS 的改善只有 %+.2f %%（要求 >= %.1f %%）：%.5f -> %.5f"
        % (100 * gain, 100 * RMS_IMPROVE_MIN, p0, p1))


def test_two_passes_are_worse_than_one(m6_curves):
    """★★ **不是越平滑越好** —— 两遍把前缘吸力峰削过头。

    实测 pooled 0.21125（一遍）-> 0.22266（两遍），而 **LE 带 +8.3 %** ——
    收益与代价落在**同一个带**上，正是 G6.3 说的「averaging smears the LE suction peak」。
    ⇒ 这一条防的是「既然一遍有用那就多来几遍」。
    """
    exp, etas, curves, _, _ = m6_curves
    p1, p2 = _pooled(curves[1], exp, etas), _pooled(curves[2], exp, etas)
    assert p2 > p1, (
        "两遍平滑的 pooled RMS（%.5f）不比一遍（%.5f）差 —— 若这是真的，"
        "本文件关于「一遍是那个点」的判断需要重测，而不是把这条删掉" % (p2, p1))


def test_the_bands_reconstruct_the_pooled_rms(m6_curves):
    """★ 分解的可信前提：LE/MID/TE 三带**按构造**必须加得回 pooled。

    `band_rms` 的 docstring 明写这是「分解不至于在量别的东西」的唯一保证 ——
    所以它该是一条断言，不是一个假定。
    """
    exp, etas, curves, _, _ = m6_curves
    for sp in (0, 1, 2):
        ss = n = 0
        for e in etas:
            br = band_rms(curves[sp], exp, e)
            for bname, _, _ in BANDS:
                for side in ("upper", "lower"):
                    v, k = br["%s_%s" % (bname, side)]
                    ss += v
                    n += k
        recon = (ss / max(n, 1)) ** 0.5
        assert abs(recon - _pooled(curves[sp], exp, etas)) < 1e-9, (
            "sp=%d 的分带加不回 pooled（%.9f vs %.9f）—— 分解在量别的东西"
            % (sp, recon, _pooled(curves[sp], exp, etas)))


def test_the_loads_path_stays_unsmoothed(m6_curves):
    """★★★ 载荷路径**必须不平滑**（G6.3，2026-07-08 实测）。

    对整个曲面积分时 ±锯齿**天然抵消**；平滑反而**抹掉前缘吸力峰**，
    在 M6 coarse 上把 `CL_p` 相对可信的 `CL_KJ` 从 2.40 % 推到 3.35 %。
    ⇒ 这一条锁的是**生产为什么这么配**，而在此之前它只活在 `surface.py` 的一句注释里。
    """
    import inspect
    assert inspect.signature(wall_force_coefficients).parameters[
        "smooth_passes"].default == 0, (
        "wall_force_coefficients 的 smooth_passes 默认不再是 0 —— "
        "载荷路径开平滑会抹掉 LE 吸力峰（G6.3：CL_p 2.40 % -> 3.35 %）")
    assert inspect.signature(section_cp_curve).parameters[
        "smooth_passes"].default == 0, (
        "section_cp_curve 的默认变了。★ 默认 0 是对的（它同时服务于诊断用途），"
        "但**截面对比与截面 cl 的调用点必须显式传 smooth_passes=1**（使用者裁决 2026-08-25）")
