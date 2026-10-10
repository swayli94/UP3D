"""A53 — 库默认值门（A 类）。**零求解**：只读签名。

★★★ **为什么需要它。** 2026-08-24 清点时实测：这些默认值决定每一个报出去的数字，
而**没有任何东西检查它们有没有被静默改掉**。它们分散在 `solve/newton.py` 的签名、
`meshgen/wing3d.py` 的关键字、以及一个模块常量里，**只活在文档的散文里**。

★★ **这道门不是「防止改动」，是「让改动必须是有意的」。** 每一条都带着它的
**来源与代价**，所以红的时候读的人立刻知道自己动了什么、以及要不要走勘误清单：

| 值 | 它是什么，代价在哪 |
|---|---|
| `ew_eta0` / `ew_eta_max` = 1e-6 | 决策记录 `phases/p2/docs/dev_phase_two/DECISION-2026-08-02-precond.md`（含被排除的选项与推翻条件）。★ 这是它当初**没能被归档**的唯一可测理由 |
| `entropy_correction` = True | ★★ **模型修正，不是稳定化技巧**。实测它让峰更**强**（M²max 6.426 ON vs 5.607 OFF）；理由是 S1b 测到 M0.80 的激波被它**移进** Euler 锚点带 |
| `tip_cap` = "round" | ★ flat 必须**按名字显式索取** —— P13/G13.3 实测 flat 的凸边在加密下**发散**（p = +0.321）⇒ 任何基于加密的结论在 flat 网格上**前提为假** |
| `kutta_estimator` = "probe" | ★ 生产的**翼身**用 `"pressure"`；R19 的坑正是改了 pressure 那一行而 NACA 配方走 probe |
| `tip_taper` = None | ★ 生产的保形翼身**显式**传 `("vanish_smooth", 0.05·b_semi)`，代价 **−1.3 % cl_p**，是**模型偏置** —— 任何翼身升力读数都带着它 |
| `m_cap` = 3.0 | 限的是**局部 Mach**，不是 M²；判据是 `q2 >= q2_at_mach(m_cap, ...)` |
| `rho_floor` = 0.05 | 限的是**人工**密度 ρ̃，库**不返回**它 |
| `_SEED_FALLBACK` = 5 | 冷启动回退，只在「无种子 ∧ 无热启动 ∧ 结果被钳制」**三条同时成立**时触发一次 |

★ **改动流程**（红了之后该怎么做，写在这里而不是让人猜）：
先确认那是有意的，再改本文件里的期望值，**并按再基线纪律 grep 被移动的数字**
（`docs/agent-rules.md` 纪律 11）—— 一个默认值动了，通常有一批已提交的数字跟着动。
"""
import inspect

from pyfp3d.meshgen.wing3d import onera_m6_wing_mesh
from pyfp3d.solve.newton import (_SEED_FALLBACK, solve_newton_lifting,
                                 solve_newton_transonic)

#: 实测于 2026-08-24（`inspect.signature`，不是从文档抄的）。
EXPECTED_LIFTING = {
    "upwind_c": 1.5,
    "m_crit": 0.95,
    "m_cap": 3.0,
    "rho_floor": 0.05,
    "n_picard_seed": 0,
    "ew_eta0": 1e-06,
    "ew_eta_max": 1e-06,
    "precond": "amg",
    "tip_taper": None,
    "gamma_target": None,
    "kutta_estimator": "probe",
    "entropy_correction": True,
}

#: 跨驱动一致性：连续路径与 Mach 爬坡路径**必须同意**这几个的取值，
#: ★ 否则「同一个求解器」在两条路上做的是两件不同的事。
EXPECTED_TRANSONIC = {
    "upwind_c": 1.5,
    "m_crit": 0.95,
}

WHY = {
    "ew_eta0": "EW forcing，决策记录在 phases/p2/.../DECISION-2026-08-02-precond.md",
    "ew_eta_max": "同上",
    "entropy_correction": "模型修正而非稳定化技巧；ON 让峰更强（6.426 vs 5.607）",
    "kutta_estimator": "生产翼身用 pressure；改错行是 R19 踩过的坑",
    "tip_taper": "生产翼身显式传 vanish_smooth 0.05·b_semi，代价 -1.3% cl_p（模型偏置）",
    "m_cap": "限的是局部 Mach 不是 M^2",
    "rho_floor": "限的是人工密度 rho_tilde，库不返回它",
    "precond": "M6 fine 的 'direct' 是 4h39m/26GB 的 splu 陷阱",
}


def _msg(fn, key, got, want):
    why = WHY.get(key, "")
    return (
        "%s 的默认 %s 变了：%r -> %r。\n"
        "  这个值是什么：%s\n"
        "  ★ 若这是有意的：改本文件的期望值，并按纪律 11 grep 被移动的数字 —— "
        "一个默认值动了，通常有一批已提交的数字跟着动。" % (fn, key, want, got, why or "见本文件 docstring 的表")
    )


def test_lifting_driver_defaults():
    """`solve_newton_lifting` 的每一个生产默认值。"""
    sig = inspect.signature(solve_newton_lifting)
    for key, want in EXPECTED_LIFTING.items():
        assert key in sig.parameters, (
            "参数 %s 从 solve_newton_lifting 消失了 —— 这不是默认值变化，是接口变化，"
            "调用方会静默拿到别的行为" % key)
        got = sig.parameters[key].default
        assert got == want, _msg("solve_newton_lifting", key, got, want)


def test_transonic_driver_agrees_with_the_lifting_driver():
    """★ 两条驱动路径对人工密度的取值必须一致。

    不一致意味着「同一个求解器」在连续路径与 Mach 爬坡路径上做的是**两件不同的事**，
    而那种差别在结果里看不出来 —— 只会表现成两条路给不同的数。
    """
    sig = inspect.signature(solve_newton_transonic)
    for key, want in EXPECTED_TRANSONIC.items():
        got = sig.parameters[key].default
        assert got == want, _msg("solve_newton_transonic", key, got, want)
        assert got == EXPECTED_LIFTING[key], (
            "transonic 与 lifting 两条驱动的 %s 默认值不一致（%r vs %r）—— "
            "同一个求解器在两条路上做的是两件事" % (key, got, EXPECTED_LIFTING[key]))


def test_seed_fallback_constant():
    """冷启动回退的种子数。★ 它不是可调参数，是一个被测量出来的恢复值。"""
    assert _SEED_FALLBACK == 5, (
        "_SEED_FALLBACK 从 5 变成 %r。它只在「无种子 ∧ 无热启动 ∧ 结果被钳制」三条"
        "同时成立时触发一次；改它会改变冷启动在超临界 M∞ 上的恢复行为" % (_SEED_FALLBACK,))


def test_wing_tip_cap_is_round_by_default():
    """★★ 圆尖是默认，flat 必须按名字显式索取。

    P13/G13.3 实测 flat 帽的凸边在加密下**发散**（峰值 Mach 指数 p = +0.321），
    ⇒ 任何**基于加密**的结论在 flat 网格上**前提为假**。默认翻成 flat 会让一整类
    结论静默失去前提。

    ★ 本条读的是 `inspect.signature`，**不是**在源码上跑正则 —— 写这道门时我第一版
    正则匹配到的是 **docstring 里的 `tip_cap="round"`**，而不是签名，于是它在基线上就红。
    「提到 ≠ 使用」在这里又出现一次，而这次是在**为了防止这类错误而写的门**里。
    """
    p = inspect.signature(onera_m6_wing_mesh).parameters["tip_cap"]
    assert p.default == "round", (
        "onera_m6_wing_mesh 的 tip_cap 默认不再是 'round'（读到 %r）。"
        "★ flat 的凸边在加密下发散（p = +0.321），默认翻成 flat 会让每一个基于加密的"
        "结论静默失去前提" % (p.default,))
