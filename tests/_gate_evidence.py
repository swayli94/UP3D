r"""C / D 类门的**证据回归锁**：新鲜计算 vs 已提交 `summary.csv`。

★★★ **为什么需要它，用测量说。** `tests/conftest.py` 的 `gate_evidence_dir` 文档串承诺
「平时跑**不写**，**断言对着已提交的 `summary.csv`** ⇒ 代码一改答案就红」。
**实测 2026-08-28 那是假的**：包住 `open` 跑全套，`cases/gates/` 下 **19 个证据 CSV 里
16 个无人读** —— 门重算一遍，然后对着测试文件里的**硬编码字面量**断言，CSV 是
**只写不读**的产物。C/D 门断言里的 81 个浮点字面量中，**71 个**在 **13 道完全不读
自己证据**的门里。

★★ 后果是**两个方向**都会静默漂移：用 `PYFP3D_GATE_FIGURES=1` 刷新了 CSV，新数进 git 而
没有任何东西比对它与门断言的那个数；反之改了门里的字面量而 CSV 没刷新，两者分家。
**与 F06 那次（字面量对字面量、参考文件从没被打开）是同一个缺陷，只是这次是整类的。**

---

## ★★★ 设计上的三个坑，本 session 各踩过一次，写在这里免得再犯

1. **不能「值 == 它自己的来源」**。若门先从 CSV 读、再断言等于 CSV，那**按构造恒真**
   —— F06 实测：把参考改成 0.64，全套依然绿。⇒ 正确形状是 **门重算 → 与 CSV 比**，
   **CSV 是实测值的唯一真值来源**。
2. **判据阈值不进 CSV**。`ORDER_MIN`、`SUP_RMS_MAX` 这些是**标定**不是测量，
   它们该留在门里写死；进了 CSV 就变成"拿门的阈值对门的阈值"。
3. ★ **刷新是两遍，不是一遍**（实测踩过两次）：同一次 `PYFP3D_GATE_FIGURES=1` 运行里，
   本锁在导图腿**写盘之前**就跑了 ⇒ 它比的是**旧** CSV，必然红一次。
   正确流程：**先带标志跑一遍刷新，再不带标志跑一遍验证**。
4. **证据要按真值精度存，不是显示精度**。D13 的样板第一次跑就红在 `asym_ibl` 的
   3.9e-05 上 —— 那不是求解器差异，是那一列按 `.4f` 写盘。⇒ **实测列一律 `.9e`**
   （见 `fmt`），否则容差得从存储格式反推，而那是把仪器的刻度当成被测量的性质。
"""
import csv
import os

MEASURED_FMT = "{:.9e}"
#: 键列（`level` / `surface` / `x_c` 这类）**不改精度**：它们是索引不是测量，
#: 改了会让查表的浮点比较失配。★ 这一条是实测逼出来的，不是洁癖。
DEFAULT_REL_TOL = 1e-6


def fmt(v):
    """实测值的写盘格式。★ 全项目统一 `.9e`，理由见模块 docstring 第 3 条。"""
    return MEASURED_FMT.format(float(v))


def assert_matches_committed(gate_dir, fresh, measured, rel_tol=DEFAULT_REL_TOL,
                             key_of=None, refresh_hint=None):
    """把**本次算出的**实测值与已提交的 `summary.csv` 逐行比对。

    Args:
        gate_dir: `cases/gates/<门号>/`（`gate_evidence_dir` fixture 给的那个）
        fresh: ``{行键: {列名: 数值}}`` —— **本次运行**算出来的，不是从 CSV 读的
        measured: 要比对的列名（**只列实测列**；`*_ref` / `*_exact` 这类另有真值来源的
            不必在这里锁，它们由门自己直接读参考文件）
        rel_tol: 相对容差。★ 解是确定性的，所以这是**机器级**余量；若它需要放宽，
            那本身就是一个要查的信号，不是可以顺手调大的旋钮。
        key_of: ``row -> 行键``；默认取 ``(level, surface, x_c)`` 里存在的那些
        refresh_hint: 红了之后给人看的刷新命令

    ★ 断言消息里带**刷新命令 + 纪律 11 提醒** —— 因为刷新证据往往同时要改门 docstring
    的表、`cases/gates/INDEX.md` 的条目和台账，漏一处就是下一次审计的发现。
    """
    path = os.path.join(str(gate_dir), "summary.csv")
    assert os.path.exists(path), (
        f"已提交的证据 {path} 不存在。生成：{refresh_hint or 'PYFP3D_GATE_FIGURES=1 pytest <本门>'}")
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert rows, f"{path} 是空的"

    if key_of is None:
        def key_of(r):
            return tuple(r[k] for k in ("level", "surface", "case", "x_c", "nx", "leg")
                         if k in r)

    #: ★★★ **重复键 = 错误，不是可以平均掉的东西**（2026-08-28 实测）：C05 的证据里
    #: 阶梯腿与反例腿都是 `nx = 200`，`(nx,)` 分不开 ⇒ 锁会拿反例行去比阶梯值。
    #: 这次两者差 50 % 所以红得很响，**但若两行值接近就会静默比错行** —— 所以守卫必须在。
    seen = {}
    for r in rows:
        k = key_of(r)
        if k in seen:
            raise AssertionError(
                f"{path} 里键 {k} 出现多次（第 {seen[k]} 行与本行）——"
                " 键列不能唯一标识一行，锁会比错行而且可能静默通过。"
                " 请给证据加一个区分列（例如 `leg`），并把它纳入 key_of。")
        seen[k] = len(seen) + 1

    n = 0
    for r in rows:
        k = key_of(r)
        if k not in fresh:
            continue
        for col in measured:
            if col not in r or r[col] in ("", "-", "nan"):
                continue
            want = float(r[col])
            got = fresh[k].get(col)
            if got is None:
                continue
            n += 1
            rel = abs(float(got) - want) / max(abs(want), 1e-30)
            assert rel <= rel_tol, (
                f"{k} 的 {col}：本次算出 {float(got):.9e}，已提交证据是 {want:.9e}"
                f"（相对差 {rel:.2e} > {rel_tol:.0e}）\n"
                f"  ★ 若这是有意的代码改动：{refresh_hint or 'PYFP3D_GATE_FIGURES=1 pytest <本门>'}\n"
                "    刷新证据，**并按纪律 11 grep 被移动的数字** —— 门 docstring 的表、"
                "`cases/gates/INDEX.md` 的条目、以及引用它们的判定文件都要一起改。")
    assert n > 0, (
        f"{path} 里没有一行与本次运行的键对上（比了 0 个数）——"
        " 键的构造变了？这种情况下断言会**静默通过**，所以这一条必须在。")
    return n
