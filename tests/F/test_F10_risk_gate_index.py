"""F10 — 风险↔门机械索引(W3.2,2026-10-12)。

**为什么有这道门**:design.md §12 列了 5 条已知风险,而「哪道门在守它」此前
**没有一条能机械地查到** —— 第 43/44 轮的线程/缓存事件就是「风险没有门」的代价。
⇒ 每条风险带一行「守它的门」指针,本测试断言**每个指针可解析**
(文件存在、测试函数存在),且「无门」的登记数有上界。

★ 反恒真设计:读的文件**路径可注入**(G-TEETH 用扰动副本喂它);
「无门」计数有上界,新增无门必须显式改锁 —— 索引不许悄悄烂。
"""
import ast
import os
import re

import pytest

from tests.conftest import REPO_ROOT

DESIGN = os.path.join(str(REPO_ROOT), "docs", "design.md")
SECTION12 = re.compile(r"## 12\..*?(?=^## )", re.S | re.M)
POINTER = re.compile(r"守它的门\*\*[:：]")
TESTREF = re.compile(r"(tests/[A-F]/[\w.]+\.py)(?:::(\w+))?")
RISK = re.compile(r"^\d+\. \*\*", re.M)

#: 登记中的「无门」数上界(2026-10-12:风险 3 强激波 —— 有运行时缓解、无门)。
#: 新增一条「无门」必须显式改这个数 —— 那是一次裁决,不是顺手。
KNOWN_UNGATED_MAX = 1


def _section12(path=None):
    with open(path or DESIGN, encoding="utf-8") as fh:
        text = fh.read()
    m = SECTION12.search(text)
    assert m, "design.md 里没有 §12"
    return m.group(0)


def _pointers(sec):
    """每条「守它的门」行 → [(file, test_or_None), ...];无门行 → None。"""
    out = []
    for line in sec.splitlines():
        if not POINTER.search(line):
            continue
        if "无门" in line:
            out.append(None)
            continue
        refs = TESTREF.findall(line)
        out.append(refs or "EMPTY")
    return out


def _assert_pointers_resolve(sec):
    for refs in _pointers(sec):
        if refs is None:
            continue
        assert refs != "EMPTY", "指针行里没有任何 tests/ 引用"
        for path, func in refs:
            full = os.path.join(str(REPO_ROOT), path)
            assert os.path.exists(full), f"指针指向不存在的文件:{path}"
            if func:
                tree = ast.parse(open(full, encoding="utf-8").read())
                names = {n.name for n in ast.walk(tree)
                         if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
                assert func in names, f"{path} 里没有 {func}"


def test_every_risk_has_a_pointer_line():
    sec = _section12()
    risks = RISK.findall(sec)
    ptrs = _pointers(sec)
    assert len(risks) == 5, f"§12 的风险数变了:{len(risks)}"
    assert len(ptrs) == len(risks), (
        f"§12 有 {len(risks)} 条风险但只有 {len(ptrs)} 条「守它的门」指针行 —— "
        "每条风险必须有一行(可以是「无门(登记)」)")


def test_every_pointer_resolves():
    _assert_pointers_resolve(_section12())


def test_the_ungated_count_is_capped():
    n = sum(1 for x in _pointers(_section12()) if x is None)
    assert n <= KNOWN_UNGATED_MAX, (
        f"「无门」登记数 {n} 超过上界 {KNOWN_UNGATED_MAX} —— "
        "新增无门是一次裁决(改 KNOWN_UNGATED_MAX 并写明登记),不是顺手")


def test_the_index_is_behaviorally_load_bearing(tmp_path):
    """G-TEETH(自身):扰动副本必须让检查红。不红 = 本测试自己是恒真。"""
    text = open(DESIGN, encoding="utf-8").read()
    bad = text.replace("test_A28_wingbody_topology.py",
                       "test_A28_DOES_NOT_EXIST.py", 1)
    assert bad != text, "扰动没打上(自检)"
    p = tmp_path / "design_bad.md"
    p.write_text(bad, encoding="utf-8")
    with pytest.raises(AssertionError, match="不存在的文件"):
        _assert_pointers_resolve(_section12(str(p)))
    # 删掉一条指针行 ⇒ 指针数 < 风险数 ⇒ 红
    lines = [ln for ln in text.splitlines() if "守它的门" not in ln or "风险 1" in ln]
    # 直接删含守门的头一行等价做法:删掉第一条指针行
    import re as _re
    text2 = _re.sub(r"^.*守它的门.*\n", "", text, count=1, flags=_re.M)
    p2 = tmp_path / "design_missing.md"
    p2.write_text(text2, encoding="utf-8")
    sec2 = _section12(str(p2))
    assert len(_pointers(sec2)) < 5
    with pytest.raises(AssertionError):
        assert len(_pointers(sec2)) == 5
