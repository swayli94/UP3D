"""F07 — 仓库不变量：**路径解析** 与 **门↔证据一致性**（F 类，默认开启）。

★★★ **这道门为什么存在，用测量说。** 2026-08-24 的重编号把测试移进 `tests/A..F/`，
一层目录的位移让 `__file__` 相对路径全部少一级。这一族错误在**一轮之内咬了 13 次**，
而每一次我都只补了当次那一种写法 —— 因为判据的模式覆盖不到同一件事的另一种合法写法：

| 写法 | 我的第几版仪器看不见它 |
|---|---|
| `parents[N]` | ✓ 看得见（照 2026-08-10 的记录标定的） |
| `.parent.parent` | ✗ 漏 20 处 |
| `os.path.join(HERE, "..", "..", "..")` | ✗ 漏 50 处 |
| `os.path.join(dirname(__file__), "..")` —— **单个 `".."` 且无逗号** | ✗ 漏，直接造成 4 failed |
| `dirname(dirname(...))` | ✗ 列了但**只在 bench 量过，没在 tests 重跑** |
| 路径**逐段拼**给 `spec_from_file_location` | ✗ 漏，造成 2 个 FileNotFoundError |

★★ 而更隐蔽的一层：`REPO_ROOT` 用在**函数体**里，**import 期不报错** ⇒
`pytest --collect-only` 全绿而全套 3 failed + 28 errors。**收集不是这一类的牙。**
★★★ 最后一次判错更值得记：我的检查器问「`REPO_ROOT` 有没有在**任何地方**被 import」，
而某文件有一个**函数内**的 import 就让它放行了 —— 可**另一个函数**里的使用点根本不在那个作用域。
**真正的问题是「在使用处是否在作用域内」，不是「是否出现过」。**

⇒ 所以这道门**不匹配任何写法**，它**实际解析**：把每个含 `__file__` 的模块级赋值 exec 出来看落点，
并按**作用域**判名字可见性。一个穷举检查，替掉一排永远补不全的正则。
"""
import ast
import glob
import os

import pytest

from tests.conftest import REPO_ROOT

TEST_FILES = sorted(glob.glob(str(REPO_ROOT / "tests" / "**" / "*.py"), recursive=True))


def _module_level_names(tree):
    """只算**模块级**可见的名字 —— 函数内的 import 不算（这正是上一版判错的地方）。"""
    out = set()
    for n in tree.body:
        if isinstance(n, (ast.Import, ast.ImportFrom)):
            out |= {(a.asname or a.name).split(".")[0] for a in n.names}
        elif isinstance(n, ast.Assign):
            for t in n.targets:
                if isinstance(t, ast.Name):
                    out.add(t.id)
                elif isinstance(t, ast.Tuple):      # ★ 元组解包：M_INF, ALPHA = ... 也曾漏掉
                    out |= {e.id for e in t.elts if isinstance(e, ast.Name)}
    return out


def test_every_file_relative_path_resolves_outside_tests():
    """每个 `__file__`-相对的模块级路径常量，**实际解析**之后不得落在 `tests/` 里面。

    ★ 判据是**解析结果**，不是写法 —— 于是它对上面表格里每一种写法都成立，
    包括还没被人写出来的那一种。
    """
    bad = []
    for p in TEST_FILES:
        src = open(p, encoding="utf-8").read()
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        for n in tree.body:
            if not isinstance(n, ast.Assign):
                continue
            seg = ast.get_source_segment(src, n) or ""
            if "__file__" not in seg:
                continue
            g = {"__file__": os.path.abspath(p), "os": os,
                 "Path": __import__("pathlib").Path}
            try:
                exec(seg, g)                                   # ★ 真的算一遍
            except Exception:
                continue
            for t in n.targets:
                if not isinstance(t, ast.Name):
                    continue
                s = str(g.get(t.id, ""))
                if not s.startswith(str(REPO_ROOT)):
                    continue
                rel = os.path.relpath(s, str(REPO_ROOT))
                #: 落在 tests/ 里、或指向 tests/cases 这种，就是少了一级
                if rel == "tests" or rel.startswith("tests" + os.sep):
                    if not rel.startswith(os.path.join("tests", "_")):
                        bad.append("%s: %s = %s" % (os.path.relpath(p, str(REPO_ROOT)), t.id, s))
    assert not bad, (
        "这些路径常量解析后落在 tests/ 内部 —— 目录深度变了而表达式没跟上：\n  "
        + "\n  ".join(bad)
        + "\n★ 修法不是补 '..'，是改用 `from tests.conftest import REPO_ROOT`（与深度无关）。")


def test_names_are_in_scope_where_they_are_used():
    """`REPO_ROOT` 在**每个使用点**都必须在作用域内。

    ★★★ 上一版检查器问的是「有没有在任何地方 import」，被一个**函数内**的 import 骗过 ——
    而另一个函数里的使用点不在那个作用域。**import 期不报错，全套才炸。**
    """
    bad = []
    for p in TEST_FILES:
        src = open(p, encoding="utf-8").read()
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        mod = _module_level_names(tree)
        for fn in [n for n in ast.walk(tree)
                   if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]:
            local = {a.arg for a in fn.args.args}
            for n in ast.walk(fn):
                if isinstance(n, (ast.Import, ast.ImportFrom)):
                    local |= {(a.asname or a.name).split(".")[0] for a in n.names}
                elif isinstance(n, ast.Assign):
                    for t in n.targets:
                        if isinstance(t, ast.Name):
                            local.add(t.id)
            for n in ast.walk(fn):
                if (isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)
                        and n.id == "REPO_ROOT"
                        and "REPO_ROOT" not in mod and "REPO_ROOT" not in local):
                    bad.append("%s:%d 在 %s()" % (
                        os.path.relpath(p, str(REPO_ROOT)), n.lineno, fn.name))
                    break
    assert not bad, (
        "REPO_ROOT 在这些使用点不在作用域内（import 期不报错，运行时才 NameError）：\n  "
        + "\n  ".join(bad))


@pytest.mark.parametrize("cls", ["C", "D"])
def test_every_cd_gate_has_an_evidence_dir_and_vice_versa(cls):
    """★★ 门 ↔ 证据目录**双向**一致 —— 按门号做键，所以可机械互查。

    ★ 占位门（`@pytest.mark.skip` 的那些）豁免：它们的门号已占住而内容未实现，
    这**正是裁决 ④「暂时没有数据的先占位」的形状**，不是缺口。
    """
    gates = []
    for p in sorted(glob.glob(str(REPO_ROOT / "tests" / cls / "test_*.py"))):
        src = open(p, encoding="utf-8").read()
        if "PLACEHOLDER" in src:            # 占位，豁免
            continue
        if "gate_evidence_dir" not in src:  # 这个门不产出图（例如纯数值判据）
            continue
        gates.append(os.path.basename(p)[len("test_"):-len(".py")])
    ev = {os.path.basename(d.rstrip("/"))
          for d in glob.glob(str(REPO_ROOT / "cases" / "gates" / "*/"))}
    missing = [g for g in gates if g not in ev]
    orphan = [e for e in ev if e.startswith(cls) and e not in gates]
    assert not missing, (
        "这些 %s 类门用了 gate_evidence_dir 但没有被提交的证据目录 "
        "（跑 PYFP3D_GATE_FIGURES=1 pytest tests/%s 生成）：%s" % (cls, cls, missing))
    assert not orphan, (
        "cases/gates/ 下这些目录没有对应的 %s 类门（门被删了？证据该跟着走）：%s" % (cls, orphan))


def test_every_tests_import_in_live_code_resolves():
    """★★★ 活代码（`bench/`、`cases/`、`pyfp3d/`）里每一条 `tests.*` import 都必须解析。

    ★ 这条是被一次**真实失败**逼出来的：2026-08-26 重生成 M3a 时
    `bench/run_m3_budget.py` 炸在 `from tests import test_p8_newton as _p8` ——
    重编号的改写**漏了 `as` 这种形式**，而它在**函数体内** ⇒ `--collect-only` 看不见；
    而 `bench/` **不在任何周期上**，所以只有真跑那个脚本才会发现。

    ⇒ 本条**不匹配写法**，它对每一条 import **实际做解析**（`find_spec` + 属性回退），
    于是 `from X import Y`、`from X import Y as Z`、`import X` 一视同仁。
    """
    import importlib
    import importlib.util
    bad = []
    files = [f for f in glob.glob(str(REPO_ROOT / "**" / "*.py"), recursive=True)
             if "/phases/" not in f and "/tests/" not in f and "/.git/" not in f]
    for f in files:
        try:
            tree = ast.parse(open(f, encoding="utf-8").read())
        except SyntaxError:
            continue
        for n in ast.walk(tree):
            pairs = []
            if isinstance(n, ast.ImportFrom) and (n.module or "").startswith("tests"):
                pairs = [(n.module, a.name) for a in n.names]
            elif isinstance(n, ast.Import):
                pairs = [(a.name, None) for a in n.names if a.name.startswith("tests")]
            for mod, name in pairs:
                try:
                    ok = importlib.util.find_spec(mod) is not None
                except (ModuleNotFoundError, ValueError):
                    ok = False
                if ok and name:
                    spec = importlib.util.find_spec(mod)
                    if spec.submodule_search_locations is not None:
                        try:
                            ok = importlib.util.find_spec("%s.%s" % (mod, name)) is not None
                        except ModuleNotFoundError:
                            ok = False
                        if not ok:
                            ok = hasattr(importlib.import_module(mod), name)
                if not ok:
                    bad.append("%s:%d -> %s%s" % (
                        os.path.relpath(f, str(REPO_ROOT)), n.lineno, mod,
                        "" if name is None else "." + name))
    assert not bad, (
        "活代码里这些 tests.* import 解析不了（bench/ 不在任何周期上，"
        "所以只有真跑才会炸）：\n  " + "\n  ".join(bad))

#: `bench/bitcheck.py` 自己写出的产物 —— 引用它们是**输出路径**，不是失效引用。
_GENERATED = {"bench/results/bit_after.npz", "bench/results/bit_before.npz"}
_PATHLIKE = __import__("re").compile(
    r"(?:bench|cases|docs|phases)/[A-Za-z0-9_/.-]+\.(?:py|md|csv|npz)")


def test_every_embedded_repo_path_in_live_code_exists():
    """★★★ 活代码里内嵌的每一条仓库路径都必须真实存在。

    ★ 这条是被**同一族错误在一天里咬两次**逼出来的（2026-08-26）：归档把文件搬进
    `phases/p*/` 之后，先是活**文档**里 158 处路径指空，扫完文档才发现**代码**里还有
    48 处 —— 因为上一版的清扫只覆盖了 `*.md`。**仪器覆盖不到的地方就是缺陷藏身的地方。**

    ★★ 两个具体后果，都不是「文档过期」那么轻：
      · `CLAUDE.md` 硬规则 #2 让人跑 `pytest tests/test_v0_freestream.py` —— 那个路径已不存在，
        **规则本身跑不动**；
      · `pyfp3d/solve/newton.py` 引用 `bench/run_capability_matrix.py`，而那是**本期改名**造成的，  (PATH-EXAMPLE)
        库文件指向一个不再存在的脚本。

    ★★★ 而修法本身也踩了一次同族的坑，记在这里：朴素子串替换
    （`docs/roadmap.md` -> `phases/p1/docs/roadmap.md`）会把**已经正确**的那些也换掉，  (PATH-EXAMPLE)
    造出 `phases/p1/phases/p1/docs/roadmap.md` —— 37 处。  (PATH-EXAMPLE)
    **替换前必须排除「已经落在正确路径里」的出现。**
    """
    bad = []
    for f in glob.glob(str(REPO_ROOT / "**" / "*.py"), recursive=True):
        rel = os.path.relpath(f, str(REPO_ROOT))
        #: 归档是历史快照、`build/` 是构建产物，两者都不在活代码里
        if rel.startswith(("phases" + os.sep, "build" + os.sep)):
            continue
        for i, line in enumerate(
                open(f, encoding="utf-8", errors="replace").read().split("\n"), 1):
            for m in _PATHLIKE.findall(line):
                #: ★★ 哨兵：讲述这类缺陷的散文必须能**举出**失效路径而不被判为犯了它。
                #: 「提到 ≠ 使用」—— 本门第一次跑就红在自己的 docstring 上，这是同一族
                #: 错误在一天里的第四次，而这次发生在为抓它而写的门里面。
                if (m in _GENERATED or "..." in m or "DELETED" in line
                        or "PATH-EXAMPLE" in line):
                    continue
                if not os.path.exists(os.path.join(str(REPO_ROOT), m)):
                    bad.append("%s:%d -> %s" % (rel, i, m))
    assert not bad, (
        "活代码里这些内嵌仓库路径不存在（搬档/改名之后没跟上）：\n  "
        + "\n  ".join(sorted(set(bad)))
        + "\n★ 若目标已归档，指向 `phases/p*/...` 的真实位置；若已删除，"
          "就地写明「已删除」并保留引文作为幸存记录 —— 不要留一个指空的指针。")
