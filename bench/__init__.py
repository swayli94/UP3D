"""`bench/` —— 产品指标门、共享库、有 cadence 的 runner 与工具。

★ 这是一个 package（2026-08-24 加），因为 `bench/recipes.py` 需要被
`tests/` 与 `cases/demo/` 以 `from bench.recipes import ...` 引用。
一次性研究脚本不在这里 —— 它们随其阶段归档进 `phases/p*/bench/`。
"""
