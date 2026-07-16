# Analysis & verification cases (Track A)

Cross-cutting *analysis* studies that measure the solver/mesh machinery rather
than add to it — profiling, method A/B comparisons, cost accounting. The
counterpart to `cases/demo/` (per-phase capability demos), but for Track A work.

Each `run_*.py` regenerates every figure/CSV in its `results/` folder and
self-checks the study's gate numbers (exit code 0 = all checks pass). Heavy 3-D
parts run only under `PYFP3D_TRANSONIC_GATES=1` and cache to gitignored `.npz`;
the committed PNG/CSV are the evidence. Roadmap gates:
[docs/roadmap/track_a.md](../../docs/roadmap/track_a.md); conclusions:
[docs/demo_report/track_a.md](../../docs/demo_report/track_a.md).

| Folder | Phase | Run | Runtime |
|---|---|---|---|
| `a1_solver_bottleneck/` | A1 conforming-vs-level-set × Picard-vs-Newton cost benchmark (GA1.1–GA1.5) | `python cases/analysis/a1_solver_bottleneck/run_a1.py` | ~5 min (+ gated 3-D `run_a1_m6.py` ~1 h) |
