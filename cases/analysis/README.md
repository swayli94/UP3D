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
| `a2_te_kutta_fidelity/` | A2 TE/Kutta fidelity attribution (GA2.1–GA2.5; S1 jitter + S2 TE Cp gap) | `python cases/analysis/a2_te_kutta_fidelity/run_a2.py` | ~90 s (+ gated `run_a2_interventions.py`) |
| `c1_ls_jacobian_fd/` | **A3 / GA3.5 + B19** — is the LS Newton Jacobian exact on a 3D mesh? (kimi code-review C1; `run_check.py` measures it before/after the B19 fix, `run_legb_probe.py` measures B19 Leg B's residual asymmetry) | `python cases/analysis/c1_ls_jacobian_fd/run_check.py` · `.../run_legb_probe.py` | ~4 min each |
| `p14_te_pressure_diag/` | P14 pressure-Kutta diagnostics (20 checks; feeds the P14 demo) | `python cases/analysis/p14_te_pressure_diag/run_diag.py` | heavy |
| `b9_fuselage_guardrail/` | B9 fuselage-Cp guardrail (GB9.6, RECORDED — no pass/fail) | `python cases/analysis/b9_fuselage_guardrail/run_guardrail.py` | heavy |

*(Rows for a2/b9/p14 added in A3 2026-07-18: they existed on disk but the
table still listed only a1. Note two rows are NOT Track A — `b9_*` and
`p14_*` are analysis studies belonging to Track B / Track P phases that live
here because they are studies, not capability demos.)*
