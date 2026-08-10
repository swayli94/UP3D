# bench/ — phase-two drift detection and A/B tools

Established by GS0.3 (see [../docs/dev_phase_two/roadmap.md](../docs/dev_phase_two/roadmap.md)).
Every development round reports against these two tools.

## 1. `run_bench.py` — the repeatable metric set

Cheap (~2 min on 16 threads, coarse meshes only) so it can be run every round.

```bash
NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 \
    python bench/run_bench.py                       # -> bench/results/bench_<stamp>.csv + diff vs baseline
python bench/run_bench.py --with-medium             # + the medium legs (~5 min more)
```

Three groups:

| group | what it prices | why it is here |
|---|---|---|
| `linalg/*` | one elliptic solve (AMG+CG), AMG setup, `splu` cost + fill, assembly | the algorithmic floor: the solver must be judged against its own ingredients, not against another language (audit §3.2) |
| `airfoil/*` | NACA0012 M0.80 α1.25 at `upwind_c` = 1.0 / 1.5 / 3.0 | product metric **M1**: shock position and cl must stop depending on the dissipation constant (audit §3.1) |
| `m1a/*` | NACA0012 **M0.72**/α1.25 (M_max ≈ 1.17) coarse + medium (+ fine under `--with-medium`) | product metric **M1a** (decision D3): the point INSIDE the measured h-convergent envelope (M_max ≲ 1.2). Unlike `airfoil/*` these rows are supposed to **stay put** — they lock what the solver can do today while GS1.6 widens the envelope |
| `wing/*` | ONERA M6 M0.8395 α3.06, `precond` direct vs amg | keeps GS3.1 honest: same answer, less time (audit §3.3) |

`bench/baseline_2026-07-28.csv` is the **frozen phase-two starting point**
(measured on DESKTOP-N6UP769, WSL2, 16 threads). `--compare` prints every metric
that moved beyond its tolerance: physics values 1e-6 relative, wall times 40 %
(shared box), GMRES counts 25 %.

⚠ The baseline is a **drift detector, not a truth reference** (roadmap principle 1).
`airfoil/*` and `wing/*` values are known to be wrong against external data —
that is exactly what S1/S2 are meant to change. When a round moves them ON
PURPOSE, the round file records the move; the baseline is then re-frozen with a
new date and the old one kept.

## 2. `bitcheck.py` — development-time bit-identity A/B

Decision D1: bit identity is meaningful only on the same machine, same
environment, before vs after an edit — so it is a tool, not a permanent test.

```bash
python bench/bitcheck.py --save bench/results/bit_before.npz    # before editing
# ... make the change ...
python bench/bitcheck.py --save bench/results/bit_after.npz
python bench/bitcheck.py --diff bench/results/bit_before.npz bench/results/bit_after.npz
```

Ten probes: residual/matrix assembly, `rho_tilde` at a **subcritical** state
(must stay bit-identical through any S1 change — ν = 0 there) and at a
**supercritical** state (what S1 changes on purpose), the lifting Laplace solve,
the coarse 2-D transonic Newton solve, and the coarse M6 wing solve.
Exit code 0 = every probe bitwise identical.

Self-test 2026-07-28: two consecutive runs on the same machine gave 10/10
bitwise identical probes, i.e. the solver is run-to-run deterministic here — so
a `bitcheck` difference means the edit did it.

`bench/results/` is gitignored (npz are ~1 MB each); commit only the baseline CSV
and quote numbers in the round file.
