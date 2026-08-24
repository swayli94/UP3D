# UP3D / pyFP3D

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in
supersonic zones, wake cut + Kutta condition for lift.

## Environment — do this exactly (it is load-bearing, not hygiene)

```bash
conda env create -f environment.yml     # or: conda env update -f environment.yml --prune
conda activate up3d
export PYTHONNOUSERSITE=1               # NOT optional -- see below
```

Three requirements, each of which has already broken this project once:

- **`PYTHONNOUSERSITE=1`** — measured: without it `import pyamg` **fails** inside
  `up3d`, because a `~/.local` pip pyamg shadows the environment's own.
  `.vscode/settings.json` and `.env` set it for IDE-launched runs; a shell run
  must export it.
- **`scipy >= 1.12`** — a *code* requirement, not a preference: `solve/linear.py`
  calls `spla.cg(..., rtol=...)` with no version shim. Pinning scipy 1.11.4
  produced 81 failures and 50 errors, every one `TypeError: unexpected keyword
  argument 'rtol'`.
- **Do not work in anaconda `base`.** It is self-inconsistent (old scipy, a
  user-site pyamg needing a newer one) and stayed silently broken for a week.

`pip install -e ".[dev]"` still works for a solver-plus-tests install, but the
conda environment above is what every committed number was produced in.

On a **headless Linux box (incl. WSL)** the `gmsh` wheel dlopens two system
libraries that are not part of the wheel:

```bash
sudo apt-get install -y libglu1-mesa libopengl0     # libGLU.so.1, libOpenGL.so.0
```

Without them `import gmsh` fails with `OSError: libGLU.so.1: cannot open shared
object file` and every mesh generator dies. If you cannot install system
packages, download the two `.deb` files, extract them locally and point
`LD_LIBRARY_PATH` at the extracted `usr/lib/x86_64-linux-gnu`.

## Meshes

The 3-D `.msh` files are **gitignored**; regenerate them with the scripts under
`cases/meshes/*/generate_*.py` (all runnable standalone, no `PYTHONPATH` needed;
the whole set takes ~3.5 min).

★ Since GS4.0 each write also emits a tracked `<name>.msh.manifest.json`
sidecar (`pyfp3d/mesh/manifest.py`) carrying the mesh's sha256 and its
node/tet/boundary counts. **Quote it whenever you quote a number.** A mesh file
is part of a result's provenance exactly like the thread count: on 2026-08-04 a
flat→round tip-cap change went through the gitignored meshes invisibly to git,
and the resulting 2.0 % discrepancy was restated as "unexplained" in eight
documents for five days.

## Running

```bash
export NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
pytest tests/                                                  # ~8 min @8 threads
PYFP3D_TRANSONIC_GATES=1 python bench/run_capability_locks.py   # fast capability tier
PYFP3D_TRANSONIC_GATES=1 pytest tests/                          # full gated set, ~1-2 h
```

- **Cap the BLAS threads too**, not just Numba: missing them costs ~33 % and
  fails a wall-clock gate.
- ★ **Wall-clock is a calibration of the machine, not a property of the solver** —
  the same answer has been measured at 1.6×, 1.9× and 5.4× spreads on this box.
  Quote any wall time together with its thread count and load average.
- ★ **Quote a gated test count with its thread count.** One leg is deliberately
  a non-strict xfail because its outcome flips between 8 and 16 threads.

Baselines (2026-08-16): ungated `479 passed + 12 skipped + 2 xfailed` @8 threads;
gated `488 passed + 1 skipped + 3 xfailed + 1 xpassed` @16 threads.

`PYFP3D_NOJIT=1` swaps `@njit` for identity so print/pdb work.

## Where to start reading

Read these in order; each is authoritative for one thing, and none of them
duplicates another (a copied entry point forks from its original — a failure
this project has logged repeatedly).

| Document | Authoritative for |
|---|---|
| [phases/p2/docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md](phases/p2/docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md) | **what the solver can do, which routes are closed, which gaps are open** |
| [phases/p3/docs/dev_phase_three/20260816-1000-gs41-initiation.md](phases/p3/docs/dev_phase_three/20260816-1000-gs41-initiation.md) | **the single entry point for phase 4** — read its §0 first |
| [phases/p2/docs/dev_phase_two/roadmap.md](phases/p2/docs/dev_phase_two/roadmap.md) | the plan: product metrics M1–M5, stages S0–S6, working principles, rulings D1–D7 |
| [phases/p4/docs/dev_phase_four/progress.md](phases/p4/docs/dev_phase_four/progress.md) | current phase, one row per round |
| [docs/design.md](docs/design.md) | theory and numerics: equations, wake/Kutta, discretization, kernel rules |
| [phases/p1/docs/inspection/](phases/p1/docs/inspection/) · [p2](phases/p2/docs/inspection/) · [p3](phases/p3/docs/inspection/) | independent audits, archived with their phase (`docs/inspection/` deleted 2026-08-24) |
| [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) | layout, per-module status, known gaps |
| [CLAUDE.md](CLAUDE.md) / [docs/agent-rules.md](docs/agent-rules.md) | working rules for the coding agent |

★ **Phase-one documents are frozen history, not plans.** `phases/p1/` and
`phases/p2/` hold the archived roadmaps, ledgers and evidence of finished
phases. Their technical facts — especially the negative results — are still
citable, but they must be cited as historical records. Anything under
`phases/p1/` that imports `pyfp3d.wake` is **not runnable**: that package was
deleted in phase 3. For a working pre-reorganisation tree use
`git worktree add ../up3d-prereorg d224223`.

`cases/reference_data/` is the only external judge in the project. **Never edit it.**

## Headless artifact convention

Visualization checks are generated by scripts (matplotlib `Agg`, PyVista
off-screen), never by GUI inspection. ★ The evidence itself is a **committed** PNG/CSV under a
tracked `results/` directory -- the old gitignored `artifacts/<gate_id>/` target was deleted
2026-08-24 precisely because it was never in HEAD. A number that exists only in prose is not
evidence.
