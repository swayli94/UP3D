# pyFP3D — Project Instructions for Claude Code

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in supersonic
zones, wake cut + Kutta condition for lift. Target: wings at M∞ 0.3–0.87,
workstation-scale (minutes for 1–3 M nodes).

## Document map (read the relevant one before coding)

Docs were split by track on 2026-07-15; the old monolith paths remain as thin
indexes, so historical references like "roadmap.md Track B" or "demo_report §P4"
resolve through one hop.

- [docs/overview.md](docs/overview.md) — human-readable snapshot: per-track
  status table, **document map** (which file is authoritative for what, and when
  to update it), regression-baseline lineage, long-standing open items.
- [docs/roadmap.md](docs/roadmap.md) — **active tracker index**: working rules,
  gate-ID/renumbering conventions, one-line status per track. The phase entries,
  gate checklists and progress ledgers live in **[docs/roadmap/](docs/roadmap/)**
  (`track_p.md` P0–P14 solver, `track_m.md` M0–M5 meshing, `track_b.md` B1–B32
  level-set wake — **B16/B17 far-field aux pin + `pin_gamma`, B18 wing-body
  transonic, B19 LS-Jacobian exactness, all ✓ CLOSED 2026-07-18; B20
  mixed-plain main-field density ADOPTED PERMANENTLY + re-baselined and B21
  N1 freeze-capture fix (restores the M6-medium M0.84 ramp; GB20.7's
  "capability loss" verdict overturned), both ✓ CLOSED 2026-07-19; B23
  junction discriminator (the pocket = the wake inboard free-edge
  singularity), B24 waterline-extension route closed (negative) and B25
  `inboard_clip` CURES the pocket (corrM 14.66→0.63, default None
  bit-identical), all ✓ CLOSED 2026-07-19; B26 post-cure LS ceiling
  re-measured = the conforming site (medium 0.7625 / coarse 0.84 reached)
  and B27 B18 demo refresh (checks 8/8 PASS, 336/336 bit-identical;
  transonic cross-model M0.65 2.4% PASS / M0.75 2.5%), both ✓ CLOSED
  2026-07-20 — the B18 "junction-limited" story is RETIRED; B28 cl_fus
  decoupling + GB9.4 re-spec (the "fuselage spurious lift" label retired;
  out-band cross-model ≤15%, medium 7.0% PASS) and B29 flat-fragment adopted
  as the wing-body LS production config, both ✓ CLOSED 2026-07-20; B30
  (b)-class ceiling attribution (conforming stall and LS+clip death = the
  SAME wing-tip P13 free-edge singularity + high-M Newton, not a wake-model
  pocket) ✓ CLOSED 2026-07-21; B31 C-class wing-tip cure (production
  pressure+taper cures the conforming 0.83 dying level via the FD-verified
  Gamma-pin row blend; LS-side C-class closed negative) and B32 conforming
  tip_taper adopted (wing-body medium ceiling M0.79 → **M0.84 reached**,
  cl_p 0.2738, 0 clamps; weld-sign per-step refresh rolled back as
  ill-posed), both ✓ CLOSED 2026-07-22; **P11
  curved wall elements ✓ CLOSED 2026-07-19 in track_p — measured NEGATIVE,
  G1.6 re-attributed to intrinsic P1 capability at h=0.08 (not the wall
  variational crime), route fork = user's call** —
  `track_v.md` V1–V4 viscous, designed-not-started, `track_a.md` A1–A3
  verification & analysis; **A3 ✓ CLOSED 2026-07-18** = the response to the
  2026-07-17 independent inspection (docs/inspection/; the 2026-07-19
  second-round inspection's N1/D1–D10 findings were executed by B21 + the
  errata wave). **Next phase = user's call.**)
  "What phase are we in" and "what gate is open" live there, nowhere
  else. Track B numerics live in a separate spec,
  [docs/design_track_b.md](docs/design_track_b.md) (it supersedes DN1).
- [docs/design.md](docs/design.md) — theory & numerics reference: equations (§2–§3),
  wake/Kutta (§4), BCs (§5), discretization (§6), Numba kernel rules (§7), solver
  strategy (§8), V0–V6 validation ladder (§10), risks/mitigations (§12); §11 is a
  pointer to roadmap.md + docs/roadmap/ since 2026-07-15.
- [docs/demo_report.md](docs/demo_report.md) — **evidence dossier index** (per-
  phase directory table); the evidence sections live in
  **[docs/demo_report/](docs/demo_report/)** (`track_p.md`, `track_m.md`,
  `track_b.md`, `track_a.md`): one self-checking demo per phase under `cases/demo/<phase>/`
  with committed figures + measured gate numbers. When a phase closes, add its
  demo section to the matching track file and a row to the index.
  **A claim without a committed artifact is not evidence.** The 2026-07-13 audit
  found the P13 M0.84 transonic result (cl_KJ 0.2866 ⇒ "the 0.019 gap is
  resolution" ⇒ "P11's lift case is refuted") existing as *prose only* — no
  script, no CSV, no cached `.npz` — after a P11 ledger status had already been
  changed on its strength. If a run is too expensive to repeat, that is the
  reason to commit its CSV, not a reason to skip it.
- [docs/analysis/](docs/analysis/) — analysis/review reports (capability
  reviews etc.), dated snapshots, non-normative. [docs/archive/](docs/archive/)
  — historical archives (e.g. the pre-2026-07-15 agent-rules narrative);
  never a coding spec (rule 11). `docs/discussion_notes/` was **DELETED
  2026-07-14** (commit 0e4895a; history via
  `git show 8aa4aee:docs/discussion_notes/<file>`). The rule stands: plan
  against **roadmap.md/roadmap-track gates + design.md numerics only**.
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) — layout, per-module status, and
  **"Known gaps"**: read it before touching the G1.6 sphere-Cp problem (formerly
  G1.2; P1 gates renumbered 2026-07-06, mapping in roadmap.md) — it is already
  root-caused, and the G1.3/G1.4 oracles ruled out boundary-data corrections
  (design.md §5.1.2); do not re-propose h-refinement, recovery tweaks, Nitsche,
  or flux-data corrections. Open route: Option C gate re-spec + curved elements.
- `cases/reference_data/` — ground truth, never edit.

## Hard rules and current phase

@docs/agent-rules.md

## Environment — run everything in `up3d` (added 2026-08-04, measured)

    conda env create -f environment.yml      # or: conda env update -f environment.yml --prune
    conda activate up3d
    export PYTHONNOUSERSITE=1                # NOT optional -- see below

`PYTHONNOUSERSITE=1` is load-bearing, not hygiene: **measured** that without it `import pyamg`
FAILS inside `up3d`, because a `~/.local` pip pyamg shadows the env's own. With it set, all nine
project imports resolve to the env. `.vscode/settings.json` and `.env` set it for IDE-launched
runs; a shell run must export it.

`scipy >= 1.12` is a **CODE REQUIREMENT**, not a preference: `pyfp3d/solve/linear.py` calls
`spla.cg(..., rtol=...)` and `spla.gmres(..., rtol=...)` with no version shim, and `tol -> rtol`
landed in scipy 1.12. Pinning scipy=1.11.4 "to keep the numerics identical" produced 81 failures
and 50 errors, every one `TypeError: cg() got an unexpected keyword argument 'rtol'`.

★ Do NOT work in anaconda `base`. It is self-inconsistent (scipy 1.11.4 whose `cg` has no `rtol`,
a user-site pyamg 5.3.0 that needs scipy >= 1.12, and a half-removed user-site scipy), and it was
inconsistent for a week without saying so -- the failure only surfaced when a code path finally
imported pyamg. Two separate environment failures in one day came out of that mixture.

## Reporting conventions fixed by measurement (2026-08-04)

- **Quote the tip taper's bias.** Production is `tip_taper=("vanish_smooth", 0.05·b_semi)` on the
  conforming wing-body and it costs **−1.3 % cl_p** (medium, measured −1.514 % at M0.50 and
  −1.31 % at M0.84). It is a MODEL bias, so any reported wing-body lift carries it. It stays at
  0.05: removing it fails M0.88 medium outright (6/6 clamps), and the cheaper r_c = 0.025 was not
  adopted.
- **The entropy correction is a model correction, NOT a stability trick.** It is default ON and
  stays on. Describing it as buying robustness is wrong twice over: measured, it makes the peak
  STRONGER (M²max 6.426 ON against 5.607 OFF), and its justification is that S1b measured the
  M0.80 shock moving INTO the Euler anchor band where isentropic sat outside it.
- **Round tip is the default; flat must be asked for by name.** `wing3d`'s `tip_cap` default is
  `"round"`, standard level names build round, and flat is only reachable through explicit
  `_flat` level names -- because P13/G13.3 measured the flat cap's sharp convex edge DIVERGING
  under refinement (peak-Mach exponent p = +0.321), so any refinement-based claim on a flat-cap
  mesh has a false premise. Measured consequence: the LE-band convergence order is 0.37 on flat
  against 0.87 on round. Only P13 and M5, which exist to MEASURE flat, use the `_flat` levels.
- **The canonical ladder is `(xcoarse, coarse, medium)`.** `fine` is out of every test and demo:
  it costs 345-1561 s per flow point, is meaningless against the 2 CPU-min target, and only ever
  added a third level where three already existed. `xcoarse` replaces it at the cheap end, which
  is where the third level was actually missing.

## Diagnosing a non-converged solve — CLASSIFY it, never report "conv=False"

Measured 2026-08-04: labelling every failure `conv=False, |R|=…, n/m clamps` and then correlating
it against a knob means a correlation analysis may be separating a MIXTURE OF DIFFERENT DISEASES
with one number. Five hypotheses died to this in one day. The modes have distinct signatures and
completely different fixes:

| mode | signature | what to fix |
|---|---|---|
| clamping | `n_limited` (m_cap) or `n_floored` (rho_floor) > 0 — **and where those cells are** | physics/geometry |
| limit cycle | residual oscillates, `descent10 ≈ 1`, residual values revisited | frozen-selection churn |
| ill-conditioning | `n_gmres_stalled` > 0, GMRES iterations climbing | preconditioner |
| line-search collapse | λ driven to ~0, steps rejected | globalisation |
| σ-transport | `accept_reason = sigma_transport_not_converged` | entropy correction |

`bench/run_le14_common_root.py::classify_failure` implements it; the solver's own `accept_reason`
outranks any inferred signature. **A "budget_limited" leg (still descending, ran out of
iterations) is NOT a failure** — mistaking one for a capability limit retracted a whole "dead
band" finding.

★ **`budget_limited` requires a MONOTONE tail, not just a descent ratio** (fixed 2026-08-05 after
the reverse error): a residual ratio over a 10-step window is an artefact of where that window
lands in a cycle. Measured — a period-3 limit cycle (9.037e-05 → 8.059e-05 → 4.026e-05,
repeating to 4-5 figures) produced descent10 = 2.0021 and so fell on the `> 2.0` side of a
threshold I had picked by hand, and that hand-picked number decided a published conclusion.
Raising n_newton_max 60 → 200 left |R| BIT-IDENTICAL, which is what exposed it. Any fixed
threshold is a CALIBRATION, not a guarantee — the same lesson as the EW forcing and the taper
r_c. Validate a classifier change on cases whose answer is already known by measurement (here:
the period-3 cycle, and LE-1's ls_naca_medium which converged in 6 more steps when uncapped).

★ **`accept_reason` strings are solver-internal, and some are ACCEPT routes, not failures.**
`assignment_cycle` means the live residual stopped improving across upwind-donor refreshes, so
the driver accepts that Mach level at the selection-discontinuity floor (`newton_ls.py:745-750`).
I read it as a failure mode from its name alone. Read the code before naming a mode.

★ **σ-transport: CLOSED 2026-08-05, and it was a false failure.** Seen three times (LE-4
`conf_wb_coarse` M0.78; p13 round-tip coarse M0.5; LE-16 r_c = 0.0325 where |R| = 2.19e-14),
always with the flow and Kutta residuals already satisfied. Root cause: a length-2 HARMLESS
donor cycle (2 of 68624 and 2 of 90099 elements, s = 1 on both) can never satisfy
`transport_sigma`'s termination test "every ancestor is a genuine root", so all 24
pointer-doubling rounds burn and the caller refuses to report convergence per the GS1.4
clamp-not-silent contract. The graph-theoretic acyclicity requirement comes ONLY from the
entropy correction — the upwind density and Jacobian Term 3 read `upstream[e]` one hop and are
indifferent to cycles — so that cycle was harmless until S1b introduced the correction.
Fix: break harmless cycles (all s exactly 1) in a LOCAL COPY of the donor map, leave shocked
cycles intact so they still collapse and are still refused, termination test unchanged.
Measured A/B: acceptance flips False -> True while sigma_min and cl_p agree to 10 digits with
identical refresh and step counts (`bench/gate_results/sigma_fix_verify{,_prefix}.csv`) —
so it changes the VERDICT, not the numbers. Note:
`docs/dev_phase_two/20260805-0200-sigma-transport-root-cause.md`.
Do not re-attribute this to the tip: the affected elements span the whole mesh (z/b 0.002-1.650).

## ★★ Prove a kernel property with an independent ORACLE, not with cases you invented

Measured 2026-08-05, at the cost of two wrong fixes in one round. Three candidate termination
criteria for `transport_sigma`, and **both wrong ones passed every case I could think of**:

| criterion | hand-built cases | random-graph oracle | outcome |
|---|---|---|---|
| "the product stopped moving" | ✗ collapse locks went red at once | — | reverted |
| "`P[A[e]] == 1`" (the ancestor contributes nothing) | **✓ all green, incl. 5 cycle types** | **✗ failed 4 of 5 graphs** | reverted |
| break harmless cycles in a copy | ✓ | ✓ | adopted |

The second one drops a shock sitting further upstream than the finite segment it inspects. It
survived the two collapse locks AND five cycle variants I built specifically to cover the case
space. What killed it was a test that walks each donor chain hop by hop and compares — **sharing
no code with the implementation** — on random graphs containing roots, long chains, a harmless
cycle and a shocked one.

So: for a kernel-level property (completeness, conservation, an accumulated path quantity),
**write the brute-force oracle**. A live solve cannot substitute — it can only show that SOME
criterion fired, never that an accumulation finished. And two traps inside the oracle itself,
both of which I first misread as kernel bugs:
- **do not demand bit-equality of a floating-point product** — the fast algorithm multiplies the
  same factors in a different ORDER and multiplication is not associative (rtol 1e-12, not
  `array_equal`);
- **honour the kernel's stated preconditions in the generated data** — my random graph gave a
  chain ROOT a shock factor < 1, which `shock_factor_sweep` cannot produce, and the documented
  consequence looked like a defect.

## Cold-start seed fallback, and the two levels of seed exposure (2026-08-05/06, measured)

`solve_newton_lifting` retries ONCE with `n_picard_seed = _SEED_FALLBACK` (5) when three
conditions hold together: no seed was requested, there is no warm start either, and the attempt
ended clamped. Each is a measurement, not a guess. The failure it recovers is NOT "seed 0 is bad" —
it is the CONJUNCTION of no seed, a cold start directly at a supercritical M∞, and a mesh fine
enough to resolve the supersonic pocket. NACA0012 M0.80: coarse converges either way, medium dies
with M_max exactly at m_cap (7265 limited / 758 floored) from Newton step 3, and the clamped cells
sit 81 % in MID / 17 % in TE / **0 % at the LE** with the peak at x/c 0.75 — the shock, not the
leading edge. A Mach ramp cures it at the same seed, because the previous level's converged
solution does the seed's job: **ramp and seed are two implementations of one function.**

- The success path is untouched BY CONSTRUCTION: GS1.4 already refuses to report a clamped state as
  converged, so the trigger can never fire on a result a caller would have used.
- The `phi_init is None` condition is what keeps the fallback OUT of a ramp's intermediate levels,
  which warm-start from the last converged level and would be made WORSE by discarding that.
- ★ **It does not fix the SOFT SHIFT.** A cold seed-0 solve can converge with zero clamps and still
  land on a different solution — M1a's fine leg moved cl 0.254830 → 0.252930 and flipped a
  criterion's sign. Nothing fires on that; it needs answer anchors, not a fallback.

★ **A non-strict xfail CANNOT detect a regression** — it is satisfied by pass and by fail alike.
The seed regression sat behind `xfail(strict=False)` from 2026-08-02 to 2026-08-05, and that mark
had been made non-strict for an unrelated reason (thread dependence). If a lock must tolerate an
environment-dependent outcome, record the outcome; do not hide the signal.

★ **Guard REPORT code the way solver code is guarded.** Three losses in one day, all in reporting
rather than numerics: a ~40-minute solve thrown away because `m_last_converged` came back None and
the reporting line did `float(None)` with nothing cached (the FOURTH instance of the caching hazard
below); a verdict loop that indexed a key legs which DIED do not have, crashing after the CSV was
already written; and two legs lost to `m_inf < m_start`, which the transonic driver raises on.
Corollary that paid off immediately: **cache before you report** — the G8.2 re-anchor's every number
was then computed from the npz with zero re-solves.

★ **`m_last_converged` / `m_final` are both None on a FAILED ramp.** The per-level `converged` flags
live in `level_results`, so "the highest level that converged" must be read from there.

## Mesh knobs are not orthogonal — measured, twice

There is **no single-variable knob in `onera_m6_wing_mesh` that leaves the rest alone**:

| knob | intended | measured contamination of the LE band |
|---|---|---|
| `h_wall` | bulk + surface | LE face count **+41.7 %** at fixed `h_edge` |
| `h_far` | far field only | wall triangles 15094 → 15142, LE spacing **0.37 %** |

So a strict 2×2 factorial on LE-vs-bulk is not constructible here. `h_far`'s contamination is ~100×
smaller and can be BOUNDED using an LE-only arm as calibration (its 26.1 % spacing change is 70.8×
larger), but a bound is not cleanliness — do not claim an arm is clean once its guard has failed.

★ **`h_edge` sizes the LE *and* the TE.** Not splitting them makes an "LE refinement" read the P13
tip free-edge singularity instead (GS2.1 addendum #1).

★ **M6 at α 3.06 is NOT subcritical at M0.70** — measured M_max 1.5358 with 214 shock cells, which
falsifies the load-bearing comment in `NEWTON_M6_RECIPE`. Measured in-window Machs where a shock
exists AND refinement stays unclamped: **0.70 and 0.75**. And refinement drives the peak into
`m_cap` at M0.75 too as soon as BOTH the LE and the far field are refined (dose-response: base
1.678 → far 1.988 → LE 2.061 → both 3.000). So "refinement hits the limiter" is not specific to
M0.8395.

## Operating hazards that have cost real time (all measured)

- **Read signatures and import paths; do not recall them.** Five wrong-from-memory calls in one
  day, one of which (`phi_init` at the top level instead of inside `newton_kw`) killed all five
  legs of a 40-minute run with TypeError.
- **`pgrep`/`pkill -f <pattern>` matches YOUR OWN command line.** `pgrep -f "pytest tests/"`
  reported "still running" for a job that had finished; `pkill -f run_le5_taper_coverage` killed
  the invoking shell (exit 144) and lost the script it was writing. **Kill by PID**, and poll a
  log's completion marker rather than the process table.
- **Cache φ AND γ AND the diagnostic history** (`residual_history`, `clamp_history`, `F_history`,
  `n_gmres_stalled`, `accept_reason`). Incomplete caching forced three re-solves of the same five
  states in one day; the third was caught only by killing a fresh run 5 minutes in.
- **Adding a mesh level means growing EVERY per-level dict**, not just `LEVELS` —
  `FUSELAGE_CREASE_MAX_DEG` and `SEAM_MAX_DEG` each failed a generation AFTER meshing.

## Workflow

1. Before coding: find the open gate in the current phase's docs/roadmap/ entry
   and plan against its acceptance criterion. Every visual gate needs a headless artifact
   (`artifacts/<gate_id>/*.png` + `summary.csv`; matplotlib `Agg`, PyVista
   off-screen — never GUI-only checks).
2. After any kernel or assembly change, run the primary regression first:
   `pytest tests/test_v0_freestream.py`
3. Full suite: `pytest tests/` — current baseline **700 passed + 28 skipped +
   2 xfailed** (2026-08-06, measured 1769.56 s @16 threads). ★ The per-round
   lineage lives in [docs/dev_phase_two/progress.md](docs/dev_phase_two/progress.md),
   not here — phase two tracks it there, and this line had been left at the
   phase-one 479 for two weeks.
   ★★ **The GATED set is a separate obligation and it is currently in debt.**
   17 test files carry `PYFP3D_TRANSONIC_GATES` markers and the full gated set has
   NOT been run since the round-tip switch of 2026-08-04. That switch regenerated
   `onera_m6` and `onera_m6_wakefree` as ROUND, while several gated files hold
   absolute anchors measured on the FLAT cap — `test_b7_onera_m6.py`,
   `test_b15_ls_newton_freeze.py`, `test_b22_ls_3d_anchors.py` (γ at rtol 1e-4,
   deliberately tight) and `test_p5_onera_m6.py`. On 2026-08-06 running just
   `test_p4_transonic` + `test_p8_newton` + `test_s1_m1a_envelope` found 2 of 2
   capability locks red. Assume the rest carry bills until measured.
   ★ Note `onera_m6_wakefree` has no `_flat` variants generated, so those tests
   cannot fall back to the flat cap without regenerating meshes.
   Superseded lineage (phase one): 479 passed + 25 skipped + 2 xfailed
   (2026-07-20, B25 inboard fragment clip, +6 passed =
   `tests/test_b1_cut_elements.py::TestInboardFragmentClip` (4) + the same
   file's foot-preference lock (1) + `tests/test_m2_wingbody.py`'s
   waterline-extension lock (1);
   measured 1100.63 s @16 threads;
   the full lineage lives in [docs/overview.md](docs/overview.md), do not
   re-grow it here). Skip
   semantics: the M6 `.msh` are gitignored — 16 M1 tests skip until
   `cases/meshes/onera_m6/generate_onera_m6.py` runs (~30 s); the wake-free
   families likewise (M3 medium ~40 s, M4 ~12 s); the heavy transonic/Newton
   gates (P4 medium + G4.3 sweep, P5, P8 G8.1/G8.2 + FD pocket, B6 M0.80
   dual-mesh + LS-Newton, B7 M6 3D) only run under `PYFP3D_TRANSONIC_GATES=1`
   and make up most of the skips. G8.3's CI reference is 301.66 s.
4. Numba debugging: `PYFP3D_NOJIT=1` swaps `@njit` for identity — print/pdb work.
5. **When a phase closes — the refresh checklist** (extended in A3 after the
   2026-07-17 audit found 17 consistency defects, most of them close-out debt):
   tick the gate in the phase's `docs/roadmap/track_*.md` entry, then update
   **all five** surfaces, because each has gone stale at least once by being
   "obvious enough to skip":
   1. that track file's **progress ledger** (bullet entry + track-status line;
      the track ledgers are wrapped bullet lists, not pipe tables, since
      2026-07-20 — append new phases as bullets),
   2. the **"Current phase"** block in docs/agent-rules.md **and its baseline line**,
   3. **docs/overview.md** (status bullet list + the regression-baseline lineage),
   4. **PROJECT_STRUCTURE.md** — the footer one-liner AND any directory tree
      the phase added files to (this is the one that silently rots),
   5. the **`cases/demo/README.md` table row / `cases/analysis/README.md` bullet**
      for the new demo or study.
   Keep the commit phase-scoped.
   ★ **Backport check.** This codebase has TWO wake paths (conforming
   `newton.py` / level-set `newton_ls.py`, and their Picard twins). When a fix
   lands on one, explicitly check whether the other needs it and record the
   answer in the phase entry — "N/A because ..." is a fine answer, silence is
   not. Two B15-era LS robustness fixes sat un-backported for three phases
   until an external review found them (A3 / kimi C2, C3).
   ★ **Re-baseline erratum checklist** (added by B22 after the 2026-07-19
   inspection: its D1/D2/D7/D8 findings were ALL products of this rule not
   existing). The five-surface list governs NEW sections; it does not catch
   OLD sections quoting numbers a re-baseline just superseded. So any commit
   that regenerates committed evidence must carry, in the phase entry, a
   checklist of every doc location that quotes the old numbers (grep the
   moved values — e.g. `grep -rn "0.2115\|2.4938" docs/`), each one either
   corrected in place or annotated "(pre-X value; superseded, see Y)". A
   number left standing silently is a future audit finding.
6. **Cost caution — do not recompute expensive artifacts casually.** Some
   evidence is committed precisely because regenerating it is slow: the P4 heavy
   demo figures (`cases/demo/p4_transonic/run_demo.py` under
   `PYFP3D_TRANSONIC_GATES=1`) are ~40 min of Picard; the medium G4.1 gate is
   ~17 min; the P5 M6 medium from-scratch continuation+polish is ~45–75 min
   (its solution npz `cases/demo/p5_onera_m6/results/medium_solution.npz` is a
   LOCAL gitignored cache like the .msh — the demo re-solves it when absent, or
   under `PYFP3D_P5_RESOLVE=1`; the committed PNG/CSV are the evidence); the
   ONERA M6 medium/fine `.msh` are minutes to regenerate. Treat the
   committed baseline as authoritative and only rerun the heavy part when a real
   solver/mesh/reference change would move those numbers AND you will commit the
   refresh. For routine edits, verify on the cheap coarse path. Prefer reading a
   committed CSV/PNG over recomputing it.

Gate IDs are `G<phase>.<n>` per Track P numbering (Track V: `GV<phase>.<n>`;
phases V1–V4 are distinct from the validation-case IDs V0–V6). Track P was
renumbered 2026-07-08 and 2026-07-11, Track B 2026-07-12 (×2) and 2026-07-13 —
docs before those dates use the then-current IDs (e.g. pre-2026-07-11 "P9
curved walls"/"P10 backlog" read as P11/P12). The one-line convention summary
and per-phase mapping notes live in docs/roadmap.md and the affected
docs/roadmap/ entries.
