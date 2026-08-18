# pyFP3D — Project Instructions for Claude Code

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in supersonic
zones, wake cut + Kutta condition for lift. Target: wings at M∞ 0.3–0.87,
workstation-scale (minutes for 1–3 M nodes).

## Document map (read the relevant one before coding)

★★ **LAYOUT CHANGED 2026-08-10 — read this before following any path below.** Finished
phases were archived into **`phases/p1/`** and **`phases/p2/`** (tracked, not gitignored, so
the evidence stays in HEAD). The link targets in this map were repointed mechanically, so a
line reading `[docs/roadmap.md](phases/p1/docs/roadmap.md)` means: **the name is historical,
the file now lives in the archive.** What stayed in place is what phase 3 still needs, and the
rule was measurable — *if a surviving test or script reads it, it stays*:

- `pyfp3d/`, `cases/reference_data/`, `cases/meshes/`, and **73 test files (every conforming
  anchor)**; the 18 pure level-set test files are archived, since deleting them is phase 3's
  first task;
- `bench/` keeps the **7-module closure** behind the five product metrics **plus
  `bench/gate_results/`** — the capability boundary cites those CSVs by path;
- `cases/demo/` keeps **9 of 33** and `cases/analysis/` **10 of 33** subdirectories, the ones
  kept tests actually read (e.g. `cases/demo/p11_curved_walls/`, which
  `tests/test_laplace_sphere.py` reads for the LIVE G1.6 Option C gate);
- `docs/` keeps overview, design, design_track_v, agent-rules, inspection/, and the **six
  phase-three-facing** files in `docs/dev_phase_two/` (roadmap, progress, the capability
  boundary, the level-set inventory, the precond decision record, the template).

**Start here for phase 3:** [docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md](docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md),
then [phases/README.md](phases/README.md) for what was archived and how to get a
guaranteed-runnable phase-1/2 tree (`git worktree add ../up3d-prereorg d224223`).

Docs were split by track on 2026-07-15; the old monolith paths remain as thin
indexes, so historical references like "roadmap.md Track B" or "demo_report §P4"
resolve through one hop.

- [docs/overview.md](docs/overview.md) — human-readable snapshot: per-track
  status table, **document map** (which file is authoritative for what, and when
  to update it), regression-baseline lineage, long-standing open items.
- [docs/roadmap.md](phases/p1/docs/roadmap.md) — **active tracker index**: working rules,
  gate-ID/renumbering conventions, one-line status per track. The phase entries,
  gate checklists and progress ledgers live in **[docs/roadmap/](phases/p1/docs/roadmap)**
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
  [docs/design_track_b.md](phases/p1/docs/design_track_b.md) (it supersedes DN1).
- [docs/design.md](docs/design.md) — theory & numerics reference: equations (§2–§3),
  wake/Kutta (§4), BCs (§5), discretization (§6), Numba kernel rules (§7), solver
  strategy (§8), V0–V6 validation ladder (§10), risks/mitigations (§12); §11 is a
  pointer to roadmap.md + docs/roadmap/ since 2026-07-15.
- [docs/demo_report.md](phases/p1/docs/demo_report.md) — **evidence dossier index** (per-
  phase directory table); the evidence sections live in
  **[docs/demo_report/](phases/p1/docs/demo_report)** (`track_p.md`, `track_m.md`,
  `track_b.md`, `track_a.md`): one self-checking demo per phase under `cases/demo/<phase>/`
  with committed figures + measured gate numbers. When a phase closes, add its
  demo section to the matching track file and a row to the index.
  **A claim without a committed artifact is not evidence.** The 2026-07-13 audit
  found the P13 M0.84 transonic result (cl_KJ 0.2866 ⇒ "the 0.019 gap is
  resolution" ⇒ "P11's lift case is refuted") existing as *prose only* — no
  script, no CSV, no cached `.npz` — after a P11 ledger status had already been
  changed on its strength. If a run is too expensive to repeat, that is the
  reason to commit its CSV, not a reason to skip it.
- [docs/analysis/](phases/p1/docs/analysis) — analysis/review reports (capability
  reviews etc.), dated snapshots, non-normative. [docs/archive/](phases/p1/docs/archive)
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

## ★★ Four criterion defects in four rounds, all the same shape (2026-08-12, measured)

Phase three's sigma rounds produced a criterion defect EVERY round, and every one was found after
the measurement rather than before it. They are one mistake wearing four hats:

| round | the defect | what it let through |
|---|---|---|
| sigma selection | **one-sided bands** -- only contemplated the control being WORSE | a reversed, family-consistent result could only land in "no direction" |
| sigma strength | **unbalanced panel** -- compared spreads across theta without fixing the converged-seed set | "non-monotone" was triggered by a set-SIZE artefact (the dying seed was the outlier) |
| soft membership | **absolute threshold** (<= 5 %) on families whose baselines differ 150x | an 8.8x DEGRADATION scored as PASS |
| magnitude-preserving | **an OUTPUT treated as an INPUT** -- "match sigma_min back to legacy" | the target was unhittable BY CONSTRUCTION, so the binding leg could only read J3 |

The common sentence: **the criterion did not cover the domain of the quantity it compares.** So
before registering a band, run it against these four questions:

1. **Is it one-sided?** State where the OPPOSITE outcome would land.
2. **Does the independent variable change which samples EXIST?** If yes, fix the sample set first --
   a spread over a shrinking set is not comparable -- and never report a spread over fewer than two
   converged legs as "small"; call it UNDEFINED.
3. **Is the threshold absolute where the baselines differ?** Make it relative to each arm's own
   baseline, or an improvement cannot be told from a tolerated regression.
4. **Is the quantity being matched an INPUT or an OUTPUT?** "Restore X to its baseline" is only
   well-posed if X is something you set. `sigma_min` is a diagnostic of the CONVERGED STATE, so
   when the state moves the target moves with it.
5. ★★ **Are the two numbers I am comparing THE SAME THING?** Added 2026-08-13 after the M6 triage
   hit this family FOUR TIMES IN ONE ROUND: cross-PIPELINE (P14's anchor against the script's own
   output), cross-LEVEL (a function's `levels=("coarse",)` DEFAULT against the evidence's medium),
   cross-PROVENANCE (an 8-thread run against a 16-thread committed row) and cross-TIME (a HEAD run
   against a reference produced before three recorded recipe changes). Before comparing, check that
   both sides come from the same pipeline, the same level, the same thread count and the same code.
   ★ All four were caught by a guard that STOPPED the round rather than by reading the numbers —
   which is the argument for writing an instrument check with a hard stop into every registration.
   ★★ **Fifth member, 2026-08-16: cross-MESH-FAMILY.** The `onera_m6` `.msh` files were regenerated
   on 2026-08-04 when the level names flipped flat → round, so every pre-08-04 anchor
   (`run_m3_budget.P14_ANCHOR` among them) is a FLAT-CAP number. Measured with the mesh file as the
   only variable: coarse cl_p 0.262123 (flat) vs 0.268115 (round) = **+2.29 %** — which is exactly
   the "two pipelines disagree on coarse by 2.0 %" debt that had stood for five days. The two sides
   were never two pipelines; they were one recipe on two mesh generations. **A mesh file is part of
   the provenance — check it like the thread count.**
   ★★★ And the process lesson, which cost more than the measurement did (the whole check was **14
   seconds** of compute): **before restating a debt for the second time, spend the ZERO compute
   needed to pin down what its referent is and where that referent lives.** That 2.0 % had been
   written down with no value and no source, and it then propagated into EIGHT files as "unexplained"
   — while the referent sat as a module constant in the very script that produced the number. A gap
   with no citable referent is not a debt, it is a sentence, and sentences replicate.

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

★ **`m_last_converged` / `m_final` say WHICH MACH THE RETURNED STATE LIVES AT** — on a failed ramp the
returned dict is the FAILED level's state, at a LOWER Mach than `m_inf`, so pairing its `cl`/`phi`
with the requested `m_inf` mislabels the row. `m_last_converged` is None when nothing converged.
★★ **ERRATUM 2026-08-16 (GS4.0), and the correction matters more than the original line**: this
paragraph used to read "both None on a FAILED ramp". They were not None — **they did not exist**, on
any surviving driver. They were set only by `newton_ls.py`, deleted with the level-set route in phase
3, and the conforming `newton.py` had never had them (discipline #9, a backport check that was never
done). So callers using `r.get("m_last_converged", <default>)` got the DEFAULT, silently, instead of
crashing: `bench/run_capability_matrix.py`'s `MACH_NOT_ATTAINED` guard compared `abs(m_att - m)` where
`m_att` was identically `m`, i.e. **a guard that could not fire**, whose own comment said it existed so
as not to trust the driver flags. ⇒ `newton.py::_ramp_honesty_fields` now provides all three
(`m_final`, `m_last_converged`, `target_reached`) on the conforming ramp, unit-tested without a solve
in `tests/test_gs40_provenance.py`. **Read them directly; do not `.get` them with a default** — a
missing key is a library regression and must be loud. `solve_newton_lifting` (single Mach) does not
have them, and there the state is at `m_inf` by construction. The per-level `converged` flags still
live in `level_results`, which is what the three fields are derived from.
★ **General lesson, logged**: when a route is deleted, the "only subtraction" ledger closed on
`tests/` (passed 457 → 457) **cannot see `bench/`** — bench scripts are in no test collection, so they
run on no cadence. **Deleting a route means grepping `bench/` and `cases/` for reads of its return
keys**, not just counting tests.

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

## A guard must cover what the CONCLUSION claims, not what is easiest to check (2026-08-09)

The LE factorial's guard G1 asked "did the far-field arm leave the rest alone?" and checked the
SURFACE mesh: wall triangle count and LE-band tangential spacing, which moved 0.343 % -- reassuringly
clean. It never checked the VOLUME. Measured afterwards, that same arm changed the median cell size
at r = 1-2 MAC by **24 %** -- the mesh hugging the wing -- and the LE band's NORMAL spacing by 1.5 %.

So `h_far` is not a far-field knob, it is a bulk-grading knob, and a conclusion of the form "the FAR
FIELD controls the LE band" could not be supported by a guard that only looked at the skin. The
criterion and every number survived; the attribution had to be rewritten to "the BULK MESH including
near-body grading". Those point at different next steps, so it is not a wording quibble.

★ It also flipped the reading of a second result: "refining the far field makes the error WORSE
(+3.4 %)" is really "COARSENING THE NEAR BODY makes it worse", which is unsurprising rather than
mysterious. A mis-scoped guard does not merely weaken a conclusion -- it can invert the sense of the
one next to it.

- **write the guard against the sentence you intend to publish.** If the claim is about the volume,
  the guard measures the volume.
- **a knob's name is not its scope.** Measure what it changed; `h_far`, `h_edge` (LE *and* TE) and
  `h_wall` have all now been caught acting outside their names in this generator.
- and the standing rule stays: there is NO clean single-variable mesh knob here, so a strict
  factorial has to come from changing the generator, not from picking a better knob.

## The gated set needs a cadence, and it now has one (2026-08-09, measured)

The 17 gated files hold the project's absolute capability anchors and the FULL gated run cost
**4 h 09** (2026-08-06, 8 threads under load), so in practice it was not run on any cadence. The
consequence, measured: the first full gated run since the 2026-08-04 round-tip switch came back
**7 failed**, and one of them (b9) had been red for **seventeen days**. Debt in the gated set is
invisible by construction — the ungated suite stays green while capability locks rot.

So there is now a FAST tier, `PYFP3D_TRANSONIC_GATES=1 python bench/run_capability_locks.py`,
**measured 2026-08-11 TWICE at the same 8 threads on the same box: 891 s and
564 s, both 5/5 green** (keep both: a 1.6x spread with an identical result is the
same wall-clock-is-a-calibration lesson as G8.2's 5.4x) (was 644 s / 7 groups: three
level-set locks left with the route in phase 3, and the conforming wing-body transonic
ceiling lock was added — a capability lock kept OUT of this tier would run only in the
2 h gated set, i.e. 6x less often than the tier it belongs to) — ★ and that cost holds ONLY with the thread caps
pinned: the same script measured **2940 s = 49 min** uncapped on this 24-core box (2026-08-10),
per group up to **20.7×**. The script now pins them itself and prints the resolved values with
the load average, because one of its locks asserts a wall-clock budget. Run it at **every close-out**; run the full gated set at
**phase boundaries**.

★ The script prints WHAT IT DOES NOT COVER every time, with each exclusion's measured cost (b22
medium ~35 min — note a strict xfail still RUNS the solve; p4 ~32 min; p5 45–75 min; b14/b15/b18
heavy ramps). A fast tier mistaken for full coverage would recreate exactly the failure it exists to
prevent, so that list is part of the output, not a comment.

★ And `PYFP3D_TRANSONIC_GATES=1` is checked at entry: without it every gated lock SKIPS and the run
reports a vacuous green.

★★ **Put a wall-clock assert LAST in its test** (2026-08-10, measured, cost ~20 min of diagnosis).
G8.2's `assert wall < 450.0` sat ABOVE its physics anchors, so when the fast tier hit a loaded machine
the leg reported red at 588 s with cl, M_max and the three shock positions **never evaluated** — the
b7/b9 trap again (first failing assert hides the rest). Measured afterwards via `_m6_case` from a bench
script (no test edit): every physics anchor reproduces the committed value to six decimals — cl
**0.268691**, M_max 1.996867, shocks 0.596316/0.540203/0.371440 — and **the same solve then took 109 s**,
against 588 s half an hour earlier on the same box (load average 22 → 15). **A 5.4× spread with a
bit-identical answer** ⇒ a fixed wall-clock bound is a CALIBRATION of the machine, not a statement about
the solver (same family as the EW forcing, the taper `r_c`, and the `descent10` threshold). Keep it as a
gate — a real 2× regression must still fail — but order it after the physics, and read a red as a timing
reading until the physics above it has passed. Evidence `bench/gate_results/g82_anchor_check.csv`,
regenerable by `bench/run_g82_anchor_check.py`.

## Retiring a criterion means grepping the TESTS, not just the demo (2026-08-09, measured)

B28 retired the "the fuselage should carry almost no lift" premise on 2026-07-20 -- its own
checks.csv says `<=5%-of-wing premise RETIRED (physical carryover; B23)` -- and replaced it with a
CROSS-MODEL gap criterion, `|conf_out - LS_out| <= 15 % |conf_out|`, which passes at 7.0 % on medium.
But `test_b9_wingbody_conforming.py` kept `abs(cl_f) < 0.15 * abs(cl_w)`: the retired premise, with
the 15 % apparently borrowed from the new criterion's cross-model number. That test is GATED, so
nobody ran it for seventeen days, and when it finally failed it looked like a regression.

- **re-thresholding a retired premise re-imports it under a new number.** The fix is to remove the
  assertion and RECORD the reading, pointing at where the live gate is.
- **the five-surface close-out ritual needs a sixth item: grep the criterion's own numbers and
  wording across `tests/`.** Demos and docs were covered; the tests were not.
- ★ and read the WHOLE test before re-specifying it. Both b7 and b9 hid later assertions behind the
  first failing one -- b7's transonic gate asserts five things and the log showed only the clamps;
  b9's second layer only surfaced after the first was fixed.

## ★★ Swapping library files for an A/B: the index is the trap, and `git status` is not the check

Two separate incidents on 2026-08-06/09, the second WORSE than the first because it reached a commit.

**Incident 1 (working tree).** A script used `cp` for backups; `$SC` was undefined inside the
heredoc so the backup failed, and the script had no `set -e` so it overwrote both library files
anyway and could not restore them. Fixed by `git checkout -- pyfp3d/`.

**Incident 2 (a COMMIT).** The "safe" replacement script used git — and still lost the work:

    git checkout <rev> -- <paths>     # writes the INDEX *and* the worktree
    git checkout -- <paths>           # restores from the INDEX -> the OLD content comes back

so the old library stayed in the index, and a later `git add -A && git commit` baked it into
`6819d38`, dropping 191 lines of committed work from HEAD. Two commits shipped with the regression
before an unrelated import error exposed it.

The rules that actually prevent this:

- **restore with `git checkout HEAD -- <paths>`**, never bare `git checkout -- <paths>`, after a
  `git checkout <rev> -- <paths>`. Or use `git stash` / `git worktree`, which never touch the index
  for this purpose.
- ★ **`git status` IS NOT THE VERIFICATION.** With index and worktree both holding the old content,
  status can look unremarkable. Verify by **importing a sentinel from the module** —
  `from pyfp3d.solve.newton import _SEED_FALLBACK` — and by running one fast test file. Do this after
  EVERY swap script and again before any `git add -A`.
- **`set -euo pipefail` plus `trap restore EXIT`** in any script that mutates tracked files — but put
  `|| true` on commands EXPECTED to fail (a pytest leg that is supposed to be red will otherwise abort
  the script under `pipefail`, which happened on the first attempt).
- **Never `git add -A` in the same session as a file-swap script** without the sentinel check first.

## Never hand-copy a file git already tracks (2026-08-06, cost: the working tree)

A background script swapped `pyfp3d/solve/newton.py` and `pyfp3d/kernels/entropy.py` to an older
revision for an A/B, using `cp` to a backup path and `cp` back afterwards. Two layers failed at
once: `$SC` was never defined INSIDE the heredoc'd script (it was an outer-shell variable and the
heredoc was quoted), so the backup `cp` wrote to `/_n.py` and got Permission denied; and that script
had no `set -e`, so it carried on and overwrote both library files, after which the restore `cp`
failed too. The working tree sat at the old revision with 190 lines of committed work missing from
it until `git checkout -- pyfp3d/` put it back.

Nothing was lost, because the changes were COMMITTED — which is the whole point:

- **the file is in git, so git is the backup.** `git checkout <rev> -- <paths>` to swap and
  `git checkout -- <paths>` to restore; or `git stash`. A hand-rolled backup dance can fail its
  first step and still proceed to its second, which is exactly what happened. Earlier A/Bs in this
  same session used `git stash` correctly; this one regressed to `cp`.
- **`set -euo pipefail` and `trap restore EXIT` in any script that mutates the tree.** With the
  trap, the tree returns to HEAD however the script exits — including when it is killed.
- **Commit before an A/B that touches library files.** That is what made this recoverable.

## Operating hazards that have cost real time (all measured)

- **Read signatures and import paths; do not recall them.** Five wrong-from-memory calls in one
  day, one of which (`phi_init` at the top level instead of inside `newton_kw`) killed all five
  legs of a 40-minute run with TypeError.
  ★★ **And a signature is not the contract** (2026-08-16, measured, cost one completed M6 medium
  solve): a dry-check that read `classify_failure`'s argument NAMES still passed lists where it
  needs arrays (it does `tail > 0`; its own call site builds them with `np.asarray`) and unpacked
  TWO returns from a function that returns FOUR. It raised in the reporting layer AFTER the solve
  finished and BEFORE the row was appended, so the solve was lost -- the same family as the
  40-minute solve destroyed by a `float(None)`. ⇒ **the dry-check must exercise RETURN ARITY and
  ARGUMENT TYPES, not just the signature** -- call the function once on a toy input before
  spending compute. Corollary already in force: put `append + write` AHEAD of any post-processing
  and wrap the post-processing, so a reporting error can only add a column.
- **`pgrep`/`pkill -f <pattern>` matches YOUR OWN command line.** `pgrep -f "pytest tests/"`
  reported "still running" for a job that had finished; `pkill -f run_le5_taper_coverage` killed
  the invoking shell (exit 144) and lost the script it was writing. **Kill by PID**, and poll a
  log's completion marker rather than the process table.
  ★★ **"Kill by PID" is not enough if you SEARCH for the PID** (2026-08-11, measured, cost one
  shell): `PID=$(ps -eo pid,args | grep "[r]un_task3_sigma_selection.py" | awk '{print $1}')`
  then `kill $PID` killed the invoking shell — because that shell's own argv contained the
  pattern (the script name appeared elsewhere in the same compound command), and the `[r]`
  bracket trick **only hides the grep process, never the parent shell**. The heredoc that was
  supposed to apply a fix never ran. ⇒ **capture the PID at LAUNCH — `cmd & echo $! > pidfile`**
  — and check liveness with `kill -0 $(cat pidfile)`. Never search for a PID by pattern.
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
3. Full suite: `pytest tests/` — current baseline **569 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-20, **measured in full @805.39 s @8 threads at
   load 17.8**, GS4.1 round 9 = the five missing turbulent terms + the lag
   equation).
   ★ +14 vs the 555 below = `TestFiveMissingTerms` (6) + `TestLagEquation` (8);
   **555 + 6 + 8 = 569 closes exactly, and skipped/xfailed did not move** — neither
   leg touched a solve.
   ★ The 568 measured 40 minutes earlier @687.82 s at load 4.8–13.2 is the SAME
   baseline one lock earlier; it was re-measured rather than incremented by
   arithmetic because the last change edited `lag_rate`'s expression (`s_eq - s`
   to `s_eq - ald*s`, bit-identical at `ald = 1.0` but still a library edit).
   ★★ And the two walls are 687.82 s at load ~5–13 against 805.39 s at load 17.8
   for the same suite on the same box — quote suite walls with their load, never
   as a cost.
   ★★ The five terms are the round's standing lesson: `c_f = max(CFT, CFL)` on
   turbulent stations, `DFAC`'s low-Hk fade, `0.995` (not 1.0) in the outer
   dissipation, the laminar-stress outer term, and the `Us` clamp were all ABSENT
   from `closures_2d.py`, and **not one of them was findable by round 8's G-USED
   guard, which walks the constants already written down** — 0.995 was a number
   never typed, DFAC a function that did not exist, the max a branch. **A guard
   over the constants you wrote cannot see the terms you did not.** Reading the
   source BLOCK whole found them; `A-WHOLE` now does that as a machine check
   (every plain assignment in XFOIL's `BLVAR` classified, unclassified = FAIL).
   ★ Measured consequence, recorded without being used to reject the fix: being
   MORE faithful made agreement with two external ZPG `c_f` correlations WORSE
   (62/69 → 43/70 stations inside), which points at those correlations'
   unestablished applicability rather than at the repair.
   ★ Fast tier at this baseline: **5/5 green, 581 s (9.7 min) @8 threads, load 4.2**.
   ★ The one located code defect still open: just behind transition
   (`x/c` 0.049–0.089) the strip disagrees with XFOIL itself — candidate
   `xblsys.f:1197 TRDIF`, which splits the transition interval into a laminar and
   a turbulent part inside ONE station interval where we switch abruptly at a
   whole station.
   Previous: **555 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-20, GS4.1 round 8 = the GCC transcription fix).
   ★ +2 vs the 553 below = `TestPostTransitionRelaxation`, added because the
   existing 48 locks ALL passed after a change that moved c_f by up to 4.5 % —
   every turbulent lock windows at Re_theta >= 800 and the change concentrates
   below it. **A lock whose window excludes the region a change acts on does not
   cover that change.** The new anchors interpolate H to FIXED Re_theta (600 /
   1000 / 3000) rather than taking a max over stations, because a max over
   stations depends on where the first station lands relative to x_tr.
   Previous: **553 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-19, GS4.1 round 6, **measured @524.70 s @8 threads at
   load 3.11**). ★ +1 vs the 552 below = the E-ATTRACT lock (two seeds 15 % apart in H
   collapse onto one curve) — that is what "equilibrium" means for a turbulent branch,
   and round 5's criterion wrongly tested it as a constant H.
   ★★ This wall also settles the two previous rounds: 1052 s and 875 s were measured at
   load 16.7 and on a busy box, and the same suite runs in ~525 s quiet. Quote suite walls
   with their load; never as a cost.
   Previous: **552 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-19, GS4.1 round 5 = the turbulent closure).
   ★ +10 vs the 542 below = `tests/test_gs41_turbulent_closure.py`. Wall 1052 s @8
   threads at load average **16.7** — the box was genuinely busy in this window, which
   also explains the 875 s reading below; quote suite walls with their load, never as a
   cost.
   ★★ The locks worth knowing: the sourced constants with their xblsys.f/xbl.f citations
   (CTCON is DERIVED from GACON/GBCON, not typed), the identity `DI = 2 c_D / H*` whose
   misreading drove a memory-written version of this closure to an unphysical H = 0.60,
   and that a ZPG turbulent plate's H genuinely DRIFTS ~5 %/decade — the last one is
   recorded because round 5's own T-EQUIL criterion wrongly demanded a constant H, which
   is registered as fix-before-reuse.
   Previous: **542 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-19, GS4.1 round 4 = the closures.py source audit;
   ★ wall 875.53 s @8 threads against 472–486 s for the same suite earlier the same day —
   a 1.8x spread with an identical result. Load average at launch was only 1.26, so the
   cause is NOT established; quote suite walls flagged and never as a cost, per the logged
   precedent of 1.6x/1.9x/5.4x spreads on bit-identical answers).
   ★ +5 vs the 537 below = `tests/test_gs41_closures_audit.py` (one finding, five
   assertions), which locks the audit's
   ROOT CAUSE (the kinetic-energy integrand is degree 21, so 8-point Gauss is short) and
   NOT the error magnitude — so it stays true whether or not the quadrature is fixed.
   ★★ The audit found `closures.py` faithful to Drela AIAA 2013-2437 on every profile and
   integral definition (1e-15), with ONE substantive divergence: `ETA_LAM` is 8 points and
   the kinetic-energy thicknesses need 11, giving phi*_1 a 4.3e-05 quadrature error. Reported,
   NOT fixed — fixing moves every committed Track V number and needs a re-baseline errata
   list. Also recorded: the library's comment justifying 8 points ("degree <= 13") is false
   for those two thicknesses.
   Previous: **537 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-19, **measured in full @485.74 s @8 threads**, GS4.1 round 3).
   ★ +16 vs the 521 below = `tests/test_gs41_closures_2d.py` (route (a2)'s correlation
   closure). **Skipped and xfailed did not move.** ★ Two of the 16 are structural rather
   than numeric: the two closure families must not import each other (`closures.py` owns the
   3-D IBL, `closures_2d.py` owns the 2-D strip — see its docstring), and the profile-path
   fixed point must still read round 1's values, so adding a route cannot move an existing
   number.
   Previous: **521 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-18, **measured in full @472.19 s @8 threads**, GS4.1 round 1).
   ★ +21 vs the 500 below = `tests/test_gs41_strip2d.py` (the 2-D strip core's locks).
   **Skipped and xfailed did not move** — the round touched no solve, and `strip2d.py` is a
   new module rather than a change to one. ★ Two of those 21 assert a **recorded FAIL's own
   numbers** (the closure family's flat-plate fixed point sits +4.52 % in H and +6.94 % in
   c_f√Re_x from Blasius): a recorded FAIL that gets silently re-baselined must be loud,
   and a non-strict xfail cannot do that job.
   Previous: **500 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-17, **measured in full @465.24 s @8 threads**, GS4.0 + the
   R1 addendum). ★ 499 + 1 = the capability-matrix stale-schema lock.
   Previous within GS4.0: **499 + 12 + 2** (@477.27 s @8 threads).
   ★ +20 vs the 479 below = `tests/test_gs40_provenance.py` (the GS4.0 instrument locks:
   the ramp honesty fields, the mesh manifest, the fast tier's node-list check). **Skipped
   and xfailed did not move**, which is the point: GS4.0 was an instrument round and had no
   licence to change a solve. ★ An independent full run of the 479 baseline was made the
   same day BEFORE any change (472.99 s @8t) and reproduced it exactly, so the +20 is a
   clean delta rather than a re-baseline.
   Previous: **479 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-12, **measured in full @494.75 s @8 threads**). It was
   first carried as 474 + 5 by arithmetic and is now measured directly; the 474 before it
   had likewise been 468 + 6. ★ Both arithmetic steps closed exactly against the later
   direct measurement, which is the only reason that bookkeeping is allowed at all — it is
   for non-interacting pure-Python asserts, never for anything that touches a solve.
   ★ Earlier the same day it went 468 → 472 → 468 and that is
   an account closing, not churn: the +4 were the TEMPORARY `sigma_scale` instrument's
   locks, and they were deleted WITH the instrument at its registered expiry
   (docs/dev_phase_three/20260812-0500-sigma-strength-verdict.md §6). A knob kept
   "in case we need it again" is how temporary knobs become permanent.
   ★★ Same disposal applied at the phase-3 close-out (2026-08-16) to `capture_select` /
   `capture_select_abs`, whose route K3 had measured HARMFUL (cl spread 1.21 % → 7.47 % on the
   balanced panel, one seed lost) — it had survived as a default-OFF option with NO registered
   expiry, which is exactly how the previous instance started. **A killed route's knob is deleted
   with the route.** Test count unmoved (it had no locks); the defect it was meant to address is
   still reported, by `sigma_freeze_report` — ★ the distinction worth keeping: **a known defect
   belongs in the report, not in an option the caller is invited to flip.**
   Same count as: **468 passed + 12 skipped +
   2 xfailed, 0 failed** (2026-08-11, measured 494.47 s @8 threads on a quiet box;
   the 930.97 s recorded below was the same 8 threads UNDER LOAD — a 1.9x spread on
   wall time with the same result, so quote suite walls flagged, never as a cost).
   ★ +11 vs the 457 below = `tests/test_meshgen_structured.py` (phase 3 task 3):
   `pyfp3d/meshgen/structured.py` had ZERO tests/ coverage — G0's bit-identical
   single-variable knobs, the one thing route (A) actually delivered, were asserted
   only by a bench script, i.e. on no cadence. Same gap round 2b closed for the
   conforming wing-body transonic lock.
   ★ The drop from 538 is **phase 3 task 1: the level-set route was DELETED**
   (ruling D5) — 9 library files / **4624 lines**, `pyfp3d/wake/` gone entirely,
   `post/unified.py` collapsed onto its conforming half. Read the numbers as an
   account that closes, not as "nothing broke":
   - deleting the 4624 library lines left **passed UNCHANGED at 457** (+1 skipped =
     the new gated wing-body lock), i.e. the deletion subtracted only;
   - ★★★ **GATED full set RE-MEASURED after GS4.0 (2026-08-18): 509 passed + 1 skipped +
     3 xfailed + 1 XPASSED, 0 failed, 1:12:38 @16 threads.** All four numbers close against
     the ungated 500 + 12 + 2: the gated run unlocks **11** skips, of which **9 became
     passed (500 + 9 = 509), 1 became xfailed (2 -> 3) and 1 became XPASSED** -- 9+1+1 = 11.
     ★ And **488 + 21 = 509 exactly**: the delta from the phase-3 handover reading is
     precisely GS4.0's instrument locks, nothing else moved. The XPASS is again the leg the
     mark itself predicts (thread-dependent) => NOT a regression.
   - ★★ **GATED full set RE-MEASURED at the phase-3 close-out (2026-08-16): 488 passed +
     1 skipped + 3 xfailed + 1 XPASSED, 0 failed, 1:11:48 @16 threads.** Read it as an
     account that closes on all four numbers: the ungated suite is 479 + 12 + 2, the gated
     run unlocks **11** skips, and **9 became passed (479 + 9 = 488), 1 became xfailed
     (2 -> 3) and 1 became XPASSED** -- 9 + 1 + 1 = 11 exactly.
     ★★★ And the xpass needs no investigation, because the mark's own reason predicts it:
     `test_p4_transonic::test_g41_transonic_medium_gate` is non-strict ON PURPOSE because the
     outcome is ENVIRONMENT-DEPENDENT -- its text records "at 16 threads this leg CONVERGES
     (|R| 2.8e-13) ... the gated suite at 8 threads then showed it NOT converging (|R|
     3.77e-05)". I ran at 16 threads and the 466-baseline ran at 8, so the ONLY non-pass/fail
     difference from that baseline is the one test whose mark says it flips with thread count.
     ⇒ no regression, and **a gated count must always be quoted with its thread count**.
   - the 2026-08-11 GATED reading it supersedes: **466 passed + 1 skipped + 4 xfailed,
     0 failed** (2:08:50 @8 threads), against 720/2/8 at 3:04:44 before: 208 archived gated items + 47 in
     three archive files that no longer collect + 5 amputated legs − 1 new lock
     = **259**, exactly the difference. The 4 strict xfails that vanished were the
     level-set ones (b7 ×3 + b22 medium), archived with their files.
   - FAST tier is **5 groups** (`bench/run_capability_locks.py`): its two level-set
     locks went with the route, and the new wing-body transonic ceiling lock was
     added — a lock outside this tier would run only in the 2 h gated set.
   ★ `pyfp3d/wake` no longer exists: `from pyfp3d.wake import ...` is a hard error,
   and that is the intended end state. Anything under `phases/p1/` still imports it
   and is NOT runnable — the archive is a historical snapshot; a working tree is
   `git worktree add ../up3d-prereorg d224223`.
4. Numba debugging: `PYFP3D_NOJIT=1` swaps `@njit` for identity — print/pdb work.
5. **When a phase closes — the refresh checklist** (extended in A3 after the
   2026-07-17 audit found 17 consistency defects, most of them close-out debt):
   ★ **Step 0, added 2026-08-09: run the FAST capability-lock tier** —
   `PYFP3D_TRANSONIC_GATES=1 python bench/run_capability_locks.py` (measured 10.7 min,
   7/7 green). The ungated suite does NOT cover the capability anchors; that is how
   seven gated failures accumulated unnoticed, one of them for seventeen days.
   ★ **Step 6, added 2026-08-09: when a criterion is RETIRED or re-specified, grep its
   numbers and wording across `tests/` too** — B28 retired the fuselage-lift premise in
   the demo and the matching assertion sat in a gated test until 2026-08-09.
   ★ **Phase THREE (opened 2026-08-11) records rounds in `docs/dev_phase_three/`**
   (its own progress.md); the PLAN and rulings D1–D5 still live in
   docs/dev_phase_two/roadmap.md, and docs/dev_phase_two/README.md holds the per-file
   checklist for when those six files may move to `phases/p2/docs/`.
   ★★★ **Phase THREE CLOSED 2026-08-16** (ruling D7; all three tasks on `main` — ① via PR #26,
   ②③ via PR #27 `d5efe53`). **PHASE 4'S SINGLE ENTRY POINT IS
   [docs/dev_phase_three/20260816-1000-gs41-initiation.md](docs/dev_phase_three/20260816-1000-gs41-initiation.md)**
   — read its §0 before anything else: it is self-contained by design (handover baselines, the items
   D7 assigned forward, the six things that must not be silently changed, and the two measured
   conclusions GS4.1 must carry — **P1 = 6.07× now but crosses the band's floor at the next level**,
   and **the crossflow a 2-D strip discards is a REAL trailing-edge structure**, J = 10× the null).
   ★ It stays where it is; phase 4's own progress.md **links** it rather than copying it — a copied
   entry point forks from its original, which is the failure this project has logged repeatedly.
   ★ **Phase TWO uses a different surface list**, because the phase-one docs below are
   frozen (docs/dev_phase_two/roadmap.md §8): (1) `docs/dev_phase_two/progress.md` — one
   row, plus **its own 阶段进度概览 and 产品指标追踪 tables** (they live in progress.md,
   NOT in the roadmap); (2) `docs/dev_phase_two/roadmap.md` when a ruling or a stage
   disposition moves; (3) the round file
   itself; (4) this file's **baseline line** when the suite count moves; (5)
   `docs/dev_phase_two/PHASE_TWO_CAPABILITY_BOUNDARY.md` when a measurement moves a
   capability claim; ★ **(6) PROJECT_STRUCTURE.md** — added 2026-08-10 after checking it
   and finding phase two absent from it ENTIRELY: neither `bench/` nor `docs/` was in the
   directory tree and the footer baseline still read the phase-one 652, so ~70 round files
   and every committed CSV were invisible to anyone navigating by that document. It is the
   surface the phase-one ritual already flags as the one that silently rots, and leaving it
   off the phase-two list is how it rotted again. Steps 0 and 6 and the erratum checklist
   apply to both lists.
   For a **phase-one** close-out:
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
   5. the **`cases/demo/README.md` table row / `bench/studies/README.md` bullet**
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
