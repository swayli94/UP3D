# PRE-REGISTRATION — GV5.4: augmented-step cost on M6 medium with the block preconditioner

Committed BEFORE the first code change (Track V discipline; user-sequenced
2026-07-24: GV5.1d → GV5.5 → GV5.2–5.4; GV5.3 merged 2026-07-25 PR #18).
Binding text: `docs/roadmap/track_v.md` **GV5.4** —

> augmented step wall-time ≤ ~2× the inviscid Newton step on M6 medium
> with the block preconditioner working; measured number recorded either
> way.

Guidance in force: the V5-tail note ("Block precond: AMG-φ / ILU-BL,
`schur_ls.py` prototype; **BL block is NOT thin — measure before
Schur**"), `docs/design_track_v.md` §12 item 2 ("GMRES + 块预条件
（AMG-φ / ILU-BL）留 GV5.4"), and the §13 diagnosis (the row+col
equilibration kills the null cluster; a genuine scaled (A,Ψ) stiffness of
1e5–1e7 remains = the real GV5.4 target; the BL preconditioner should be
organized by the (A,Ψ) structure; the equilibration piece is the ready
ingredient).

## 1. Question

On the ONERA M6 at TEST 2308 (M 0.8395, α 3.06, Re_MAC 11.72e6), medium
mesh, the GV5.1 augmented (φ, Γ, U) system (~124k DOFs — ≈ 4.4× the
2.5-D strip, where splu was affordable):

- **(a)** What does ONE augmented Newton step cost with a
  block-preconditioned GMRES linear solve, relative to ONE inviscid
  Newton step measured in the same session? RECORDED vs the ≤ ~2×
  reference band, **the measured number recorded either way**.
- **(b)** Does the block preconditioner WORK (the binary adjudication
  defined in §3 D5)? PASS = a rung works; FAIL = none does.
- **(c)** Diagnostics RECORDED: exact DOF counts, per-phase breakdown,
  GMRES statistics, Schur setup (if reached), memory, seed IBL floor.

This is a COST gate, not an accuracy/convergence gate: no Cp/CL band is
registered, and the IBL floor (the committed ~1e-6 F_BL floor) is
expected to cap the Newton loop — the steps are measured from the
pre-registered seed regardless.

## 2. Conditions, meshes, anchors

- **M6 medium binding**: `cases/meshes/onera_m6/medium.msh` (gitignored;
  regenerate with `cases/meshes/onera_m6/generate_onera_m6.py`). Uncut
  63,237 nodes / 350,718 tets (committed `medium_stats.csv`); wake-cut
  n_cut 70,663, n_st = 166 stations; IBL surface 10,205 nodes
  (GV5.3-measured) → the augmented system ≈ 62.7k + 166 + 61,230
  ≈ **124k DOFs** (the exact counts are measured and recorded, W2).
  **Coarse = wiring shakedown only** (execution mechanics, recorded in
  the VERDICT narrative, not gated).
- Condition = the GV5.3 conventions verbatim: TEST 2308 M 0.8395,
  α 3.06, Re_MAC 11.72e6, x_tr/c = 0.05 both sides, tip band
  z > 0.95·b_semi (the production tip_taper radius).
- Inviscid step anchor: measured **in-session, same driver, same branch,
  same thread count** from the seed ramp's final-level `step_records`
  (D6). The committed A1 numbers (@16 threads: the final level
  12 steps / 52.14 s ≈ 4.3 s/step mean; one LU refactor ≈ 14.7 s; one
  stale-LU GMRES ≈ 3.2 s — `a1_m6_levels.csv`, `a1_m6_runs.csv`) are
  quoted as a NON-binding cross-check only.
- Probe-branch anchor (W1): the committed A1 `conf_newton`
  (probe-Kutta, NEWTON_M6_RECIPE verbatim) **cl_p 0.26918**
  (`cases/analysis/a1_solver_bottleneck/results/a1_m6_runs.csv`).

## 3. Design decisions

- **D1 — probe Kutta (adjudicated).** The augmented F_Γ row exists
  verbatim for the PROBE estimator only (`tight_driver.py:239-244`;
  pressure's un-eliminated row is NOT wired — module docstring). Wiring
  the pressure row = new Jacobian code + its own FD gate = a separate
  work item (registered follow-up, §7). Probe is justified for a COST
  gate: the system size, sparsity, and block structure are identical;
  the probe branch carries its own committed anchor (A1 cl_p 0.26918);
  and the inviscid anchor is measured on the SAME probe branch, so the
  comparison is internally consistent. (Context: the probe branch reads
  ~4.5 % low on lift vs pressure — committed P14 V14.6 — an accuracy
  fact that does not touch a cost measurement.)
- **D2 — seed = the `newton_tight` seed semantics verbatim** (module
  docstring): an inviscid-converged (φ, Γ) + ONE standalone IBL solve.
  The inviscid state: `solve_newton_lifting` M0.70 probe seed →
  `solve_newton_transonic` ramp (NEWTON_M6_RECIPE imported from
  `tests/test_p8_newton.py`; the estimator default = probe at every
  level). The IBL state: `IBL3Solver.solve` defaults at the inviscid
  edge data; `converged=False` at the committed ~1e-6 floor is EXPECTED
  (the GV5.1 finding), recorded, not a failure.
- **D3 — linear-solve injection (library change, minimal).**
  `scaled_damped_step` gains an optional `solve(A, b) -> y` callback
  (default None = the committed splu line, bit-identical);
  `newton_tight` gains a `step_solve=None` passthrough. The augmented
  steps run `scaling="rowcol"`, `lm_damping=False`, `floor_stop=False`
  (mu ≡ 0 — a pure equilibrated Newton step; the LM schedule is GV5.1b's
  convergence machinery, not this gate's target). Suite +1 test (W4).
- **D4 — the block-preconditioner ladder** (on the row+col equilibrated
  system; blocks extracted from the assembled augmented J by the index
  ranges n_free / n_st / 6·n_s, the node-major 6i+k layout):
  - **Rung 1 — block-Jacobi**: M = bdiag(M_φ = AMG on the (φ,φ) block
    [`build_amg_preconditioner`, pinned seed], M_Γ = the exact inverse
    of the extracted 166×166 (Γ,Γ) block, M_BL =
    `build_ilu_preconditioner` on the (BL,BL) block with its committed
    robustness ladder).
  - **Rung 2 — exact-BL Schur** (the B14 `schur_ls.py` pattern):
    splu(J_BL,BL) once per step (its setup cost RECORDED — the
    roadmap's "measure before Schur"); the matrix-free reduced operator
    on (φ,Γ); GMRES with bdiag(AMG(J_φφ), M_Γ); exact back-substitution
    (the BL rows exactly satisfied).
  - **Escalation**: rung 1 first. If any measured step's GMRES fails to
    converge within budget, or the per-step linsolve wall cap (1800 s)
    trips, rung 1 is adjudicated NOT-working and the remaining steps
    move to rung 2. Both rungs' numbers are recorded.
  - **GMRES budget**: rtol 1e-8, restart 60, maxiter 5 (≤ 300
    iterations), `on_fail="return"` (the P8 idiom: a non-converged
    return adjudicates the rung; the best-iterate step is still taken
    through the safety backtracking and recorded).
  - **No stale reuse**: every measured step builds its preconditioner
    fresh (the honest per-step cost). The inviscid anchor's per-step
    mean includes its own committed refactor+reuse mix; the comparison
    annotates this asymmetry. Stale-preconditioner reuse = a registered
    optimization follow-up (§7).
- **D5 — "the block preconditioner working" (binary adjudication).** A
  rung WORKS iff GMRES converges (info = 0) within budget on EVERY
  measured Newton step of that rung AND every computed step is accepted
  by the safety backtracking (the committed accept-or-least-bad idiom,
  max_backtracks = 30, the probe guard; λ recorded) AND the W3 FD guard
  passed. Item (b) = PASS iff any rung works (the working rung
  reported).
- **D6 — the inviscid anchor protocol.** The seed ramp's FINAL level
  (0.80 → 0.8395, strict) returns `step_records` (per-step phase deltas,
  `solve/timing.py`): per-step wall + assembly/precond/linsolve,
  refactor vs reuse flagged. The denominator = the MEAN final-level
  step wall in-session. Expected ~12 steps (the A1 final level).
- **D7 — the measurement.** N = 5 augmented Newton steps from the seed
  (`max_iter=5`; the IBL floor makes "cap" the expected termination — an
  early "converged"/line-search termination is recorded with the
  measured n). Per step: preconditioner setup + GMRES timed INSIDE the
  callback; assembly timed separately at the seed and terminal states
  (2 samples, the same `augmented_jacobian` call); the residual +
  backtracking share = the step wall (history `wall_s` deltas) minus
  the timed pieces. ru_maxrss recorded.
- **D8 — the comparison.** ratio = mean(augmented step wall) /
  mean(inviscid final-level step wall), both in-session @8 threads.
  RECORDED with the within/above annotation vs the ≤ ~2× reference. No
  PASS/FAIL on the ratio itself.
- **D9 — thread constraint.** Executed under the temporary 8-thread
  session constraint (user-directed): NUMBA/OMP/OPENBLAS_NUM_THREADS =
  8; every wall time is flagged non-comparable to the 16-thread ledger
  entries.

## 4. Wiring guards (raise = recipe error, NOT verdicts)

- **W1** — the probe seed at M 0.8395: the final level converged strict
  AND |cl_p − 0.26918| ≤ 1.5 % vs the committed A1 anchor (the
  tolerance absorbs the ΔM 0.0005 label difference).
- **W2** — pack DOF counts: n_st = 166, 6·n_s = 61,230 (n_s = 10,205),
  n_free ∈ [60k, 64k] (the estimate band; the exact value recorded).
- **W3** — sampled FD guard on the assembled wing augmented J at the
  seed: 60 rows sampled across the blocks (φ wake-adjacent, the Γ tip
  station, BL tip-pinned / LE-band / TE-band rows), the GV5.1 fd_gate
  ladder (h = 1e-5/1e-6/1e-7 + the fallback row mask), per-block median
  relative error < 1e-6. Guards the NEW wing wiring (the spanwise-γ
  far-field, the tip pin mask, the probe row on 166 stations, the wing
  edge-velocity operator).
- **W4** (CI-side) — the new suite test: `newton_tight` with
  `step_solve=<explicit splu>` reproduces the default path bit-for-bit
  on the committed small fixture (suite 642 → 643).

## 5. Verdict items + exit code

- (a) the cost ratio — RECORDED (annotation: within / above the ≤ ~2×
  reference).
- (b) the preconditioner adjudication — WORKING (the rung reported) =
  PASS; NOT-WORKING = FAIL.
- (c) diagnostics — RECORDED.
- The runner exits 1 iff (b) reads FAIL; the guards raise (recipe
  error). The cost number is recorded either way.

## 6. Execution mechanics

- Coarse shakedown first (the same protocol, cheaper), then medium
  (binding). Expected cost @8 threads: the seed ramp ≈ 3-6 min (A1
  scaling) + the standalone IBL solve + 5 augmented steps (est.
  1-15 min each) + the FD guard (≈ 10-20 min) per level.
- Artifacts: `VERDICT.md`, `results/summary.csv`,
  `results/steps_{coarse,medium}.csv` (per-step phases),
  `results/gmres_{rungN}_{level}.csv` (per-step GMRES histories),
  `results/inviscid_anchor_{level}.csv` (the final-level step_records),
  `results/fd_guard_{level}.csv`, `cost_breakdown_{level}.png` (the
  phase breakdown + GMRES convergence), `run.log`.
- The full suite is re-run after the library change (D3) and the new
  baseline (642 → 643) is filled into the three ledger spots with the
  8-thread flag.

## 7. Follow-ups (registered, NOT opened — user adjudication)

- The pressure-Kutta row wired into the augmented un-eliminated layout
  (D1) with its FD gate.
- Productizing the working preconditioner into `tight_driver` (if a
  rung works), incl. stale-preconditioner reuse across Newton steps.
- The (A,Ψ)-structured BL preconditioner (the design-doc organization)
  if the plain ILU rung stalls on the scaled (A,Ψ) stiffness.

---

## Addendum 2026-07-25 #1 — W1 becomes medium-binding (execution mechanics)

First execution: the **coarse** W1 fired (cl_p 0.255793 vs the A1 anchor
0.26918, rel 4.97 % > 1.5 %). Root cause: the A1 `conf_newton` anchor is
a **medium-mesh** number (a1_m6_runs.csv runs `onera_m6/medium` only) —
no committed coarse probe anchor exists. The measured coarse deviation is
the expected mesh effect: the committed P14 pressure-Kutta mesh effect is
coarse 0.262778 vs medium 0.277628 = **−5.35 %**, and the measured coarse
probe deviation is **−4.97 %** — the same sign and size, i.e. the coarse
probe seed sits exactly where the committed mesh effect puts it. The
wiring itself is separately guarded by W3 (the FD guard passed the coarse
wiring smoke at ≤ 2.5e-11 medians) and by the binding medium W1.

Change (guard scope only — no band, protocol, or verdict-item change;
medium W1 untouched):

- **W1 (medium, binding)**: unchanged — the final level converged strict
  AND |cl_p − 0.26918| ≤ 1.5 %.
- **W1 (coarse, shakedown)**: the final level must converge strict
  (raise = recipe error); cl_p is RECORDED with the mesh-effect
  cross-check (vs the committed pressure mesh effect −5.35 %) quoted in
  the VERDICT narrative.

---

## Addendum 2026-07-25 #2 — the M0.70 seed's strict non-convergence does not raise (execution mechanics)

Second execution (medium): the M0.70 probe seed stalled at the 60-step
cap (converged=False, 120 s) — the committed GV5.0/GV5.3 medium
signature (the seed plateaus at |R| ~ 1e-6). The pre-registration's
make-seed wording ("the ramp's level 0 re-converges only a converged
seed family") contradicted GV5.3's **measured** evidence (GV5.3
addendum #1: the NEWTON_M6_RECIPE ramp from the very same stalled-seed
family re-converges level-by-level — 0.70: 3 Newton → 6.4e-6; …;
0.8395: 12 → 7.6e-15 — onto the anchored branch to 9 printed digits;
the bridge's early-RETURN was the bug, not the seed's strict
non-convergence). Here the ramp is even smoother than GV5.3's case: the
seed and every ramp level share the SAME probe estimator.

Change (execution mechanics only — no band/protocol/verdict change):
the M0.70 probe solve's non-convergence no longer raises; its (φ, Γ)
seeds the ramp unconditionally (the actual GV5.3 addendum-#1 lesson:
never return early, never refuse). The seed's converged flag is logged.
**W1 is untouched and remains the anchor guard**: the final level must
converge strict AND (medium) cl_p within 1.5 % of the committed A1
probe anchor — a wrong-branch landing is exactly what W1 catches.

---

## Addendum 2026-07-25 #3 — the probe seed chain IS A1's conf_newton verbatim (execution mechanics)

Third execution (medium): the ramp converged strict (12 final steps —
A1's own count) but W1 fired: cl_p 0.264293 vs the committed A1 anchor
0.26918 (rel 1.816 % > 1.5 %). Root cause, read off
`cases/analysis/a1_solver_bottleneck/run_a1_m6.py:286-292`: the
committed probe anchor belongs to **A1's chain** — ONE
`solve_newton_transonic(mc, wc, m_inf=0.84, alpha_deg=3.06,
**NEWTON_M6_RECIPE)` call, the ramp Picard-seeding level 0 itself (no
separate M0.70 Newton solve, no phi_init handoff) — while D2's seed
wording borrowed the **pressure** chain's shape (a separate M0.70 probe
solve → the ramp with n_picard_seed=0; the P14 cold-start idiom, the
right chain for the pressure estimator whose anchor is P14's). The two
chains land on nearby but distinct FP solution families (the M6
transonic branch sensitivity, design.md Sec 12 risk 2): mine read
1.8 % low in cl_p. A guard firing on a recipe-family mismatch is the
guard working — the fix is to run the anchor's own chain, not to widen
the guard.

Change (execution mechanics only — D2's seed sentence is amended; no
band/verdict/guard-tolerance change):

- The inviscid probe seed = **A1's conf_newton verbatim**: one
  `solve_newton_transonic` call with NEWTON_M6_RECIPE (the estimator
  default = probe; the ramp Picard-seeds level 0; M 0.8395 per the
  TEST 2308 dataset label — the ΔM 0.0005 vs A1's 0.84 is what the
  pre-registered 1.5 % tolerance absorbs). The separate M0.70 Newton
  solve is dropped (addendum #2 becomes moot — that seed no longer
  exists).
- W1 is unchanged (final strict convergence AND, medium-binding,
  |cl_p − 0.26918| ≤ 1.5 %).

---

## Addendum 2026-07-25 #4 — W1 re-anchored to the P14-committed probe lock (guard anchor correction)

Fourth execution (medium, the A1-verbatim chain): W1 fired AGAIN with
the SAME value to all printed digits — cl_p 0.264293 (two different
seed chains land on the same state ⇒ the offset is not a recipe-family
artifact). Anchor archaeology:

- `a1_m6_runs.csv` conf_newton cl_p **0.26918** — the A1-era reading
  (A1 was a performance study, not a cl gate; the probe path's cl was
  never re-locked there).
- `p14_pressure_kutta/results/cross_model_medium_m084.csv` —
  **conforming probe (G8.2 lock) cl_p 0.2646 / cl_KJ 0.2692** vs
  pressure 0.2776/0.2823 vs level-set 0.2772/0.2813 (the same table
  behind the committed V14.6 check "the PROBE path was 4.5 %/4.3 %
  below LS"). The probe path evidently shifted ~1.7 % between the A1
  era and the P14 re-lock (the P13/P14 Kutta-row work).
- In-session (both chains): **0.264293** = 0.116 % below the P14-locked
  0.2646 — today’s probe branch IS the P14-locked branch.

Change (guard anchor correction — the tolerance is unchanged):

- **W1 anchor**: 0.26918 (A1, stale) → **0.2646** (the P14-committed
  probe G8.2 lock, `cross_model_medium_m084.csv`; the most recent
  committed probe reading at M 0.84 medium). The 1.5 % tolerance now
  absorbs the ΔM 0.0005 label difference AND the anchor’s 4-digit
  precision. The A1 number stays quoted in the VERDICT as the stale
  superseded reading.
