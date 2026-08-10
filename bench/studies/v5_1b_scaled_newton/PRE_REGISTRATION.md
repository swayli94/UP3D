# PRE-REGISTRATION — GV5.1b: scaled + damped augmented Newton, the pre-floor quadratic window

Committed BEFORE the first execution, per discipline. Gate:
`docs/roadmap/track_v.md` V5 ("next = GV5.1b design" after the
IBL-floor diagnosis). Design inputs: the committed diagnosis
`../v5_ibl_floor/results/findings.md` (= `docs/design_track_v.md`
§13): (i) the raw cond 4e10–4e13 of J_BL,BL is mostly a scaling
artifact — one row+column 2-norm equilibration leaves cond
2e4/7e5/1e7 with 0/0/2 sub-1e-6 singular values; (ii) a genuine
turbulent (A, Ψ) stiffness of 1e5–1e7 remains after scaling;
(iii) the floor residual lives in J's range, so the quadratic
window exists BEFORE the floor, not beyond it; (iv) globalization
alone cannot pass the floor (Q7: the pseudo-time controller bottoms
out) — damping must act on the Jacobian, not only on the step
length.

## Architecture (solver-internal only; assembly untouched)

`newton_tight` gains a scaled + damped linear step. The assembled
residual F and Jacobian J are **bit-identical to GV5.1** — the
GV5.1 FD verdicts (2e-8..5e-9) stand; band (a) asserts this by
regression. Per Newton iteration:

1. Equilibration: diagonal R, C from the current full augmented J
   (row / column 2-norms, zero-safe: norm 0 → 1). Scaled system
   J̃ = R·J·C, F̃ = R·F; solve for δy, unscale δx = C·δy.
2. Levenberg-style damping in the scaled space:
   (J̃ + μI) δy = −F̃ via splu (sparsity preserved, no normal
   equations). Deterministic schedule: μ₀ = 1e-6; on a line-search
   rejection μ ← min(10·μ, 1e2) and the step is retried; on an
   accepted step μ ← max(μ/3, 1e-12). The P8/P14 backtracking with
   the committed probe guard applies unchanged on each damped
   trial step.
3. Termination: GV5.1's converged criterion kept; ADD the
   pre-registered **floor-reached** stop = merit relative decrease
   < 1e-4 over 3 consecutive accepted steps → terminate and read
   the bands on the history (replaces GV5.1's λ-collapse crawl to
   the iteration cap).

## Gate bands

- **(a) implementation exactness (PASS/FAIL).** Algebraic identity
  tests at machine precision: R·J·C equals explicit diagonal
  scaling; the μ = 0 damped step equals the undamped step;
  δx = C·δy round-trip; μ-schedule transitions as specified. The
  committed tight regression (20 tests) stays green; the FD
  assembly verdicts are unchanged (assembly untouched — recorded).
  New tests `tests/test_v5_tight_scaled.py` (~6); the k=1 smoke
  follows the NOJIT-skip precedent (`test_v3_coupling.py:204`).
- **(b) pre-floor quadratic window (medium binding, coarse
  recorded).** Protocol = the GV5.1 amended protocol verbatim
  (seed = HEAD-regenerated loose-converged state, committed recipe;
  wiring guard |dcl_k0| ≤ 1e-8 per leg; augmented Newton as
  polish). The **floor band** is pre-registered at 10× the
  diagnosed standalone floors: coarse 3.16e-5, medium 1.71e-5.
  PASS = ≥ 3 consecutive contractions with median p ∈ [1.5, 2.5],
  p computed on the F_BL max-norm sequence (the GV5.1 read),
  measured on iterates with F_BL above the floor band; if the run
  converges outright the band is trivially PASS. Merit-sequence p
  recorded alongside. Honestly recorded: termination class
  (converged / floor-reached / cap), final F_BL vs the diagnosed
  floor, final merit vs the GV5.1 committed final merits (read
  from `../v5_tight_coupling/results/`, never recomputed).
- **(c) counts (recorded; aspirational band).** N_polish to enter
  the floor band (or converge) and N_total = N_loose + N_polish vs
  the committed loose outer counts (4 coarse / 5 medium).
  Aspirational PASS: N_polish ≤ 2× loose; 3+ recorded with user
  adjudication (the GV5.1 honesty note carried over). Standalone
  k=1-seed retry: recorded only, no band.

## Fallbacks and aborts (pre-registered)

- If no slope-2 window appears even before the floor: record the
  full merit/F_BL/μ/λ history; the partial-success read is whether
  the scaled+damped run descends below the GV5.1 committed final
  merit and whether termination is floor-reached rather than
  λ-collapse. Both outcomes are evidence for the next formulation
  step (TE-band (B, δ) truncation), not a gate crash.
- If the wiring guard fails on either leg: abort that leg, record,
  user adjudication (the GV5.1 pattern).
- V4-reopen trigger: not invoked by this gate's failure (the
  diagnosis already localizes the stall inside the BL block).

## Out of scope

Closure/formulation edits; TE-band discretization changes; GMRES
and block preconditioners (GV5.4); loose-loop criterion changes;
3-D M6 tight solves (the 2.5-D strip is the GV5.1x testbed).

## Artifacts

`run.py` (regenerates everything; seeds by the committed recipe,
loose numbers read from `../v3_loose_coupling/results/`, GV5.1
histories read from `../v5_tight_coupling/results/`) ·
`results/newton_history_{coarse,medium}.csv` (F blocks, merit, μ,
λ, p, termination class) · `results/compare.csv` (vs GV5.1 and
loose committed) · `results/summary.csv` (one row per band) ·
VERDICT at wrap-up.
