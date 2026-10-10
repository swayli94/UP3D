# GV6.1 VERDICT — conforming wake-sheet δ* source

- Date: 2026-07-25 · Branch: `kimi/track-v6-gv6-1`
- Binding text: `phases/p1/docs/roadmap/track_v.md` GV6.1
- Pre-registration: [PRE_REGISTRATION.md](PRE_REGISTRATION.md) **including
  the 2026-07-25 addendum** (the per-face ½ṁ_wake recipe; committed BEFORE
  the first code change, commit `f4e90de`, per the §7 amendment
  discipline — a design-review catch, not an execution patch)
- Evidence: `results/summary.csv`, `results/mms_probes.csv`,
  `results/smoke_history_{on,off}.csv`, `results/wake_field.csv`,
  `results/te_cp.csv`, `results/gv6_1_smoke.png`,
  `results/ab_cache_mode_isolation.csv` — regenerate with
  `python bench/studies/v6_1_wake_sheet/run.py` (exit 1 on any honest
  FAIL)
- **VERDICT: GV6.1 PASS — 6 PASS / 0 FAIL / 7 RECORDED.** Every binding
  band ((a)(i)/(a)(ii) bit-identity, (b) sign-pin MMS, (c) TE-continuity)
  PASSED on first execution of the final harness. Test lane:
  `tests/test_v6_wake_sheet.py` 7/7 (JIT) + 4 pass / 3 skip (NOJIT);
  regression `test_v0_freestream` + `test_v2_*` + `test_v3_coupling`
  30/30. The default-OFF flag is bit-inert; V6's conforming sheet-source
  channel is closable, the measured-effect question moves to GV6.2.

## Result table

| gate | metric (band) | measured | verdict |
|---|---|---|---|
| (a)(i) | flag-ON prescribed-zero δ*_wake vs flag-OFF, coarse 3 outer (bit-identical) | phi_equal=True gamma_equal=True | PASS |
| (a)(ii) | flag-OFF vs gate-free library @13916b5, fresh-compile worktree legs (bit-identical) | phi_equal=True gamma_equal=True | PASS |
| (b) | ejects-away sign at every probe (v+ > 0, v− < 0) | sign_ok=True | PASS |
| (b) | antisymmetry max\|v+ + v−\|/m0 (≤ 5%) | 0.0081 | PASS |
| (b) | jump max\|(v+ − v−) − m0\|/m0 ≤ 5% — also pins the addendum's per-face ½ factor | 0.0044 | PASS |
| (c) | W2 TE-continuity identity (1e-12 rel) held at every outer | 3 outers clean | PASS |
| (d) | cl flag OFF / ON (final) | 0.28246 / 0.28261 | RECORDED |
| (d) | on/off Δ-cl (in-session A/B; scatter caveat) | +0.00015 | RECORDED |
| (d) | outer count ON vs OFF | 3 vs 3 | RECORDED |
| (d) | ṁ_wake scale max over outers | 1.800e-02 | RECORDED |
| (d) | δ*_wake,TE / θ_wake,TE (final) | 6.61795e-03 / 4.47577e-03 | RECORDED |
| (d) | TE-region (x/c > 0.9) max \|ΔCp\| on/off | 0.00250 | RECORDED |

## Per-gate analysis

### (a) δ*_wake = 0 bit-identity — both legs

(a)(i): the flag-ON loop with a prescribed ZERO δ*_wake field assembles the
exact zero RHS (m_wake ≡ 0 → b_wake ≡ 0 by the GV2.1(b) discipline) and
reproduces the flag-OFF loop bit-for-bit, in-process, coarse 3 outer.

(a)(ii): the flag-OFF loop of THIS tree (HEAD + working-tree delta
overlaid) is bit-identical to the gate-free library at the pinned baseline
`13916b5`. **Both legs run as subprocesses on FRESH worktrees** (fresh
numba compile) — see the harness-fix section below for why this
discipline is load-bearing.

### (b) sign-pin MMS — the recipe's sign AND the addendum's ½ factor

Dead-air coarse strip (m_inf = 0, the manually-driven Laplace reduction
through the production Tᵀ route — at m_inf = 0 the Picard driver IS this
incompressible solve, the G3.3 equivalence), uniform ṁ0 = 0.01 through the
PRODUCTION assembly (`assemble_wake_sheet_rhs`, ½ṁ per coincident face
copy). The probes over the strip's middle third: v+ > 0 and v− < 0 at
EVERY probe (a thickening wake ejects fluid away from the sheet on both
sides — the V2 blowing-positive convention), antisymmetry
max|v+ + v−|/m0 = 0.0081, and the jump max|(v+ − v−) − m0|/m0 = 0.0044.
**The jump clause empirically pins the 2026-07-25 addendum's per-face ½
factor**: without it the Tᵀ fold would realize [ρv_n] = 2ṁ0 and the jump
would measure ≈ 2m0 — a 100% error, two decades outside the 5% lock.

### (c) TE-continuity (W2) — runtime-asserted, every outer

The producer's construction identity (δ*_wake(0) = δ*_upper(TE) +
δ*_lower(TE), 1e-12 rel) is asserted inside every
`wake_transpiration_source` call; it never fired over the band-(d) ON
run's 3 outer iterations. TE-copy wiring is correct by construction check,
not by inspection.

### (d) GV3.1 smoke — RECORDED (non-binding), the GV6.2 input

Medium mesh, the committed GV3.1 recipe verbatim (M 0.5 / α 2° / Re 3e6 /
xtr 0.05, loose Picard leg), flag ON vs OFF in-session on the same seed:

- Δ-cl = **+0.00015** on cl ≈ 0.2825 (≈ +0.05%), both legs converged in
  3 outers (tol_ds 1e-3).
- TE-region max |ΔCp| = **0.00250** (x/c > 0.9; the ON leg's Cp sits
  slightly higher over most of the region, see `te_cp.csv` and the right
  panel of `gv6_1_smoke.png`).
- Producer state: δ*_TE = 6.62e-3 c, θ_TE = 4.48e-3 c (H_TE ≈ 1.48),
  ṁ_wake max 1.8e-2 concentrated at the TE and decaying over ~2 c
  (L_rel = 1.0 c pinned MODEL CHOICE; δ*_wake relaxes δ*_TE → θ_TE, the
  H → 1 far-wake limit, middle panel).

These numbers are the input to GV6.2's significance question (the on/off
effect vs the A4 input band), which GV6.1 deliberately does NOT judge.
Caveat carried from GV3.1: the committed medium fixed point is not
cross-run reproducible (the v5_tight_coupling scatter caveat); only the
in-session ON/OFF pairing is meaningful.

## Implementation fixes made during gate execution

No solver-code or recipe change was made after execution began. Two
HARNESS fixes were required (bands unchanged; evidence
`results/ab_cache_mode_isolation.csv`):

1. **(a)(ii) cache-mode discipline (the substantive one).** The first
   (a)(ii) attempt compared the in-process flag-OFF leg against a
   fresh-worktree baseline leg and FAILED at max|Δφ| = 4.2e-5 — despite
   the flag-OFF code path being statement-identical to the baseline
   (coupling.py/wake_sheet.py contain zero numba and are bit-inert under
   the flag). A five-leg bisection (isolate3) proved the cause is NOT
   the new code: with fresh caches everywhere, baseline ≡ HEAD ≡
   HEAD+our-code EXACTLY (D24 = D15 = 0), while the SAME baseline sources
   diverge at 6.5e-6–8.3e-6 in φ between a fresh-compile run and a
   cache-load run (D12), entering at outer k ≥ 1 (k = 0 inviscid exact).
   A follow-up leg (isolate4: warm tree, ONLY `pyfp3d/viscous/__pycache__`
   wiped) recovered the all-fresh result EXACTLY — **numba cache-load is
   not bit-faithful to fresh-compile, and the infidelity lives entirely
   in `pyfp3d/viscous/`** (cache-hit code for some viscous-chain njit
   function(s) — closures/ibl3 family with many cache generations —
   produces different FP results than compiling the same source fresh;
   numba 0.62.1 / py3.13.9). Cache-load IS deterministic run-to-run
   (D23 = 0). The (a)(ii) harness — in the test and in this runner —
   therefore runs BOTH legs as fresh-compile worktree subprocesses; the
   band itself (bit-identity) was NOT relaxed. Consequence for the repo:
   cross-mode absolute numbers can differ at the ~1e-5-in-φ level, so any
   future bit-identity A/B must pin one cache mode on both legs, and
   cross-run comparisons of committed viscous-run numbers already carry
   this floor (consistent with the v5_tight_coupling scatter caveat).
2. **Test-helper parsing fix (mechanical).** `git status --porcelain`
   output was read through a `.strip()`ing helper, which ate the leading
   status column of the first line (" M pyfp3d/…" → path mangled); the
   overlay now reads stdout unstripped. Test infrastructure only.

## Numerical settings (as run)

- Thread cap **8** (NUMBA/OMP/OPENBLAS; temporary user-directed session
  constraint, PRE_REGISTRATION §7 — wall times non-comparable).
- Meshes: committed `cases/meshes/naca0012_2.5d/coarse.msh` ((a)(b);
  496 wake nodes / 248 stations / s_max 14.5 c) and `medium.msh` ((c)(d);
  986 / 493 / 14.5 c).
- (a): coarse 3 outer, Picard lifting driver, Re 3e6 / M 0.5 / α 2°.
- (b): dead-air Laplace through `PicardOperator` + `WakeConstraint` +
  `farfield_dirichlet`(u_inf = 0) + `reduced_rhs`, CG rtol 1e-11 with the
  AMG preconditioner; probes in the strip's middle third.
- (d): CouplingConfig defaults (n_outer_max 10, tol_ds 1e-3, omega 1.0,
  xtr 0.05); both legs converged at k = 3.
- Producer constants: L_rel = 1.0 c (pinned MODEL CHOICE, RECORDED — the
  sensitivity sweep is a GV6.2 item, not tuned here); W2 rtol 1e-12.
