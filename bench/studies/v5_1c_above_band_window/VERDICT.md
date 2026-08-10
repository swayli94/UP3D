# GV5.1c VERDICT — the above-band seed: the pre-floor slope-2 window measured — NO quadratic regime above the floor

Gate: `docs/roadmap/track_v.md` V5 follow-up **GV5.1c** (2026-07-24,
user-directed). Pre-registration: `PRE_REGISTRATION.md` (committed
1e90d59 BEFORE the first execution). Runner: `run.py` (regenerates
every artifact in `results/`; the GV5.1b runner's seed/protocol
functions imported, not mirrored; wiring guard PASS both legs —
coarse 1.56e-12, medium 1.54e-9). Executed 2026-07-24 on
`kimi/track-v5-gv5-1c` under the **temporary 8-thread constraint**
(this session only, user-directed; runner default 16 — wall times NOT
comparable to the 16-thread ledger entries). Row-level record:
`results/summary.csv` — **2 PASS / 1 FAIL / 7 RECORDED**.

**Answer to the gate question, up front:** the pre-floor slope-2
window is now **measured, not just reframed — and the answer is NO**.
With genuinely above-band seeds delivered as pre-registered (F_BL
3.219e-1 coarse / 1.819e-1 medium ≈ 1e4× the floor band, vs the
GV5.1b amended seeds' 1.00× floor), the scaled + damped augmented
Newton shows **no quadratic contraction regime anywhere above the
floor band**: the clean-descent steps are line-search-capped halvings
(λ = 0.5, contraction exactly 0.30 dex → p = 1.00 **by construction**,
the backtracking cap — not Newton asymptotics), and the trajectory
then **stalls mid-range** (F_BL ~ 3e-2 → 1.3e-2 coarse / 2.2e-2
medium over 10 iterations), never reaching the floor band (still
4262× / 12867× the floor at the cap). Binding medium read: 7
above-band triples, median p = 0.56 ∉ [1.5, 2.5] → **FAIL** (honest);
coarse median p = 1.00 (RECORDED). The obstacle to the tight Newton
is therefore not only the formulation floor at ~3e-6: there is a
**mid-range descent barrier 3–4 decades above it**. Whether a
quadratic basin exists ADJACENT to the floor is not probed by these
far seeds — that is the near-band-seed follow-up question (candidate
GV5.1d, user adjudication). The damping arm is inert once more (μ
rejection-retries = 0, μ decay 1e-6 → 5e-11; the line search does
all the globalization). Band (a) PASS both levels: the machinery is
exact, the e2 tolerance pre-registered cond-aware this time (no
adjudication needed).

## 1. Configuration as executed

Per the pre-registration, no deviations. Base seed = the GV5.1
amended protocol verbatim (HEAD-regenerated loose-converged state,
committed GV3.1 recipe, wiring guard |dcl_k0| ≤ 1e-8; the GV5.1b
runner's `build_loose_state` imported): coarse converged n_outer=4
cl=0.26791639 (bit-identical to the committed 16-thread regen),
medium converged n_outer=3 cl=0.28245999 (**a new fixed point — the
trajectory scatter again, §4**). Perturbation: U[:, 0] (δ) × (1+ε)
at the free BL nodes (inflow band untouched), ε = 1.0e4 both levels
by the deterministic calibration (2 evaluations each: the ε = 1e4
bracket edge landed in the T1 window [5e-2, 5e-1] directly — the
F_BL response to δ-scaling saturates, §3.4). Solve:
`newton_tight(scaling="rowcol", lm_damping=True, floor_stop=True;
max_iter=10, tol=1e-8, tol_abs=1e-10)` — the GV5.1b machinery
verbatim, assembly untouched. The T2 escalation did NOT fire (≥ 3
above-band triples on both T1 trajectories).

## 2. Gate bands

| band | criterion | coarse | medium (binding) | verdict |
|---|---|---|---|---|
| (a) suite + live identities on the perturbed-seed J | suite green; e1 ≤ 1e-12, e2 ≤ max(1e-10, 10·κ₁·eps) PRE-REGISTERED, e3 ≤ 1e-6 | 37 green (fleet 28 + 9 new); e1 2.6e-16 / e2 2.06e-9 (tol 3.9e-2) / e3 3.99e-10 | 37 green; e1 2.3e-16 / e2 2.40e-9 (tol 5.2e-2) / e3 2.57e-10 | **PASS** (both) |
| (b) pre-floor slope-2 window | ≥ 3 above-band triples, median p ∈ [1.5, 2.5] | 9 triples, median p = 1.00; regression slope 0.75 | 7 triples, median p = 0.56; regression slope 0.62 | **RECORDED** (coarse) / **FAIL** (medium) |
| (b) floor-band entry | recorded | never (final F_BL 1.344e-2 = 4262× floor) | never (final F_BL 2.203e-2 = 12867× floor) | RECORDED |
| (c) counts | N_polish ≤ 2× loose aspirational | 10 vs 8 — NOT met (far seed; traversal count) | 10 vs 10 — met (degenerate) | RECORDED |
| calib | seed F_BL in T1 window [5e-2, 5e-1] | ε 1.0e4 → 3.219e-1 (1.02e4× band) | ε 1.0e4 → 1.819e-1 (1.06e4× band) | RECORDED |

Wall: coarse leg 189 s, medium leg 170 s (8 threads, machine
unloaded; NOT comparable to the GV5.1b 16-thread/loaded entries).

## 3. The descent anatomy — λ-capped halving, then a mid-range stall

F_BL sequences (max-norm; floor bands coarse 3.16e-5 / medium 1.71e-5):

- coarse: 3.219e-1 → 1.608e-1 → 8.025e-2 → 3.999e-2 → 3.500e-2 →
  3.492e-2 → 1.728e-2 → 1.620e-2 → 1.519e-2 → 1.512e-2 → 1.344e-2
  (cap, iter 10)
- medium: 1.819e-1 → 9.083e-2 → 4.535e-2 → 3.084e-2 → 3.081e-2 →
  3.062e-2 → 3.680e-2 (bounce) → 2.357e-2 → 2.214e-2 → 2.211e-2 →
  2.203e-2 (cap, iter 10)

1. **The clean-descent phase is the line search, not Newton.** The
   first 3 (coarse) / 2 (medium) accepted steps all sit at λ = 0.50
   with contraction exactly ~0.30 dex — the full Newton step overshoots,
   the first halving is accepted. The resulting p reads are
   log(0.5)/log(0.5) = 1.00 **by construction**: they measure the
   backtracking cap, and even taken at face value they are NOT in
   [1.5, 2.5].
2. **The mid-range stall.** From iteration ~4 the trajectory enters a
   slow crawl: λ collapses to 1e-3–0.1, F_BL wanders (medium even
   bounces 3.06e-2 → 3.68e-2 at iter 6 while the merit still
   decreases — the F_φ block trades against F_BL), and the
   stall-phase p values (0.00–340) are noise-amplified ratios of
   near-zero contractions, not physics. The merit keeps creeping
   (hence termination = cap, not floor_reached), but 10 iterations
   only buy F_BL ~ 3e-2 → 1.3e-2/2.2e-2 — the trajectory never
   comes within 4 decades of the floor band.
3. **μ inert from the far seed too.** Zero rejection-retries on both
   legs; μ decays 1e-6 → 5.1e-11 monotone. The Levenberg arm never
   engages; the P8/P14 line search carries all the globalization.
   The row/column scaling remains the active solver ingredient
   (band (a) exact).
4. **The F_BL response to δ-scaling saturates.** δ × 10001 lifts the
   residual max-norm only to ~0.2–0.3 (from ~3e-6/1.8e-6): the
   max-norm residual is not linear in the δ perturbation at wild
   states (recorded; not chased — the calibration is deterministic
   and landed in-window at the bracket edge, 2 evals per level).

## 4. Honesty notes

- **The medium fixed point scattered again** — a fourth fixed point
  (cl 0.28245999, n_outer = 3, unperturbed F_BL 1.824e-6 = 1.07× the
  diagnosed floor 1.712e-6), vs the GV5.1b 16-thread regen
  (0.28137437, 1.712e-6) and the committed GV3.1 history (0.27190697,
  1.160e-6). Mechanism = the committed GV5.1 §4 trajectory scatter on
  the near-null manifold, here triggered by the 8-thread environment
  (parallel-reduction order). The wiring guard (|dcl_k0| ≤ 1e-8,
  1.54e-9) validates the recipe, not the fixed point — by design
  (GV5.1 Addendum). The band-(b) read is unaffected: the seed is
  above-band by construction on ANY trajectory. Coarse reproduces
  bit-for-bit at 8 threads (cl 0.26791639, F_BL 3.153842e-6).
- A `RuntimeWarning: invalid value encountered in sqrt` fired in the
  edge-data chain at the wild perturbed state (mach² < 0 at extreme
  q²); the probe guard / inf-safe paths handled it — no non-finite
  state was ever accepted (recorded; expected at δ × 1e4).
- The seeds are deliberately extreme (δ × 10001, iter-1 ds_change ~
  56–65): the first iterations are the line search undoing the
  perturbation. The measured absence of slope-2 is the machinery's
  true behavior from these seeds; it does NOT exclude a local
  quadratic basin adjacent to the floor (not probed here).
- Wall times (189 s / 170 s at 8 threads, unloaded machine) are NOT
  comparable to the 16-thread ledger entries (GV5.1b: 317 s /
  1629 s under external load).

## 5. Consequences (recorded, not decided here)

1. **The far-seed window hypothesis is measured negative.** "The
   augmented Newton converges quadratically before the floor" is no
   longer untested — from ~1e4×-band seeds the descent is
   globalization-limited (p = 1.00 by construction) then stalls
   mid-range. The gate FAIL is the honest binding read.
2. **The new target is the mid-range stall / the basin question.** A
   near-band seed protocol (seed at 10–100× floor, e.g. a small-ε
   perturbation or a truncated IBL pseudo-time state) would probe
   whether a quadratic basin exists adjacent to the floor —
   candidate **GV5.1d**, user adjudication, against the alternative
   that the whole above-floor landscape is non-quadratic and only
   formulation work (GV5.5) moves it.
3. **GV5.5 (TE-band (B, δ) formulation) motivation strengthened:**
   the tight-Newton obstacle is bigger than the floor itself.
4. The scaled + damped machinery stays delivered and exact (band (a)
   PASS, cond-aware tolerance now pre-registered); the μ arm remains
   inert even from far seeds — recorded for GV5.4's preconditioner
   design.
5. The V4-reopen trigger stays parked (per the pre-registration, not
   invoked by this gate's failure). GV5.2/GV5.3/GV5.4/GV5.5
   sequencing = the user's call.

## 6. Artifact index

- `results/summary.csv` — 2 PASS / 1 FAIL / 7 RECORDED + wall rows.
- `results/seed_calibration_{coarse,medium}.csv` — the calibration
  trace (2 evals each, bracket-hi in window).
- `results/newton_history_{coarse,medium}.csv` — the T1 polish
  histories (F blocks, merit, μ, λ, p_fbl, p_merit, ds_change,
  termination).
- `results/compare.csv` — vs GV5.1b / GV5.1 / loose committed (read,
  never recomputed): seeds, identities, p sequences, contraction
  factors, regression slopes, counts, μ schedule.
- Tests: `tests/test_v5_above_band_seed.py` (9, synthetic maps).
  Suite pre/post execution: 37 passed (tight fleet 28 + 9 new),
  twice (the band-(a) discipline); full-suite baseline re-run
  alongside (see `docs/overview.md`).
