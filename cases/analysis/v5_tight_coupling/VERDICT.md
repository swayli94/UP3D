# GV5.1 VERDICT — augmented (φ, Γ, U) Newton: exactness PASS, quadratic tail FAIL (IBL floor)

Gate: `docs/roadmap/track_v.md` **GV5.1** (2026-07-23). Pre-registration:
`PRE_REGISTRATION.md` (committed c56046d before the first execution;
addendum 98f9d81 pre-execution; addendum 2, 606b149, user-adjudicated
seed amendment, committed before the amended execution). Runner:
`run.py` (regenerates every artifact in `results/`). Executed
2026-07-23 on `kimi/track-v5-gv5-1`, M6 coarse + medium.
Full row-level record: `results/summary.csv` — **9 PASS / 1 FAIL / 36
RECORDED** across both protocols (k=1-seed and amended).

**Answer to the gate question, up front:** the augmented Newton's
Jacobian is **exact** — band (a) FD exactness PASSES at both levels, at
the seed AND at the last Newton iterate (worst sweet-spot 5.07e-9
medium / 2.25e-8 coarse vs tol 1e-5, 0 masked rows, frozen-veps
omission ~2–3e-8 scaled as pre-registered). The quadratic tail —
band (b), binding on medium — is an **HONEST FAIL**: the polish runs
10 iterations without converging and shows no slope-2 regime
(p = 0.02 / 0.50 / 16.07). The blocker is measured and it is NOT the
coupling: F_BL sits from iteration 0 at 1.708e-06, the medium
loose-final state's own IBL floor — the steady BL block has an
intrinsic residual floor on a near-null manifold
(cond(J_BL,BL) ~ 4e10) that the standalone pseudo-time IBL solve also
cannot converge below. The FP/coupling side is healthy throughout
(F_φ at 1.16e-7 by iteration 1). Tight coupling therefore does not
yet beat the loose loop on iteration count (band (c) not met:
N_polish = 10, N_total = 14/13 vs loose 4/5), and the follow-up that
matters is the IBL floor, not the Jacobian.

## 1. Configuration as executed

Pre-registered system, built and verified in three stages before the
gate runs:

- Monolithic augmented state (φ_free, Γ, U) with the un-eliminated
  [J_ff B; K −I] pressure-Kutta form; splu direct solve at 2.5-D scale;
  P8/P14 line search (probe guarded: ArithmeticError + numba-boundary
  SystemError → merit = inf, keep halving — the medium iter-6 crash
  class). GMRES + block preconditioner deferred to GV5.4 as
  pre-registered.
- J_φ,BL = −(Tᵀ W S P L D)[free] from the shipped linear pieces
  (Stage 1, `pyfp3d/viscous/tight.py`); J_BL,φ = J_e·D_ue·G with the
  closure-packet edge-derivative extension douts_e (30,2) =
  d/d(re_d, e_prime) threaded through an 8-wide derivative stack, state
  columns bit-identical (Stage 2: `closures.py`, `ibl3.py`,
  `surface_mesh.build_edge_jacobian_pattern`, isentropic field
  derivatives); J_φφ augmentation dṁ/dφ through ρ_e·u_e (Stage 3,
  `tight_driver.py`). Edge basis pinned to 7 per-node scalars
  (q, ρ, μ, M, û); global diffusion scales veps = eps·max(q) frozen
  within a Newton step (decision 5, omission measured per direction).
- FD protocol as pre-registered: zones state-independent (fixed
  geometry), discrete switch kept, central FD, kink-row masking
  ≤ 2 % of F_BL rows — 0 rows masked in the event, both levels.
- The FD gates arbitrated three real assembly bugs during Stage 2:
  drhom wa-weight (ibl3), the s1e two-factor diffusion chain (ibl3),
  and closures dR range(6)→range(7) (poisoned re_d column in the
  turbulent branch). All fixed before the gate runs.

Seed protocol, as amended by user adjudication (Addendum 2): the
original pre-registered k=1 seed gave a measured basin failure (Newton
crawls on the stalled k=1 IBL floor, coarse RECORDED); the amended
protocol seeds from the loose loop's converged state (committed GV3.1
recipe, regenerated) and reads band (b) on the polish history, with
standalone N_aug ≤ 2 recorded NOT met and no further standalone
retries. The V4-reopen trigger was considered and NOT invoked — the
stall mechanism is understood; globalization is recorded as a
follow-up.

The medium seed required a second user adjudication (§4): the
committed GV3.1 medium fixed point is not reproducible from any
current checkout (environment drift amplified at the IBL floor). The
HEAD-regenerated seed was accepted under a wiring guard
(|dcl_k0| ≤ 1e-8: PASS, 1.309e-9), with pointwise offsets RECORDED
under the scatter caveat.

## 2. Gate bands

| band | criterion | coarse | medium (binding) | verdict |
|---|---|---|---|---|
| (a) FD exactness, seed | worst sweet-spot < 1e-5 | 2.246e-8 | 5.074e-9 | **PASS** |
| (a) FD exactness, endpoint | worst sweet-spot < 1e-5 | 2.244e-8 | 5.074e-9 | **PASS** |
| (a) masked rows | ≤ 2 % of F_BL rows | 0 / 1236 | 0 / 2460 | **PASS** |
| (a) veps frozen omission | recorded (decision 5) | 3.0e-8 scaled | 2.0e-8 scaled | RECORDED |
| (b) quadratic tail | converged, median(last 3 p) ∈ [1.5, 2.5] | converged=False, p = 0.98/3.68/0.57 | converged=False, p = 0.02/0.50/16.07 | **FAIL** (medium binding; coarse recorded) |
| (c) N_aug ≤ 2 | standalone from k=1 seed | not met (crawl on the k=1 floor) | — | NOT met (recorded, pre-registered honesty note) |
| (c) polish counts | recorded honest counts | N_polish = 10, N_total = 14 vs loose 4 | N_polish = 10, N_total = 13 vs loose 5 | RECORDED |
| seed cross-check | coarse: converged, \|dcl_k0\| ≤ 1e-8, \|dcl_final\| ≤ 1e-3, \|dds_max\| ≤ 1e-2 rel, n_outer ±1 | dcl_k0 = 1.6e-12, dcl_final = 1.5e-4, dds rel 2.7e-3, n_outer 4 vs 4 | — | **PASS** |
| seed wiring guard | medium: \|dcl_k0\| ≤ 1e-8 (user-adjudicated) | — | dcl_k0 = 1.309e-9 | **PASS** |
| seed pointwise offsets | recorded (scatter caveat, §4) | dcl_final = 1.5e-4 | dcl_final = 9.5e-3, dds_max = 3.4e-3 (rel 0.50), n_outer 3 vs 5 | RECORDED |

Wall times (idle machine): coarse leg 31 s total (loose 28 s, Newton
1 s / 10 iters); medium leg 76 s (loose 66 s, Newton 5 s / 10 iters).
The augmented machinery is NOT the cost driver — the loose seed
regeneration is.

## 3. band (b): the quadratic tail, measured

Medium polish history (`gv5_1_newton_history_medium.csv`): F_BL enters
at 1.709994e-06 — the loose-final state's own IBL floor — and moves
only in the fourth significant digit over 10 iterations
(→ 1.708240e-06). The line search collapses immediately
(λ = 9.8e-4 → 2.4e-4 → 1.5e-5 → … → 0), so the iterate freezes
(ds_change_rel at the last step = 5.5e-8, four decades below
tol_ds = 1e-3). The local-contraction ratio p never enters a slope-2
regime: 0.02 / 0.50 / 16.07 over the last three. Meanwhile the
coupling side is exact and done: F_φ absorbs its transpiration kick
within one iteration (3.2e-9 → 1.16e-7 scaled block, i.e. the kick
resolved at iteration 1) and F_Γ sits at 1.5e-9 throughout.

Coarse, same mechanism (RECORDED): F_BL pinned at 3.11e-06 = the
coarse loose-final IBL floor; p = 0.98/3.68/0.57.

Mechanism, as pre-registered as the basin risk and now measured at
the converged state too: cond(J_BL,BL) ~ 4e10 with 501/1236 singular
values below 1e-6·max at the k=1 state — the steady IBL block has a
near-null manifold along which the residual cannot be reduced below
its ~1e-6 truncation floor. The evidence that this is intrinsic to
the BL block rather than a tight-coupling defect: the standalone
pseudo-time IBL solve from the same loose-final state also stops
non-converged at the floor (medium 1.712e-06, coarse 3.154e-06;
committed loose finals 1.160e-06 / 3.113e-06 — same class). The
augmented Newton inherits the floor verbatim; exactness of the
coupling blocks does not help because the direction it computes is
dominated by the BL block's null structure, and the line search
correctly refuses to move along it.

## 4. Finding: the loose medium trajectory is chaotic at the IBL floor

Surfaced by the seed-regeneration cross-check and diagnosed in
`results/gv5_1_medium_seed_diagnosis.md` (user-adjudicated): the
committed GV3.1 medium fixed point is one sample of an
IBL-floor-trajectory scatter, not a reproducible point. Three
code/environment trajectories converge to three different fixed points:

- committed GV3.1: cl = 0.2719, ds_max = 6.84e-3 (n_outer = 5)
- HEAD regen: cl = 0.2814, ds_max = 3.45e-3 (n_outer = 3)
- pre-Stage-2 checkout (c2dc325, same env): cl = 0.2217, ds_max =
  9.73e-3 (n_outer = 6)

Mechanism: 1e-12-level perturbations (douts derivative columns,
environment drift) amplified on the cond ~ 4e10 near-null manifold at
the 100-iteration-capped IBL inner solves. Controlled exclusions: the
recipe diffs CLEAN line-by-line, the mesh md5 is unchanged, k = 0
bit-matches to 1.3e-9, and HEAD is bit-identical run-to-run. Coarse
is well-conditioned by comparison (0.14 % k=1 ds scatter).

Consequences already applied: all GV5.1 medium seed-vs-committed
comparisons carry the scatter caveat; the wiring guard
(|dcl_k0| ≤ 1e-8) passed at 1.309e-9. Standing caveat for the
project: **GV3.1-era committed medium numbers (cl, ds) are
trajectory samples with a demonstrated ±0.03 cl scatter** — quote
them with that band until the IBL floor is fixed.

## 5. Consequences (recorded, not decided here)

1. **The IBL floor is the real blocker for tight coupling** — not the
   Jacobian, not the coupling blocks, not the line search. The
   follow-up that unblocks band (b): diagnose the J_BL,BL near-null
   structure at the converged state (which rows/variables carry the
   4e10 conditioning — closure DELTA_MIN floor region? TE band?
   re_d column scaling?), then a formulation-level regularization or
   a projected/templated solve for the BL block. Until then, Newton
   convergence below ~1e-6 F_BL is not reachable at any seed.
2. **Tight coupling does not yet beat loose on iteration count**
   (N_total 14/13 vs 4/5). The loose loop's ds-based stopping
   criterion (tol_ds = 1e-3) converges long before the F_BL floor
   matters — the two loops are not stopping at the same tolerance on
   the same functional. Any future N_aug comparison should be read
   against an F-norm-matched loose criterion.
3. **V4-reopen trigger remains parked** (considered in Addendum 2,
   not invoked): the stall mechanism is understood and is not the
   transpiration/Kutta closure class V4 covers. Globalization of the
   augmented Newton (trust region / pseudo-transient continuation) is
   recorded as a secondary follow-up behind the IBL floor work.
4. **GV5.2 / GV5.3 / GV5.4 sequencing is open for adjudication.** The
   exact Jacobian machinery (Stage 1–3 deliverables, FD-verified) is
   reusable regardless: GV5.4's GMRES + block preconditioner work now
   has an exact operator to precondition. Whether the IBL-floor
   follow-up enters as GV5.1b or a new gate, and whether GV5.2/5.3
   proceed on the loose loop meanwhile, is a user decision.
5. The coarse seed cross-check band (environment-drift class,
   ~1.5e-4 cl) and its documentation in run.py's docstring stand as
   the template for future seed-regeneration checks: bit-exact k=0,
   band-checked final, wiring-guarded.

## 6. Artifact index

- `results/summary.csv` — 46 rows, 9 PASS / 1 FAIL / 36 RECORDED,
  both protocols.
- `results/gv5_1_fd_report.csv` — per-scope FD ladder (seed and
  endpoint, both levels): sweet-spot errors, masked-row counts, veps
  omission per direction.
- `results/gv5_1_newton_history_{coarse,medium}.csv` — amended-seed
  polish histories (block residuals, λ, ds_change, per-iter wall).
- `results/gv5_1_newton_history_coarse_k1seed.csv` — the original
  k=1-seed protocol's crawl (basin-risk evidence).
- `results/gv5_1_seed_reproducibility_{coarse,medium}.csv` — loose
  recipe regeneration vs committed GV3.1 (cross-check bands, wiring
  guard).
- `results/gv5_1_medium_seed_diagnosis.md` — the §4 scatter finding,
  full exclusion chain.
- `results/gv5_1_compare.csv` — augmented vs committed loose:
  iteration counts, ds_max, cl_p/cl_kj, IBL floors, with the scatter
  caveat on all medium pointwise rows.
- Implementation (FD-verified, regression-covered):
  `pyfp3d/viscous/tight.py`, `pyfp3d/viscous/tight_driver.py`,
  Stage-2 derivative stack in `closures.py` / `ibl3.py` /
  `surface_mesh.py` / isentropic helpers; tests
  `tests/test_v5_tight_jacobian.py` (8), `tests/test_v5_tight_edge.py`
  (7), `tests/test_v5_tight_system.py` (5), shared fixture
  `tests/v5_state.py`.
