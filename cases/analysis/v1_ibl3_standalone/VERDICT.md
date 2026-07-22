# GV1.1 VERDICT — standalone IBL3 verification

Date: 2026-07-22 · Branch: `kimi/track-v1-ibl3-core` · Binding text:
`docs/roadmap/track_v.md` GV1.1(a)–(e) (2026-07-22 re-spec) ·
Pre-registration: [PRE_REGISTRATION.md](PRE_REGISTRATION.md) (committed
before the first execution; user directive: **no gate re-spec**) ·
Evidence: `results/*.csv`, `results/gv1_1_panels.png`, `results/summary.csv`
(regenerate: `python cases/analysis/v1_ibl3_standalone/run.py`).

**Current run (stabilized scheme, ε_s = 0.02): 9 PASS / 2 FAIL.**
First execution (isotropic-only diffusion): 8 PASS / 3 FAIL — the (e) FAIL
exposed the streamwise grid mode and was fixed by the D-HB tensor
follow-up; see "Addendum: (e) stabilization" below.

## Result table

| gate | metric (band) | measured | verdict |
|---|---|---|---|
| (a) | max interior \|H−2.59\|/2.59 (≤ 2 %) | 3.76 % (H ∈ [2.6036, 2.6875] on 200×32) | **FAIL** |
| (a) | δ*(x) power-law exponent ([0.48, 0.52]) | 0.5288 | **FAIL** |
| (a) | δ* exponent, march reference (same fit) | 0.5149 | RECORDED |
| (a) | δ* exponent, FE downstream half (x > 1.2) | 0.5101 | RECORDED |
| (b) | max \|cf_FE − cf_march\|/cf_march at matched Re_θ (≤ 5 %) | 0.07 % (Re_θ ∈ [161, 653]) | PASS |
| (c)/P1 | H monotone non-decreasing (m=−0.0904) | rise 0.5042, worst dip 3.3e-3 | PASS |
| (c)/P2 | H rise ratio ([1.05, 4.0]) | 1.191 (H 2.642 → 3.147) | PASS |
| (c)/P3 | m=−0.05 rise < m=−0.0904 rise | 1.108 < 1.191 | PASS |
| (d) | quasi-2-D lock, laminar (\|B\|,\|Ψ\|<1e-10, Cτ2=CTAU_LAM) | 0 / 0 / pinned | PASS |
| (d) | quasi-2-D lock, turbulent (\|B\|,\|Ψ\|,\|Cτ2\|<1e-10) | 0 / 0 / 0 | PASS |
| (e) | refinement error strictly decreasing; order RECORDED | errH 4.31e-4→2.21e-4→1.12e-4 (orders 0.96, 0.99); errds 1.02e-5→6.81e-6→4.39e-6 (orders 0.59, 0.63) | PASS |

Gate summary: **9 PASS / 2 FAIL** (+2 RECORDED diagnostics).
`run.py` exits 1 (the two (a) rows are honest FAIL).

## Analysis

**(a) H band — FAIL as pre-registered expectation.** The D13 laminar profile
family's self-consistent flat-plate fixed point is H* ≈ 2.7083 (+4.54 % vs
Blasius 2.59), established in Stage 2 by three independent constructions and
recorded as a known risk in the pre-registration. The measured interior H
rises from the Blasius Dirichlet inflow (2.5906) to 2.6875 at the outflow,
still approaching the asymptote; the max deviation 3.76 % sits at the
outflow stations. This is closure-family physics, not an implementation
error: the FE solution agrees with the closure's own 2-D ODE march to
2.2e-4 / 1.1e-4 in H on the two finer meshes (see (e)), i.e. the solver
faithfully integrates the family it was given. Verdict stands as written;
no gate re-spec (user directive).

**(a) δ* exponent — FAIL, same root cause.** The measured window
x ∈ [0.22, 2.18] is dominated by the family fixed-point adjustment
transient (H: 2.59 → 2.69), during which δ* grows faster than √x. The
RECORDED diagnostics make the cause quantitative: the march reference
fitted over the same stations gives 0.5149 (itself near the band edge),
and the FE fit restricted to the downstream half (x > 1.2, where H has
nearly equilibrated) gives 0.5101, inside the band. The 0.5288 full-window
value is the transient's fingerprint, not discretization error.

**(e) refinement — PASS after the stabilization fix (addendum below).**
errH 4.31e-4 → 2.21e-4 → 1.12e-4 (orders 0.96, 0.99 — the O(ε_s·h)
artificial-diffusion error the pre-registration anticipated) and errds
1.02e-5 → 6.81e-6 → 4.39e-6 (orders 0.59, 0.63) strictly decrease on both
measures. The first execution (isotropic-only diffusion) FAILED this gate:
the strict decrease broke 100×16 → 200×32 at the outflow strip. Root cause
and fix: see the addendum.

## Addendum: (e) stabilization (2026-07-22, post-first-execution)

First-execution finding (8 PASS / 3 FAIL run, git history): with the
corrected reference (implementation fix 2/3 below), the FE-vs-march errors
split into two mechanisms:

- Near inflow (Dirichlet adjustment strip), diffusion-dominated error
  ∝ ε·h, decreasing with refinement at order ≈ 1 — proper convergence.
- At the outflow strip, a 2h grid-scale mode: alternating-sign
  station-to-station error, seeded near the inflow, amplified downstream
  with growth rate ∝ 1/h (at ε=0: −2.4e-6 at x=1.0 → −4.4e-5 at x=2.1 on
  200×32; the same run on 100×16 stays ≤ 4e-6 everywhere). The isotropic
  D-HB diffusion (ν = ε·q·h̄) damps it per cell at rate ∝ ε/h and loses to
  the ∝1/h growth everywhere inside the design knob band (ε=0.005, and
  ε=0.01 at the band edge, both leave the 200×32 outflow error above the
  100×16 level).

Fix (user-directed, same day): the anisotropic streamwise tensor listed as
the D-HB follow-up — add ε_s·max(q)·h̄·(s1·∇G)(s1·∇N) to the diffusion
term, with s1 the element-averaged edge-velocity direction (edge data, not
state ⇒ Jacobian stays first-order exact; FD both lanes green). Calibration
sweep on the (a)/(e) family (ε_s ∈ {0.01, 0.02, 0.03, 0.05}, all strictly
decreasing on both measures):

| ε_s | errH@200×32 | errds@200×32 | orders H | orders ds |
|---|---|---|---|---|
| 0.01 | 8.77e-5 | 3.55e-6 | 0.97, 0.60 | 0.43, 0.53 |
| **0.02** | **1.12e-4** | **4.39e-6** | **0.96, 0.99** | **0.59, 0.63** |
| 0.03 | 1.56e-4 | 5.16e-6 | 0.96, 0.98 | 0.67, 0.69 |
| 0.05 | 2.43e-4 | 6.69e-6 | 0.95, 0.98 | 0.76, 0.77 |

Adopted: **ε_s = 0.02** — the knee: 0.01 damps only marginally (ds order
~0.5), 0.02 restores H order ≈ 1.0 with real margin at ~1.3× the minimal
absolute error; the O(ε_s·h) modeling error grows linearly beyond. The
solver default `eps_diff_s = 0.02`; 0.0 recovers the isotropic-only scheme.
SUPG (consistent, zero modeling error) remains the longer-term upgrade
route for V3+ if the O(ε_s·h) error becomes limiting; it needs the dropped
second-derivative terms or a non-exact Jacobian, out of V1 patch scope.

**(b), (c), (d) — PASS.** (b) cf(Re_θ) matches the closure's own turbulent
2-D march to 0.07 % over Re_θ ∈ [161, 653] (the 1/7th power-law
correlation is RECORDED in `gv1_1b_cf_retheta.csv`, not graded). (c) both
FS branches converge (14 Newton iterations each after the PTC fix below)
and the H-rise indicator behaves: monotone within tolerance for m=−0.0904,
rise ratio 1.191 in band, m=−0.05 rise strictly smaller (1.108). (d) the
structural quasi-2-D lock holds to machine zero on both regimes, as
designed (crossflow equations homogeneous in (B, Ψ, Cτ2); laminar Cτ
pinned at CTAU_LAM by D-TR).

## Implementation fixes made during gate execution

All three are recorded in `docs/design_track_v.md` §9; the first two are
bugs caught by the gate, the third is a measurement refinement. No gate
band was touched.

1. **PTC backtracking merit** (`pyfp3d/viscous/ibl3.py::solve`). The line
   search judged steps on the pure steady residual while the step is the
   Newton step of the pseudo-time residual F_pt = R + w(G−G_old). On the
   FS decelerating branch every step was rejected from a near-solution
   seed and CFL collapsed to its floor. Merit changed to F_pt (the
   function the step actually linearizes); FS branches then converge in 14
   iterations, flat-plate gates unchanged (`test_v1_ibl3.py` both lanes
   green).
2. **2-D reference march start-station teleport** (`run.py::march_2d`).
   The march initialized its integration at the first *recording* station
   (x0+0.1) with the inflow seed instead of marching from the inflow x0.
   The ODE is x-autonomous, so the whole reference trajectory was shifted
   by 0.1, and its first record was the un-marched seed. This produced a
   mesh- and ε-independent error plateau (3.3e-2 in H) that masked the
   true (e) behavior. Fixed by integrating from the physical inflow;
   corrected march agrees with the FE at the first comparison station to
   four significant digits.
3. **Exact-landing reference evaluation** (`run.py::march_2d`). The march
   now integrates each inter-station segment with an integer number of
   RK4 substeps and is invoked per mesh with that mesh's centerline
   stations, removing the ~1e-6 linear-interpolation noise floor that had
   hidden the 100×16 → 200×32 increase.

Also diagnosed and *rejected* as causes: ε_diff (sweep 0.005 → 0 showed
the plateau unchanged — diffusion was not the plateau), triangle diagonal
orientation (alternating-diagonal mesh: outflow mode unchanged), analytic
Jacobian under variable u_e (FD agreement 1.1e-9), and the 2-D-reduction
flux identities (Jx_x−u·M_x = ρq²θ and E_x−q²M_x = ρq³θ_s verified to
machine precision against the FE flux tables, derivatives included).

## Numerical settings (as run)

- eps_diff = 0.005 (D-HB, mid-band; band [0.001, 0.01]), V_eps = eps·max(q),
  h̄ = √(2A_tri); **eps_diff_s = 0.02** (D-HB streamwise tensor follow-up,
  calibration table in the addendum); c_l = 0.09 (D-CT-2).
- Newton/PTC: tol = 1e-9 relative on the steady residual inf-norm, CFL ramp
  1 → 1e8 ×2, halving backtracking on the F_pt merit with cfl/8 retry on
  rejection (≤ 10), scipy spsolve.
- Meshes: (a)/(e) 50×8, 100×16, 200×32 on x∈[0.2, 2.2], z∈[−0.2, 0.2];
  (b) 100×16 same domain; (c) 80×8 on x∈[0.4, 1.4], m ∈ {−0.05, −0.0904}.
- 2-D reference marches: RK4, 4000 substeps per span, exact station
  landing, started at the physical inflow.
- Thread cap 16 (NUMBA/OMP/OPENBLAS).
