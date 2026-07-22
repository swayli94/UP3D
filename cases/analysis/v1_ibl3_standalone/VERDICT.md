# GV1.1 VERDICT — standalone IBL3 verification

Date: 2026-07-22 · Branch: `kimi/track-v1-ibl3-core` · Binding text:
`docs/roadmap/track_v.md` GV1.1(a)–(e) (2026-07-22 re-spec) ·
Pre-registration: [PRE_REGISTRATION.md](PRE_REGISTRATION.md) (committed
before the first execution; user directive: **no gate re-spec**) ·
Evidence: `results/*.csv`, `results/gv1_1_panels.png`, `results/summary.csv`
(regenerate: `python cases/analysis/v1_ibl3_standalone/run.py`).

## Result table

| gate | metric (band) | measured | verdict |
|---|---|---|---|
| (a) | max interior \|H−2.59\|/2.59 (≤ 2 %) | 3.77 % (H ∈ [2.6036, 2.6875] on 200×32) | **FAIL** |
| (a) | δ*(x) power-law exponent ([0.48, 0.52]) | 0.5287 | **FAIL** |
| (a) | δ* exponent, march reference (same fit) | 0.5149 | RECORDED |
| (a) | δ* exponent, FE downstream half (x > 1.2) | 0.5100 | RECORDED |
| (b) | max \|cf_FE − cf_march\|/cf_march at matched Re_θ (≤ 5 %) | 0.07 % (Re_θ ∈ [161, 653]) | PASS |
| (c)/P1 | H monotone non-decreasing (m=−0.0904) | rise 0.5038, worst dip 3.6e-3 | PASS |
| (c)/P2 | H rise ratio ([1.05, 4.0]) | 1.191 (H 2.642 → 3.146) | PASS |
| (c)/P3 | m=−0.05 rise < m=−0.0904 rise | 1.108 < 1.191 | PASS |
| (d) | quasi-2-D lock, laminar (\|B\|,\|Ψ\|<1e-10, Cτ2=CTAU_LAM) | 0 / 0 / pinned | PASS |
| (d) | quasi-2-D lock, turbulent (\|B\|,\|Ψ\|,\|Cτ2\|<1e-10) | 5e-19 / 7e-18 / 0 | PASS |
| (e) | refinement error strictly decreasing; order RECORDED | errH 9.0e-5→4.5e-5→6.3e-5; errds 2.8e-6→1.5e-6→2.5e-6 | **FAIL** |

Gate summary: **8 PASS / 3 FAIL** (+2 RECORDED diagnostics).
`run.py` exits 1 (honest FAIL).

## Analysis

**(a) H band — FAIL as pre-registered expectation.** The D13 laminar profile
family's self-consistent flat-plate fixed point is H* ≈ 2.7083 (+4.54 % vs
Blasius 2.59), established in Stage 2 by three independent constructions and
recorded as a known risk in the pre-registration. The measured interior H
rises from the Blasius Dirichlet inflow (2.5906) to 2.6875 at the outflow,
still approaching the asymptote; the max deviation 3.77 % sits at the
outflow stations. This is closure-family physics, not an implementation
error: the FE solution agrees with the closure's own 2-D ODE march to
8.9e-5 / 6.3e-5 in H (see (e)), i.e. the solver faithfully integrates the
family it was given. Verdict stands as written; no gate re-spec (user
directive).

**(a) δ* exponent — FAIL, same root cause.** The measured window
x ∈ [0.22, 2.18] is dominated by the family fixed-point adjustment
transient (H: 2.59 → 2.69), during which δ* grows faster than √x. The
RECORDED diagnostics make the cause quantitative: the march reference
fitted over the same stations gives 0.5149 (itself near the band edge),
and the FE fit restricted to the downstream half (x > 1.2, where H has
nearly equilibrated) gives 0.5100, inside the band. The 0.5287 full-window
value is the transient's fingerprint, not discretization error.

**(e) refinement — FAIL on a genuine numerical finding.** With the
corrected reference (see "Implementation fixes" below), the FE-vs-march
errors are clean enough to resolve the discretization's behavior, and they
split into two mechanisms:

- Near inflow (where the Dirichlet adjustment lives), the error is
  diffusion-dominated, ∝ ε·h: 4.5e-5 → 1.6e-5 in H for 100×16 → 200×32 —
  proper convergence, order ≈ 1 as the pre-registration anticipated.
- At the outflow strip a 2h grid-scale mode dominates: alternating-sign
  station-to-station error, seeded near the inflow, amplified downstream
  with growth rate ∝ 1/h (at ε=0: −2.4e-6 at x=1.0 → −4.4e-5 at x=2.1 on
  200×32; the same run on 100×16 stays ≤ 4e-6 everywhere). The D-HB
  isotropic diffusion (ν = ε·q·h̄) damps it per-cell at rate ∝ ε/h, which
  loses to the ∝1/h growth inside the design knob band: ε=0.005 and
  ε=0.01 (band edge) both leave the 200×32 outflow error above the 100×16
  level (ε=0.01: errH strictly decreasing but errds 2.64e-6 → 3.07e-6).

Verdict: the strong-form Galerkin discretization with the D-HB diffusion
has an under-damped streamwise grid mode on fine meshes — exactly the kind
of deficiency GV1.1(e) exists to catch. Recorded as V1 limitation with
follow-up: streamwise-upwind / SUPG-type stabilization of the defect
convection (Drela's production IBL codes upwind these terms), or the
anisotropic diffusion tensor already listed in D-HB. Absolute level of the
contamination at production settings: ≤ 6.3e-5 in H, ≤ 2.5e-6 in δ* on the
finest gate mesh — benign for V1's purposes but structural for V3+.

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
  h̄ = √(2A_tri); c_l = 0.09 (D-CT-2).
- Newton/PTC: tol = 1e-9 relative on the steady residual inf-norm, CFL ramp
  1 → 1e8 ×2, halving backtracking on the F_pt merit with cfl/8 retry on
  rejection (≤ 10), scipy spsolve.
- Meshes: (a)/(e) 50×8, 100×16, 200×32 on x∈[0.2, 2.2], z∈[−0.2, 0.2];
  (b) 100×16 same domain; (c) 80×8 on x∈[0.4, 1.4], m ∈ {−0.05, −0.0904}.
- 2-D reference marches: RK4, 4000 substeps per span, exact station
  landing, started at the physical inflow.
- Thread cap 16 (NUMBA/OMP/OPENBLAS).
