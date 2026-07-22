# GV1.1 Pre-registration (written before any gate execution)

Binding gate text: `docs/roadmap/track_v.md` GV1.1(a)–(e) (2026-07-22 re-spec).
User directive 2026-07-22: **no gate re-spec** — bands below are the
execution-level operationalization the gate text itself requires ("band
pre-registered"), not a redefinition of the gate.

## Gate (a) — laminar flat plate → Blasius

- Case: prescribed `u_e = (1, 0, 0)`, rho = 1, mu = 1e-5, domain
  x in [0.2, 2.2], z in [-0.2, 0.2]; inflow Dirichlet = `blasius_seed(x0)`
  (H = 2.5906, theta matched exactly by construction).
- **Band (binding)**: every interior node (all boundary strips excluded)
  must have H within ±2 % of 2.59, i.e. H in [2.5382, 2.6418].
- **Power law (operationalized)**: least-squares fit of delta*(x) = C x^p
  over interior stations must give p in [0.48, 0.52] (gate text: delta* ∝ sqrt(x)).
- Known risk, recorded up front (family consistency analysis, Stage 2
  evidence): the D13 laminar profile family's self-consistent flat-plate
  fixed point is H* ≈ 2.7083 (+4.54 % vs Blasius), i.e. the closure
  equilibrates to its own fixed point within ~2 length units of the inflow.
  If the measured H(x) leaves the band downstream, the verdict is FAIL as
  written, with the fixed-point analysis recorded as the cause.

## Gate (b) — turbulent flat plate C_f(Re_theta)

- Case: prescribed `u_e = (1, 0, 0)`, rho = 1, mu = 1e-5, domain
  x in [0.2, 2.2]; forced-turbulent everywhere (D-TR); inflow = power-law
  turbulent seed at Ctau equilibrium.
- Reference (per design doc §3.2): the **same closure's** 2-D ODE march
  (von Kármán + kinetic-energy + stress-lag), generated inside the gate
  script; no external correlation is pass/fail.
- **Band (binding)**: C_f(x) from the FE solution within ±5 % of the ODE
  reference at every interior station.

## Gate (c) — decelerating u_e separation indicator

- Case: `u_e = u0 * x^m`, laminar; band of interest m in [-0.12, -0.05]
  (Falkner-Skan decelerating branch; FS separation at m ≈ -0.0904).
  Runs: m = -0.05 and m = -0.0904.
- **Pre-registered assertions** (indicator = H rise):
  - P1: for m = -0.0904, interior H(x) is monotone non-decreasing in x
    (tolerance: local scatter < 0.5 % of the total rise);
  - P2: H(last interior station) / H(first interior station) in [1.05, 4.0]
    for m = -0.0904 (clear rise, no Goldstein crossing inside the domain);
  - P3: the m = -0.05 run shows a strictly smaller rise ratio than
    m = -0.0904 (trend with m inside the band).
  - H values themselves: RECORDED (no magnitude band).

## Gate (d) — quasi-2-D structural lock

- On the (a) and (b) solutions: max |B| and max |Psi| over all nodes
  < 1e-10; max |C_tau2| < 1e-10 on turbulent nodes; laminar nodes have
  C_tau2 pinned at CTAU_LAM = 1e-8 exactly (D-TR pin, by design).

## Gate (e) — refinement order

- The (a) case run on three meshes (h, h/2, h/4 in both directions).
- Error measure: max-norm over interior of |H_FE - H_march(x)| and
  |delta*_FE - delta*_march(x)| vs the closure's own laminar 2-D ODE march
  (self-consistent reference; the family fixed point H* ≈ 2.7083 is the
  march's asymptote).
- **Assertion**: errors strictly decrease with refinement; the measured
  order p = log2(e_h/e_{h/2}) is RECORDED (no order band pre-claimed —
  artificial diffusion O(h) expected to dominate).

## Numerical settings (recorded with the artifacts)

- eps_diff = 0.005 (D-HB knob, mid-band), V_eps = eps * max(q).
- Newton: tol = 1e-9 relative on the steady residual inf-norm, CFL ramp
  1 → 1e8 geometric ×2, halving backtracking; c_l = 0.09 (D-CT-2).

## Amendment (2026-07-22, post-first-execution)

Written after the first execution; **no gate band above was changed**.
Three execution-level corrections, all in the reference/measurement code
and the solver globalization — the gates, bands, meshes, and eps_diff
setting are as pre-registered:

1. `run.py::march_2d` start-station teleport bug: the 2-D reference march
   now starts at the physical inflow (was: first recording station), and
   lands exactly on each station (integer substeps per segment; was:
   linear recording interpolation, ~1e-6 noise floor). This is the same
   pre-registered reference ("the closure's own 2-D ODE march"), evaluated
   correctly. The per-mesh (e) comparison evaluates the march at each
   mesh's own centerline stations (same pre-registered error measure,
   exact evaluation points).
2. `pyfp3d/viscous/ibl3.py::solve` backtracking merit changed from the
   pure steady residual to the pseudo-time residual F_pt (the function the
   Newton step linearizes). Required for GV1.1(c) convergence at all;
   convergence criterion unchanged (pure steady residual, tol 1e-9).
3. Two RECORDED diagnostics added under gate (a) (march-reference delta*
   exponent, downstream-half FE exponent) to quantify the fixed-point
   transient's effect on the power-law fit. They carry no pass/fail.

First-execution superseded numbers are not listed here; the git history of
this directory shows the sequence. Final verdicts: VERDICT.md.
