# PRE-REGISTRATION — GV5.5: TE-band (B, δ) formulation — breaking the IBL floor

Committed BEFORE the first code change, per discipline (the GV5.5
registration's own clause: "success criterion sketch — to be
pre-registered BEFORE code per discipline when the item opens").
Gate: `docs/roadmap/track_v.md` V5 **GV5.5** (registered 2026-07-24,
user-directed; opened 2026-07-24, user-directed). Design inputs: the
committed IBL-floor diagnosis (`../v5_ibl_floor/results/findings.md` =
`docs/design_track_v.md` §13 item 3 — the floor lives in the TE-band
(B, δ) equations essentially entirely inside J's range; the
pseudo-time controller bottoms out there; globalization alone cannot
pass it, GV5.1b) and the basin-hunting verdicts (GV5.1c/1d: no
quadratic regime anywhere from 1e4× down to 24× the floor — the floor
and its flat neighborhood are ONE formulation-level obstacle).

**Naming caveat (binding):** the "(B, δ) equations" are the diagnosis's
index labels — mechanically **row 0 = x-momentum** (the δ/θ carrier)
and **row 2 = kinetic-energy** (the shape/A carrier). On the quasi-2-D
testbed the crossflow state B is machine-zero (the quasi-2-D lock);
nothing in this item touches crossflow physics.

## Route choice (made at opening, per the registration)

The registration offered two candidates: (a) TE natural-outflow
discretization work on the (B, δ) equations; (b) closure
regularization in the TE band. **Chosen first: (a), as a row-level
TE-outflow treatment.** Rationale: Q4 killed the closure-floor
hypothesis (no clamp active), Q6 showed the floor is not
diffusion-tunable, and the residual localizes exactly at the
truncated-tent outflow rows (x/c 0.96–1.00) — the discretization,
not the closures, is the prime suspect. Route (b) and the heavier
(a)-variant (upwind boundary-flux, needs the boundary-edge machinery
that `sm.boundary_edges` only dormantly provides) are the
pre-registered escalation ladder if the chosen variant RECORDED-fails.

**Variant V1 — TE-outflow row replacement (first-order
extrapolation).** At each TE outflow node `i` (with its upstream chain
neighbor `up`), replace the two Galerkin rows by algebraic outflow
conditions:

- row 0: `R[i,0] = U[i,0] − U[up,0]` (δ, zero streamwise gradient —
  the momentum/θ carrier);
- row 2: `R[i,2] = U[i,1] − U[up,1]` (A, zero streamwise gradient —
  the energy/shape carrier).

Carrier pairing (recorded design decision): in 2-D the
momentum+energy pair determines (δ, A); extrapolating (δ, A) keeps
the constraint structure on the same variables and does NOT strand A
at the TE (a δ/B pairing would leave A determined only through
closure couplings — a new near-null direction by construction). The
quasi-2-D lock is untouched (B's own rows are not replaced). The
replacement needs only the CSR entries (i,i) and (i,up) — in-pattern
by construction (the TE node shares its last element with `up`);
a construction-time guard raises if a pair is out-of-pattern. J rows
zeroed with diag 1 and the −1 off-diagonal; the J_e (edge-Jacobian)
rows zeroed (the condition has zero edge derivative), exactly the
`_apply_rows` / `_apply_rows_edge` precedent. **Default OFF: flag
`te_extrapolate` default False = legacy path bit-identical** (the
established flag precedent; the committed GV1.1/V3/GV5.x evidence
base is untouched by construction). The (i, up) pairs are supplied by
the case layer as frozen data (`te_pairs`), derived from the strip
chain / TE station table — recorded in the runner.

**V0 — control**: the same standalone re-solve with the flag OFF must
reproduce the committed floors to roundoff (guard against seed drift).

## Testbed, seeds, solvers

The 2.5-D NACA0012 strip at the GV3.1 recipe (the GV5.1x testbed),
coarse + medium. Seeds = the GV5.1 amended protocol verbatim
(HEAD-regen loose-converged states; wiring guard converged +
|dcl_k0| ≤ 1e-8, read — never recomputed; abort + record on failure).
Solvers: (1) the standalone strip re-solve = the diagnosis Q7
pseudo-time protocol (imported from the committed v5_ibl_floor
runner, not mirrored), run to its stall with the flag OFF (V0) and ON
(V1); (2) secondary read: the GV5.1b tight-Newton polish from the
amended seed with the flag ON, the floor_stop terminal F_BL quoted
against the GV5.1b committed finals (read, never recomputed).

## Metrics (binding definitions)

The TE-row residuals are zeroed by construction under V1, so the
floor is read TWO ways to keep the metric honest:

- **m1 = the variant-system floor**: max-norm steady residual of the
  V1 system at its own converged/stalled state (the solver's metric).
- **m2 = the original-system residual at the V1 state** (the physical
  read): the flag-OFF residual evaluated at V1's terminal state — not
  gameable by row replacement. **m2 is the binding metric.**

Reference floors (committed, read — never recomputed): coarse
3.154e-6 / medium 1.710e-6 (the loose-final floors; the diagnosis Q5).

## Gate bands

- **(a) implementation exactness (PASS/FAIL).** The committed suite
  green with the flag default OFF (bit-identity by construction) +
  new synthetic tests `tests/test_v5_te_outflow.py`: row-replacement
  structure (replaced rows exact, other rows untouched), FD-vs-analytic
  J with the flag ON on a small strip, J_e row-zeroing, the
  out-of-pattern guard, the quasi-2-D lock retained with the flag ON,
  default-off bit-identity. Live FD check with the flag ON at the
  strip seed (the established FD discipline, tolerance as committed).
- **(b) floor descent (medium binding, coarse recorded).** **PASS =
  m2 ≤ 0.5× the committed floor on BOTH levels** (coarse < 1.577e-6,
  medium < 8.55e-7) **with all (c) guards green.** 0.5–0.9× = RECORDED
  (partial move; mechanism discussed; user adjudication). ≥ 0.9× =
  RECORDED (no move — the row-treatment route is dead; escalation =
  the user's call among the ladder: upwind boundary-flux (a)-variant /
  closure regularization (b) / accept the floor). m2 worse than the
  committed floor (> 1×) = RECORDED + the variant stays default-OFF
  (nothing rolls back — the flag is never on by default). The factor
  0.5× is chosen against Q6's wobble band (eps_diff ×4 moved the floor
  ≤ 24 %): 0.5× is a ≥ 2σ-class move, honestly distinguishable from
  knob-jiggling. m1 recorded alongside (not binding); the per-node
  residual anatomy vs the diagnosis's `residual_anatomy_*.csv` quoted.
- **(c) evidence-base guards (the registration's re-check clause;
  bands re-quoted with the flag ON, not silently inherited).**
  (i) plate H smoke with the flag ON: the GV1.1 laminar/turbulent
  plate cases re-run with `te_extrapolate=True`; H inside the
  re-quoted bands (laminar 2.55–2.75 / turbulent 1.2–2.0 — the
  committed GV1.1 bands, quoted, re-measured); (ii) loose-loop smoke
  with the flag ON: the strip loose loop converges ≤ 10 outer with cl
  within the A4 input band (±2.5 %) of the committed cl; (iii) the
  tight fleet + full suite green (flag OFF). A guard breach =
  RECORDED, variant stays default-OFF, user adjudication.

## Fallbacks and aborts (pre-registered)

- Wiring guard failure: abort, record, user adjudication (the GV5.1
  pattern).
- V0 control off the committed floors by more than roundoff
  (> 1 %): the seed/state drifted — abort the (b) read, RECORDED,
  user adjudication (the 8-thread medium scatter is EXPECTED to
  change the medium seed: the 4th fixed point cl 0.28245999; the
  floors are then quoted against the seed's OWN flag-OFF floor, the
  relative read, with the scatter flagged).
- The V1 solve diverges / stalls worse than V0 on either level:
  RECORDED with the trajectory, the variant stays default-OFF.
- The V1 tight polish (secondary read) behaves pathologically
  (non-finite step, rejection storm): recorded, not a crash.
- This item does NOT close V5: GV5.2/5.3/5.4 stay open regardless;
  the V4-reopen trigger stays parked.

## Out of scope

Upwind boundary-flux discretization (the (a)-ladder variant; needs
the boundary-edge table work); closure regularization (route (b));
wake coupling (V6); 3-D M6 application of the variant (a recorded
follow-up once the 2.5-D floor responds); GMRES/block preconditioners
(GV5.4). No change to any default code path: the flag is default-OFF
everywhere, including production case builders.

## Environment note (temporary, this session only)

Executed under a **temporary 8-thread constraint** (user-directed
2026-07-24, this session only): `NUMBA_NUM_THREADS=8 /
OMP_NUM_THREADS=8 / OPENBLAS_NUM_THREADS=8` vs the ledger-standard 16,
applied via the environment at execution time (the runner keeps the
16 defaults untouched) and recorded in summary.csv / VERDICT; wall
times are flagged NOT comparable to the 16-thread ledger entries. The
medium loose regen scatters across thread counts (the GV5.1 §4
mechanism) — see the V0 fallback clause.

## Artifacts

`pyfp3d/viscous/ibl3.py` (the flag + row machinery, default OFF) ·
case-layer `te_pairs` plumbing (frozen data) ·
`tests/test_v5_te_outflow.py` (synthetic + FD) ·
`cases/analysis/v5_5_te_floor/run.py` (regenerates everything) ·
`results/floor_probe_{coarse,medium}.csv` (V0/V1 trajectories) ·
`results/residual_anatomy_{coarse,medium}.csv` (per-node, vs the
diagnosis's committed anatomy) · `results/guards.csv` (plate H smoke,
loose smoke, FD live) · `results/summary.csv` (one row per band) ·
VERDICT at wrap-up.
