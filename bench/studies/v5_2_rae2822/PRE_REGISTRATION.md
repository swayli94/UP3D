# PRE-REGISTRATION — GV5.2: RAE2822 transonic VII vs committed experiment

Committed BEFORE the first code change, per discipline. Gate:
`phases/p1/docs/roadmap/track_v.md` V5 **GV5.2** (registered in the V5 design;
user-sequenced 2026-07-24: GV5.1d → GV5.5 ✓ → **GV5.2** → GV5.3 →
GV5.4). Executed under the **temporary 8-thread session constraint**
(runner default 16; wall times flagged non-comparable).

## 1. Question and scope

Does the loose VII loop — the committed GV3.1 recipe with the GV3.2
transonic-point Newton-driver protocol — reproduce the **committed**
RAE2822 experimental Cp (`cases/reference_data/rae2822_experiment/`,
ground truth, never modified): shock location inside a pre-registered
band (PASS/FAIL, medium binding), Cp RMS RECORDED against the A4 input
band?

Explicitly OUT of scope: tight-coupling transonic (the IBL-floor work
lives in GV5.1b/5.5, closed); a fine mesh level; any CL metric (none
registered — the roadmap quotes shock + Cp RMS only); wall-interference
corrections beyond the dataset-labeled conditions; wing-body.

## 2. Geometry (new machinery: the point-set airfoil path)

- **Source**: NASA GRC NPARC validation archive `geom.txt` (Cook,
  McDonald & Firmin, AGARD-AR-138, Table 6.1 measured ordinates; 65
  stations; x/c, z/c lower, z/c upper; sharp TE at (1, 0)),
  `https://www.grc.nasa.gov/www/wind/valid/raetaf/geom.txt`, fetched
  2026-07-24. Committed as `cases/meshes/rae2822_2.5d/rae2822.dat`
  (geometry INPUT — the nl7301 `.dat` precedent, NOT
  `cases/reference_data/`); source URL + sha256 recorded in the
  generator header.
- **New point-set path** in `pyfp3d/meshgen/planar.py`: load the
  ordinates, monotone-PCHIP resample, cosine-clustered `n_half` points
  per side (mirroring `naca0012_coordinates`' clustering), TE forced to
  the single sharp point (1, 0) — the ordinate set already closes there.
- **Mesh family** `cases/meshes/rae2822_2.5d/generate_rae2822.py`
  mirroring `generate_naca0012.py`: embedded wake sheet (the M0-style
  family the V-gates use), R_FAR = 15 chords, the same
  Distance+Threshold size-field recipe, `extrude_single_layer`, the same
  tags (`wall` / `farfield` / `symmetry` / `wake`). Levels: **coarse
  h_wall = 0.020 / medium 0.010 ONLY** (fine registered, not built —
  two levels suffice for the gate).

## 3. Band (a): the A4 TE-wedge pre-check (RECORDED-with-fallback)

- Measure the TE wedge angle on BOTH built meshes: the
  `post/surface.py::_wall_vertex_normals` guard ratio
  (|Σ area·n̂|/Σ area → degrees) plus a direct geometric measure on the
  resampled contour (spline tangents at the TE, upper/lower).
- Fallback clause (pre-registered): if the guard ratio < 0.05 (wedge
  ≲ 6°), the TE-band u_e recovery falls back to linear+smoothed (the
  existing per-zone machinery), recorded in the gate.
- Prior estimate from the ordinate set: ~13° total (upper ~11.4° /
  lower ~1.5° over the last 4.8 %c) → the guard is EXPECTED not to
  fire; recorded either way.

## 4. VII protocol

- The committed **GV3.1 loose recipe verbatim**:
  `run_loose_coupling(driver, case, cfg, probe)`, ω = 1.0,
  n_outer_max = 10, tol_ds = 1e-3 — with the driver =
  **`make_newton_lifting_driver` + NEWTON_ARGS** (upwind_c 1.5,
  m_crit 0.95, m_cap 3.0, rho_floor 0.05, tol_residual 1e-10), the
  GV3.2 transonic-point protocol (the only committed VII transonic
  precedent; IMPORTED from the GV3 runner, not re-invented).
- Conditions = the committed datasets' labeled (corrected) conditions
  **verbatim**, no further correction:
  - **P1** = (M 0.725, α 2.55°, Re 6.5e6) — `ExpCase7_RAE2822_…dat`
  - **P2** = (M 0.73, α 3.19°, Re 6.5e6) — `Expe_RAE2822_…dat`
- Transition: **forced x_tr/c = 0.03 upper/lower** (AGARD practice for
  these cases; `CouplingConfig.x_tr_upper/lower`), binding-recorded.
- Levels: coarse (recorded / crash-stop) + **medium (binding)**.
- Validity envelope (the V5 scope guards): attached / mildly-shocked,
  M_shock ≲ 1.3 — the computed pre-shock peak Mach is RECORDED; a point
  exceeding it reads as outside-envelope RECORDED, not FAIL.

## 5. Metrics and bands

- **Band (b) shock location** (PASS/FAIL, medium binding, coarse
  recorded): x_shock = the x/c of max |dCp/dx| on the upper-surface
  computed wall Cp at the loose-final state (the wall-Cp machinery on
  the final phi). **Addendum 2026-07-24 (pre-execution clarification):
  the search is restricted to the compression branch (positive dCp/dx)
  inside the mid-chord window x/c ∈ [0.2, 0.9]** — the unwindowed max
  |dCp/dx| picks the LE suction spike (experimentally ~0.3 at
  x/c = 0.0002), not the shock; the window covers both committed
  brackets with margin and excludes the LE/TE recoveries. Experimental
  brackets from the committed data (the
  two stations flanking the shock jump): P1 [0.525, 0.55], P2 [0.55,
  0.575]. **Pre-registered acceptance = the bracket widened ±0.03c**
  (the G4.1 inviscid band): P1 **[0.495, 0.580]**, P2 **[0.520,
  0.605]**. PASS = the computed x_shock inside at BOTH points, medium.
- **(c) Cp RMS** (RECORDED, per the roadmap): RMS(Cp_cfd − Cp_exp) at
  the experimental stations (computed Cp linearly interpolated), per
  point per side, quoted with the A4 medium input band (~2.5 %
  peak-rel u_e; LE band 4–7 %) annotated separately from viscous-model
  error. Both committed layouts parsed (P1 two-zone upper/lower; P2
  wrapped single zone).
- **(d) convergence/guards** (RECORDED): loose converged ≤ 10 outer per
  point per level; final IBL floors; cl/cd_p histories; computed peak
  Mach; per-leg wall times.

## 6. Failure clauses (pre-registered)

- Shock outside the band on medium at either point → honest FAIL (not
  a crash).
- Loose non-convergence at a point → that point RECORDED (the
  recipe-limit clause, GV3.3 precedent); the other point still read.
- TE-wedge guard firing → the §3 fallback, recorded, the gate proceeds.
- Computed M_peak > 1.3 at a point → outside-envelope RECORDED for
  that point.

## 7. Tests

- `tests/test_meshgen_rae2822.py` (pure-python, cheap): ordinate
  load/resample (monotone, bounded by the data range, TE closed at
  (1, 0), cosine clustering); the wedge measure on the analytic
  contour; the Cp-compare helpers on BOTH committed layouts. Mesh
  generation itself is validated by band (a) + the committed stats.
- Tight fleet + full suite green; the new tests add to the baseline
  (three ledger baseline lines fill on suite completion).

## 8. Artifacts

- `cases/meshes/rae2822_2.5d/`: `rae2822.dat`, `generate_rae2822.py`,
  `coarse.msh` / `medium.msh` + stats (the NACA family commits its
  `.msh`; follow).
- `bench/studies/v5_2_rae2822/`: PRE_REGISTRATION (this file), run.py,
  VERDICT.md, results: `te_wedge.csv`, `shock_{level}.csv`,
  `cp_compare_{point}_{level}.csv`, `convergence_{point}_{level}.csv`,
  `summary.csv`.
- Seven-surface ledger sweep post-execution (track_v / roadmap /
  overview / agent-rules / PROJECT_STRUCTURE / cases-analysis README /
  design_track_v.md §18).

## 9. Runtime estimate (8 threads)

Meshgen minutes; VII legs: coarse ~5–10 min, medium ~20–40 min per
point (Newton-driver transonic outer = FP Newton + IBL solve), total
~1–1.5 h for the four legs. A leg exceeding 2 h reads as a recipe
problem — stop and adjudicate, do not silently grind.

## Addendum 2026-07-24 (execution mechanics — no band/protocol change)

Found and fixed between the first implementation commit and execution;
none of these touches the metrics, bands, or the VII recipe:

1. **Lower-surface ordinate sign.** Cook Table 6.1's `z/c lower` column
   is POSITIVE-DOWN (distance below the chord line; it changes sign
   where the aft lower surface rises above the chord). Physical
   `z_lower = −(tabulated)`; verified against the published signatures
   (max thickness 12.1 % @ 37.9 %c; camber 1.3 % @ 75.7 %c) and a
   segment-intersection sweep (0 self-intersections after the fix; 2
   before). Generator + `rae2822.dat` header updated; tests lock the
   signatures and the no-self-intersection property.
2. **`cut_wake` Kutta-probe fallback** (`pyfp3d/mesh/wake_cut.py`).
   RAE2822's reflex camber places BOTH TE flank neighbours on the +y
   side of the TE node, so the global-hint sign rule found no lower
   probe and raised. New fallback: when both strict passes come up
   one-sided, re-classify with the local TE-wedge bisector normal
   (aligned with the global hint). The fallback fires ONLY where the
   old code would have raised — previously working meshes are
   bit-identical. Regression: `test_kutta_probes_cambered_te`; the
   hard-rule-7 sweep (`test_topology_asserts_all_wake_meshes`) now also
   covers both new RAE meshes.
3. **Runner Cp side split** (`_wall_cp_sides`, `_peak_mach`): switched
   from centroid-y to the outward-normal y sign (the D11 idiom,
   `post/surface_ls.py`) — a centroid-y split mislabels the reflex band
   (9 aft lower triangles above y = 0 on coarse) as upper, which would
   contaminate the band-(b) shock window and the per-side RMS.
4. **Band (a) measured on the UNCUT mesh** (the A4 runner's method):
   on the cut mesh the TE strips no longer share the TE edge, so the
   max-crease measure does not see the TE wedge. Values recorded:
   coarse 9.46° / medium 9.92° mesh-crease vs 12.91° ordinate fit
   (x ≥ 0.95 secant sum; the gap = reflex-camber curvature over the fit
   window). Quadratic recovery available on both levels; the ≈6° guard
   clears by both measures → the §3 fallback does NOT fire.

## Addendum 2026-07-24 #2 (FP-driver stall rescue — pre-execution)

The first execution attempt crashed at coarse P1 outer iter 1 on the
GV3.3 loud-fail guard (`FP driver did not converge`). Diagnosis run
BEFORE any remedy (standalone, coarse mesh, inviscid k=0 call):

- M 0.725 / α 2.55 single-shot Newton (the GV3.2 protocol verbatim):
  NOT converged in 30 Newton iterations — the residual drops to
  ~2.7e-6 and then sits on a 4-iteration identical plateau (a
  shock-cell limit cycle), tol_residual = 1e-10 unreachable.
- M 0.72 / α 2.55 same protocol: converged in 7 iterations to 7.7e-13.
  The mesh and α are fine; the stall is specific to M 0.725+.

Pre-registered rescue (this addendum, committed before re-execution):
the FP driver keeps the GV3.2 single-shot as first choice for every
warm-started (k ≥ 1) solve, and falls back to the library's DESIGNATED
transonic path — `solve/newton.py::solve_newton_transonic`, upward Mach
continuation from m_start = 0.70 with warm starts from converged levels
only and a STRICT final level at NEWTON_ARGS' tol_residual = 1e-10
(design.md Sec 8.1; the solve_newton_lifting docstring itself routes
transonic levels through it) — for (i) the k=0 cold start and (ii) any
in-loop single-shot that reports converged=False. Validated standalone:
M 0.70 level 5 Newton iters to |R| = 3.8e-13, M 0.725 final level 8
iters to |R| = 4.0e-11 (accept=tol, strict). If BOTH paths fail the
GV3.3 guard still raises — no silent acceptance anywhere.

Unchanged (still binding): NEWTON_ARGS verbatim, the loose recipe
(ω = 1.0, ≤ 10 outer, tol_ds = 1e-3), the conditions, all bands and
metrics. What changes is ONLY how each inner FP solve reaches its
converged state; the outer fixed point and every metric are computed
from converged states exactly as before. The per-call driver path
(single_shot / continuation) is counted in summary.csv.

## Addendum 2026-07-24 #3 (stall-acceptance tier — pre-execution)

Addendum #2's chain (single-shot strict → continuation strict) was
validated end-to-end and FAILED at outer iter 2 — a DIFFERENT failure
mode, measured (coarse P1, k=2 FP call, mdot_max = 1.30e-2):

- warm single-shot: 30 Newton iters, residual frozen at ~1.0e-9 for the
  last 6+ iterations (a hard shock-cell plateau, NOT slow convergence —
  raising n_newton_max to 60 changes nothing);
- strict continuation (M 0.70 → 0.725, with dm-halving retry via
  M 0.7125): final level plateaus at 5.2e-10 — also short of 1e-10.

The plateau state is converged to ~7 orders of magnitude of residual
reduction — 6+ orders tighter than the OUTER fixed-point tolerance
(tol_ds = 1e-3) the solution feeds into; tol_residual = 1e-10 is
inherited subsonic practice (cheap under quadratic convergence),
unreachable at a transonic shock-cell plateau. The library already
carries the honesty-guarded acceptance for exactly this state:
`solve_newton_lifting(accept_on_stall=True)` accepts a live plateau
(accept_reason "stall") ONLY when the Kutta constraint is converged
(f_norm < tol_gamma), no upwind limiter/floor activity is active, and
the plateau detector (live_stalled) fires — otherwise it keeps
iterating/reporting non-convergence.

Pre-registered FP-call chain (this addendum; supersedes #2's), ordered
cheap → deep, first success wins, every attempt's
(path, accept_reason, converged) logged:

- warm-started (k ≥ 1): single strict → single stall-accept →
  continuation strict → continuation stall-accept → the GV3.3 loud
  raise;
- cold start (k = 0): continuation strict → continuation
  stall-accept → raise.

Strict 1e-10 remains the FIRST choice at every call; a stall acceptance
is recorded per leg (summary.csv: fp_calls / fp_continuation /
fp_stall_accepted) and reported in the VERDICT. No silent acceptance:
if every tier fails, the guard still raises and the point reads
RECORDED per §6. All bands, metrics, conditions, NEWTON_ARGS and the
loose recipe remain binding and unchanged.
