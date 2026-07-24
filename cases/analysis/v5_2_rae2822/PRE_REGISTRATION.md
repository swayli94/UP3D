# PRE-REGISTRATION — GV5.2: RAE2822 transonic VII vs committed experiment

Committed BEFORE the first code change, per discipline. Gate:
`docs/roadmap/track_v.md` V5 **GV5.2** (registered in the V5 design;
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
  the final phi). Experimental brackets from the committed data (the
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
- `cases/analysis/v5_2_rae2822/`: PRE_REGISTRATION (this file), run.py,
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
