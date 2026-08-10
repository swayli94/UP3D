# GV3.1 / GV3.2 Pre-registration (written before any gate execution)

Binding gate text: `docs/roadmap/track_v.md` GV3.1/GV3.2 (2026-07-22
re-spec, V3 — Loose coupling). User directive 2026-07-22 (Track V
opening): **no gate re-spec** — the bands below are the execution-level
operationalization the gate text itself requires, not a redefinition.

The pytest suite (`tests/test_v3_coupling.py`) covers the wiring and a
coarse-mesh smoke of the same machinery; THIS gate is the
committed-evidence execution (CSV/PNG artifacts + VERDICT) per the project
evidence discipline.

## Case conditions (binding, both gates)

- Geometry/mesh: `cases/meshes/naca0012_2.5d/{medium,coarse}.msh`,
  wake-cut (`cut_wake(read_mesh(...))`); the IBL surface is built on the
  cut wall. **medium is the primary level** (A4 input band ≈2.5 %);
  coarse is the cross-check (RECORDED, not banded).
- Conditions: M∞ = 0.5, α = 2°, Re = 3.0e6 per chord, forced transition
  x_tr/c = 0.05 upper AND lower (XFOIL side: `VPAR XTR 0.05 0.05`).
- Driver chain: compressible Picard `solve_subsonic_lifting` with the V2
  `body_source_rhs` channel (loose loop of
  `pyfp3d/viscous/coupling.py::run_loose_coupling`).
- IBL configuration (binding): eps_diff = 0.005, eps_diff_s = 0.02 (the
  V1 calibration, unchanged); transition flags per side from x_tr; u_e
  per the A4 per-zone discipline (LE band x/c ≤ 0.05 linear +
  crease-gated smoothed, elsewhere quadratic); inflow Dirichlet over the
  stagnation **band** x/c ≤ 0.02 with per-node laminar (Blasius-matched)
  seed states frozen at outer iteration k = 1 (the single-station pinning
  variant was debugged 2026-07-22 into the near-singular near-separation
  basin even at α = 0 — recorded in VERDICT; the band restores the V1
  x0-line discipline).
- Reference: `cases/reference_data/naca0012_viscous_xfoil/` (generation
  script + committed CSVs; XFOIL 6.99, 280 panels, Ncrit default):
  `delta_star_cf_alpha2_m05_xtr005.csv` (PRIMARY), `polar_summary.csv`,
  `inviscid_summary.csv`. `xtr030` is a RECORDED sensitivity only (its
  upper surface transitions naturally at x/c = 0.2668 ahead of the 0.30
  trip — see the reference README).

## GV3.1 — δ*/c and C_f bands vs XFOIL + viscous Δcl

**Profile metric (binding).** Per-side station profiles: our station
table groups the extruded strip by (x, y); the two TE copies are split by
`side_node`. For each side (upper/lower) and each station, δ*/c
(`OUT_DS1`) and C_f (`OUT_CF1`) are averaged over the station's span
copies of that side. XFOIL's CSV is linearly interpolated onto our
station x/c per side. Relative error per station:
`err(x) = (ours(x) − xfoil(x)) / xfoil(x)`.

**Windows (binding):**

- **Banded window**: x/c ∈ (0.05, 0.95] — outside the LE/stagnation zone
  per the gate text.
- **LE zone** x/c ≤ 0.05: RECORDED, not gated (input-limited 4–7 % per
  A4).
- **Near-TE zone** x/c > 0.95: RECORDED, not gated. Rationale (declared
  before execution): our BL domain truncates at the TE (outflow), while
  XFOIL continues the BL into the near-wake — the edge-velocity models
  differ structurally there, independent of discretization quality.

**Bands (binding), medium level**: at EVERY station of the banded window,
|err| ≤ **25 %** for δ*/c and ≤ **15 %** for C_f. Quoted alongside, never
summed: the A4 inviscid-input band ≈2.5 % (medium). These are
method-to-method bands (D13-implementation vs XFOIL's Drela–Giles
closure), declared wide enough to absorb closure-family differences; the
VERDICT quotes the per-station errors in full.

**Δcl (binding)**: cl from `post/forces.py::wall_force_coefficients`
(s_ref = the strip's z thickness; probe-logged at every outer iteration).
Viscous decrement Δcl = cl_inviscid − cl_coupled > 0 (direction), AND its
magnitude matches XFOIL's own decrement (inviscid cl from
`inviscid_summary.csv` minus xtr005 viscous cl from `polar_summary.csv`)
within a factor of **2** (band: 0.5 ≤ Δcl_ours/Δcl_xfoil ≤ 2.0).

**RECORDED (no band)**: coarse-level profile errors in the same windows;
LE-zone and near-TE values; the k=1 inflow q_ref; per-iteration IBL
residual floors (the weakly-observable-direction floor ~3e-6, documented
in the VERDICT); xtr030 sensitivity; cl per outer iteration; ṁ stats.

## GV3.2 — loose-loop convergence + transonic recorded point

- **Binding**: on the GV3.1 medium case, ‖Δδ*‖/‖δ*‖ < 1e-3 in ≤ 10 outer
  iterations, under-relaxation factor recorded honestly (executed with
  ω = 1.0 unless instability forces otherwise — then the used ω is the
  recorded one).
- **RECORDED (no pass/fail)**: one transonic-attached point, M = 0.72,
  α = 2°, Re = 3e6, x_tr = 0.05/0.05, coarse mesh, conforming-Newton
  driver (`solve_newton_lifting` + `external_rhs`): outer-iteration count
  (or divergence), relaxation, IBL residual floor. Feeds the V4 skip
  decision; the DN6-predicted near-separation risk is measured, not
  chased.

## Numerical settings (recorded with the artifacts)

- Thread cap 16: NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
  OPENBLAS_NUM_THREADS=16.
- Picard driver defaults (subsonic lifting, compressible); Newton leg:
  upwind_c = 1.5, m_crit = 0.95, m_cap = 3.0, rho_floor = 0.05,
  tol_residual = 1e-10.
- IBL: tol 1e-9 (steady max|R|), pseudo-time cfl0 = 1, growth 2,
  max_iter 100 per outer iteration, warm-started.
- Loose loop: ω = 1.0, n_outer_max = 10, tol_ds = 1e-3 on the raw
  successive δ* fields.

## Addendum 2026-07-22 (execution-discovered, disclosed before rerun)

**C_f reference-frame fix.** The first execution surfaced that the
pre-registered metric compared two DIFFERENT cf normalizations:

- XFOIL's DUMP column `Cf` is **freestream-normalized**,
  `CF = TAU/(0.5*QINF**2)` (XFOIL 6.99 source, `src/xoper.f:1970,2178`;
  the local-frame variant `TAU/QLOC` in `blplot.f` is the interactive
  plotter only, not the dump).
- Our closure's `OUT_CF1` is **local-edge-normalized** (D13 (61),
  `cf = 2 τ/(ρ_e u_e²)` — `pyfp3d/viscous/closures.py`; the
  Blasius/power-law seeds confirm the local frame).

The banded cf metric below is therefore executed on the **freestream
frame** — our per-node cf converted before station-averaging,
`cf_fs = cf_local · ρ_e |u_e|²` (FP units u∞ = ρ∞ = 1) — which is the
like-with-like quantity the band was declared for. This is an
execution-level operationalization fix, NOT a re-spec: the windows,
bands (25 % / 15 %), every-station rule, and Δcl criterion are
unchanged. The as-registered raw local-frame comparison is still
computed and kept as RECORDED rows in summary.csv for transparency.
δ* is a pure length ratio (both sides δ*/c) and is unaffected.

**Pinned-band labeling (same discovery).** Stations inside the Dirichlet
inflow band (x/c ≤ inflow_band_x = 0.02) are prescribed per-node Blasius
boundary data, not solution — comparing them against XFOIL is meaningless
by construction (their closure outputs are the seed values, e.g. a
nonphysical cf at the exact stagnation station). They are labeled
`pinned` in the profiles CSV and excluded from every error statistic,
INCLUDING the RECORDED LE-zone max|err| (which now covers (0.02, 0.05]).
No banded-window station is affected; no gate metric changes.
