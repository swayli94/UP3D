# CFL3D reference data — Euler and RANS, 2-D (gates D05 / D06 / D08 / D09)

Answers the reference-data request in
[docs/dev_phase_six/20260824-0400-gate-taxonomy-analysis.md](../../../docs/dev_phase_six/20260824-0400-gate-taxonomy-analysis.md)
§4 (Euler, cases 2D-1…2D-6) and §4.7 (RANS, cases R-1…R-6), under
**裁决三 — three-layer reference, aligned by model level**:

| our mode | reference | this directory | gate |
|---|---|---|---|
| full potential (inviscid) | **Euler** | `euler_naca0012/`, `euler_rae2822/` | D05, D06 |
| full potential + IBL (viscous-coupled) | **RANS** | `rans_naca0012/`, `rans_rae2822/` | D08, D09 |

The 3-D ONERA M6 rows of the request (3D-1, 3D-2, R-7 ⇒ gates D07 / D10) are
**NOT in this directory** — see “What is missing” at the end.

---

## ★★★ Read this before writing a criterion against these files

**A CFL3D solution is another NUMERICAL solution, not truth.** Every criterion
built on it must be phrased as *“difference from a recognised solver’s
Euler/RANS solution at a stated grid and convergence caliber”*, never as
*“error”*. That is why every row carries its grid level, its residual drop, its
wall spacing and (for RANS) its y+, and why every condition is run on a
**three-rung ladder**: a single grid is one number with no error bar, and a
reference without an uncertainty cannot carry a gate.

Two quantities in these files are the **uncertainty**, and a criterion tighter
than them is UNDEFINED, not failed:

- `grid_convergence.csv` → `delta_L2_L3` / `rel_delta_L2_L3` — the
  discretisation error bar of the reference itself;
- `turbulence_spread.csv` → `spread` / `rel_spread` (RANS only) — **|SST − SA|**,
  the resolution noise floor that comes from the turbulence model being a free
  choice. Same logic as the A4 2.5 % input band, and the reason GS4.1 round 15
  demoted E-CF from a gate to RECORDED: *a reference whose uncertainty has no
  boundary cannot be used to attribute anything.*

**Criterion shape is already ruled (裁决一, 2026-08-24)** and these files are
laid out for it: gate on **(a) shock position**, **(b) pre-shock Cp**,
**(c) cl**; record **(d) post-shock Cp** without a gate. Full potential is
isentropic and irrotational, so a real model difference **behind** the shock is
expected, and writing the criterion as “seven-station pointwise Cp difference”
would score a known model difference as a defect.

---

## 1. Provenance

### Solver

- **Code**: CFL3D v6.7, Aerolab modified build (NASA CFL3D HEAD 776d516,
  2022-08-12, plus the Aerolab additions), repo
  `CFL3D_6.7_THU_Aerolab` at commit `9fe2acb` (2026-08-03).
- **Binary**: `tools/cfl3d_seq` — an **external tool, NOT a repo dependency**
  (`tools/` is gitignored, exactly like `tools/xfoil/`). Build it with
  `cd build && make -f makefile_linux_gfortran cfl3d_seq` (gfortran, double
  precision `-fdefault-real-8 -fdefault-double-8`), then drop the executable at
  `tools/cfl3d_seq`. This build opens `cfl3d.inp` in the working directory
  explicitly — there is no stdin redirection — and reads the grid as **stream**
  binary PLOT3D (no Fortran record markers).
- **Deck**: written by the generator (see below); one deck per run is kept in
  the work directory, which is gitignored — the committed CSVs are the evidence.

### Grid generator — VENDORED

`cgrid_gmsh.py` in this directory is a **verbatim copy** of
`verification/src/cgrid_gmsh.py` from that same CFL3D repo/commit
(Apache-2.0, same licence as CFL3D itself):

    md5  1a1bffc2422435aa67a28f5a91dcb09e

It is vendored rather than located under `tools/` because without it the
dataset cannot be regenerated from this repository, and the project rule is
that committed evidence must be regenerable. The md5 above is how the copy is
checked against its source; **do not edit the file** — the deck extensions this
dataset needs are applied by `cfl3d_runner.py::_patch_deck` afterwards, so the
vendored copy stays byte-comparable to upstream.

Single-block C-grid, sharp trailing edge only. `j` runs lower far field → lower
TE → LE → upper TE → upper far field, `k` from the wall to a far field at
**80 chords**, `i` is the two symmetry planes of the 2-D case. Wall BC is
**1005** (inviscid flow tangency) for Euler and **2004** (adiabatic no-slip)
for RANS; far field 1000/1001; the wake cut is 1-to-1 blocking.

### ★ Geometry — from OUR generators, not from the CFL3D repo

The airfoil coordinates are written by
`cfl3d_runner.py::write_geometry` from **`pyfp3d.meshgen.planar`**:

| section | source | why |
|---|---|---|
| NACA0012 | `naca0012_coordinates` — analytic, **closed-TE** coefficient set (−0.1036), sharp point at (1, 0) | the same section every pyFP3D NACA mesh is built from |
| RAE2822 | `pointset_airfoil_coordinates` over `cases/meshes/rae2822_2.5d/rae2822.dat` (Cook Table 6.1 ordinates, PCHIP) | the same section `cases/meshes/rae2822_2.5d/` is built from |

The contour is emitted at `n_half = 1001` because the generator re-splines it
in arc length and then asks for ~9e-6 chord spacing at the nose; a source
coarser than that interpolates the nose shape.

★ **Why this matters**: the CFL3D repo ships its own `NACA0012.dat` /
`RAE2822.dat` (129 points each). Using them would put a **geometry difference
inside every comparison** — the “are the two numbers the same thing?” family of
error this project has logged six times (cross-pipeline, cross-level,
cross-provenance, cross-time, cross-thread, cross-mesh-generation). Here the
reference and the pyFP3D solution stand on a bit-identical section by
construction.

### ★ Angle of attack — experimental, uncorrected

**User ruling 2026-08-24**: “与实验数据的对比中，使用实验攻角还是使用修正攻角很难说，
我们目前先使用实验攻角。” Every case here uses the **experimental α with no
wind-tunnel correction**, the same caliber the pyFP3D side uses, so the
reference introduces **no correction constant of its own**. The published
corrected conditions (AGARD-AR-138, NASA TMR) differ mainly in α; the usual
consequence of leaving them uncorrected is a computed cl above the measured
one.

### Freestream / equation set

- γ = 1.4, Pr = 0.72, Pr_t = 0.9, T∞ = 460 °R (CFL3D defaults / deck value).
- **Euler**: `IVISC(I) = IVISC(J) = IVISC(K) = 0` — no viscous term and no
  turbulence model, wall BC 1005. The deck still carries `REUE,MIL`; it is
  set to the matching experiment’s Re only so the row reads sensibly.
  ★ **Its inertness is measured, not asserted**: the same case (M 0.80,
  α 1.25°, Euler L1) at `re_mil` 6.0 and 12.0 gives **bit-identical** CL, CD,
  CMZ (to 12 digits) and a bit-identical wall Cp over all 201 surface points.
  A claim about the solver gets an A/B, not a sentence.
- **RANS**: `IVISC = 7` (k-ω SST) and `IVISC = 6` (Spalart-Allmaras), **both
  models for every condition**, per §4.7, plus the keyword `edvislim 1.e05`
  on every RANS row — see §2.1, it is measured rather than assumed.

---

## 2. The grid ladder

Refinement is applied in **every** direction, not one at a time, so the ladder
is a refinement sequence rather than a single-knob sweep.

The **tangential and wake ladder is shared by both equation sets**
(`n_foil` 101 → 141 → 201, `n_wake` 41 → 57 → 81), so the Euler and RANS grids
differ **only** in the wall-normal direction. That is deliberate: an
Euler-vs-RANS read at the same condition then differs in the wall treatment and
not in the surface discretisation, which is what makes the
`FP − Euler` / `(FP+IBL) − RANS` error decomposition clean.

### Euler ladder (`EULER_LEVELS`)

| level | n_foil | n_wake | n_grow | block `nj × nk` | first cell height |
|---|---|---|---|---|---|
| L1 | 101 | 41 | 49 | 281 × 49 | 2.000e-3 c |
| L2 | 141 | 57 | 69 | 393 × 69 | 1.414e-3 c |
| L3 | 201 | 81 | 97 | 561 × 97 | 1.000e-3 c |

### RANS ladder (`RANS_LEVELS`)

| level | n_foil | n_wake | n_grow | block `nj × nk` | y+ target |
|---|---|---|---|---|---|
| L1 | 101 | 41 | 97 | 281 × 97 | 1.0 |
| L2 | 141 | 57 | 129 | 393 × 129 | 1.0 |
| L3 | 201 | 81 | 161 | 561 × 161 | 1.0 |

### ★ Why `n_wake` is 41/57/81 and not 29/41/57 — a measurement

The wake cut grows geometrically from the trailing-edge spacing out to 80
chords, so **too few wake points make its growth ratio large enough to
destabilise the transonic RANS cold start**. Isolated at L1 on
M 0.778 / α 2.03° / Re 6e6 / SST, one variable at a time:

| leg | outcome |
|---|---|
| `n_wake` 29, CFL 1 | **diverged** — NaN, block 3 cycle 8 |
| `n_wake` 29, CFL 1, `NITFO` 500 | **diverged** — NaN, cycle 10 |
| `n_wake` 29, CFL 0.5 / ramp 10 / `NITFO` 500 | ok, cl 0.333826 |
| **`n_wake` 41, CFL 1** | **ok**, cl 0.332273 |
| `n_wake` 29, CFL 1, **the CFL3D repo’s own NACA0012 ordinates** | **diverged** — NaN, cycle 17 |

★★ The last leg is the one that matters for provenance: swapping in the
upstream section changes nothing, so the failure was the **wake grading** and
**not our geometry**. Every Euler case had converged happily on the 29-point
wake — an equation set that tolerates a grid is not evidence that the grid is
good, which is why this was isolated instead of patched with a lower CFL.

★ y+ is held at **1 on every rung**, not refined: §4.7 item 3 — *a RANS result
whose boundary layer is not resolved is not a reference*. The first cell height
is sized from **each case’s own chord Reynolds number**
(`cf = 0.0576 Re_c^-0.2`, `y/c = y+/(Re_c √(cf/2))`), so it differs per case;
the measured value is in `forces.csv` (`h1_wall`), and the achieved y+ in
`yplus_avg` / `yplus_max`. The normal ladder therefore refines the **outer**
boundary layer and wake and brings the wall-normal growth ratio down
(≈1.19 → 1.14 → 1.11). Starting it as low as the Euler ladder’s 49 points
would put that ratio at ≈1.42, which is not a RANS grid.

**Multigrid constraint**: every block edge and every BC-segment end (both
trailing-edge j indices) satisfies `2² m + 1`, matching `ncg = 2` / `mseq = 3`
in the deck. The generator refuses to emit a block that does not.

---

### 2.1 ★★★ The eddy-viscosity cap `edvislim 1.e05`, and why it is on *every* rung

CFL3D's default is `edvislim = 1e10`, i.e. effectively unbounded. This dataset
caps the eddy viscosity at 10⁵·μ on **every** RANS row. Both halves of that
sentence are measurements.

**Why a cap at all.** Refining the RANS ladder made SST diverge (NaN on the
coarsest mesh-sequence level) while SA was untouched. A factorial on
M 0.778 / α 2.03° / Re 6e6, one count at a time:

| grid | model | outcome |
|---|---|---|
| 281 × 97 (L1) | SST | ok, 3.42 decades |
| 481 × 97 (`n_foil` only) | SST | ok, but only **0.79** decades |
| 361 × 97 (`n_wake` only) | SST | ok, 3.39 decades |
| **281 × 161 (`n_grow` only)** | SST | **diverged, cycle 9** |
| 561 × 161 (L3) | **SA** | **ok, 3.66 decades** |
| 561 × 161 (L3) | SST | diverged, cycle 11 |

⇒ the trigger is the **wall-normal count**, and it is **SST-specific**: SA
converges on the identical grid. With the cap, that same SST L3 case reaches
`|R|` 2.54e-09 (**3.44 decades**) and continues its own ladder —
cl **0.332325 → 0.339302 → 0.342969**, monotone, ~3.4 decades at every rung.

★ Four things were tried first and are recorded as measured negatives, so
nobody spends the time again: less mesh sequencing (`mseq` 3→2→1 — `mseq=1`
diverges at cycle 9 too, so it is not the sequencing); the leading-edge cell
aspect ratio (`le_factor` 0.06→0.24, `dx_le/h1` 0.65→10.3 — still diverges,
and this was the most plausible hypothesis since `dx_le/h1` degrades
monotonically along the ladder); `nsubturb 3` (**worse** — NaN moves from
cycle 11 to cycle 7); and `nfreeze` (refused outright by `precfl3d`:
*“not allowed to freeze turb model immediately upon mesh-sequencing-up to next
level”*, i.e. incompatible with `mseq > 1`).

**Why on every rung.** A limiter switched on only for the level that needed it
would make the three rungs **three different models** — precisely the
cross-provenance error this dataset exists to avoid. That obliges measuring
what the cap does where the un-capped run already converged, so it was A/B'd
with the cap as the only variable (M 0.778, `n_foil`/`n_wake`/`n_grow` fixed):

| rung | model | `edvislim` | cl | cd | cd_friction | \|R\| |
|---|---|---|---|---|---|---|
| L1 | SST | 1e10 (default) | 0.332325 | 0.025884 | 0.005371 | 6.714991e-09 |
| L1 | SST | **1.e05** | 0.332325 | 0.025884 | 0.005371 | 6.649566e-09 |
| L1 | SA | 1e10 | 0.413929 | 0.030476 | 0.006301 | 1.011592e-08 |
| L1 | SA | **1.e05** | 0.413929 | 0.030476 | 0.006301 | 1.011592e-08 |
| L2 | SST | 1e10 | 0.339302 | 0.026098 | 0.005385 | 1.895687e-08 |
| L2 | SST | **1.e05** | 0.339299 | 0.026100 | 0.005385 | **4.465264e-09** |
| L2 | SA | 1e10 | 0.422484 | 0.030850 | 0.006343 | 4.763216e-09 |
| L2 | SA | **1.e05** | 0.422484 | 0.030850 | 0.006343 | 4.763216e-09 |

★★ **SST**: the answer is unchanged — six digits at L1, five at L2 — while
`|R|` moves, and at L2 it moves the right way (2.75 → **3.38** decades). So the
cap **acted during the transient and left the converged solution alone**, which
is what makes it admissible rather than a thumb on the scale.

★★ **SA**: bit-identical including `|R|`, and that is structural, not lucky —
`edvislim` appears in `twoeqn.F` / `threeeqn.F` / `foureqn.F` and **not** in
`spalart.F`, so it cannot touch the SA arm at all. The SA rows are therefore
CFL3D's untouched default behaviour.

**The alternative that was rejected.** A gentler cold start also gets the
failing legs to *run*:

| `CFL_START` × `CFL_RAMP` (terminal CFL) | M0.778 L3 | M0.5 xtr005 L2 |
|---|---|---|
| 0.1 × 10 (1.0) | diverged, cycle 19 | diverged, cycle 15 |
| 0.1 × 20 (2.0) | diverged, cycle 19 | diverged, cycle 16 |
| **0.05 × 20 (1.0)** | ran | ran |
| **0.02 × 50 (1.0)** | ran | ran |

but *running* is not *converging*, and the caliber says so: all four reach only
**~1.1–1.6 decades**, and the two recipes that work disagree by **4.9 %** in cl
at M0.778 L3 (0.373755 vs 0.356217) while agreeing to **0.05 %** at
M0.5 xtr005 L2 (0.248883 vs 0.249002). ⇒ a startup-rescued L3 row would be a
**snapshot whose value depends on how it was started** — disqualifying for a
reference. It is kept only as `RANS_FALLBACK`, a single retry for a leg that
diverged anyway, and any row it produces is identifiable by its `recipe` column
and its low `resid_decades`.

★★★ **Read `resid_decades` before quoting any RANS row.** The A/B above fixes
what it is worth: rows at **≳ 3 decades** were startup-independent to 5–6
digits, and a row at **~1 decade** was startup-dependent at 4.9 %. A criterion
built on this dataset should require the former.

## 3. Conditions

### Euler — `euler_naca0012/` (D05), `euler_rae2822/` (D06)

| request | case tag | M | α | pairs with |
|---|---|---|---|---|
| 2D-1 | `n0012_m0500_a2.00` | 0.500 | 2.00° | `naca0012_m05` (an independent second opinion on our own panel-method reference) |
| **2D-2** | `n0012_m0800_a1.25` | 0.800 | 1.25° | ★★★ **the M1 target itself** — `naca0012_m080/shock_reference.csv` |
| 2D-3 | `n0012_m0720_a0.00` | 0.720 | 0.00° | subcritical → supercritical transition; α = 0 removes lift and wake coupling ⇒ single-variable read on artificial dissipation |
| 2D-3 | `n0012_m0750_a0.00` | 0.750 | 0.00° | as above |

#### ★★ 2D-3: use `cp_min` and `cd`, NOT the shock position

Measured: M 0.720 at α = 0 is **genuinely subcritical** (`cp_min` −0.6687 vs
`cp_critical` −0.6996 ⇒ `has_shock = 0` on both surfaces at every rung), and
M 0.750 carries only a **weak supersonic bubble**, not a shock. The evidence is
in the dataset’s own `n_cells` column: at L3 the detected “jump” spans
**83 stations**, i.e. a smooth recompression rather than a 1–2 cell shock, and
the Cp* crossing then moves **6.1 %** between L2 and L3 (0.2778 → 0.2958).

⇒ For 2D-3 the gate-worthy quantities are **`cp_min`** — which is converged to
4e-05 across the whole ladder (−0.778977 / −0.779017 / −0.778999) — and **`cd`**.
The shock position at these two conditions is **RECORDED only**. This does not
weaken 2D-3’s purpose: it was requested as a single-variable read on artificial
dissipation, and a converged suction peak at α = 0 with no lift or wake coupling
is exactly that.

★ These two conditions also carry the dataset’s strongest self-checks, which
need no reference at all: `cl` is **exactly** 0.000000 at all three rungs, and
the upper- and lower-surface shock positions come out **bit-identical**
(0.295782 both) — the symmetry the geometry and boundary conditions demand.
| 2D-4 | `n0012_m0778_a2.03` | 0.778 | 2.03° | `naca0012_experiment` M0.778 — same condition ⇒ model-form and discretisation error can be separated |
| 2D-4 | `n0012_m0803_am0.10` | 0.803 | −0.10° | `naca0012_experiment` M0.803 |
| 2D-5 | `rae2822_m0725_a2.55` | 0.725 | 2.55° | `rae2822_experiment` ExpCase7 |
| 2D-6 | `rae2822_m0730_a3.19` | 0.730 | 3.19° | `rae2822_experiment` M0.73 |

### RANS — `rans_naca0012/` (D08), `rans_rae2822/` (D09)

Every row exists **twice**, once per turbulence model (`sst`, `sa`).

| request | case tag | M | α | Re_c | transition |
|---|---|---|---|---|---|
| R-1 | `n0012_m0500_a2.00_xtr005` | 0.500 | 2.00° | 3.0e6 | trip **x/c = 0.05** |
| R-1 | `n0012_m0500_a2.00_xtr030` | 0.500 | 2.00° | 3.0e6 | trip **x/c = 0.30** |
| R-2 | `n0012_m0778_a2.03` | 0.778 | 2.03° | 6.0e6 | fully turbulent |
| R-3 | `n0012_m0803_am0.10` | 0.803 | −0.10° | 6.5e6 | fully turbulent |
| R-4 | `n0012_m0352_a12.86` | 0.352 | 12.86° | 3.0e6 | fully turbulent — ★★ **near stall** |
| R-5 | `rae2822_m0725_a2.55` | 0.725 | 2.55° | 6.5e6 | trip **x/c = 0.03** |
| R-6 | `rae2822_m0730_a3.19` | 0.730 | 3.19° | 6.5e6 | trip **x/c = 0.03** |

#### ★ Correction to the request: R-1 is at Re 3.0e6, not 6e6

The request table says “R-1 NACA0012 M 0.5, α 2°, **Re 6e6**, x_tr 0.05 与
0.30 — 对上 `naca0012_viscous_xfoil` 的两套”. The committed XFOIL set it is
meant to pair with is at **Re = 3.0e6** (`polar_summary.csv`, and the README
header of that directory). Run at 6e6 the three-way comparison
(XFOIL / RANS / us) would have stood on **two different Reynolds numbers**, so
R-1 here is at **3.0e6** — the value that makes the comparison well-posed.
Recorded rather than silently changed, because it is a departure from the
written request.

#### Transition treatment (§4.7 item 2)

Track V specifies transition **explicitly** (`x_tr/c` = 0.03 for RAE2822, 0.05
for M6; the XFOIL set has 0.05 and 0.30), so the tripped rows here are aligned
with it rather than being fully turbulent. Mechanism: CFL3D’s per-block
laminar-region card `ILAMLO ILAMHI JLAMLO JLAMHI KLAMLO KLAMHI`. Inside that
index box the turbulence **production** term is zeroed (`cutoff = 0`,
`twoeqn.F:2411` for k-ω/SST, `spalart.F` for SA), so no turbulence is generated
ahead of the trip and μ_t stays at freestream level. In the C-grid’s index
space “x/c ≤ x_tr on **both** surfaces” is one contiguous `j` interval
straddling the leading edge; CFL3D halves the index pair itself for each
coarser mesh-sequence level (`global.F:1361`), so the card is written for the
finest level only.

★★ **`i_lam_forcezero 1` is NOT used, and that is a measurement.** The keyword
additionally forces `vist3d = 0` inside the box, which looks like the stronger
and more faithful way to impose a trip. With it on, **every** tripped case NaNs
at cycle 4 on the coarsest mesh-sequence level; with it off and everything else
identical the same case converges. Isolated one variable at a time at L1
(M 0.5, α 2°, Re 3e6, SST):

| leg | outcome |
|---|---|
| no trip | ok — cl 0.248504, \|R\| 4.65e-09 |
| **trip 0.05, no `forcezero`** | **ok — cl 0.250708, \|R\| 4.65e-09** |
| trip 0.05, `forcezero`, k box 40 cells | **diverged** — NaN, block 3 cycle 4 |
| trip 0.05, `forcezero`, full k box | **diverged** — NaN, block 3 cycle 4 |

The k extent of the box makes no difference, so it is the hard zeroing itself
and not the region size. `forcezero_echo` therefore reads **0 on every row**:
the column records that the flag was not used, rather than being dropped and
leaving the question open.

★ **The trip is verified, not trusted.** The solver echoes the laminar region
it actually used, and the generator reads that echo back into the
`laminar_echo` / `forcezero_echo` columns and **raises** if a tripped case
reports no laminar region (or a fully-turbulent case reports one). A silently
unpatched deck would run fully turbulent while the CSV claimed a trip — the
“mentioned ≠ used” failure this project found by runtime probe in F06.

★ **R-4 is deliberately not gate material.** 12.86° on a NACA0012 sits at the
edge of stall; that is where the two turbulence models disagree most, which is
exactly what makes it the best measurement of the noise floor and the worst
choice for a tight gate. Its startup is also special
(`CFL_START = 0.05`, `CFL_RAMP = 20`, `NITFO = 1500`, `NCYC = 2000`) — measured
by the upstream verification harness, where a cold start at CFL 1 diverges in
~15 cycles — and the answer retains a residual CFL sensitivity. Read it as a
spread measurement, with its startup quoted.

---

## 4. Files and columns

Column names are written in every file (§4.4 item 6).

### `forces.csv` — one row per (case, turbulence model, grid level)

| column | meaning |
|---|---|
| `request` | data-request id (2D-1 … R-6) |
| `case`, `geometry`, `model`, `turb_model`, `level` | identity |
| `mach`, `alpha_deg`, `re_chord` | freestream; α is the **experimental, uncorrected** value |
| `x_tr`, `x_tr_actual` | requested trip and the last wall station actually inside the laminar box |
| `nj`, `nk`, `n_surface_cells` | block size and airfoil surface cell count |
| `h1_wall` | first cell height, chords |
| `yplus_avg`, `yplus_max` | achieved y+ (blank for Euler) |
| `cl`, `cd`, `cd_pressure`, `cd_friction`, `cm_quarter_chord` | integrals; **reference chord and area = 1, moment centre x/c = 0.25** (`SREF = CREF = BREF = 1`, `XMC = 0.25`) |
| `resid_final`, `resid_decades`, `ncyc_total` | convergence caliber on the finest mesh-sequence level. ★ Read this before quoting a RANS row — see §2.1 |
| `recipe` | the startup actually used, e.g. `cfl1x5_nitfo500_ncyc2000_mseq3`. A row rescued by `RANS_FALLBACK` carries a **`_fallback`** suffix — needed because R-4's designed near-stall override has the same knob values, so the knobs alone cannot tell a chosen startup from a retried one. Measured: the fallback fired **zero** times on this dataset, and the two `cfl0.05x20_nitfo1500_ncyc2000_mseq3` groups are R-4's override |
| `status`, `wall_s` | `ok` / `cached` / `diverged` / `timeout`, and wall time (**not a cost** — the box was loaded; quote it with its load or not at all) |
| `keywords` | the keyword block **echoed back by the solver**, not the one requested. A requested keyword missing from the echo raises |
| `laminar_echo`, `forcezero_echo` | what the solver said it did about transition. `forcezero_echo` is 0 on every row by design (§3) |
| `note` | why the condition is in the set |

★ **`cd_friction` is the new axis.** Drag had **no reference of any kind** in
this project before this dataset, and `cd_pressure` / `cd_friction` are
separated here because full potential + IBL can produce both (§4.7 item 5).
On the Euler rows `cd_friction ≡ 0` and `cd` is wave drag plus numerical
dissipation — which doubles as a self-check: a shock-free subsonic Euler case
must give ≈ 0 drag by d’Alembert.

### `cp_<case>_<turb_model>_<level>.csv` — the primary data

`x_c, y_c, cp, mach_local, surface` (`surface` ∈ {`upper`, `lower`}), one row
per surface **cell centre** on the `k = 1` line, x ascending. Euler rows use
`turb_model = none`.

### `shock.csv` — DERIVED

`has_shock, x_shock, n_cells, monotone, cp_min, cp_pre_shock, cp_post_shock,
cp_critical` per (case, level, surface).

★ **Computed with OUR OWN operator**, `pyfp3d.post.shock.shock_metrics` — the
same function the pyFP3D side is read with — because using CFL3D’s shock
definition on one side and ours on the other would compare two different
quantities. `cp_critical` is the exact isentropic sonic Cp*, and the shock
position is the last supersonic→subsonic Cp* crossing, linearly interpolated.
`cp_pre_shock` / `cp_post_shock` are sampled **0.05 c** either side.

⇒ These columns **move if that operator changes**. The primary data are the
`cp_*.csv` curves; `--derive-only` rebuilds this table from them with **no
solver run**, which is also how a criterion change is re-costed for free.

### `grid_convergence.csv` — the error bar

Per (case, turbulence model, quantity): the value on `L1`, `L2`, `L3` and
`delta_L2_L3` / `rel_delta_L2_L3`. Quantities are cl, cd, cd_pressure,
cd_friction, cm and the per-surface shock position.

★ Reported as an **interval, not an extrapolation**: three rungs at ratio ≈√2
do not by themselves establish an asymptotic order on a shocked solution, and
claiming one would be the same over-reach as the retired
“the 0.019 gap is resolution”.

★★ **`delta_L2_L3` (absolute) is the error bar; `rel_delta_L2_L3` is a
convenience and is MEANINGLESS where the quantity passes through zero.** Named
instances in this dataset, so nobody quotes them as regressions:

| row | `rel_delta` | why it is not a defect |
|---|---|---|
| `n0012_m0720_a0.00` / `n0012_m0750_a0.00`, `cd` | −70 % / −68 % | α = 0 inviscid ⇒ cd → 0 by d’Alembert. The absolute deltas are −1.2e-04 and −1.4e-04, i.e. the reference is **converging to zero as it should**. |
| `n0012_m0803_am0.10`, `cl` | ≈ −3 % | cl ≈ −0.029, near zero at α = −0.10°. Absolute delta 1.0e-03. |
| `n0012_m0500_a2.00`, `cd` | −47 % | subsonic inviscid ⇒ cd → 0. Absolute delta −1.2e-04. |

This is the mirror image of criterion defect #3 in the phase-three list (an
absolute threshold where the baselines differ 150×): here a **relative**
measure over a baseline that is legitimately ~0.

### `turbulence_spread.csv` — the noise floor (RANS only)

Per (case, level, quantity): `sst`, `sa`, `spread = |sst − sa|`, `rel_spread`.
**A criterion tighter than this is UNDEFINED, not failed.**

---

## 5. Regenerating

    conda activate up3d && export PYTHONNOUSERSITE=1
    cd cases/reference_data/cfl3d
    python generate_cfl3d_reference.py --set euler --jobs 8
    python generate_cfl3d_reference.py --set rans  --jobs 8
    python generate_cfl3d_reference.py --derive-only     # no solver run

Needs `tools/cfl3d_seq` (see §1) and the `up3d` env (gmsh ≥ 4.11, numpy ≥ 2,
scipy). Run directories go to `tools/cfl3d_work/` (gitignored); they are
resumable — a run whose `cfl3d.out`, `cfl3d.prt` and `grid_record.csv` are
present is reported `cached` and not repeated unless `--rebuild` is given.

★ Grids are built **serially on the main thread** and only the solver runs are
parallel: gmsh installs a SIGINT handler in `initialize()`, and
`signal.signal` raises outside the main thread. Each `cfl3d_seq` is a
single-threaded process, so `--jobs N` uses N cores.

---

## 6. The 3-D ONERA M6 grid (gates D07 / D10)

`wing3d_otip.py` builds a **multiblock C-H grid with an O-block over the tip**
for the M6, plus the matching 7-block `cfl3d.inp`.

★ **It is NOT in this commit.** The 2-D datasets above are complete and
verified, and are banked here on their own; the 3-D generator lands next,
after a known defect is fixed — it hard-codes one wall spacing for both
equation sets (see “Wall spacing” below), which over-refines the Euler grid by
~600x and leaves the RANS y+ uncomputed rather than derived from y+ and Re.

    python wing3d_otip.py --check-2d --ref ref_2d.npz   # stage 1, the section
    python wing3d_otip.py --out ./m6 --model euler      # grid + deck

    python wing3d_otip.py --check-2d --ref ref_2d.npz   # stage 1, the section
    python wing3d_otip.py --out ./m6 --model euler      # grid + deck

Read the module docstring before changing anything — the topology is the whole
difficulty.

### ★★★ Why not a single block: a blunt tip has nowhere to go

The wing surface is the `k = 1` coordinate surface, and the tip cross-section is
the *area enclosed by* the section contour — in index space that contour is the
LINE `k = 1, j ∈ [jte1, jte2]`. An area bounded by a grid line is not a face of
the block. The ways out are (a) a separate tip block, (b) a collapsed/singular
grid line, or (c) closing the section to zero thickness.

★★ A first attempt took route (c) and was **deleted**. It is worth stating why,
because the failure is instructive: closing the section deforms **the very
geometry the reference exists to measure**, and it is the same mistake this
directory refuses elsewhere, where it declines the CFL3D repo's own airfoil
ordinates precisely so that no geometry difference enters the comparison.
User ruling 2026-09-03: *“M6 机翼不能改变机翼形状, 不可以厚度收到零, 本来就不
应该使用单块 C-H 网格来生成”*. **The M6 tip must not be reshaped.**

### The topology, and the fact that makes it cheap

Route (a), following the two reference implementations this file was written
against — the user's own CFL3D-oriented generator `tools/cgrid/`
(`examples/wing-simple-OTip/wing-onera-m6.py`, `cgrid/wing.py::gen_grid_o_tip`)
and [pyHyp's M6 example](https://github.com/mdolab/pyhyp/tree/main/examples/m6):

| block | dims (M6 reference parameters) | role |
|---|---|---|
| 1 | 61 × 441 × 81 | main C-H, root → tip |
| 2 | 61 × 61 × 17 | tail block — the blunt TE |
| 3 | 49 × 153 × 17 | **O-block over the tip** |
| 4 | 49 × 61 × 17 | O-block in the tip wake void |
| 5 | 61 × 441 × 33 | main block → spanwise far field |
| 6 | 61 × 153 × 17 | tip O-block → far field |
| 7 | 61 × 61 × 17 | wake O-block → far field |

`NBLI = 18`, grid normalised to **root chord = 1**. The outermost main block is
pushed outward by one tip thickness to open a *basin*, and the two O-blocks
Coons-fill it; the wing itself is never touched.

★★ **A blunt TE does not need a blunt-TE mesher.** `cgrid` strips the TE
thickness, meshes the **sharp** section, then puts the thickness back and fills
the gap with one extra block. So the sharp-TE gmsh C-grid already vendored here
(`cgrid_gmsh.CGridGmsh`) is a drop-in for the hard part. The index conventions
already agree — both use `nj = 2(n_foil+n_wake) − 3`, `jTE0 = n_wake`,
`jTE1 = n_wake + 2 n_foil − 2` — which is why the port is clean.

★ **The blunt TE is load-bearing, not cosmetic**: it supplies one of the tip
O-block's four Coons edges (`n_tail` points across the TE base). A sharp TE
collapses that edge. Geometry basis is therefore the real M6 (CST section,
`rel_thick` 0.09779, **`rel_tail` 0.00141**, flat tip) — user ruling
2026-09-03. Our own pyfp3d M6 (sharp TE, `tip_cap="round"`) differs, and that
difference is **ours** and already documented, not the reference's.

★ gmsh does the **meshing** (the 2-D transfinite fill, where its value is);
the 3-D assembly is transfinite interpolation and distribution algebra in
numpy, the same split `cgrid_gmsh.py` already uses here.

### ★★★ Verified against an independent implementation

`tools/cgrid` is not just a model — it is the **verification baseline**, being a
separate implementation of the same topology. Agreement is therefore evidence,
not self-consistency. Measured:

| check | result |
|---|---|
| block count / dims / total points | **identical** — 7 blocks, 3,530,151 points |
| all six BC tables (I0…KDIM), field by field | **identical**, 0 mismatches |
| 1-to-1 blocking, 18 pairs / 36 entries | **identical** |
| wing half-thickness at all 7 experimental stations | ≤ **8.1e-07** |
| **tip half-thickness** | **0.027478** = full (0.09779/2 × 0.4529/0.8059) ⇒ **the wing is not deformed** |
| **wetted area, same deck on both grids** | 2.3874353 vs 2.3874508 = **6.5e-06 relative** |
| 2-D section: wall curve vs `cgrid`'s | **9.19e-06** chord (different point distributions) |
| negative cells | 0 |

★★ And a stronger check than intended: running **your grid with my deck** and
**your grid with your deck** gave **bit-identical** output
(`|R| 0.13133101E-01`, `cl 0.86412949E+02`) ⇒ the two decks are functionally the
same file, verified by the solver rather than by diffing.

### Wall spacing: sized per equation set

★★★ **Euler must NOT have y+ clustering; RANS with SA or SST must, in 2-D and
in 3-D alike.**  How each side is built:

| | 2-D (committed) | 3-D |
|---|---|---|
| Euler | `h1` set DIRECTLY: 2.0e-3 / 1.414e-3 / 1.0e-3 chord per rung. No y+ columns are even written | ★ **defective** — see above |
| RANS | `y_plus = 1.0` on every rung, `h1` computed from **each case's own** chord Reynolds number (`cf = 0.0576 Re^-0.2`, `y/c = y+/(Re √(cf/2))`) ⇒ h1 4.35e-06 … 8.73e-06 as Re varies 3e6 → 6.5e6. **Achieved y+_avg 0.90–1.14**, y+_max 1.27–2.70 | ★ **defective** — see above |

★ The 2-D numbers above are read back out of the committed `forces.csv`
(`h1_wall`, `yplus_avg`, `yplus_max`), not asserted: the boundary layer is
resolved on every RANS rung, and the Euler rungs carry no wall clustering at
all.

### Convergence caliber — Euler runs, RANS does not

Measured on the 0.59 M coarse level, M 0.8395 / α 3.06 (Euler unless stated):

| startup | outcome |
|---|---|
| CFL 1.0 → 5.0, `NITFO` 0 | diverged |
| CFL 0.8 → 5.0, `NITFO` 0 | diverged |
| **CFL 0.5 → 2.0, `NITFO` 500** | **ok**, cl 0.341995, \|R\| 1.4e-07 |
| **CFL 0.2 → 5.0, `NITFO` 1000** | **ok**, cl 0.318199, \|R\| 8.3e-11, drop 4.33 |
| **CFL 0.1 → 10, `NITFO` 1000** | **ok**, cl 0.317064, \|R\| 1.6e-09 |
| RANS (`ivisc` 7), any of the above | **diverged in ~1.2 s**, \|R\| → 1.9e+17 |

★ Euler is **startup-independent** (0.318199 vs 0.317064 = 0.35 %), which is the
property a reference needs.

★★★ **RANS diverges on the reference grid with the reference's own deck too**, so
it is not a regression introduced here. Note what that revealed: the reference
example *generates* its grid and deck and never runs the solver, so “the
reference deck is known to run” was an **unverified assumption of mine** — the
deck is verified as a faithful port, not as a converging RANS setup.
⇒ **D07 (Euler) is unblocked; D10 (RANS) needs the RANS setup diagnosed first.**

