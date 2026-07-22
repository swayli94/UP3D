# GV3.3 Pre-registration — fuselage body-of-revolution smoke

Written 2026-07-22 before any GV3.3 execution. Binding gate text:
docs/roadmap/track_v.md GV3.3 (2026-07-22, user-directed). User
directive: no gate re-spec — the bands below are the execution-level
operationalization the gate text itself requires.

## Case conditions (binding)

- Geometry/mesh: `cases/meshes/fuselage_bor/coarse.msh` (PRIMARY),
  FuselageParams() defaults (L = 4.0295, r_f = 0.15, 2-diameter ellipsoid
  nose, cone afterbody, tail sphere cap), full-2π revolve in a full
  R_FAR = 25 MAC sphere, groups exactly {wall, farfield}. medium level
  RECORDED if generated.
- Conditions: M∞ = 0.3, α = 0 (non-lifting), Re = 3.0e6 per body length
  (re_chord = 3.0e6 / 4.0295 = 7.443e5 per unit length; RECORDED).
- Forced transition at x_tr = 0.05 of body length (D-TR);
  `build_closed_body_case` with stag_band_frac = 0.05 (A4 LE-band u_e
  discipline at BOTH stagnation bands — item (c), RECORDED).
- Driver chain: `solve_subsonic` (compressible Picard, non-lifting, no
  wake/Kutta anywhere) + V2 `body_source_rhs`, via
  `coupling.py::run_loose_coupling` with
  `make_picard_nonlifting_driver`.
- Numerical settings as in the GV3.1/GV3.2 pre-registration (IBL tol
  1e-9, max_iter 100, warm-started; loose loop ω = 1.0, n_outer_max 10,
  tol_ds 1e-3; eps_diff 0.005 / eps_diff_s 0.02; thread cap 16).

## Assertions (binding bands)

(a) **Azimuthal δ* collapse.** Wall nodes are binned into x-stations
(bin width ≈ 2·h_body over the body). Window: station centers with
(x − x_nose_tip)/L ∈ [0.20, 0.95] (outside both stagnation bands, inside
the tail-cap zone boundary). At EVERY window station:
azimuthal coefficient of variation σ/μ of δ* ≤ **0.15**.
The full (max−min)/μ distribution is RECORDED, not banded.

(b) **Crossflow ≈ 0** (axisymmetric flow on a genuinely 3-D surface):
on the final state, max|B| ≤ 0.05·max|A| AND
max|C_τ2| ≤ 0.05·max|C_τ1| (U-columns 2 and 5 vs 1 and 4).
RECORDED: max|Ψ|, max|DS2|/max|DS1|, max|CF2|/max|CF1|.

## RECORDED (no pass/fail)

(c) Nose/tail stagnation-band seeding: inflow pinned node count, LE-band
node count, q_ref at the pin.
(d) Transpiration on/off effect: max/mean |ΔCp| on the wall between the
k = 0 (inviscid) and final probes; meridian-averaged H(x) profile —
the tail-cone adverse-gradient H rise is expected; an indicated tail
separation is recorded and MASKED (excluded from (a)'s window only if
δ* changes sign there — declared in the VERDICT if used), not chased.
Loose-loop iteration count and ω; IBL residual floor.
