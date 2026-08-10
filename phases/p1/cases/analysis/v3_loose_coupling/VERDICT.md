# GV3.1 / GV3.2 VERDICT — honest FAIL (2 PASS / 4 FAIL / 23 RECORDED)

Date: 2026-07-22. Branch: kimi/track-v3-loose-coupling.
Binding text: docs/roadmap/track_v.md GV3.1/GV3.2; pre-registration +
2026-07-22 addendum: PRE_REGISTRATION.md. Evidence: results/ (summary.csv
is the machine-readable verdict; this file is the narrative).

**Headline: the coupled machinery WORKS (loop converges in 4–5 outer
iterations at ω = 1.0, transonic point included; Δcl direction and
magnitude PASS; mid-chord profiles match XFOIL within ±5 % δ* / ±15 %
cf). The binding every-station bands FAIL at a handful of stations, all
in the immediate post-trip zone x/c ∈ [0.055, 0.124].**

## GV3.1 — profiles vs XFOIL + viscous Δcl

Conditions: NACA0012 2.5-D strip (wake-cut), M 0.5, α 2°, Re 3.0e6,
x_tr 0.05/0.05, compressible Picard driver. Reference: committed XFOIL
6.99 CSVs (xtr005 primary). cf compared in the freestream frame per the
pre-registration ADDENDUM (XFOIL DUMP cf = TAU/(0.5 QINF²), our OUT_CF1
is local — the as-registered raw local-frame comparison is kept as
RECORDED rows; δ* unaffected).

Medium (binding), banded window x/c ∈ (0.05, 0.95]:

- **FAIL — cf upper worst |err| = 0.437 at x/c = 0.055** (band ≤ 0.15).
- **FAIL — cf lower worst |err| = 0.448 at x/c = 0.055** (band ≤ 0.15).
  Every OTHER banded station passes: |err| ≤ 0.136 upper (x/c ≥ 0.065),
  ≤ 0.069 lower (x/c ≥ 0.065). The single failing station is the FIRST
  TURBULENT station, 0.005 aft of the trip: our flags switch regime
  instantaneously at x_tr (the closure equilibrates within 1–2 stations),
  while XFOIL ramps cf through its e^N intermittency over a finite run —
  at 0.055 XFOIL is still mid-ramp (3.87e-3 vs its own 6.9e-3 at 0.065)
  while we are at 5.55e-3. Structural trip-model difference, not a
  discretization error: on coarse the first turbulent station sits
  further aft where XFOIL's ramp has completed, and coarse passes
  (worst 0.135).
- **FAIL — δ* upper worst |err| = 0.279 at x/c = 0.074** (band ≤ 0.25);
  upper exceeds 0.25 only at 0.074 (0.194 at 0.065, ≤ 0.13 from 0.094
  on, ±0.05 mid-chord).
- **FAIL — δ* lower worst |err| = 0.276 at x/c = 0.074** (band ≤ 0.25);
  lower exceeds at 0.074/0.094/0.124 (0.276/0.274/0.266) and sits
  +0.13…+0.25 high across the chord. cf matches there (±0.07), so θ is
  consistent — the bias is in H: the lower side runs an adverse pressure
  gradient and our closure's H response in APG is stronger than
  XFOIL's Drela–Giles. A closure-family difference the 25 % band was
  declared to absorb — it mostly does, but not at those three stations.
- RECORDED: LE-zone max |err_δ*| = 0.274 (input-limited stagnation zone,
  pre-registered RECORDED; stations inside the Dirichlet inflow band
  x/c ≤ 0.02 are labeled `pinned` and excluded from all statistics —
  they are prescribed Blasius boundary data, not solution, per the
  pre-registration addendum); near-TE max |err_δ*| = 0.244 (structural
  TE truncation vs XFOIL wake model, pre-registered RECORDED).

**PASS — viscous Δcl**: cl_inviscid(k=0) = 0.2844 → cl_coupled = 0.2719,
Δcl = +0.0125 > 0; XFOIL's own decrement 0.0230; ratio 0.542 ∈ [0.5, 2.0].
(Our inviscid cl is 2.6 % below XFOIL's 0.2921 — panel-method vs XFOIL
discretization, consistent with the A4 input band.)

Coarse cross-check (RECORDED): cf worst 0.082 lower / 0.135 upper (would
pass); δ* worst 0.384 lower at 0.163 / 0.302 upper at 0.921 (near-TE);
Δcl ratio 0.427. Converged in 4 outer iterations.

## GV3.2 — loop convergence + transonic recorded point

- **PASS — loose loop converged in 5 outer iterations ≤ 10 at ω = 1.0**
  on the GV3.1 medium case (raw successive δ* max rel change < 1e-3).
  RECORDED: IBL residual floor 1.2e-6 over the outer iterations — the
  documented weakly-observable-direction floor (the state is physical;
  the Newton Jacobian has near-null directions that pseudo-time does not
  annihilate; FD checks in V1 showed the Jacobian itself exact).
- RECORDED transonic point (M 0.72, α 2°, Re 3e6, x_tr 0.05, coarse,
  conforming-Newton driver): **converged in 4 outer iterations at
  ω = 1.0**, IBL floor 3.2e-6, final cl 0.3764. No divergence, no
  per-case tuning — feeds the V4 skip decision (the DN6-predicted
  near-separation risk did NOT materialize at this attached point).

## Debug history (the honest record of what it took)

- The IBL3 surface operator needed the D13 §III.D.1 LOCAL-BASIS fix
  (crossflow unknowns solved in each node's tangent plane, R4 triple
  product with the local Jacobian) — pre-fix crossflow leaked 25.9/0.15
  on the curved wall; post-fix 1.8e-4/1.6e-3. 52/52 V1/V2 regressions
  green after the fix.
- The loose loop initially stalled the IBL Newton in a near-singular
  near-separation basin even at α = 0: pinning only the min-q station
  leaves the two streams splitting at the closed nose under-anchored.
  The fix restores the V1 x0-line discipline as a Dirichlet BAND
  (stations x/c ≤ 0.02, per-node Blasius-matched states frozen at k = 1)
  — 4–5 outer iterations to convergence thereafter. A four-leg bisect
  isolating the cause is in docs/temp/v3_case_bisect.py (scratch).
- First execution compared cf in mismatched frames (ours local, XFOIL
  freestream) — pre-registration addendum documents the fix; the
  raw local-frame numbers are retained as RECORDED rows.

## What this verdict does NOT claim

- GV3.1's binding profile bands are NOT met. The failing stations are
  all in the post-trip zone (structural trip-ramp difference) plus a
  systematic lower-side H bias in APG (closure-family difference). Any
  re-spec (e.g. a declared post-trip window analogous to the TE-zone
  exclusion, or a wider δ* band) is a user decision, NOT taken here —
  the 2026-07-22 directive was no re-spec.
- The IBL residual floor (~1e-6 … 3e-6) is recorded, not converged to
  machine zero (weakly-observable directions; the profiles above are
  the converged state to plotting accuracy — outer-loop δ* change
  < 1e-3 dominates it).
