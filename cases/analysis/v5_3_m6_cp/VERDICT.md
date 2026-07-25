# VERDICT — GV5.3: M6 wing direction+magnitude check vs committed Cp

**Outcome: band (b) honest FAIL** (medium binding: 1/5 unmasked stations
improved + pooled RMS 0.1288 → 0.1299 INCREASED) — **band (a) RECORDED
input-limited** (Δcl_KJ −2.20 %, direction DOWN both estimators but under
the A4 2.5 % floor). 0 PASS / 1 FAIL / 17 RECORDED.

Executed 2026-07-25 under the temporary 8-thread session constraint (wall
times flagged non-comparable). Protocol, bands, failure clauses:
[PRE_REGISTRATION.md](PRE_REGISTRATION.md) (committed ba636d9 BEFORE
code; addendum 2026-07-25 #1 committed 36b80f9 BEFORE the re-execution).

## 1. Question and protocol (as pre-registered)

On the ONERA M6 at the committed experimental condition (TEST 2308:
M 0.8395, α 3.06, Re_MAC 11.72e6, forced x_tr/c 0.05), does the loose
VII loop (GV3.1 recipe verbatim ω = 1.0, ≤ 10 outer, tol_ds = 1e-3;
the GV5.0 wing case, tip band z > 0.95·b_semi pinned + ṁ-masked; the P14
transonic FP recipe verbatim, NEWTON_M6_RECIPE imported) move (a) the CL
DOWN from the same-mesh inviscid baseline beyond the A4 floor, and (b)
the wall Cp CLOSER to the committed 7-station experiment than the
same-mesh k = 0 inviscid baseline? Coarse recorded, medium binding.

## 2. Execution narrative and wiring guards

Two executions. The FIRST (this date, same protocol) ran coarse clean
and had the medium k = 0 cold start land on the HALF-CONVERGED M0.70
seed state (cl 0.2157/0.2260): the driver's cold start carried the GV5.0
bridge's M0.5 short-circuit (return the probe seed early when it does
not converge to 1e-10), and on medium @8 threads the M0.70 probe solve
stalled at |R| ~ 1e-6 so the Mach ramp never ran. **W1 fired exactly as
designed** and the level read RECORDED per the recipe-limit clause; the
loose loop itself CONVERGED in 9 outer from that poisoned seed (the k =
1 warm solve jumped back to cl 0.277). Root-caused by a measured
diagnostic (same worktree, same thread count): the P14 ramp from the
very same failed seed converges level-by-level onto the anchored branch
(cl_p 0.27726 / cl_KJ 0.28188, identical to the LOCAL-P14-cache-seeded
state to 9 printed digits) — NOT an 8-thread branch scatter. Addendum
#1 removed the short-circuit (execution mechanics only; the ramp now
runs unconditionally). Coarse results from the first execution were
bit-identical to the re-execution's (the short-circuit never fired on
coarse); the committed results below are the re-execution's, both
levels from one process.

Re-execution wiring guards: **W1 PASS both levels** (k = 0 cl_p 0.2625 /
cl_kj 0.2685 coarse vs anchors 0.2628/0.2688; cl_p 0.2773 / cl_kj 0.2819
medium vs 0.2776/0.2823 — the ΔM = 0.0005 label difference well inside
the 1 % guard). **W2 PASS both levels** (every zone's max-Cp point at
x/c < 0.05 = the LE stagnation point; pooled k = 0 RMS with the chosen
Z/L-sign side mapping 0.95 / 0.63 vs the flipped mapping 2.65 / 2.68
coarse/medium). W3 never fired (all 7 stations extract 37–58 points per
side at both levels).

## 3. Band (a) — ΔCL direction gate (medium binding): RECORDED input-limited

`results/history_{level}.csv`, cl read at k = 0 (inviscid) and the
terminal outer:

| level  | cl_KJ k=0 → terminal | Δcl_KJ (rel)    | Δcl_p (rel)     | read |
|--------|----------------------|-----------------|-----------------|------|
| coarse | 0.2685 → 0.2658      | −0.0028 (−1.03 %) | −0.0036 (−1.35 %) | input-limited RECORDED |
| medium | 0.2819 → 0.2757      | −0.0062 (−2.20 %) | −0.0067 (−2.40 %) | input-limited RECORDED |

The direction is DOWN at both levels with BOTH estimators (the
physically expected viscous decambering), and the two estimators agree
within 0.2–0.3 % (no estimator-disagreement note). But the medium
binding move −2.20 % stays under the A4 2.5 % input floor ⇒ the
pre-registered "smaller move" clause: RECORDED, flagged input-limited —
NOT a PASS, and (per the registration) a smaller move is not a FAIL
either. Note the medium move concentrates in the LAST two outer
iterations (§5): k = 1–8 sit at Δcl ≈ −0.3 %, k = 9–10 carry the rest.

## 4. Band (b) — per-station Cp RMS direction gate (medium binding): **FAIL**

`results/cp_rms_{level}.csv` (RMS of computed-minus-experiment Cp over
each zone's points, both sides pooled; computed Cp from
`section_cp_curve` at the run Mach, interpolated linearly per side):

| eta   | coarse inv → visc     | medium inv → visc     | medium better |
|-------|-----------------------|-----------------------|---------------|
| 0.20  | 0.1402 → 0.1428 (+0.0026) | 0.0988 → 0.0969 (−0.0019) | YES |
| 0.44  | 0.1725 → 0.1750 (+0.0026) | 0.1159 → 0.1160 (+0.0001) | no  |
| 0.65  | 0.2060 → 0.2087 (+0.0027) | 0.1400 → 0.1409 (+0.0009) | no  |
| 0.80  | 0.1822 → 0.1839 (+0.0017) | 0.1174 → 0.1180 (+0.0006) | no  |
| 0.90  | 0.2503 → 0.2530 (+0.0027) | 0.1549 → 0.1583 (+0.0034) | no  |
| pooled (unmasked) | +0.0024         | 0.1288 → 0.1299 (+0.0010) | —   |

**FAIL**: 1/5 unmasked stations improved (band requires ≥ 4/5) and the
pooled RMS INCREASED at the binding level (coarse likewise 0/5 with a
pooled increase). Every per-station |ΔRMS| < 0.05 ⇒ all flagged
input-limited per the pre-registered A4 Cp-scale annotation — the FAIL
is a DIRECTION verdict: the viscous correction does not move the Cp
toward the experiment anywhere except (marginally) the root station;
the magnitude of the adverse shift is small. The tip-masked stations
behave as constructed (η = 0.96: +0.0022/+0.0020, η = 0.99:
+0.0009/+0.0029 coarse/medium — no large-shift anomaly).

The overlays (`results/cp_overlay_{level}.png`) show why: the inviscid
baseline already mismatches the experiment in the inviscid-family way —
LE suction peak too shallow (exp ≈ −1.0…−1.2 at x/c ≈ 0.02–0.04 vs
computed ≈ −0.8) and the shock aft — and the loose viscous correction
at this strength barely perturbs the upper-surface curves (dashed on
solid), the small deviations it adds landing on the wrong side of the
measurement at 4 of 5 stations.

## 5. (c) Convergence and guards (RECORDED)

| level  | converged | outer | wall (8 thr) | fp calls (cont / stall-acc) | IBL floor | mdot 1st→last |
|--------|-----------|-------|--------------|------------------------------|-----------|----------------|
| coarse | no (cap)  | 10    | 1721 s       | 21 (1 / 10)                  | 5.37e-6   | 0.006 → 0.289  |
| medium | no (cap)  | 10    | 12479 s      | 22 (1 / 10)                  | 1.88e-6   | 0.005 → 0.106  |

- Neither level contracts to tol_ds = 1e-3 in ≤ 10 outer (the GV5.0/GV5.2
  loose-loop signature): coarse ds_change_rel oscillates 0.09–0.67 with
  mdot climbing to 0.33; medium sits in an extremely tight limit cycle
  k = 1–8 (ds_change_rel 0.010–0.025, cl_kj pinned to 4 digits at
  0.2810) then **a late separation-patch event at k = 9** (ds_max
  0.0021 → 0.0039, ds_change_rel 0.448, mdot 0.007 → 0.085 → 0.106,
  cl_kj 0.2810 → 0.2757) — bounded, no runaway (the loud-fail guard
  never fired), the terminal state finite and physical.
- The IBL hits its 100-iteration cap at every outer at both levels; the
  IBL residual floor is small (1.9e-6 medium).
- The pre-registered FP rescue chain was load-bearing at both levels:
  10/21 (coarse) and 10/22 (medium) FP calls stall-accepted at the
  honesty-guarded plateaus (one continuation call each; strict remained
  the first choice at every call; per-attempt path + accept_reason in
  run.log). The medium k = 0 cold start: the M0.70 probe seed
  non-converged (logged), the unconditional ramp converged strictly
  (addendum #1 path).
- Crossflow: max|B|/max|A| = 0.026 (coarse) / 0.055 (medium) — the
  GV5.0-class crossflow magnitudes, no 3-D blow-up.

## 6. Reading

1. **Band (b): honest FAIL** — the loose VII loop does NOT move the M6
   transonic Cp toward the committed experiment at the unmasked stations
   (1/5 medium, 0/5 coarse; pooled worse at both levels). This is the
   3-D counterpart of the GV5.2 RAE2822 finding: at transonic conditions
   the loose displacement-thickness feedback is too weak to repair the
   inviscid-family mismatch (shallow LE suction, aft shock), and what
   little it adds lands on the wrong side of the measurement.
2. **Band (a): RECORDED input-limited** — the CL direction is
   physically correct (DOWN, both estimators, both levels) but the
   binding move −2.20 % sits under the A4 input floor; most of it
   arrives with the k = 9–10 separation-patch event, i.e., the move is
   still growing when the 10-outer recipe cap stops the loop.
3. The 3-D loose loop does not converge at M 0.8395 either (both levels
   capped), extending the GV5.0 (M 0.5) non-convergence to the
   transonic condition — though the medium limit cycle is tight and the
   late separation event bounded, so the terminal states are usable
   reads, not garbage.
4. Programmatic: the W1 wiring guard earned its keep (caught a real
   driver bug at the first execution; root-caused and fixed under
   addendum #1 with the diagnostic committed to the record). The gate
   data argues the same conclusion as GV5.2: further reads belong to
   the tight/augmented coupling, not to more loose-loop tuning. No
   follow-up opened by this gate (user adjudication).

## 7. Artifacts

- `cases/analysis/v5_3_m6_cp/`: PRE_REGISTRATION.md (+ addendum #1),
  run.py, VERDICT.md, `results/summary.csv`,
  `results/history_{coarse,medium}.csv`,
  `results/cp_stations_{k0,final}_{coarse,medium}.csv`,
  `results/cp_rms_{coarse,medium}.csv`,
  `results/cp_overlay_{coarse,medium}.png`,
  `results/cl_history_{coarse,medium}.png`. (run.log stays in the
  working tree — `*.log` is gitignored; the per-attempt fp-chain
  narrative it carries is summarized in §5 and the counts in
  `results/summary.csv`.)
- No library or test changes in this gate (the runner + the committed
  recipes only); the full-suite baseline therefore re-quotes the GV5.2
  count (642 + 25 + 2), re-measured at the ledger commit.
