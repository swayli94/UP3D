# VERDICT — GV5.5: TE-band (B, δ) formulation — breaking the IBL floor

**2 PASS / 1 FAIL / 9 RECORDED — the V1 TE-outflow row replacement does NOT
break the floor: the binding m2 lands 4–5 orders of magnitude ABOVE it (the
pre-registered "worse" clause); the flag stays default-OFF.**

Binding text: `PRE_REGISTRATION.md` (committed `db2af1e`, before the first
code change). Implementation commit `0ea1b11` (rebased onto `origin/main`
157d7c6). Executed 2026-07-24 under the **temporary 8-thread session
constraint** (runner default 16; wall times non-comparable).

## 1. Question and protocol

Question (pre-registered): does the TE-band (B, δ) formulation-level
treatment — route (a), variant **V1 = TE-outflow row replacement**
(first-order extrapolation: row 6i+0 `R = δ_i − δ_up`, row 6i+2
`R = H_i − H_up`, exact Jacobian rows, CSR in-pattern by construction,
default-OFF flag `te_extrapolate`) — break the IBL floor that basin hunting
(GV5.1b/1c/1d) exhausted?

Protocol per level (coarse = crash-stop, medium = binding): rebuild the
loose state and amended GV5.1 seeds verbatim (wiring guard
`|dcl_k0| ≤ 1e-8`, PASS both legs); V0 control (flag OFF, diagnostic Q7
pseudo-time protocol) vs the committed floors; band (a) live FD on the
variant system at the seed; V1 solve (flag ON) from the same seed;
GV5.1b tight polish on the variant system as the secondary read; guards
(plate H bands flag-ON, loose-loop smoke flag-ON, tight fleet + full suite
flag-OFF — the latter run separately, §5).

Metrics (pre-registered): **m1** = the variant system's own floor
(recorded); **m2** = the ORIGINAL-system residual at the V1 terminal state
(BINDING — not gameable by row replacement). PASS = m2 ≤ 0.5× the
committed floor on BOTH levels with guards green; 0.5–0.9× partial;
≥ 0.9× no-move; worse = RECORDED, variant stays default-OFF.

## 2. Seeds and V0 control

| level | committed floor | V0 floor | control_rel | clause |
|---|---|---|---|---|
| coarse | 3.154e-6 | 3.1537e-6 | 1.06e-4 | OK (roundoff) |
| medium | 1.712e-6 | 1.8238e-6 | 6.53e-2 | **scatter clause fired** |

The medium 8-thread scatter struck again: today's loose regen landed on the
**same 4th fixed point** as the GV5.1c/1d 8-thread runs (cl 0.28245999;
V0 floor 1.8238e-6, matching GV5.1c's unperturbed F_BL 1.824e-6 = 1.07×
floor). Per the pre-registered clause, the medium (b) read uses the seed's
own flag-OFF floor as `floor_ref = 1.8238e-6`; coarse is bit-close to
committed. Wiring guards PASS both legs (|dcl_k0| ≤ 1e-8).

## 3. Band (a): implementation exactness, live on the variant system

Jacobian-action FD (central, eps 1e-6, 4 random directions, at the amended
seed, flag ON): max rel err **1.79e-7 (coarse)** / **1.09e-8 (medium)**,
both < 1e-5 → **PASS ×2**. The replaced rows are analytically exact
(tests/test_v5_te_outflow.py gates the same on the unit fleet, 9 tests).

## 4. Band (b): m1 / m2 — the binding read

| level | variant residual at seed | m1 (variant floor) | m2 (orig. system @ V1) | ratio vs floor_ref | clause |
|---|---|---|---|---|---|
| coarse | 9.821 | 4.703 (stall) | 1.752e-2 | **5554×** | worse |
| medium | 4.846 | 9.863e-1 (stall) | 4.487e-1 | **245998×** | worse |

Both legs: the pseudo-time protocol (the diagnostic Q7 recipe that takes
the flag-OFF system to its floor in 11–21 iterations from the same seed)
descends the variant residual only ~2–5×, then **stalls with every step
rejected** (cfl driven to the 1e-3 floor; coarse frozen at it=11–17,
medium grinds 75 iterations to 0.986 and freezes). Neither V1 solve
converged.

Secondary read (GV5.1b tight polish on the variant system, from the seed):
coarse terminal F_max 7.32e-5 (23.8× the committed GV5.1b final 3.07e-6);
medium F_max **3.98** (the polish diverges; committed 1.708e-6). No floor
break through the polish either.

**Anatomy** (`results/residual_anatomy_{level}.csv`, the diagnosis's column
naming): the m2 peak is NOT at the TE — it sits at **x_c ≈ 0.027 (the LE
suction zone), rows F_B/F_Psi** (coarse node 53: F_B −1.75e-2; medium node
310: F_Psi −0.449, F_B −0.271, both sides symmetric). The TE rows
themselves are quiet in the top ranks. The V1 terminal state is globally
corrupted, worst at the LE — consistent with a variant system that has no
solution anywhere near the amended seed.

Reading (recorded, not a new gate): the seed is O(5–10) away from the
variant system in residual norm precisely because the replaced rows measure
the natural TE jump (δ_TE − δ_up, H_TE − H_up ≠ 0 at the loose-converged
state). The hard first-order extrapolation constraint is inconsistent with
the coupled momentum/closure system at strip resolution — the pseudo-time
finds no descent path and the damage migrates to the stiffest region (LE).

**Band (b) verdict: FAIL — the pre-registered "worse" clause, both levels
(medium binding).** The V1 row replacement does not break the floor; it
destroys solvability 4–5 decades above it.

## 5. Guards (c)

| guard | coarse | medium |
|---|---|---|
| plate H flag-ON, lam ∈ [2.55, 2.75] | [2.606, 2.687] PASS | (same run) PASS |
| plate H flag-ON, turb ∈ [1.2, 2.0] | [1.509, 1.872] PASS | (same run) PASS |
| loose smoke flag-ON | **FAIL** — hit the 10-outer cap, converged=False, cl_rel = 2.62% (> 2.5%) | PASS — 3 outer, converged, cl_rel = 2.49% (marginal, 0.01 pt under the band edge) |
| band (a) FD | PASS | PASS |

The flag-ON loose loop on the 4th-fixed-point medium seed converges in 3
outer with cl shifted 2.49% from the flag-OFF cl — at the A4 band edge;
the coarse leg fails the same band at 2.62% and does not converge within
the cap. The TE treatment moves the loose-loop fixed point by ~2.5–2.6%
in cl — physically noticeable, one more reason the variant stays OFF.

Tight fleet + full suite flag-OFF: `tests/test_v5_te_outflow.py` (9) PASS;
`tight fleet + V1/V3/V5 regression 71 passed`; `test_v0_freestream` 6
passed; full-suite baseline re-measured post-execution (§7) — flag OFF is
bit-identical to legacy by construction (test 1 of the 9) and by the V0
control above.

## 6. Interpretation and consequences

1. **The floor stands.** Basin hunting (GV5.1b/1c/1d) showed no quadratic
   basin above or adjacent to the floor; GV5.5 now shows the first
   formulation-level TE treatment does not move the floor either — it
   breaks the system's solvability instead. The (B, δ) TE-band stiffness
   diagnosed in `v5_ibl_floor` is not bypassed by the natural-outflow row
   replacement at this resolution.
2. **Where the damage goes:** not the TE (the replaced rows are satisfied
   by construction at convergence attempts) but the LE suction zone,
   rows F_B/F_Psi — the stiffest coupled block. Any follow-up formulation
   work must watch the LE, not only the TE band.
3. **Escalation ladder stays CLOSED** (pre-registered): the upwind
   boundary-flux (a)-variant and closure regularization (b) remain
   registered-not-opened; opening either is the user's adjudication call.
4. The medium 8-thread scatter (4th fixed point, cl 0.28245999) reproduced
   for the third time (GV5.1c/1d/GV5.5) — the GV5.1 §4 caveat stands;
   the scatter clause handled it exactly as pre-registered.
5. The flag `te_extrapolate` stays **default OFF** everywhere
   (`CouplingConfig`, `IBL3Solver`, v1 runner); legacy paths are
   bit-identical.

## 7. Artifacts

- `PRE_REGISTRATION.md` (committed `db2af1e`, before first code change)
- `run.py` — the protocol runner (this directory)
- `results/summary.csv` — V0 control + V1 m1/m2 + tight reads, both levels
- `results/floor_probe_{coarse,medium}.csv` — V0/V1 pseudo-time residual
  trajectories
- `results/residual_anatomy_{coarse,medium}.csv` — original-system per-node
  residuals at the V1 terminal (diagnosis column naming)
- `results/guards.csv` — plate H / loose smoke / FD verdicts
- `run.log` — the full execution log (V0/V1/tight verbose; uncommitted,
  regenerable by re-running `run.py`)
- Code: `pyfp3d/viscous/ibl3.py` (flag + row replacement + guards),
  `pyfp3d/viscous/coupling.py` (`CouplingConfig.te_extrapolate`,
  `te_outflow_pairs`), `bench/studies/v1_ibl3_standalone/run.py`
  (`run_fe(..., te=)` plate helper), `tests/test_v5_te_outflow.py` (9 tests)
- Suite baseline: filled post-execution in `docs/agent-rules.md` Baseline /
  `docs/overview.md` regression-baseline / `PROJECT_STRUCTURE.md` Default
  suite (8-thread constraint quoted in all three)

Executed: 2026-07-24, 8 threads (temporary session constraint), legs
~6 min (coarse) / ~12 min (medium) including the per-leg loose regen,
V0/V1 solves, tight polish and guards.
