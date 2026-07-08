# P5 medium gate — root cause RE-DIAGNOSIS (T1–T4): single-station Kutta closure

**Date:** 2026-07-08 (second re-diagnosis) **Verdict:** ✅ the medium-gate
outboard-TE M>2 cluster is a **single-station Kutta-closure failure**, NOT a
trailing-edge P1 discretization singularity. The earlier "TE discretization
singularity (G1.6 family), NOT a wake/Kutta change" conclusion (see
[INVESTIGATION_gamma_smoothing.md](INVESTIGATION_gamma_smoothing.md) and the
prior roadmap/demo_report/agent-rules text) was **wrong** and is superseded by
this record.

This is the lab notebook for the T1–T4 tests. One-line conclusions live in
[docs/roadmap.md](../../../docs/roadmap.md) P5, [docs/demo_report.md](../../../docs/demo_report.md)
P5, and [docs/agent-rules.md](../../../docs/agent-rules.md).

## The failure

Medium (350.7k) FAILS `physical`: M_max 5.204, 8 floored / 4 limited, V6 1.73%.
Localization (unchanged from the first re-diagnosis, `diagnose_medium.py`):
26 M>2 cells = **18 on the wing at the outboard TE** (z/b 0.80–0.81) + **8 at
the far-field sphere** (beyond the tip). The single M_max=5.20 cell is far-field.

## Why the earlier "TE singularity" reading was believable but wrong

The cluster co-locates with the steepest spanwise Γ roll-off, sits on
well-shaped cells, and *sharpens* under refinement (coarse 1.47 → medium 5.20).
Those facts are all real — but they are equally consistent with a
**circulation deficit at one station** (which also sharpens with h, because a
finer mesh resolves a sharper overspeed around the same Γ error). The
refutation below shows Γ, not h, is the lever.

## Tests (T1–T4) — all reproducible from the local solution cache + mesh

| # | Test | Result | Reads on the hypotheses |
|---|---|---|---|
| **T1** | straddle census: any tet with both a master and a slave wake node? (cache-only) | **0** of 350 718 (and 0/55 531 coarse) | Rules out "Γ jump compressed into one P1 cell". Cut topology clean. |
| **T2** | fixed-Γ (cached Γ) warm-start re-solve, omega_rho 1.0 vs 0.5, 1500 iters (medium) | omega_rho=0.5 cut the eval drho ~4× (0.18→0.047) but the defect band stayed **bit-identical** (M_max 3.099, 18 cells, 8/4 floored/limited). coarse: cached is already an exact fixed point. | **Refutes H3** (under-convergence / frozen transient): better density convergence at the *same* Γ does nothing. The defect is a genuine fixed point of that Γ. |
| **T3** | single-station scan: cached Γ everywhere but **st133 (z/b 0.801)** set to {0.050, 0.059(=its own Kutta target), 0.065}, fixed-Γ warm-start, omega_rho=0.5 | band M_max **3.099 → 2.574 → 1.155 → 1.188**; band M>2 cells **18 → 2 → 0 → 0**; floored/limited **8/4 → 1/1 → 0/1 → 0/1**; monotone in Γ₁₃₃ | **Decisive. Confirms H7, refutes H1.** Same mesh/h/TE-elements — only ONE station's Γ changed — collapses the cluster. A 1/h discretization singularity cannot depend on the circulation value. |
| **T3b** | set ALL stations to their cached-measured Kutta targets (changes only st133: 0.0431→0.0592), fixed-Γ warm-start | wing M>2 **18 → 1**, floored **8 → 0**, max\|F\| 0.016 → 0.0028 | The targets read off the (hot) cached state are already ~right; only station 133's circulation lagged. 165/166 stations were already correct. |
| **T4** | **full fix** warm-start: far-field taper ON + st133 = target | M_max **5.204 → 2.048**, floored/limited **8/4 → 0/0**, wing M>2 **18 → 1**, far-field M>2 **8 → 0** | Both clusters heal together; `physical` gate passes. The 1 residual wing cell (M≈2.0, tip TE corner) is the mild bounded H1 contribution, within the gate; it clears as Γ₁₃₃ iterates to its self-consistent ~0.063. |

## Root cause (mechanism)

The P5 continuation's frozen-Γ eval density solve runs at omega_rho=1.0, where
it does **not** converge — it floors at a drho≈0.18 limit cycle (T2). Each
per-station Kutta target is therefore read off a noisy, unconverged density
state. The per-station secant closes fine at the 165 benign stations, but
**stalls at the single steepest-Γ-roll-off station** (st133, z/b 0.801): that
station has the largest Γ change to make and the noisiest local field, and the
`max_gamma_evals=10` budget runs out with it still **32% under-circulated**
(Γ 0.0431 vs its own target 0.0592, self-consistent ≈0.063; |F|/Γ = 37% there
vs ≤5% everywhere else). The circulation deficit drives a real overspeed around
the sharp TE at that station; the artificial-density floor + speed limiter then
freeze it into a spurious M≈3 fixed point (the mechanism
`physics/isentropic.py::q2_at_mach` documents: "positivity guards stabilize
spurious converged states"). The far-field 8-cell cluster is an independent,
already-confirmed artifact of the span-uniform 2D-vortex far field prescribing
a jump across its branch ray *beyond the tip*, where no wake cut exists.

The probe-assignment degeneracy on the swept unstructured TE (adjacent stations
share Kutta probe nodes / probes fall off the station's own z-plane — 35 upper +
41 lower shared of 166 stations; st133/134 share their upper probe) is the
**latent reason st133 specifically is the one that stalls** (its target is
coupled to a neighbor, so the independent per-station secant is slower there).
It is NOT the direct cause — the target it reads is still ~right — so it was
left as a known-robustness item rather than fixed for this gate.

## The fix (recipe-only, P4 2.5D path bit-identical)

**First attempt — REFUTED by the from-scratch run.** "omega_rho=0.5 inside
the continuation's eval solves + max_gamma_evals=16" looked sufficient from
the warm-start tests, but a from-scratch medium continuation with it DIVERGED
(M=0.84 level: 16 evals, M_max 29.3, 33/40 floored/limited — worse than the
10-eval baseline's bounded-but-wrong M 5.2). Lesson: the active per-station
secant at the top Mach level is itself unstable on the 3D mesh; its 10-eval
budget was *early-stopping regularization*. Under-relaxing the density inside
the secant loop, or giving the secant more evals, feeds the overshoot.

**Landed fix** (what the T3b/T4 operation actually was — closure WITHOUT the
secant):

1. `farfield_spanwise_gamma=True` — 3D spanwise-tapered vortex far field
   (`constraints/dirichlet.py::farfield_dirichlet`), Γ(z) per station, 0 at/
   beyond the sheet tip. Removes the 8 far-field cells.
2. `n_kutta_polish=4` (`solve/continuation.py::solve_transonic_lifting`) —
   a fixed-Γ Kutta-closure POLISH after the Mach continuation: set Γ to the
   measured `kutta_targets(phi)`, re-solve fixed-Γ with under-relaxed density
   (`omega_rho_polish=0.5`), repeat. Secant-free (a damped fixed-point map:
   no slope estimate, no overshoot), contractive — measured |F| halves per
   step (1.9e-3 → 5.8e-4 over 4 steps) with floored/limited 0/0 throughout.

Both parameters default off (0 / False), so every P2/P3/P4 and 2.5D path is
bit-identical (full default suite re-verified: 140 + 4 + 2). From-scratch
medium verification with the landed fix: M_max 1.995, 0 floored / 0 limited,
Kutta |F| 5.8e-4; demo 16/16 PASS.

## Contrast with the refuted gamma-smoothing route

The gamma-smoothing route (INVESTIGATION_gamma_smoothing.md) was correctly
refuted: OVERRIDING the Γ profile with an externally-imposed smooth curve did
not heal the defect (fixed-Γ, |dΓ/dz| flattened 7×, still M~5.3). That is a
DIFFERENT operation from what fixes it here: driving the ONE stuck station to
**its own Kutta-consistent target** (the value the solver itself measures),
i.e. *closing the constraint*, not smoothing over it. Smoothing moved Γ away
from the self-consistent value; closure moves it toward it. This is why the
earlier "dead route" verdict on smoothing was right, yet the cause still turned
out to be the Kutta condition.

## Meta-lessons

1. The first re-diagnosis over-trusted three circumstantial signals
   (co-location, good cell shape, refinement-sharpening) and reached for the
   exotic explanation (TE singularity) before running the cheap decisive test.
   T3 — vary the one suspected input (Γ at the hot station) and watch the
   output — is the falsification that should have come first. "Refinement
   sharpens it" does not imply "h is the cause": a fixed Γ error also sharpens
   under refinement.
2. Warm-start evidence does not validate a from-scratch recipe. Every T-test
   ran fixed-Γ from the cached state (bypassing the secant); the first fix
   attempt generalized "omega_rho=0.5 helps" into the ACTIVE secant loop and
   diverged from scratch. The landed polish reproduces the *exact operation*
   the successful tests performed (fixed-Γ target application), not a
   paraphrase of it.
