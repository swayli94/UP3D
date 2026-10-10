"""Is a converged solve's ANSWER usable? -- the check the project did not have.

Binding text: phases/p5/docs/dev_phase_five/20260824-0500-r23-prereg.md.

★★★ Why this exists. Measured 2026-08-24 (R22): at medium / alpha 1.25 / seed 0,
`upwind_c = 1.10` converged to |R| = 2.28e-13 with ZERO clamps, and its Gamma was
0.0227 against the C = 1.5 reference's 0.1738 -- a factor 7.6 -- with cl_p 0.0403
against 0.3413 and x_shock 0.657 outside the committed 0.62 +- 0.03 band. It is a
genuine root of the discrete system with almost no lift, and it passed EVERY usability
check the project had: `converged`, zero clamps, residual at 1e-13. GS1.4's
clamp-not-silent contract cannot see this class, and the capability boundary's anomaly 1
said so in advance -- "it reports no error ... what is needed is an answer anchor".

★★ WHAT THIS IS AND IS NOT.
  IS:     outlier detection against the candidate set's own robust consensus.
  IS NOT: correctness certification. The consensus is not a truth; if every leg in a
          set were wrong together, nothing here would notice. A spurious root whose
          cl_p happens to land near the consensus also passes.

★★ TWO LIMITATIONS measured on first use (R23), stated here rather than discovered later:
  1. With only TWO legal legs the median IS the midpoint, so a spurious root drags its own
     reference halfway. Measured: on R22's set the consensus came out 0.190803, pulled
     from 0.341340 by the very leg being tested, and the good leg's own ratio rose to
     1.79 against a 3.0 threshold. It still caught the outlier, but the margin is thin
     and shrinks as the legal set shrinks. With 1 legal leg it cannot work at all.
  2. `assess_set`'s per-axis spread POOLS every other axis. Measured: "seed 5 spread
     42.61%" mixes coarse and medium legs, so it is NOT a (c) reading. A spread by one
     axis is only meaningful with the others held fixed -- the caller must do that.

★ ANCHOR_RATIO_MAX is a CALIBRATION, not a guarantee -- the same status this project
records for the EW forcing, the taper r_c and the descent10 threshold. Derived, not
picked: across the nine converged, zero-clamp legs in `bench/gate_results/
m1_gate_default.csv`, cl_p spans 0.289614 .. 0.424730, a maximum LEGITIMATE ratio of
1.47 (1.18 within medium alone), while the spurious root above sits at 8.48x the medium
median. The threshold below sits in that 5.8x gap. Re-derive it if that evidence moves.
"""
from __future__ import annotations

import numpy as np

#: ★ a CALIBRATION (see the module docstring for its derivation and its evidence)
ANCHOR_RATIO_MAX = 3.0
#: the two clamp counters GS1.4 keys on
_CLAMP_KEYS = ("n_limited", "n_floored")


def clamped_ever(clamp_history) -> bool:
    """★ R15: the scalars n_limited / n_floored are FINAL-STEP counts, not cumulative,
    so "0/0" does not mean "never clamped". Measured: a leg reporting 0/0 had clamped on
    steps 9 to 26. This reads the history instead."""
    if clamp_history is None:
        return False
    h = np.atleast_2d(np.asarray(clamp_history, float))
    return bool(h.size and (h > 0).any())


def assess(cl_p, *, converged, n_limited=0, n_floored=0, clamp_history=None,
           consensus=None, ratio_max=ANCHOR_RATIO_MAX):
    """One leg. Returns a dict; `usable` is True only if all three checks pass.

    `consensus` is the candidate set's robust reference (see assess_set). With None the
    anchor check is SKIPPED and `anchor_ratio` is None -- and `reason` says so, because
    silently passing an unanchored leg is how C = 1.10 got in.
    """
    clamped = bool(n_limited) or bool(n_floored)
    rec = dict(cl_p=(float(cl_p) if cl_p is not None else None),
               converged=bool(converged), clamped=clamped,
               n_limited=int(n_limited), n_floored=int(n_floored),
               clamped_ever=clamped_ever(clamp_history),
               anchor_ratio=None, usable=False, reason="")
    why = []
    if not converged:
        why.append("not converged")
    if clamped:
        why.append(f"clamped {n_limited}/{n_floored} (GS1.4)")
    if consensus is None or cl_p is None:
        why.append("NO ANCHOR: consensus not supplied -- the anchor check did NOT run")
    else:
        a, b = abs(float(cl_p)), abs(float(consensus))
        rec["anchor_ratio"] = (max(a, b) / min(a, b)) if min(a, b) > 0 else float("inf")
        if rec["anchor_ratio"] > ratio_max:
            why.append(f"anchor_ratio {rec['anchor_ratio']:.2f} > {ratio_max} "
                       f"(cl_p {a:.6g} vs consensus {b:.6g}) -- OUTLIER, not 'wrong'")
    rec["usable"] = not why
    rec["reason"] = "ok" if not why else "; ".join(why)
    return rec


def assess_set(legs, *, axes, ratio_max=ANCHOR_RATIO_MAX):
    """A whole candidate set.

    `legs`: iterable of dicts with cl_p / converged / n_limited / n_floored and whatever
    axis keys you name. `axes`: the axis names the GATE ITSELF samples -- required, not
    optional. ★★ R22's erratum: I asserted "only one usable point at medium" from the
    seed-0 column while the gate samples two seeds, with the seed sitting in the same CSV
    I had read twice. Naming the axes makes that omission impossible to make silently.

    The consensus is the MEDIAN cl_p over converged, unclamped legs -- median so one
    outlier cannot drag it.

    Per-axis spread is UNDEFINED with fewer than two usable legs (R13: never report a
    spread over fewer than two converged legs as "small"; call it UNDEFINED).
    """
    if not axes:
        raise ValueError("assess_set: `axes` must name the axes the gate samples "
                         "(R22 erratum) -- pass them explicitly, e.g. ('seed',)")
    legs = list(legs)
    base = [float(l["cl_p"]) for l in legs
            if l.get("converged") and not l.get("n_limited") and not l.get("n_floored")
            and l.get("cl_p") is not None]
    consensus = float(np.median(base)) if base else None
    out = []
    for l in legs:
        r = assess(l.get("cl_p"), converged=l.get("converged", False),
                   n_limited=l.get("n_limited", 0), n_floored=l.get("n_floored", 0),
                   clamp_history=l.get("clamp_history"), consensus=consensus,
                   ratio_max=ratio_max)
        r.update({a: l.get(a) for a in axes})
        out.append(r)
    spreads = {}
    for a in axes:
        for v in sorted({r[a] for r in out}, key=lambda x: (x is None, x)):
            u = [r["cl_p"] for r in out if r[a] == v and r["usable"]]
            spreads[(a, v)] = (
                dict(n=len(u), spread="UNDEFINED (<2 usable legs)") if len(u) < 2 else
                dict(n=len(u), rel_min=(max(u) - min(u)) / abs(min(u)),
                     rel_max=(max(u) - min(u)) / abs(max(u)),
                     rel_mean=(max(u) - min(u)) / abs(np.mean(u))))
    return dict(consensus=consensus, legs=out, spreads=spreads,
                n_usable=sum(r["usable"] for r in out))
