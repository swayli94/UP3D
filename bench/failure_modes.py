"""LE-14: one quantity for three failures -- what actually separates converging from failing?

Three independent things measured today all trade tip-loading freedom against high-Mach
convergence:

  removing the taper        M0.88 wing-body medium fails (6/6 clamps)
  the r_c dead band         0.025 OK, 0.030 FAIL, 0.0375 FAIL, 0.045 OK, 0.05 OK
  the round tip             LS M6 medium envelope M0.84 -> M0.70

★ NO NEW HYPOTHESIS. Five mechanism guesses of mine were refuted by my own measurements today
(dissipation strength, the entropy correction weakening the shock, supersonic connectivity, the
root/symmetry termination, taper-as-the-recipe-gap). The two rounds that produced anything --
LE-11 and LE-12 -- localised and measured instead of testing a guess. So this round asks which
MEASURABLE separates the known success/failure pattern, and lets that pick the mechanism.

★★ FAILURE-MODE CLASSIFICATION FIRST -- a method gap the user identified, and it may invalidate
the separation question itself. Every failure today has been recorded as "conv=False, |R| = ...,
n/m clamps" and then correlated against a knob, without ever characterising WHAT fails. Five
distinguishable modes have distinct signatures and completely different fixes:

  local Mach too high / clamping   n_limited (m_cap) or n_floored (rho_floor) > 0, AND where
                                   those cells sit -- tip, LE or TE
  limit cycle                      residual oscillates with no trend (descent10 ~ 1, values
                                   revisited)
  ill-conditioning                 GMRES iteration counts climb, n_gmres_stalled > 0
  line-search collapse             lambda driven to ~0, steps rejected
  sigma-transport failure          accept_reason = sigma_transport_not_converged

Lumping all five under one "not converged" label means any correlation analysis may be trying to
separate a MIXTURE OF DIFFERENT DISEASES with one number -- which would explain why five
mechanism hypotheses failed today. So each leg is classified before the separation question is
asked, and if the three failures turn out to be different modes, "one common root cause" is the
wrong question and that is the round's result.

★ THE DEAD BAND IS THE HYPOTHESIS-KILLER, and that makes the criterion sharp: any quantity
MONOTONE in r_c cannot explain 0.025 OK / 0.030 FAIL / 0.0375 FAIL / 0.045 OK. A quantity only
"explains" the pattern if it separates all three successes from all three failures with NO
overlap -- across a non-monotone knob. That is a hard test to pass by accident.

The physically motivated candidate comes from tip_taper_factors' own derivation rather than from
me: the tip-edge singularity is driven by the TRAILING vorticity gamma = -dGamma/dz, not by bound
Gamma. All three knobs change it -- the taper directly, r_c through its profile, the round tip by
freeing the tip loading. So the primary measurable is |dGamma/dz| at the outermost stations, with
three others recorded alongside so the round is not a single-candidate bet:

  A  |dGamma/dz| over the outermost 1 / 2 / 3 station intervals   (the derivation's own quantity)
  B  Gamma at the outermost station                                (bound, not trailing -- the
     derivation says this is NOT the driver, so it is the control: if B separates and A does not,
     the derivation is wrong)
  C  the tip residual-density peak                                 (LE-11's quantity, which
     separated 5-15x on the entropy/taper legs)
  D  M^2 max within 10 % of the tip

One geometry, one knob, six configurations whose outcomes are already known -- so the pattern is
not being discovered here, only explained. gamma IS cached this time; LE-10 stored only phi and
LE-11 could not reconstruct the full residual because of it.

Outputs (TRACKED):    bench/gate_results/le14_common_root.csv
Outputs (gitignored): bench/gate_results/le14_cache/<r_c>.npz
"""

import csv
import math
import os
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import capability_matrix as cap                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry          # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared            # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

OUT = os.path.join(HERE, "gate_results")
CACHE = os.path.join(OUT, "le14_cache")
os.makedirs(CACHE, exist_ok=True)
CSV = os.path.join(OUT, "le14_common_root.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_TARGET, ALPHA = 0.88, 3.06
#: r_c with its KNOWN outcome from LE-5/LE-7/LE-8, so the pattern is given, not discovered
R_CS = [(0.0, False), (0.025, True), (0.030, False), (0.0375, False),
        (0.045, True), (0.05, True)]
KEYS = ["r_c", "known_outcome", "converged", "failure_mode", "mode_evidence",
        "n_limited", "n_floored", "clamp_site", "res_final", "descent10",
        "res_revisits", "n_gmres_stalled", "accept_reason", "f_final",
        "A_dgdz_1", "A_dgdz_2", "A_dgdz_3", "B_gamma_tip",
        "C_tip_res_peak", "D_m2_tip_max", "wall_s", "note"]


def classify_failure(hist, clamps, fhist, gstall, areason, nlim, nflr):
    """Which of the five failure modes is present? Signatures, not guesses.

    Returns (mode, evidence). Modes are not exclusive in principle, so the first matching
    signature in order of specificity wins and the evidence string carries the rest -- a leg
    that both clamps and limit-cycles must not be silently filed as one of them.
    """
    ev = []
    d10 = float("nan")
    revisits = 0
    if len(hist) >= 11 and hist[-1] > 0:
        d10 = float(hist[-11] / hist[-1])
    monotone_tail = True
    n_up = 0
    if len(hist) >= 6:
        #: ★ REVISIT TOLERANCE WIDENED and MONOTONICITY ADDED, 2026-08-05, because the
        #: original pair of tests misclassified a textbook limit cycle as budget_limited.
        #: MEASURED on the wing-body r_c = 0.0375 leg: the tail is a period-3 cycle,
        #: 9.037e-05 -> 8.059e-05 -> 4.026e-05 repeating to 4-5 significant figures, and
        #: raising n_newton_max 60 -> 200 left |R| and descent10 BIT-IDENTICAL. Two faults:
        #:   (1) the 1e-12 relative revisit tolerance is far tighter than a cycle reproduces
        #:       to in floating point over many steps, so `revisits` counted 0;
        #:   (2) descent10 = 2.0021 is just "the cycle's high point is 2x its low point" --
        #:       an artefact of WHERE the 10-step window lands in the period, with no
        #:       bearing on convergence -- and it fell on the >2.0 side of a threshold I
        #:       picked by hand. A hand-picked threshold decided a published conclusion.
        #: So descent10 alone can no longer assign budget_limited: the tail must ALSO be
        #: monotone. A cycle always fails that, whatever its window ratio.
        tail = hist[-12:]
        pos = tail[tail > 0]
        revisits = int(sum(1 for i in range(len(pos))
                           for j in range(i + 1, len(pos))
                           if abs(pos[i] - pos[j]) <= 1e-4 * pos[i]))
        d = np.diff(tail)
        n_up = int(np.sum(d > 0))
        monotone_tail = bool(np.all(d < 0))
    if areason and areason not in ("None", "?", "nan"):
        ev.append(f"accept_reason={areason}")
    if nlim or nflr:
        ev.append(f"clamps {nlim}/{nflr}")
    if gstall:
        ev.append(f"gmres_stalled={gstall}")
    if d10 == d10:
        ev.append(f"descent10={d10:.3f}")
    if revisits:
        ev.append(f"res_revisits={revisits}")
    if len(hist) >= 6:
        ev.append(f"tail_monotone={monotone_tail} (n_up={n_up})")
    if len(fhist):
        ev.append(f"F_final={float(fhist[-1]):.2e}")
    #: order matters: a named accept_reason is the most specific statement the solver
    #: itself makes, so it wins over inferred signatures
    if areason and "sigma_transport" in areason:
        mode = "sigma_transport"
    elif gstall:
        mode = "ill_conditioning"
    elif nlim or nflr:
        mode = "clamping"
    elif revisits >= 2 or not monotone_tail:
        #: a non-monotone tail IS a cycle, whatever descent10 says about it
        mode = "limit_cycle"
    elif d10 == d10 and d10 > 2.0 and monotone_tail:
        #: budget_limited now requires BOTH: still descending by a margin AND descending
        #: monotonically. Validation case -- LE-1's ls_naca_medium was confirmed genuinely
        #: budget-limited by measurement (raising the cap converged it in 6 more steps) and
        #: its tail is 5.6e-07 -> 1.58e-10 strictly decreasing, so it still classifies here.
        mode = "budget_limited"
    else:
        mode = "unclassified"
    return mode, "; ".join(ev), d10, revisits


def main():
    print(f"wing-body medium, M{M_TARGET} / alpha {ALPHA} -- one knob (r_c), "
          f"six known outcomes")
    print("criterion: a quantity explains the pattern only if it separates the "
          "3 successes from\nthe 3 failures with NO overlap -- across a "
          "NON-MONOTONE knob\n")
    rows = []
    for rc, known in R_CS:
        npz = os.path.join(CACHE, f"rc{rc}.npz")
        mc, wc = cut_wake(read_mesh(MP))
        t = (None if rc == 0.0
             else tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                                    rc * B_SEMI))
        t0 = time.perf_counter()
        if os.path.exists(npz):
            d = np.load(npz)
            phi, gam, conv = d["phi"], d["gamma"], bool(d["conv"])
            nlim, nflr, res = int(d["nlim"]), int(d["nflr"]), float(d["res"])
            hist = d["hist"] if "hist" in d else np.array([])
            clamps = d["clamps"] if "clamps" in d else np.array([])
            fhist = d["fhist"] if "fhist" in d else np.array([])
            gstall = int(d["gstall"]) if "gstall" in d else 0
            areason = str(d["areason"]) if "areason" in d else "?"
            wall = 0.0
        else:
            try:
                seed = solve_newton_lifting(mc, wc, m_inf=cap.WB_MSTART,
                                            alpha_deg=ALPHA, **cap.CONF_SEED_KW)
                nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                          phi_init=seed["phi"], gamma_init=seed["gamma"],
                          n_picard_seed=0)
                if t is not None:
                    nk["tip_taper"] = t
                r = solve_newton_transonic(mc, wc, m_inf=M_TARGET,
                                           alpha_deg=ALPHA,
                                           m_start=cap.WB_MSTART, dm=cap.DM,
                                           dm_min=0.01, freeze_tol=1e-5,
                                           intermediate_tol=1e-4, newton_kw=nk)
            except Exception as exc:                               # noqa: BLE001
                print(f"  r_c {rc:<7} DIED {type(exc).__name__}", flush=True)
                rows.append(dict(r_c=rc, known_outcome=known, converged=False,
                                 note=f"{type(exc).__name__}: {exc}"))
                continue
            phi = np.asarray(r["phi"]); gam = np.asarray(r["gamma"])
            #: ★ the diagnostic history, cached so the failure MODE can be classified without
            #: a third re-solve. LE-10 cached only phi and LE-11 then could not reconstruct the
            #: residual; LE-14's first launch repeated the omission for the history and was
            #: killed 5 minutes in rather than paying for a fourth pass over these states.
            hist = np.asarray(r.get("residual_history", []), dtype=float)
            clamps = np.asarray(r.get("clamp_history", []), dtype=float)
            fhist = np.asarray(r.get("F_history", []), dtype=float)
            gstall = int(r.get("n_gmres_stalled") or 0)
            areason = str(r.get("accept_reason"))
            #: ★ GS4.0 (2026-08-16): `m_final` / `target_reached` existed only on the
            #: DELETED level-set driver, so this `.get` chain had been falling through
            #: to M_TARGET since phase 3 -- making the second conjunct below identically
            #: True and silently degrading `conv` to the bare driver flag. Read the key
            #: with no default now: a missing key is a library regression and must be
            #: loud. (phases/p3/docs/inspection/20260816-2200-independent-audit-zh.md §7.1)
            m_att = r["m_final"]
            conv = bool(r["target_reached"])
            nlim = int(r.get("n_limited") or 0); nflr = int(r.get("n_floored") or 0)
            res = float(r["residual_history"][-1])
            wall = time.perf_counter() - t0
            #: ★ gamma cached too -- LE-10 stored only phi and LE-11 then could not
            #: reconstruct the full residual. Same states, second re-solve avoided.
            np.savez_compressed(npz, phi=phi, gamma=gam, conv=conv,
                                nlim=nlim, nflr=nflr, res=res, hist=hist,
                                clamps=clamps, fhist=fhist, gstall=gstall,
                                areason=areason)

        #: --- A: trailing vorticity |dGamma/dz| at the outermost stations -------
        o = np.argsort(wc.station_z)
        z, g = wc.station_z[o], gam[o]
        dg = np.abs(np.diff(g) / np.diff(z))
        A = [float(np.max(dg[-k:])) if len(dg) >= k else float("nan")
             for k in (1, 2, 3)]
        B = float(g[-1])
        #: --- C/D: tip residual-density peak and tip M^2, from the same field ----
        Bm, V = precompute_element_geometry(mc.nodes, mc.elements)
        grad = np.einsum("eaj,ea->ej", Bm, phi[mc.elements])
        q2 = np.einsum("ej,ej->e", grad, grad)
        m2 = mach_number_squared(q2, M_TARGET)
        cen = mc.nodes[mc.elements].mean(axis=1)
        tipm = np.abs(cen[:, 2] - B_SEMI) < 0.10 * B_SEMI
        D = float(m2[tipm].max()) if np.any(tipm) else float("nan")
        C = float(np.abs(np.einsum("eaj,ej->ea", Bm, grad)).sum(axis=1)[tipm].max()
                  ) if np.any(tipm) else float("nan")
        mode, ev, d10, revisits = classify_failure(hist, clamps, fhist, gstall,
                                                   areason, nlim, nflr)
        #: WHERE the clamped cells sit -- the mode alone does not say tip / LE / TE
        site = "n/a"
        if nlim or nflr:
            hot = tipm & (m2 >= np.quantile(m2[tipm], 0.999)) if np.any(tipm) else tipm
            xcs = cen[:, 0]
            site = ("tip" if np.any(tipm) and m2[tipm].max() >= m2.max() * 0.999
                    else f"elsewhere (argmax x={xcs[int(np.argmax(m2))]:.3f}, "
                         f"z/b={cen[int(np.argmax(m2)), 2] / B_SEMI:.3f})")
        print(f"  r_c {rc:<7} known={known!s:5s} got={conv!s:5s} "
              f"MODE={mode:16s} site={site:12s}", flush=True)
        print(f"           {ev}", flush=True)
        print(f"           A1={A[0]:.5f} A2={A[1]:.5f} A3={A[2]:.5f}  "
              f"B={B:+.6f}  C={C:.4e}  D={D:.3f}"
              f"{'' if wall == 0 else f'  ({wall:.0f}s)'}", flush=True)
        rows.append(dict(r_c=rc, known_outcome=known, converged=conv,
                         failure_mode=mode, mode_evidence=ev, clamp_site=site,
                         descent10=(round(d10, 4) if d10 == d10 else None),
                         res_revisits=revisits, n_gmres_stalled=gstall,
                         accept_reason=areason,
                         f_final=(float(fhist[-1]) if len(fhist) else None),
                         n_limited=nlim, n_floored=nflr, res_final=res,
                         A_dgdz_1=round(A[0], 8), A_dgdz_2=round(A[1], 8),
                         A_dgdz_3=round(A[2], 8), B_gamma_tip=round(B, 8),
                         C_tip_res_peak=C, D_m2_tip_max=round(D, 6),
                         wall_s=round(wall, 1), note=""))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== failure modes first: is 'one common root cause' even the right question? ===")
    modes = {r["r_c"]: r.get("failure_mode") for r in rows if not r["converged"]}
    print(f"  failing legs: {modes}")
    if len(set(modes.values())) > 1:
        print("  ★ the failures are DIFFERENT MODES -- 'one common root cause' is the")
        print("  wrong question, and any single quantity separating them would be")
        print("  separating a mixture of diseases. Read the per-mode evidence instead.")
    else:
        print("  all failing legs share one mode, so a common quantity is worth seeking.")

    print("\n=== which quantity separates the pattern? ===")
    good = [r for r in rows if r["converged"]]
    bad = [r for r in rows if not r["converged"] and r.get("A_dgdz_1") is not None]
    if not good or not bad:
        print("  incomplete -- cannot read"); return 0
    for key, label in (("A_dgdz_1", "A |dG/dz| outermost 1"),
                       ("A_dgdz_2", "A |dG/dz| outermost 2"),
                       ("A_dgdz_3", "A |dG/dz| outermost 3"),
                       ("B_gamma_tip", "B Gamma at tip (control)"),
                       ("C_tip_res_peak", "C tip residual peak"),
                       ("D_m2_tip_max", "D M^2 max near tip")):
        gv = [r[key] for r in good if r.get(key) is not None]
        bv = [r[key] for r in bad if r.get(key) is not None]
        if not gv or not bv:
            continue
        sep = (min(gv) > max(bv)) or (max(gv) < min(bv))
        print(f"  {label:26s} ok {min(gv):.5g}..{max(gv):.5g}   "
              f"fail {min(bv):.5g}..{max(bv):.5g}   "
              f"{'★ SEPARATES' if sep else 'overlaps'}")
    print("\n  A separating and B not  => the derivation is right: TRAILING vorticity")
    print("  drives it, bound Gamma does not, and all three of today's trades are")
    print("  one mechanism.")
    print("  B separating and A not  => the derivation is wrong and bound loading is")
    print("  the variable.")
    print("  none separating         => the three failures are NOT one mechanism, and")
    print("  the round-tip envelope loss has to be accepted as its own cost.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
