"""Q1 — does a structured body-fitted grid make the entropy correction's UPSTREAM DONOR
easier to determine? (phase 3 task 2, route A)

The user's FIRST concern, pre-registered in
docs/dev_phase_three/20260811-0300-hex-mesh-prereg.md §4 Q1.

★★ Measured on NACA0012 at the M1 condition (M0.80 / alpha 1.25), NOT on the cylinder: the
entropy correction only acts where the flow is supersonic, so at M = 0 it is inert and a
cylinder reading would say nothing about it. (The cylinder was Q2's case, and only because it
has an exact solution.) The user made this point explicitly, and it also drives the second
design decision below.

★ Two facts verified by READING the code rather than recalling it, because getting either wrong
would mean measuring the wrong object:
  1. `UpwindOperator.__init__` defaults to `weighted=False` -- the integer WALK -- and
     `solve/newton.py` constructs it that way explicitly. So `upstream_elements` IS the donor
     map the entropy correction consumes. (The nearby comment reads as if the weighted P6 flux
     were the default; the signature and the call site say otherwise.)
  2. `upstream_map` RAISES for `weighted=True`: "the kernel-mode flux has a dense neighbourhood
     dependence and builds no single donor map". So the question "is the donor easier to
     determine" is only meaningful on the walk path -- which is the production path.

★★ And the diagnostics are reported BOTH globally and restricted to the SUPERSONIC ZONE,
because the entropy correction only reads donors where it acts; a global average would dilute
exactly the answer being asked for.

Bands (fixed before running, pre-registration §4 Q1):
    D1  median ambiguity margin: hex >= 2x unstructured        -> PASS
    D2  self-donor rate u(e) = e: hex <= 1/2 unstructured      -> PASS
    D3  donor cycles: hex = 0                                  -> PASS
    D4  hop-saturation share                                   RECORDED, no criterion

Premise (measured, never assumed): both meshes must actually carry a supersonic zone at the
condition solved. The project has a documented case of assuming subcriticality wrongly --
NEWTON_M6_RECIPE claimed M0.70 was subcritical and it measures M_max 1.5358 with 214 shock
cells -- so `n_supersonic > 0` is asserted, not hoped for.

Outputs (TRACKED): bench/gate_results/hex_q1_donor.csv
"""

import csv
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)

from pyfp3d.kernels.upwind import UpwindOperator                 # noqa: E402
from pyfp3d.mesh.reader import read_mesh                         # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                        # noqa: E402
from pyfp3d.meshgen.extrude import extrude_single_layer          # noqa: E402
from pyfp3d.meshgen.planar import naca0012_coordinates          # noqa: E402
from pyfp3d.meshgen.structured import airfoil_o_grid_2d         # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared       # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting            # noqa: E402

CSV = os.path.join(HERE, "gate_results", "hex_q1_donor.csv")
#: ★★ M0.72, NOT M0.80. The first run used M0.80 (the M1 condition) and 3 of 4 legs came back
#: conv=False, one with M_max 344.77 -- a diverged field, not a solution -- so its donor
#: diagnostics were not readings at all (GS1.4: a non-converged/clamped state is not a
#: solution). The root cause was choosing the worst available condition: M1's own measurement
#: records NACA0012 M0.80 medium as NOT converging at seed 0 (7265 limited / 758 floored),
#: rescued only by the thread-dependent seed fallback.
#: M0.72 is M1a's IN-ENVELOPE condition: it converges cleanly at every level with 0 clamps, and
#: it is still transonic (M_max ~ 1.17 > 1), so the entropy correction is LIVE and the donor
#: question is meaningful. Donor determinacy does not need to double as a transonic-robustness
#: test -- mixing the two is what made the first run unreadable.
M_INF = 0.72
#: ★★ alpha = 0 DROPPED, and the reason is a process error worth recording. The registration
#: justified it as "still transonic on NACA0012 at M0.80" -- true there (the first run measured
#: M_max 1.2465 at alpha = 0) -- but when the condition moved to M0.72 I carried that claim
#: across WITHOUT re-measuring it. The premise assert fired: at M0.72 alpha = 0 has NO
#: supersonic cell, so the entropy correction is inert and a donor reading there would not be
#: about it. Exactly the documented trap of taking a previous round's attribution into a new
#: measurement as if it were a result.
#: Cost of dropping it, stated rather than glossed: the alignment control is gone. It is
#: limited -- alpha = 1.25 is already within 1.25 deg of the chord, so the flow is nearly as
#: grid-aligned as it can get -- but "nearly" is an argument, not a measurement, and the clean
#: alignment control would need its own supercritical-and-converging condition.
ALPHAS = (1.25,)
SOLVE_KW = dict(upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6, freeze_refresh_max=8,
                precond="direct", direct_refactor_every=4, n_newton_max=80)
D1_FACTOR, D2_FACTOR = 2.0, 0.5


def donor_diagnostics(mesh, grad, supersonic):
    """D1-D4 on the donor map the entropy correction consumes.

    D1 recomputes the FIRST HOP of `upstream_elements` -- among the 4 face neighbours, the one
    with the most negative centroid displacement projected on grad[e] -- and reports how
    clearly that winner is separated from the runner-up. A structured grid whose lines follow
    the flow should separate them cleanly; a tet split introduces diagonal faces whose offsets
    can be nearly perpendicular, which is registered risk 1.
    """
    el = np.ascontiguousarray(mesh.elements)
    upw = UpwindOperator(mesh.nodes, el, weighted=False)
    u = upw.upstream_map(grad).copy()
    nb, cen = upw.face_neighbors, upw.centroids

    #: projected displacement to each of the 4 face neighbours
    disp = np.full((len(el), 4), np.nan, dtype=np.float64)
    for f in range(4):
        j = nb[:, f]
        ok = j >= 0
        d = np.einsum("ij,ij->i", cen[j[ok]] - cen[ok], grad[ok])
        disp[ok, f] = d
    neg = np.where(disp < 0.0, disp, np.nan)          # only upstream neighbours qualify
    n_neg = np.sum(~np.isnan(neg), axis=1)
    srt = np.sort(np.where(np.isnan(neg), 0.0, neg), axis=1)   # most negative first
    best, second = srt[:, 0], srt[:, 1]
    #: margin = how much more upstream the winner is than the runner-up, relative to the winner
    with np.errstate(invalid="ignore", divide="ignore"):
        margin = np.where(best < 0.0, (second - best) / np.abs(best), np.nan)
    #: elements with no upstream neighbour at all have no margin to speak of
    margin[n_neg == 0] = np.nan

    def cycles_in(donor):
        """Count cycles in the donor map (self-donors excluded -- those are D2)."""
        n = len(donor); state = np.zeros(n, dtype=np.int8); n_cyc = 0
        for s in range(n):
            if state[s]:
                continue
            path, e = [], s
            while state[e] == 0:
                state[e] = 1; path.append(e)
                nxt = donor[e]
                if nxt == e:
                    break
                e = nxt
            if state[e] == 1 and donor[e] != e and e in path:
                n_cyc += 1
            for q in path:
                state[q] = 2
        return n_cyc

    out = {}
    for tag, mask in (("all", np.ones(len(el), dtype=bool)), ("supersonic", supersonic)):
        m = margin[mask]
        m = m[~np.isnan(m)]
        out[f"d1_margin_median_{tag}"] = float(np.median(m)) if len(m) else float("nan")
        out[f"d2_self_donor_frac_{tag}"] = float(np.mean(u[mask] == np.flatnonzero(mask))) \
            if mask.any() else float("nan")
        out[f"d_zone_cells_{tag}"] = int(mask.sum())
    out["d3_cycles"] = cycles_in(u)
    out["d4_no_upstream_frac"] = float(np.mean(n_neg == 0))
    return out


def solve_and_diagnose(mesh_obj, label, alpha):
    mc, wc = cut_wake(mesh_obj)
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=alpha, **SOLVE_KW)
    grad = np.ascontiguousarray(r["grad"]) if "grad" in r else None
    if grad is None:                                    # recompute if the driver did not keep it
        from pyfp3d.mesh.metrics import precompute_element_geometry
        B, _ = precompute_element_geometry(mc.nodes, mc.elements)
        grad = np.ascontiguousarray(np.einsum("eaj,ea->ej", B, r["phi"][mc.elements]))
    q2 = np.einsum("ej,ej->e", grad, grad)
    m2 = mach_number_squared(q2, M_INF)
    supersonic = m2 > 1.0
    #: ★★ premise, measured not assumed, and it took TWO clauses. The first version only
    #: checked "is there a supersonic zone" and let three non-converged legs through -- a guard
    #: that did not cover what the conclusion claims, the same shape of error as phase two's
    #: surface-only G1. A donor map read off a diverged field describes garbage gradients.
    assert supersonic.any(), (f"{label} alpha={alpha}: NO supersonic cell at M{M_INF} -- the "
                             f"entropy correction is inert, so a donor reading here would not "
                             f"be about it")
    assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0, (
        f"{label} alpha={alpha}: not a solution (converged={r['converged']}, "
        f"{r['n_limited']} limited / {r['n_floored']} floored, "
        f"|R|={r['residual_history'][-1]:.2e}, M_max={float(np.sqrt(m2.max())):.4f}) -- "
        f"GS1.4: donor diagnostics on a non-converged state are not readings")
    d = donor_diagnostics(mc, grad, supersonic)
    #: ★ built by update(), not dict(**a, **b): this is the THIRD "dict() got multiple values"
    #: of the session (the zone counts collide with n_supersonic), and the incremental form with
    #: prefixed diagnostic keys is the fix that prevents the class, not just this instance.
    row = dict(case=label, alpha=alpha)
    row.update(converged=bool(r["converged"]),
               res_final=float(r["residual_history"][-1]),
               n_limited=int(r["n_limited"]), n_floored=int(r["n_floored"]),
               m_max=float(np.sqrt(m2.max())), n_supersonic=int(supersonic.sum()),
               n_tets=len(mc.elements))
    row.update(d)
    return row


def main():
    rows = []
    unst_path = os.path.join(REPO, "cases", "meshes", "naca0012_2.5d", "medium.msh")
    surf = naca0012_coordinates(n_half=81)[:-1]
    print(f"Q1: donor determinacy on NACA0012 at M{M_INF} (the M1 condition), walk flux\n")
    for alpha in ALPHAS:
        mu = read_mesh(unst_path)
        rows.append(solve_and_diagnose(mu, "unstructured", alpha))
        #: ★ read the COMMITTED family rather than building inline: the user pointed out the
        #: first version's hex mesh existed only inside this script, so the reading rested on a
        #: mesh nobody could inspect or regenerate. Now it is
        #: cases/meshes/naca0012_hex_2.5d/, with stats CSV and an inspection PNG like every
        #: other family.
        mh = read_mesh(os.path.join(REPO, "cases", "meshes", "naca0012_hex_2.5d",
                                    "medium.msh"))
        rows.append(solve_and_diagnose(mh, "structured_hex", alpha))
        for r in rows[-2:]:
            print(f"  a={alpha:<5} {r['case']:15} conv={str(r['converged']):5} "
                  f"M_max {r['m_max']:.4f} sup {r['n_supersonic']:>6}/{r['n_tets']:<7} "
                  f"D1 {r['d1_margin_median_supersonic']:.4f} "
                  f"D2 {r['d2_self_donor_frac_supersonic']:.4f} D3 {r['d3_cycles']}",
                  flush=True)

    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")

    print("\n=== reading (bands fixed in the pre-registration §4 Q1) ===")
    ok_any = False
    for alpha in ALPHAS:
        ru = next(r for r in rows if r["case"] == "unstructured" and r["alpha"] == alpha)
        rh = next(r for r in rows if r["case"] == "structured_hex" and r["alpha"] == alpha)
        for zone in ("supersonic", "all"):
            d1u, d1h = ru[f"d1_margin_median_{zone}"], rh[f"d1_margin_median_{zone}"]
            d2u, d2h = ru[f"d2_self_donor_frac_{zone}"], rh[f"d2_self_donor_frac_{zone}"]
            p1 = d1h >= D1_FACTOR * d1u
            #: ★★ VACUITY GUARD. The first version read `d2h == 0.0` as a PASS when the
            #: unstructured side was ALREADY 0 -- a criterion satisfied by "there was no problem
            #: to begin with", which is evidence of nothing. Same family as the project's
            #: recorded lesson that a non-strict xfail cannot detect a regression. A band only
            #: counts when the baseline had something to improve.
            vacuous = (d2u == 0.0)
            p2 = (not vacuous) and d2h <= D2_FACTOR * d2u
            print(f"  a={alpha:<5} [{zone:10}] D1 {d1u:.4f} -> {d1h:.4f} "
                  f"({d1h / d1u if d1u else float('nan'):.2f}x, need >= {D1_FACTOR}) "
                  f"{'PASS' if p1 else 'FAIL'}   "
                  f"D2 {d2u:.4f} -> {d2h:.4f} (need <= {D2_FACTOR}x) "
                  f"{'PASS' if p2 else ('VACUOUS (baseline already 0)' if vacuous else 'FAIL')}")
            if zone == "supersonic":
                ok_any |= (p1 or p2)
        d3_vac = ru["d3_cycles"] == 0
        print(f"  a={alpha:<5} D3 cycles {ru['d3_cycles']} -> {rh['d3_cycles']} (need 0) "
              f"{'VACUOUS (baseline already 0)' if d3_vac else ('PASS' if rh['d3_cycles'] == 0 else 'FAIL')}")
        print(f"  a={alpha:<5} D4 RECORDED no-upstream frac "
              f"{ru['d4_no_upstream_frac']:.4f} -> {rh['d4_no_upstream_frac']:.4f}")
    print("\n  ⇒ " + ("Q1 has a band PASS in the supersonic zone: route (A) buys donor "
                      "determinacy where the correction acts."
                      if ok_any else
                      "Q1 FAILS in the supersonic zone. With Q2 already failing both bands, "
                      "kill criterion 1 fires and (B) becomes the candidate -- but per the Q2 "
                      "verdict (B) may NOT be justified by wall accuracy."))
    return 0


if __name__ == "__main__":
    sys.exit(main())
