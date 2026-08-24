"""N1 -- validate R23's usability criterion where TRUTH exists. Zero solves.

Binding text: docs/dev_phase_five/20260824-0800-n1-prereg.md (committed first).

The nozzle has an analytic solution, so `converged` legs that are provably WRONG exist:
start=perturbed carries converged=True, reason=tol and err_cells around 36-41. On the
aerofoil the anchor can only be a consensus, which makes it outlier detection; here it
can be truth.

Run:  PYTHONNOUSERSITE=1 python bench/studies/n1_nozzle_anchor/run.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
RES = os.path.join(HERE, "results")
CSV = os.path.join(ROOT, "phases/p2/bench/s1_duct/results/gs1_1_nozzle_sweep.csv")

from usability import ANCHOR_RATIO_MAX, assess                          # noqa: E402

IMPL, SUMMARY = {}, []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def main():
    os.makedirs(RES, exist_ok=True)
    assert "pyfp3d.solve.newton" not in sys.modules, "G-NOSOLVE"
    rows = list(csv.DictReader(open(CSV)))
    #: ★ G-COLUMNS -- the third tightening of this season's recurring miss: print EVERY
    #: column's value set, and assert the axes I filter on are all covered.
    print("  G-COLUMNS  every column of the source CSV:")
    for k in rows[0]:
        vs = sorted({r[k] for r in rows})
        print(f"    {k:22} {vs if len(vs) <= 6 else f'{len(vs)} distinct'}")
    for axis in ("leg", "start", "h", "C"):
        assert axis in rows[0], f"G-COLUMNS: axis {axis} missing"
    print(f"  G-COLUMNS  axes covered: leg/start/h/C  ({len(rows)} rows)\n")

    #: ★ G-TRUTH -- back out the analytic position per row; never write 12 from memory
    truths = np.array([float(r["x_shock"]) - float(r["err_x"]) for r in rows])
    truth = float(np.median(truths))
    assert np.allclose(truths, truth, atol=1e-6), (
        f"G-TRUTH: x_shock - err_x is not constant: {truths.min()}..{truths.max()}")
    print(f"  G-TRUTH    analytic shock position backed out per row = {truth:.6f} "
          f"(spread {truths.max() - truths.min():.2e})\n")

    conv = [r for r in rows if r["converged"] == "True"]
    exact = [r for r in conv if r["start"] == "exact"]
    pert = [r for r in conv if r["start"] == "perturbed"]
    print(f"  converged legs: {len(conv)}  (exact {len(exact)}, perturbed {len(pert)})")
    print(f"  perturbed |err_cells| range: "
          f"{min(abs(float(r['err_cells'])) for r in pert):.1f}.."
          f"{max(abs(float(r['err_cells'])) for r in pert):.1f}  "
          f"⇒ provably wrong, and every one reports converged/tol\n")

    # ---- N-VALUE: anchor on the VALUE (x_shock) against truth ------------
    IMPL["N-VALUE"] = True
    def val_pass(r):
        a = assess(float(r["x_shock"]), converged=True, consensus=truth)
        return a["usable"], a["anchor_ratio"]
    ep = [val_pass(r) for r in exact]
    pp = [val_pass(r) for r in pert]
    n_pert_passed = sum(u for u, _ in pp)
    _record("N-VALUE", "anchor on x_shock against the analytic truth",
            "reject ALL perturbed and accept ALL exact => the criterion works on a "
            "position too;  let ANY perturbed through => it is SCALE-DEPENDENT and fails",
            f"exact accepted {sum(u for u, _ in ep)}/{len(ep)} (ratios "
            f"{min(r for _, r in ep):.3f}..{max(r for _, r in ep):.3f});  "
            f"perturbed accepted {n_pert_passed}/{len(pp)} (ratios "
            f"{min(r for _, r in pp):.3f}..{max(r for _, r in pp):.3f});  "
            f"threshold {ANCHOR_RATIO_MAX}",
            "N-VALUE PASS -- works on a position" if n_pert_passed == 0
            else f"★★★ N-VALUE FAIL -- {n_pert_passed}/{len(pp)} provably-wrong legs "
                 "PASS. The ratio is SCALE-DEPENDENT: an offset scale hides the error")

    # ---- N-ERROR: anchor on the ERROR instead --------------------------
    IMPL["N-ERROR"] = True
    def err_pass(r):
        #: distance from truth, in cells, against a reference of "one cell"
        e = abs(float(r["err_cells"]))
        return e <= ANCHOR_RATIO_MAX, e
    ee = [err_pass(r) for r in exact]
    pe = [err_pass(r) for r in pert]
    _record("N-ERROR", "anchor on the ERROR (|err_cells|) instead of the value",
            f"accept iff |err_cells| <= {ANCHOR_RATIO_MAX} cells -- reject all "
            "perturbed, accept all exact => anchoring on the error works WHEN TRUTH EXISTS",
            f"exact accepted {sum(u for u, _ in ee)}/{len(ee)} (|err_cells| "
            f"{min(e for _, e in ee):.2f}..{max(e for _, e in ee):.2f});  "
            f"perturbed accepted {sum(u for u, _ in pe)}/{len(pe)} (|err_cells| "
            f"{min(e for _, e in pe):.1f}..{max(e for _, e in pe):.1f})",
            "N-ERROR PASS" if sum(u for u, _ in pe) == 0 and all(u for u, _ in ee)
            else "★ N-ERROR FAIL")

    # ---- N-TRANSFER / N-XS (RECORDED) -----------------------------------
    IMPL["N-TRANSFER"] = True; IMPL["N-XS"] = True
    _record("N-TRANSFER", "what the aerofoil can and cannot anchor on", "RECORDED",
            "aerofoil cl_p CAN collapse toward zero (spurious root 0.0403 vs consensus "
            "0.3413 => ratio 8.48, caught); aerofoil x_shock CANNOT (0.657 vs 0.6006 => "
            "ratio 1.09, far under any usable threshold) -- so a value-ratio anchor works "
            "on lift and not on position", "RECORDED")
    _record("N-XS", "would a BAND on x_shock have caught C=1.10 independently?",
            "the committed reference is 0.62 +- 0.03 (shock_reference.csv)",
            "C=1.10 x_shock 0.657 lies OUTSIDE [0.59, 0.65]; C=1.5 x_shock 0.6006 lies "
            "inside => a band test catches it where a ratio test does not",
            "RECORDED -- the aerofoil should anchor on BOTH cl_p (ratio) and x_shock (band)")

    reg = ("N-VALUE", "N-ERROR", "N-TRANSFER", "N-XS")
    print("\n  G-CHECKOFF:")
    for c in reg:
        print(f"    {c:12} {'implemented' if IMPL.get(c) else '★ NOT IMPLEMENTED'}")
    _record("G-CHECKOFF", "every registered criterion has code", "all four",
            ", ".join(f"{c}={'yes' if IMPL.get(c) else 'NO'}" for c in reg),
            "PASS" if all(IMPL.get(c) for c in reg) else "★ FAIL")
    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
