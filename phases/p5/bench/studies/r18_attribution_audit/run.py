"""R18 -- audit my own attribution across every cached leg. Zero solves.

Binding text: phases/p5/docs/dev_phase_five/20260823-1500-r18-prereg.md (committed first).

R14 concluded "the TE blob is created by C=1.0" from a control that changed upwind_c AND
convergence status together -- the defect I used to void R16. This pools every cached leg
and asks which candidate actually separates TE-clean from TE-supersonic.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r18_attribution_audit/run.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "bench"))
RES = os.path.join(HERE, "results")

from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared, q2_at_mach   # noqa: E402

TE_LO, TE_HI = 0.80, 1.20
CANON = ("xcoarse", "coarse", "medium")
S = "bench/studies"
#: (leg, npz, level, m_inf, alpha, upwind_c)  -- G-PROV prints all of it
LEGS = [
    ("R12 xcoarse C1.5 a1.25", f"{S}/r12_h_pricing/results/xcoarse.npz", "xcoarse", 0.80, 1.25, 1.5),
    ("R12 coarse  C1.5 a1.25", f"{S}/r12_h_pricing/results/coarse.npz", "coarse", 0.80, 1.25, 1.5),
    ("R12 medium  C1.5 a1.25", f"{S}/r12_h_pricing/results/medium.npz", "medium", 0.80, 1.25, 1.5),
    ("R12 fine    C1.5 a1.25", f"{S}/r12_h_pricing/results/fine.npz", "fine", 0.80, 1.25, 1.5),
    ("R14 medium  C1.0 a1.25 s0", f"{S}/r14_medium_coverage/results/medium_c10_s0.npz", "medium", 0.80, 1.25, 1.0),
    ("R15 medium  C3.0 a1.25 s5", f"{S}/r15_modes/results/C3.0_s5_control.npz", "medium", 0.80, 1.25, 3.0),
    ("R15 medium  C3.0 a1.25 s0", f"{S}/r15_modes/results/C3.0_s0_target.npz", "medium", 0.80, 1.25, 3.0),
    ("R15 medium  C1.0 a1.25 s5", f"{S}/r15_modes/results/C1.0_s5_target.npz", "medium", 0.80, 1.25, 1.0),
    ("R16 medium  C1.0 a0.0", f"{S}/r16_alpha_dose/results/alpha0.0.npz", "medium", 0.80, 0.0, 1.0),
    ("R16 medium  C1.0 a0.5", f"{S}/r16_alpha_dose/results/alpha0.5.npz", "medium", 0.80, 0.5, 1.0),
    ("R17 medium  C1.0 a0 M.82", f"{S}/r17_mach_alpha0/results/M0.82.npz", "medium", 0.82, 0.0, 1.0),
    ("R17 medium  C1.0 a0 M.84", f"{S}/r17_mach_alpha0/results/M0.84.npz", "medium", 0.84, 0.0, 1.0),
    ("R17 medium  C1.0 a0 M.86", f"{S}/r17_mach_alpha0/results/M0.86.npz", "medium", 0.86, 0.0, 1.0),
]
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def separation(rows, key, higher_is_dirty=True):
    """a = extreme of the CLEAN group, b = extreme of the DIRTY group; separated iff
    they do not overlap. No threshold is chosen -- the two numbers are the result."""
    clean = [r[key] for r in rows if r["te_sup"] == 0]
    dirty = [r[key] for r in rows if r["te_sup"] > 0]
    if not clean or not dirty:
        return None, None, None
    if higher_is_dirty:
        a, b = max(clean), min(dirty)
        return a, b, a < b
    a, b = min(clean), max(dirty)
    return a, b, a > b


def main():
    os.makedirs(RES, exist_ok=True)
    assert "pyfp3d.solve.newton" not in sys.modules, "G-NOSOLVE"
    print("  G-NOSOLVE  no solver entry point imported; cached npz only")
    meshes, rows, dropped = {}, [], []
    for name, rel, lv, m_inf, alpha, C in LEGS:
        if lv not in CANON:
            dropped.append(name); continue          # G-SCOPE
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print(f"  [{name}] no cache"); continue
        if lv not in meshes:                        # G-MESH: its own mesh
            mc, _ = cut_wake(read_mesh(f"{ROOT}/cases/meshes/naca0012_2.5d/{lv}.msh"))
            B, _ = precompute_element_geometry(mc.nodes, mc.elements)
            meshes[lv] = (mc, B, mc.nodes[mc.elements].mean(axis=1))
        mc, B, cent = meshes[lv]
        d = np.load(p)
        phi = np.asarray(d["phi"], float)
        g = np.empty((len(mc.elements), 3)); q2 = np.empty(len(mc.elements))
        element_velocity_q2(mc.elements, B, phi, g, q2)
        #: G-MINF: the leg's OWN m_inf
        M = np.sqrt(np.maximum(mach_number_squared(q2, m_inf, 1.4), 0.0))
        up = cent[:, 0] < TE_LO
        te = (cent[:, 0] >= TE_LO) & (cent[:, 0] <= TE_HI)
        gam = np.asarray(d["gamma"], float)
        rows.append(dict(leg=name, level=lv, m_inf=m_inf, alpha=alpha, upwind_c=C,
                         conv=int(bool(d["conv"])),
                         n_limited=int(d["nlim"]), n_floored=int(d["nflr"]),
                         upM=round(float(M[up].max()), 4),
                         te_maxM=round(float(M[te].max()), 4),
                         te_sup=int((M[te] > 1).sum()),
                         gamma=round(float(gam[0]), 6) if gam.size else None,
                         cap_ratio=round(float(q2.max() / q2_at_mach(3.0, m_inf, 1.4)), 4)))
    print(f"  G-SCOPE    excluded {len(dropped)} non-canonical leg(s): {dropped}")
    print(f"  G-PROV     {len(rows)} legs pooled\n")
    rows.sort(key=lambda r: r["upM"])
    print(f"  {'upM':>9}{'TEsup':>7}{'TEmaxM':>9}{'conv':>6}{'C':>6}{'alpha':>7}"
          f"{'M_inf':>7}{'Gamma':>11}  leg")
    for r in rows:
        print(f"  {r['upM']:9.4f}{r['te_sup']:7d}{r['te_maxM']:9.4f}{r['conv']:6d}"
              f"{r['upwind_c']:6.1f}{r['alpha']:7.2f}{r['m_inf']:7.2f}"
              f"{(r['gamma'] if r['gamma'] is not None else float('nan')):11.5f}  {r['leg']}")
    with open(os.path.join(RES, "pooled.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print()
    a, b, sep = separation(rows, "upM")
    _record("A-SEP", "does upM separate TE-clean from TE-supersonic legs",
            "separated iff max(upM | clean) < min(upM | dirty); the two numbers ARE the "
            "result, no threshold is chosen",
            f"max(clean) {a:.4f}  <  min(dirty) {b:.4f}  ->  gap band ({a:.4f}, {b:.4f})"
            if sep else f"max(clean) {a}  vs  min(dirty) {b}  -> OVERLAP",
            "★★★ A-SEP: upM SEPARATES" if sep else "A-SEP: upM does NOT separate")

    riv = {}
    for key, hid in (("upwind_c", False), ("alpha", True), ("m_inf", True), ("conv", False)):
        aa, bb, ss = separation(rows, key, higher_is_dirty=hid)
        riv[key] = (aa, bb, ss)
    _record("A-RIVAL", "the same separation test on the four rival candidates",
            "any candidate that also separates must be reported ALONGSIDE upM -- with "
            "two separating candidates this round cannot say which is the cause",
            "; ".join(f"{k}: clean-extreme {v[0]} vs dirty-extreme {v[1]} -> "
                      f"{'SEPARATES' if v[2] else 'overlaps'}" for k, v in riv.items()),
            "A-RIVAL: " + ", ".join(k for k, v in riv.items() if v[2]) + " also separate"
            if any(v[2] for v in riv.values()) else "A-RIVAL: none of the rivals separate")

    c_sep = riv["upwind_c"][2]
    _record("A-R14", "does R14's 'C=1.0 creates the TE blob' still stand",
            "upM separates and upwind_c does not => R14's wording must become "
            "'C=1.0 is a ROUTE to the condition, not the cause'",
            f"upM separates: {sep};  upwind_c separates: {c_sep}",
            "★★ A-R14: R14's wording must be corrected" if sep and not c_sep else
            "A-R14: R14's wording survives" if c_sep else "A-R14: UNDEFINED")

    dirty = [r for r in rows if r["te_sup"] > 0]
    clean = [r for r in rows if r["te_sup"] == 0]
    _record("A-GAMMA", "Gamma against TE involvement, pooled",
            "RECORDED -- how R17's premise failure looks across the whole sample",
            f"dirty legs |Gamma| "
            f"{[abs(r['gamma']) for r in dirty if r['gamma'] is not None]};  "
            f"clean legs |Gamma| "
            f"{[abs(r['gamma']) for r in clean if r['gamma'] is not None]}", "RECORDED")

    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
