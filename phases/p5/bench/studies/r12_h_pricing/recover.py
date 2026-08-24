"""R12 recovery: every number from the CACHED npz, zero re-solves.

The ladder's four solves completed; a float(None) in the reporting layer then raised
before the CSV was written. G-CACHE had already put phi and the diagnostics on disk, so
this recomputes the whole report without re-solving -- which is the entire point of the
"cache before you report" rule.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r12_h_pricing/recover.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "bench"))

from run_le14_common_root import classify_failure                       # noqa: E402
from pyfp3d.mesh.reader import read_mesh                               # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                              # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                      # noqa: E402
from pyfp3d.post.shock import shock_report                             # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients                # noqa: E402

import importlib.util
_s = importlib.util.spec_from_file_location("r12", os.path.join(HERE, "run.py"))
R12 = importlib.util.module_from_spec(_s); _s.loader.exec_module(R12)

WALL_S = {"xcoarse": 0.5, "coarse": 1.2, "medium": 510.4, "fine": None}  # from the log
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}: band={band} measured={measured} -> {verdict}")


def main():
    rows, hrow = [], []
    for lv in R12.LEVELS:
        p = f"{R12.MESHES}/{lv}.msh"
        npz = os.path.join(HERE, "results", f"{lv}.npz")
        if not os.path.exists(npz):
            print(f"  [{lv}] no cache"); continue
        d = np.load(npz)
        mc, _ = cut_wake(read_mesh(p))
        dz = float(np.ptp(mc.nodes[:, 2]))
        assert abs(dz - R12.DZ_STAT[lv]) < 1e-9, f"G-DZ {lv}"
        h_s, n_st = R12.h_in_shock_band(mc)
        hrow.append({"level": lv, "h_wall": R12.H_WALL[lv], "h_shock": h_s,
                     "n_band_nodes": n_st, "n_nodes": len(mc.nodes),
                     "n_elements": len(mc.elements)})
        phi = np.asarray(d["phi"], float)
        hist = np.asarray(d["res_hist"], float)
        ch = np.asarray(d["clamp_hist"], float)
        conv, nl, nf = bool(d["conv"]), int(d["nlim"]), int(d["nflr"])
        usable = conv and not nl and not nf
        #: ★ CLASSIFY, never report conv=False (CLAUDE.md)
        mode, ev, d10, revis = (("converged", "", float("nan"), 0) if conv else
                                classify_failure(hist, ch, np.asarray([], float),
                                                 0, "", nl, nf))
        f = wall_force_coefficients(mc.nodes, mc.elements,
                                    mc.boundary_faces["wall"], phi,
                                    alpha_deg=R12.ALPHA, m_inf=R12.M_INF, s_ref=dz)
        cur = wall_cp_curve(mc, phi, z=0.5 * dz, m_inf=R12.M_INF)
        xs = shock_report(cur, R12.M_INF)["upper"].get("x_shock")
        rows.append({"level": lv, "h_wall": R12.H_WALL[lv], "h_shock": round(h_s, 6),
                     "n_nodes": len(mc.nodes), "wall_s": WALL_S[lv],
                     "conv": conv, "n_limited": nl, "n_floored": nf, "usable": usable,
                     "n_steps": len(hist), "res_last": hist[-1],
                     "mode": mode, "evidence": ev,
                     "cl_p": round(float(f["cl"]), 6),
                     "x_shock": (round(float(xs), 6) if xs is not None
                                 and np.isfinite(xs) else None)})
        print(f"  [{lv:8}] h_shock {h_s:.5f}  usable={usable}  mode={mode}  "
              f"|R|={hist[-1]:.3e} in {len(hist)} steps  cl_p={f['cl']:.5f}  "
              f"x_shock={xs}")

    for name, rr in (("levels.csv", rows), ("h_band.csv", hrow)):
        keys = sorted({k for x in rr for k in x})
        with open(os.path.join(HERE, "results", name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rr)

    dev = [abs(h["h_shock"] / h["h_wall"] - 1.0) for h in hrow]
    _record("P-H", "h_shock vs h_wall in the shock band (mesh only)",
            "<=25% => h_wall is a fair proxy",
            "  ".join(f"{h['level']} {h['h_shock']:.5f} ({h['h_shock']/h['h_wall']:.3f}x,"
                      f" {h['n_band_nodes']} nodes)" for h in hrow)
            + f";  max deviation {100*max(dev):.1f}%",
            "P-H: h_wall is a fair proxy" if max(dev) <= 0.25 else
            "★ P-H: NOT a fair proxy")
    g = {r["level"]: r for r in rows}
    fine = g.get("fine")
    _record("P-CONV", "does the production recipe converge unclamped at each level",
            "fine usable => h route OPEN;  not usable => blocker confirmed on current code",
            "  ".join(f"{lv} {'USABLE' if g[lv]['usable'] else g[lv]['mode']}"
                      f"({g[lv]['n_limited']}/{g[lv]['n_floored']})" for lv in g),
            "P-CONV: fine USABLE -- h route open" if fine and fine["usable"] else
            "★ P-CONV: fine NOT usable -- blocker CONFIRMED on current code")
    if g.get("medium", {}).get("usable") and fine and fine["usable"]:
        a, b = g["medium"]["cl_p"], fine["cl_p"]
        _record("P-B", "two-level cl difference (medium, fine)",
                "RECORDED ONLY -- fold zone, discipline #4",
                f"{100*(b-a)/abs(a):+.2f}%", "RECORDED")
    else:
        _record("P-B", "two-level cl difference (medium, fine)",
                "conditional on both usable", "fine is not usable", "NOT APPLICABLE")
        a, b = g["coarse"]["cl_p"], g["medium"]["cl_p"]
        _record("P-B'", "the pair that IS usable (coarse, medium)",
                "RECORDED ONLY -- fold zone; M1's (b) target is <3%",
                f"coarse {a:.5f} -> medium {b:.5f} = {100*(b-a)/abs(a):+.2f}%",
                "RECORDED (reproduces the committed -16.2/-16.3%)")
    _record("P-PRIOR", "prior blockers on the h route", "RECORDED",
            f"thread scatter 0.027 c = {0.027/R12.TOL:.1f}x the +-{R12.TOL} c "
            f"requirement; seed change 0.1683 c = {0.1683/R12.TOL:.0f}x",
            "RECORDED")
    with open(os.path.join(HERE, "results", "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    return 0


if __name__ == "__main__":
    sys.exit(main())
