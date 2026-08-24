"""R14 -- medium's convergence coverage: WHERE the clamped cells are.

Binding text: phases/p5/docs/dev_phase_five/20260823-0700-r14-prereg.md (committed first).

The modes are already classified and committed (clamping at C=1.0 on both seeds,
limit_cycle at C=3.0/seed 0, clean at C=1.5). CLAUDE.md's clamping signature is the
counts "and where those cells are", and only the counts exist -- so this locates them.

Run:  PYTHONNOUSERSITE=1 python bench/studies/r14_medium_coverage/run.py
"""
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "bench"))
RES = os.path.join(HERE, "results")

from run_gs40h_strength_and_2p5d import NACA_KW                         # noqa: E402
from pyfp3d.kernels.gradient import element_velocity_q2                 # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry             # noqa: E402
from pyfp3d.mesh.reader import read_mesh                                # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                               # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared, q2_at_mach   # noqa: E402
from pyfp3d.post.section_cut import wall_cp_curve                       # noqa: E402
from pyfp3d.post.shock import shock_report                              # noqa: E402
from pyfp3d.post.surface import wall_force_coefficients                 # noqa: E402
from pyfp3d.solve.newton import solve_newton_lifting                    # noqa: E402

LEVEL, C_ARM, SEED = "medium", 1.0, 0        # the cheaper of the two clamping legs
CANON = ("xcoarse", "coarse", "medium")      # G-SCOPE
MESH = f"{ROOT}/cases/meshes/naca0012_2.5d/{LEVEL}.msh"
M_INF, ALPHA, GAMMA = 0.80, 1.25, 1.4
ANCHOR_NLIM = 1135                            # G-REPRO, committed m1_gate_default.csv
BANDS = (("LE", 0.0, 0.1), ("MID", 0.1, 0.5), ("SHOCK", 0.5, 0.7), ("TE", 0.7, 1.01))
SUMMARY = []


def _record(tag, metric, band, measured, verdict):
    SUMMARY.append((tag, metric, band, measured, verdict))
    print(f"  [{tag}] {metric}:\n        band={band}\n        measured={measured}\n"
          f"        -> {verdict}", flush=True)


def main():
    os.makedirs(RES, exist_ok=True)
    t0 = time.perf_counter()
    assert LEVEL in CANON, "G-SCOPE"
    for v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        print(f"  G-THREADS  {v} = {os.environ.get(v, '<unset>')}")
    print(f"  G-THREADS  load average {os.getloadavg()}")
    kw = dict(NACA_KW); kw["upwind_c"] = C_ARM; kw["n_picard_seed"] = SEED
    print(f"  G-RECIPE   NACA_KW + upwind_c={C_ARM} n_picard_seed={SEED}")
    print(f"  G-SCOPE    level={LEVEL} (canonical ladder; fine is NOT touched)",
          flush=True)

    mc, wc = cut_wake(read_mesh(MESH))
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA, **kw)
    wall = time.perf_counter() - t0
    phi = np.asarray(r["phi"], float)
    nlim, nflr = int(r["n_limited"]), int(r["n_floored"])
    #: G-CACHE -- everything on disk BEFORE any number is reported
    np.savez_compressed(
        os.path.join(RES, "medium_c10_s0.npz"), phi=phi,
        conv=bool(r["converged"]), nlim=nlim, nflr=nflr,
        gamma=np.asarray(r.get("gamma", []), float),
        res_hist=np.asarray(r["residual_history"], float),
        clamp_hist=np.asarray(r["clamp_history"], float),
        limited_mask=np.asarray(r.get("limited_mask", []), float))
    print(f"  solve {wall:.1f}s  conv={r['converged']} clamps={nlim}/{nflr} "
          f"steps={r['n_newton']}", flush=True)

    _record("G-REPRO", "n_limited reproduces the committed run",
            f"{ANCHOR_NLIM}", f"{nlim}",
            "G-REPRO PASS" if nlim == ANCHOR_NLIM else
            f"★ G-REPRO FAIL ({nlim} vs {ANCHOR_NLIM}) -- kill criterion 6")

    # ---- V-WHERE ----------------------------------------------------------
    nodes, elements = mc.nodes, mc.elements
    B, _ = precompute_element_geometry(nodes, elements)
    g = np.empty((len(elements), 3)); q2 = np.empty(len(elements))
    element_velocity_q2(elements, B, phi, g, q2)
    m2 = mach_number_squared(q2, M_INF, GAMMA)
    cent = nodes[elements].mean(axis=1)
    #: ★★ Identify limiter-active elements with the limiter's OWN comparison, read from
    #: pyfp3d/physics/isentropic.py::limit_q2_field: m_cap is a cap on the local MACH,
    #: and the test is q2 >= q2_at_mach(m_cap, ...). My first draft wrote
    #: `mach_number_squared(q2) >= m_cap`, which compares M^2 against a cap on M -- wrong
    #: quantity -- and read a key name (`m_cap_value`) that does not exist. Either would
    #: have produced a plausible-looking but meaningless band histogram.
    import inspect
    _mcap_def = inspect.signature(solve_newton_lifting).parameters["m_cap"].default
    m_cap = float(kw.get("m_cap", _mcap_def))
    cap_q2 = float(q2_at_mach(m_cap, M_INF, GAMMA))
    lim = q2 >= cap_q2 * (1.0 - 1e-12)
    print(f"  m_cap={m_cap} (local MACH) -> q2 cap {cap_q2:.6f}; elements at/over: "
          f"{int(lim.sum())}  (solver counted {nlim})", flush=True)
    src, n_src = (lim, int(lim.sum()))
    rows = []
    if n_src:
        x = cent[src, 0]
        for nm, lo, hi in BANDS:
            k = int(((x >= lo) & (x < hi)).sum())
            rows.append({"band": nm, "x_lo": lo, "x_hi": hi, "n": k,
                         "frac": round(k / n_src, 4)})
        peak = float(x[np.argmax(m2[src])])
        shock_mid = sum(r_["frac"] for r_ in rows if r_["band"] in ("SHOCK", "MID"))
        le = next(r_["frac"] for r_ in rows if r_["band"] == "LE")
        for r_ in rows:
            print(f"    {r_['band']:6} x/c [{r_['x_lo']:.2f},{r_['x_hi']:.2f})  "
                  f"n={r_['n']:5d}  {100*r_['frac']:5.1f} %")
        _record("V-WHERE", "where the limiter-active cells sit, by x/c band",
                "SHOCK+MID >= 70% and LE ~ 0 => C=1.0 is insufficient upwinding;  "
                "LE >= 50% => a geometric hot spot",
                "  ".join(f"{r_['band']} {100*r_['frac']:.1f}%" for r_ in rows)
                + f";  peak M^2 at x/c {peak:.4f}",
                "V-WHERE: insufficient upwinding at C=1.0" if shock_mid >= 0.70
                and le < 0.05 else
                "★ V-WHERE: LE-dominated => geometric hot spot" if le >= 0.50 else
                "V-WHERE: UNDEFINED -- neither pattern dominates (recorded, not attributed)")
    else:
        _record("V-WHERE", "limiter-active cells", "n/a", "none found",
                "★ V-WHERE: could not locate -- m_cap identification failed")
    with open(os.path.join(RES, "clamp_bands.csv"), "w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # ---- V-CL (RECORDED, must not be used) --------------------------------
    dz = float(np.ptp(nodes[:, 2]))
    fo = wall_force_coefficients(nodes, elements, mc.boundary_faces["wall"], phi,
                                 alpha_deg=ALPHA, m_inf=M_INF, s_ref=dz)
    cur = wall_cp_curve(mc, phi, z=0.5 * dz, m_inf=M_INF)
    xs = shock_report(cur, M_INF)["upper"].get("x_shock")
    _record("V-CL", "the clamped state's own numbers",
            "RECORDED ONLY -- GS1.4: a clamped state's numbers are NOT results",
            f"cl_p {float(fo['cl']):.5f}  x_shock {xs}", "RECORDED (must not be used)")

    _record("V-READ", "the mode decomposition already committed",
            "zero-compute; declared in advance to be nearly a deduction",
            "C=1.0 clamping on BOTH seeds (1135/0, 2587/391); C=3.0/seed0 limit_cycle "
            "(0/0); C=1.5 clean on both => (c) blocked at BOTH ENDS of its sweep; "
            "seed 0 has 1/3 usable C values, seed 5 has 2/3 => neither can populate (c)",
            "V-READ: confirmed from bench/gate_results/m1_gate_default.csv")

    with open(os.path.join(RES, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tag", "metric", "band", "measured", "verdict"])
        w.writerows(SUMMARY)
    print(f"\n  {time.perf_counter() - t0:.1f} s", flush=True)
    print("\n★ kill criterion 5: (c) is NOT re-specified here. This round reports that "
          "its lower end is unsupported at medium; the specification is the user's.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
