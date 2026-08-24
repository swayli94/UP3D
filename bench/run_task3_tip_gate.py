"""Tip-gate redesign: sweep the taper RADIUS and read the LE-RMS SPANWISE PROFILE.

Pre-registered in phases/p3/docs/dev_phase_three/20260813-2100-tip-gate-redesign-prereg.md.

★ Why this gate and not the previous one: GS2.1's tip gate used peak LOCATION as a proxy for the tip
free-edge singularity, and its own author recorded why that failed -- the peak MOVED to the outboard
LEADING edge, a different object, which the criterion still counted as "in the tip band". So this gate
uses no proxy: it does a dose-response on the quantity being explained (the LE-RMS spanwise profile)
using an EXISTING knob (the taper radius), with a dimensionless self-normalising criterion (the
outboard/inboard ratio) instead of a hand-picked threshold like M_max < 1.6.

★★ r_c = 0 is deliberately NOT in the scan: that configuration is already measured not to converge,
and a non-solution is not a reading.

Verdict: phases/p3/docs/dev_phase_three/20260813-2300-tip-gate-verdict.md -- W2, the tip is not the driver
(delta-R = 0.0027 while cl_p moves 3.13 %, so the knob demonstrably reaches the solve).

Outputs (TRACKED): bench/gate_results/task3_tip_gate.csv
"""
import os, sys, time, csv
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.surface import planform_area, wall_force_coefficients
from run_m3_budget import (ALPHA, B_SEMI, ETAS, M_INF, N_UNMASKED, parse_experiment, solve)
from run_task3_le_registration import le_band_rms

RCS = (0.025, 0.05, 0.10)          # declared in advance
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "gate_results", "task3_tip_gate.csv")
print("resolved threads: " + ", ".join(f"{k}={os.environ.get(k)}" for k in
      ("NUMBA_NUM_THREADS","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS")))
print(f"load average: {os.getloadavg()}\n")
exp = parse_experiment()
mc, wc = cut_wake(read_mesh(os.path.join(HERE, "..", "cases", "meshes",
                                        "onera_m6", "medium.msh")))
s_ref = planform_area(mc.nodes, mc.boundary_faces["wall"])
etas = ETAS[:N_UNMASKED]
rows=[]
for rc in RCS:
    t0=time.perf_counter()
    r = solve(mc, wc, entropy=True, kutta="pressure", taper=True, probe_seed=0, taper_rc=rc)
    wall=time.perf_counter()-t0
    conv=bool(r.get("converged")); res=r["residual_history"][-1]
    rec=dict(taper_rc=rc, converged=conv, res_final=res,
             n_limited=int(r.get("n_limited") or 0), n_floored=int(r.get("n_floored") or 0),
             wall_s=round(wall,1))
    if conv:
        phi=np.asarray(r["phi"])
        curves={e: section_cp_curve(mc, phi, eta=e, b_semi=B_SEMI, m_inf=M_INF) for e in ETAS}
        for e in etas:
            v,_,ok = le_band_rms(curves, exp, [e], 0.0)
            rec[f"le_{e:.2f}"]=v
        for band in ("LE","MID","TE"):
            v,_,_ = le_band_rms(curves, exp, etas, 0.0, band=band)
            rec[f"rms_{band}_upper"]=v
        f=wall_force_coefficients(mc.nodes, mc.elements, mc.boundary_faces["wall"], phi,
                                  alpha_deg=ALPHA, s_ref=s_ref, m_inf=M_INF)
        rec["cl_p"]=float(f["cl"])
        rec["ratio"]=rec["le_0.90"]/rec["le_0.20"]
    rows.append(rec)
    print(f"  r_c={rc:.3f}  conv={conv!s:5} |R|={res:.2e} lim/flr={rec['n_limited']}/{rec['n_floored']}"
          + (f"  R={rec['ratio']:.4f}  LE={rec['rms_LE_upper']:.6f}  cl_p={rec['cl_p']:.6f}" if conv else "")
          + f"  ({wall:.0f}s)", flush=True)
keys=sorted({k for r in rows for k in r})
with open(OUT,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=keys); w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT}")
