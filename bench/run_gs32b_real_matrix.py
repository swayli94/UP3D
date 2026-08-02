"""GS3.2 (b) verification: does the strength-measure finding carry to the REAL matrix?

Result (2026-08-02): YES, and the numbers are nearly identical to the pure-Laplacian
sweep -- NACA 2.5-D medium 278 -> 8 iterations (4.20x total time), M6 medium 13 -> 10 at
0.47x, i.e. 2.1x SLOWER. Which additionally shows the anisotropy that matters is
GEOMETRIC (the mesh), not the density weighting.

⚠ Both states come from a single solve_newton_lifting without the Mach ramp and report
conv=False. That is fine for this test and stated rather than hidden: what it needs is a
REPRESENTATIVE A_ff, and a non-converged transonic iterate is exactly what the
preconditioner faces during the ramp.

Original docstring follows.

Does the strength-measure finding carry to the REAL preconditioned matrix?

The sweep measured a pure Laplacian (rho = 1). Production preconditions with the
DENSITY-WEIGHTED Picard matrix A_ff assembled at the current transonic state, which
carries its own anisotropy from the artificial density. Direction should carry; the
numbers need not. This extracts A_ff exactly as newton.py builds it.
"""
import os, sys, time
os.environ.setdefault("NUMBA_NUM_THREADS","16"); os.environ.setdefault("OMP_NUM_THREADS","16")
os.environ.setdefault("OPENBLAS_NUM_THREADS","16")
sys.path.insert(0, "/home/lrz/codes/UP3D"); sys.path.insert(0, "/home/lrz/codes/UP3D/bench")
import numpy as np
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.newton import NewtonWorkspace, solve_newton_lifting
from run_gs32b_amg_sweep import CONFIGS, run_one

CASES = [("naca0012_2.5d", "medium", 0.80, 1.25),
         ("onera_m6", "medium", 0.84, 3.06)]
for fam, lvl, m_inf, alpha in CASES:
    mc, wc = cut_wake(read_mesh(f"/home/lrz/codes/UP3D/cases/meshes/{fam}/{lvl}.msh"))
    kw = dict(m_inf=m_inf, alpha_deg=alpha, precond="amg", n_newton_max=60)
    if fam == "onera_m6":
        kw["farfield_spanwise_gamma"] = True
    r = solve_newton_lifting(mc, wc, **kw)
    ws = NewtonWorkspace(mc, wc, alpha_deg=alpha); ws.set_mach(m_inf)
    phi = np.asarray(r["phi"])[:ws.n_red][ws.free].copy()
    g = np.atleast_1d(np.asarray(r["gamma"]))
    st = ws.eval_residual(phi, g, 1.5, 0.95, 3.0, 0.05)[2]
    # exactly newton.py's amg branch
    A_pic = ws.op.assemble_matrix(st["rho_t"])
    A_red = (ws.con.T.T @ (A_pic @ ws.con.T)).tocsr()
    A = A_red[ws.free][:, ws.free].tocsr()
    print(f"\n=== {fam}/{lvl} M{m_inf}: REAL A_ff, {A.shape[0]} dofs, "
          f"conv={r['converged']} ===", flush=True)
    base = None
    for tag, kwargs in CONFIGS:
        if tag in ("theta0.0", "affinity", "lloyd", "blockgs", "gs2sweep"):
            continue
        res = run_one(A, tag, kwargs)
        if res.get("note"): print(f"  {tag:18s} {res['note']}"); continue
        if tag == "default": base = res
        sp = base["total_s"]/res["total_s"] if base else float("nan")
        print(f"  {tag:18s} {res['cg_iters']:5d} it  setup {res['amg_setup_s']:6.3f}s"
              f"  cg {res['cg_s']:6.3f}s  TOTAL {res['total_s']:6.3f}s  ({sp:.2f}x)",
              flush=True)
