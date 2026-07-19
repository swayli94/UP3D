"""B24 shared harness: extended-TE-polyline LS solves + the pre-registered
E1 metrics (PRE_REGISTRATION.md).

A-side (control) = the CURRENT 2-point polyline; all A solves are committed
already (B23 D1 caches) -- this module LOADS them read-only and measures
them with the SAME code as the B side for a same-extractor comparison.
B-side (treatment) = te_polyline(extend="waterline") (B1; delta>0 = B3).

Pre-registered invariants asserted here:
  * n_te_nodes stays at the M2-locked counts (76 coarse / 150 medium) --
    the wall_nodes filter must keep the waterline nodes out of the TE set
    (R3, recorded as excluded pre-run).
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B23 = REPO_ROOT / "cases/analysis/b23_junction_discriminator"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(B23))

from wb_common import FUS, LS_MESH_DIR, SOLVE_KW, build_ls, load_mesh  # noqa: E402
from run_w2_decomp import (BW_LADDER, X_BAND0,  # noqa: E402
                           fuselage_cl_parts)
from pyfp3d.meshgen.wingbody import te_polyline  # noqa: E402
from pyfp3d.post.surface import planform_area  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

M_INF = 0.5
X_TAIL = float(FUS.x_tail_start)              # 2.3878: waterline -> strip join
Z_JUNC = float(FUS.r_f)
TE_NODES_BASELINE = {"coarse": 76, "medium": 150}
ALPHA_REF = 3.06


def x_far_of(mesh) -> float:
    """Far-field x extent of THIS mesh (the strip's end station)."""
    ff = np.unique(mesh.boundary_faces["farfield"])
    return float(mesh.nodes[ff][:, 0].max())


def build_ls_ext(mesh, alpha_deg: float, mode: str = "waterline",
                 delta: float = 0.0):
    """WakeLevelSet on the EXTENDED te_polyline + CutElementMap + MvOp."""
    a = np.radians(alpha_deg)
    te = te_polyline(FUS, extend=mode, delta=delta, x_far=x_far_of(mesh))
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return wls, cm, mvop


def solve_side(mesh, level: str, side: str, alpha_deg: float,
               mode: str = "waterline", delta: float = 0.0, m_inf: float = M_INF):
    """Cached solve for one E1 leg. side "A" loads the committed B23 D1 cache
    (read-only; no re-solve); side "B1"/"B3" solves with the extended sheet.
    Returns (rec, mvop, wls) with rec carrying phi_ext + convergence meta."""
    if side == "A":
        cache = B23 / "results" / f"d1_{level}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
        if not cache.exists():
            raise FileNotFoundError(f"A-side cache missing: {cache}")
        d = np.load(cache, allow_pickle=True)
        _, _, mvop = build_ls(mesh, alpha_deg)
        wls = mvop.levelset
        rec = dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                   wall_s=float(d["wall"]), n=int(d["n"]), conv=bool(d["conv"]),
                   gamma=float(np.mean(mvop.te_jump(d["phi_ext"]))),
                   nlim=-1, nflr=-1, cached=True)
        print(f"  [solve A {level} a={alpha_deg}] CACHED (b23 d1)", flush=True)
        return rec, mvop, wls

    tag = f"e1_{side.lower()}_{level}"
    cache = OUT / f"{tag}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
    wls, cm, mvop = build_ls_ext(mesh, alpha_deg, mode=mode, delta=delta)
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        rec = dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                   wall_s=float(d["wall"]), n=int(d["n"]), conv=bool(d["conv"]),
                   gamma=float(d["gamma"]), nlim=int(d["nlim"]),
                   nflr=int(d["nflr"]), cached=True)
        print(f"  [solve {tag} a={alpha_deg}] CACHED", flush=True)
        return rec, mvop, wls
    t0 = time.perf_counter()
    r = solve_multivalued_lifting(mvop, mesh, m_inf, alpha_deg=alpha_deg,
                                  **SOLVE_KW)
    wall_s = time.perf_counter() - t0
    rec = dict(phi_ext=r["phi_ext"], res=list(r["residual_history"]),
               wall_s=wall_s, n=r["n_outer"], conv=bool(r["converged"]),
               gamma=float(r["gamma"]), nlim=int(r["n_limited"]),
               nflr=int(r["n_floored"]), cached=False)
    np.savez(cache, phi_ext=r["phi_ext"],
             res=np.asarray(r["residual_history"]), wall=wall_s,
             n=r["n_outer"], conv=r["converged"], gamma=r["gamma"],
             nlim=r["n_limited"], nflr=r["n_floored"])
    print(f"  [solve {tag} a={alpha_deg}] conv={r['converged']} "
          f"n_outer={r['n_outer']} res={r['residual_history'][-1]:.2e} "
          f"gamma={r['gamma']:.4f} nlim={r['n_limited']} nflr={r['n_floored']} "
          f"({wall_s:.0f}s)", flush=True)
    return rec, mvop, wls


def measure_e1(mesh, mvop, wls, phi_ext, alpha_deg, level: str,
               m_inf: float = M_INF, topn: int = 200):
    """The pre-registered E1 metric set (b23-consistent corridor defs)."""
    m2 = np.asarray(mvop.element_mach2(phi_ext, m_inf, 1.4, 1.0))
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    x, z = cents[:, 0], cents[:, 2]
    s, _, q = wls.evaluate(cents)

    def region(mask):
        if not mask.any():
            return np.nan, np.nan, np.nan, 0
        m2m = np.where(mask, m2, -1.0)
        i0 = int(np.argmax(m2m))
        return float(np.sqrt(m2[i0])), float(x[i0]), float(z[i0]), 0

    corr = (z < 0.5) & (x > 0.8)
    sheet = corr & (np.abs(s) < 0.03)
    tip = z > 0.8
    c_mmax, c_x, c_z, _ = region(corr)
    s_mmax, s_x, s_z, _ = region(sheet)
    t_mmax, t_x, t_z, _ = region(tip)

    order = np.argsort(m2)[::-1][:topn]
    tc = cents[order]
    ss, _, qq = wls.evaluate(tc)

    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    from pyfp3d.post.unified import wall_forces
    cl_p = float(wall_forces(mesh, mvop=mvop, phi_ext=phi_ext,
                             alpha_deg=alpha_deg, s_ref=s_ref, m_inf=m_inf,
                             wall_tag="wall")["cl"])
    w2 = fuselage_cl_parts(mesh, mvop, phi_ext, alpha_deg, s_ref)

    ff = np.unique(mesh.boundary_faces["farfield"])
    n_aux_ff = int((mvop.cm.ext_dof_of_node[ff] >= 0).sum())
    n_te = int(len(mvop.cm.te_nodes))
    base = TE_NODES_BASELINE.get(level)
    if base is not None:
        assert n_te == base, (f"TE-node drift: {n_te} != {base} ({level}) -- "
                              f"the extension must not touch the TE set (R3)")

    rec = dict(
        n_te_nodes=n_te, n_aux_farfield=n_aux_ff,
        n_cut=int(len(mvop.cm.cut_elems)),
        corr_mmax=c_mmax, corr_x=c_x, corr_z=c_z,
        sheet_mmax=s_mmax, sheet_x=s_x, sheet_z=s_z,
        tip_mmax=t_mmax, tip_x=t_x, tip_z=t_z,
        n_sup=int((m2 > 1.0).sum()), n_sup_corr=int((m2[corr] > 1.0).sum()),
        top_med_abs_s=float(np.median(np.abs(ss))),
        top_med_q=float(np.median(qq)),
        top_med_x=float(np.median(tc[:, 0])),
        pocket_peak_x=c_x,
        pocket_past_tail=bool(c_x > X_TAIL),
        cl_p=cl_p,
        cl_fus=w2["all"], cl_fus_band=w2["band_bw0.06"],
        cl_fus_out=w2["out_bw0.06"], cl_fus_poles=w2["poles"],
    )
    return rec
