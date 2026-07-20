"""B26 shared harness: pocket-healed LS transonic ceiling re-measurement
(PRE_REGISTRATION.md, 2026-07-19).

A-side = inboard_clip=None (bit-identical == the B18 recipe re-run on the
CURRENT code; the primary control per the P14 same-code discipline, with the
committed B18/GB20.5 numbers as the historical anchor). C-side = the B25
inboard fragment clip (conforming sheet topology). The B18 ramp recipe is
frozen verbatim (pre-reg 2.1):

    LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
                      freeze_tol=1e-4, freeze_max_clamped=8,
                      intermediate_tol=1e-3, n_seed=30,
                      direct_refactor_every=1000, n_newton_max=80)

alpha = 3.06 fixed, m_start=0.50 -> m_target=0.84, dm=0.05, dm_min=0.01
(default), upwind defaults (c=1.5, m_crit=0.95, m_cap=3.0, rho_floor=0.05).

Zero library changes: build_ls (A) is wb_common's default clip; build_ls_clip
(C) is B25's (CutElementMap(inboard_clip=make_inboard_clip(FUS))). Ramp
results cache to results/g1_*.npz (gitignored); committed evidence =
g1_summary.csv / g1_levels.csv / g1_peaks.csv / g1_ceiling.png.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B23 = REPO_ROOT / "cases/analysis/b23_junction_discriminator"
B24 = REPO_ROOT / "cases/analysis/b24_wake_inboard_end"
B25 = REPO_ROOT / "cases/analysis/b25_inboard_fragment_clip"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(B23))
sys.path.insert(0, str(B24))
sys.path.insert(0, str(B25))

from wb_common import (FUS, LS_MESH_DIR, TE_JUNC, Z_JUNC, build_ls,  # noqa: E402
                       load_mesh, measure)
from wb25 import build_ls_clip, measure_f1  # noqa: E402
from pyfp3d.meshgen.fuselage import radius_at  # noqa: E402
from pyfp3d.solve.continuation import mach_schedule  # noqa: E402
from pyfp3d.solve.newton_ls import solve_multivalued_newton_transonic  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

ALPHA = 3.06
M_START = 0.50
M_TARGET = 0.84
DM = 0.05
LS_RAMP_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
                  freeze_tol=1e-4, freeze_max_clamped=8,
                  intermediate_tol=1e-3, n_seed=30,
                  direct_refactor_every=1000, n_newton_max=80)
NOMINAL = [float(m) for m in mach_schedule(M_TARGET, m_start=M_START, dm=DM)]

# conforming cl_p anchors (B18 committed cl_vs_mach.csv) for the trend check
CONF_CL = {0.50: 0.2173, 0.65: 0.2321, 0.79: 0.2579, 0.84: 0.2617}
# B18 committed ceiling anchors (checks.csv): death Mach / Mmax per level
B18_DEATH = {"coarse": 0.55, "medium": 0.50}
B18_DEATH_MMAX = {"coarse": 1.31, "medium": 5.22}
# GB20.5 same-code A-medium anchor (c1 legb_b18_hypothesis.csv, main row):
# m_final=0.5, m_last=None, Mmax=5.2201, nlim/nflr 3/3, res 1.14e-13, ~807 s


def build_side(mesh, side: str):
    """A = default q>=0 clip (bit-identical B18); C = B25 fragment clip.
    Returns (wls, cm, mvop)."""
    if side == "A":
        return build_ls(mesh, ALPHA)
    return build_ls_clip(mesh, ALPHA)


_LV_KEYS = ("m_inf", "tag", "gamma", "cl_kj", "mach_max", "n_newton",
            "converged", "accept_reason", "residual_norm", "n_limited",
            "n_floored", "froze", "n_freeze_refresh", "n_freeze_reverts",
            "n_refactor", "n_schur_fallback", "wall_s", "n_lin_iters",
            "n_lin_solves")


def _slim_levels(levels):
    """Per-level scalar records (the heavy timings/residual_history/
    step_records blobs are dropped; the classification needs the scalars)."""
    return [{k: lv[k] for k in _LV_KEYS} for lv in levels]


def run_ramp(mesh, level: str, side: str):
    """The frozen B18 freeze-ramp on one leg, cached to
    results/g1_{side}_{level}.npz. Returns (rec, wls, mvop); rec holds the
    honesty fields (target_reached / m_last / m_final), the slimmed per-level
    records and the FINAL level's phi_ext (living at m_final -- the dying
    level's state when the ramp fails, which is exactly the state carrying
    the nlim/nflr that killed it)."""
    wls, cm, mvop = build_side(mesh, side)
    cache = OUT / f"g1_{side}_{level}.npz"
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        rec = dict(phi_ext=d["phi_ext"],
                   levels=json.loads(str(d["levels_json"])),
                   mach_schedule=[float(x) for x in d["mach_schedule"]],
                   target_reached=bool(d["target_reached"]),
                   m_last=(float(d["m_last"])
                           if np.isfinite(d["m_last"]) else None),
                   m_final=float(d["m_final"]), wall_s=float(d["wall_s"]),
                   cached=True)
        print(f"  [ramp {side} {level}] CACHED", flush=True)
        return rec, wls, mvop
    t0 = time.perf_counter()
    r = solve_multivalued_newton_transonic(
        mvop, mesh, M_TARGET, alpha_deg=ALPHA, m_start=M_START, dm=DM,
        verbose=True, **LS_RAMP_KW)
    wall_s = time.perf_counter() - t0
    rec = dict(phi_ext=r["phi_ext"], levels=_slim_levels(r["levels"]),
               mach_schedule=[float(x) for x in r["mach_schedule"]],
               target_reached=bool(r["target_reached"]),
               m_last=(float(r["m_last_converged"])
                       if r["m_last_converged"] is not None else None),
               m_final=float(r["m_final"]), wall_s=wall_s, cached=False)
    np.savez(cache, phi_ext=r["phi_ext"],
             levels_json=json.dumps(rec["levels"]),
             mach_schedule=np.asarray(rec["mach_schedule"]),
             target_reached=rec["target_reached"],
             m_last=(rec["m_last"] if rec["m_last"] is not None else np.nan),
             m_final=rec["m_final"], wall_s=wall_s)
    print(f"  [ramp {side} {level}] reached={rec['target_reached']} "
          f"m_last={rec['m_last']} m_final={rec['m_final']} "
          f"({wall_s:.0f}s)", flush=True)
    return rec, wls, mvop


def measure_leg(mesh, wls, mvop, phi_ext, m_inf: float, level: str):
    """Final-state metrics at m_inf: the B25 F1 set (corridor corrM/n_sup/
    peak (x,z,q), strip aux-jump anchor watch, slivers, census, the 76/150
    n_te assertion) + the b23 topk machine's GLOBAL peak attribution
    (dist_fus_surface / z-z_junc / dist_te_junc)."""
    rec = measure_f1(mesh, mvop, wls, phi_ext, ALPHA, level, m_inf=m_inf)
    pk, _top = measure(mesh, mvop, phi_ext, ALPHA, m_inf=m_inf)
    rec.update(pk_mmax=pk["mmax"], pk_x=pk["argmax_x"], pk_y=pk["argmax_y"],
               pk_z=pk["argmax_z"], pk_dist_fus=pk["argmax_dist_fus"],
               pk_z_minus_zjunc=pk["argmax_z_minus_zjunc"],
               pk_dist_te_junc=pk["argmax_dist_te_junc"])
    return rec


def corridor_peaks(mesh, wls, mvop, phi_ext, m_inf: float, topk: int = 5):
    """b23 topk machine restricted to the junction corridor (z<0.5, x>0.8):
    separates a junction-strip peak from a legitimate wing-shock peak (pre-reg
    metric 4: at transonic the corridor legitimately contains the wing shock,
    so judge by SITE, not by presence)."""
    m2 = np.asarray(mvop.element_mach2(phi_ext, m_inf, 1.4, 1.0))
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    corr = (cents[:, 2] < 0.5) & (cents[:, 0] > 0.8)
    idx = np.where(corr)[0]
    idx = idx[np.argsort(m2[idx])[::-1][:topk]]
    rows = []
    if not len(idx):
        return rows
    q = wls.evaluate(cents[idx])[2]
    for rank, (i, qi) in enumerate(zip(idx, q), start=1):
        c = cents[i]
        rr = float(np.hypot(c[1], c[2]))
        rows.append(dict(rank=rank, mach=float(np.sqrt(m2[i])),
                         x=float(c[0]), y=float(c[1]), z=float(c[2]),
                         q=float(qi),
                         dist_fus_surface=abs(rr - radius_at(FUS, float(c[0]))),
                         z_minus_zjunc=float(c[2] - Z_JUNC),
                         dist_te_junc=float(np.linalg.norm(c - TE_JUNC))))
    return rows


def classify(rec):
    """Pre-registered failure taxonomy (2.3 item 9):
      (a) strict-gate rejection -- Newton residual tight (< 1e-6) but
          nlim/nflr > 0 (the B18 signature);
      (b) Newton non-convergence (residual stall / freeze reverts exhausted);
      (c) the dm-halving cascade fired -- suffix "+dm" (death after the
          cascade = dm_min exhausted)."""
    if rec["target_reached"]:
        return "reached"
    lv = rec["levels"][-1]
    halved = len(rec["mach_schedule"]) > len(NOMINAL)
    tight = lv["residual_norm"] < 1e-6
    clamps = lv["n_limited"] + lv["n_floored"]
    cls = "a" if (tight and clamps > 0) else "b"
    return cls + ("+dm" if halved else "")
