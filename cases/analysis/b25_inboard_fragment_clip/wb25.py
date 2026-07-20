"""B25 shared harness: inboard-fragment-clip LS solves + the pre-registered
F1 metrics (PRE_REGISTRATION.md).

A-side (control) = the CURRENT q >= 0 clip; the committed b23 D1 caches are
loaded read-only and measured with the SAME code as the C side. The b23 set
has NO coarse alpha=2 leg: that single control leg is solved with the
unchanged default path (bit-identical code, pre-reg S6) and cached HERE --
the committed b23 artifacts stay untouched.

C-side (treatment) = CutElementMap(inboard_clip=make_inboard_clip(FUS)):
the sheet runs inboard until it hits the fuselage surface (trace on the
wall) or, aft of the body, the z = 0 symmetry plane -- the conforming
fragment topology. No knob (a clean A/C pair).

Pre-registered invariants asserted here (via wb24.measure_e1):
  * n_te_nodes stays at the M2-locked counts (76 coarse / 150 medium) --
    the clip must not touch the TE set.
"""

import itertools
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B23 = REPO_ROOT / "cases/analysis/b23_junction_discriminator"
B24 = REPO_ROOT / "cases/analysis/b24_wake_inboard_end"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(B23))
sys.path.insert(0, str(B24))

from wb_common import FUS, LS_MESH_DIR, SOLVE_KW, build_ls, load_mesh  # noqa: E402
from wb24 import ALPHA_REF, M_INF, measure_e1  # noqa: E402
from pyfp3d.meshgen.fuselage import make_inboard_clip  # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI  # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

X_TAIL_TIP = float(FUS.x_tail_tip)      # 2.4178: body end -> symmetry strip
Z_JUNC = float(FUS.r_f)
ROOT_FRAC = 0.25     # inboard quarter of the span = "root" for the profile
                     # distortion metric (pre-reg 2.3.3)


def build_ls_clip(mesh, alpha_deg: float):
    """WakeLevelSet on the DEFAULT 2-point TE polyline (the ruled sheet is
    unchanged; only the cut classification moves) + CutElementMap with the
    fragment clip + MultivaluedOperator."""
    a = np.radians(alpha_deg)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]),
                       inboard_clip=make_inboard_clip(FUS))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return wls, cm, mvop


def _solve_and_cache(mesh, mvop, cache: Path, tag: str, alpha_deg: float,
                     m_inf: float):
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        rec = dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                   wall_s=float(d["wall"]), n=int(d["n"]), conv=bool(d["conv"]),
                   gamma=float(d["gamma"]), nlim=int(d["nlim"]),
                   nflr=int(d["nflr"]), cached=True)
        print(f"  [solve {tag} a={alpha_deg}] CACHED", flush=True)
        return rec
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
    return rec


def solve_side(mesh, level: str, side: str, alpha_deg: float,
               m_inf: float = M_INF):
    """Cached solve for one F1 leg. side "A" = q>=0 control (b23 D1 cache,
    read-only; the missing coarse alpha=2 control is solved with the
    unchanged default path into b25's own cache). side "C" = fragment clip.
    Returns (rec, mvop, wls)."""
    if side == "A":
        cache_b23 = B23 / "results" / f"d1_{level}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
        _, _, mvop = build_ls(mesh, alpha_deg)
        wls = mvop.levelset
        if cache_b23.exists():
            d = np.load(cache_b23, allow_pickle=True)
            rec = dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                       wall_s=float(d["wall"]), n=int(d["n"]),
                       conv=bool(d["conv"]),
                       gamma=float(np.mean(mvop.te_jump(d["phi_ext"]))),
                       nlim=-1, nflr=-1, cached=True)
            print(f"  [solve A {level} a={alpha_deg}] CACHED (b23 d1)",
                  flush=True)
            return rec, mvop, wls
        cache = OUT / f"f1_a_{level}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
        print(f"  [solve A {level} a={alpha_deg}] no b23 cache -- solving "
              f"the default-path control", flush=True)
        rec = _solve_and_cache(mesh, mvop, cache, f"A {level}", alpha_deg,
                               m_inf)
        return rec, mvop, wls

    wls, cm, mvop = build_ls_clip(mesh, alpha_deg)
    cache = OUT / f"f1_c_{level}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
    rec = _solve_and_cache(mesh, mvop, cache, f"C {level}", alpha_deg, m_inf)
    return rec, mvop, wls


def _min_dihedral_deg(p: np.ndarray) -> float:
    """Smallest interior dihedral angle (deg) of tet p (4, 3). Face-normal
    dot products are invariant under a global orientation flip, so no
    winding assumption is needed."""
    faces = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
    n = np.array([np.cross(p[b] - p[a], p[c] - p[a]) for a, b, c in faces])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    best = np.pi
    for i, j in itertools.combinations(range(4), 2):
        cos_ij = float(np.clip(n[i] @ n[j], -1.0, 1.0))
        best = min(best, np.pi - np.arccos(cos_ij))
    return float(np.degrees(best))


def measure_f1(mesh, mvop, wls, phi_ext, alpha_deg, level: str,
               m_inf: float = M_INF):
    """The pre-registered F1 metric set = b24 measure_e1 (b23-consistent
    corridor definitions, n_te assertion) + the B25 additions
    (PRE_REGISTRATION.md 2.3, items 2/6/7)."""
    rec = measure_e1(mesh, mvop, wls, phi_ext, alpha_deg, level, m_inf=m_inf)
    el = mesh.elements.astype(np.int64)
    cents = mesh.nodes[el].mean(axis=1)

    # (2) q station of the corridor peak element (vs q = 0 / the body trace)
    m2 = np.asarray(mvop.element_mach2(phi_ext, m_inf, 1.4, 1.0))
    corr = (cents[:, 2] < 0.5) & (cents[:, 0] > 0.8)
    if corr.any():
        i0 = int(np.argmax(np.where(corr, m2, -1.0)))
        rec["corr_q"] = float(wls.evaluate(cents[i0:i0 + 1])[2][0])
    else:
        rec["corr_q"] = np.nan

    # (6) census: cut set split into the wing band and the NEW inboard
    # strips (beside the body / on the symmetry strip aft of it), plus the
    # aux DOFs the new strips put on the symmetry boundary (risk S2)
    cen_cut = cents[mvop.cm.cut_elems]
    q_cut = wls.evaluate(cen_cut)[2] if len(cen_cut) else np.empty(0)
    inb = q_cut < 0.0
    body_side = cen_cut[:, 0] <= X_TAIL_TIP if len(cen_cut) else inb
    rec["n_cut_wing"] = int((~inb).sum())
    rec["n_cut_inboard_body"] = int((inb & body_side).sum())
    rec["n_cut_inboard_sym"] = int((inb & ~body_side).sum())
    sym = mesh.boundary_faces.get("symmetry")
    rec["n_aux_symmetry"] = (
        int((mvop.cm.ext_dof_of_node[np.unique(sym)] >= 0).sum())
        if sym is not None else -1)

    # (7) sliver watch on the inboard strip: the sheet grazes the body
    # surface at the trace (risk S1) -- min/5th-percentile dihedral and
    # volume of the strip's cut tets
    strip = mvop.cm.cut_elems[inb]
    if len(strip):
        pts = mesh.nodes[el[strip]]
        dih = np.array([_min_dihedral_deg(p) for p in pts])
        mat = np.stack([pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0],
                        pts[:, 3] - pts[:, 0]], axis=1)
        vol = np.abs(np.linalg.det(mat)) / 6.0
        rec["sliver_dih_min_deg"] = float(dih.min())
        rec["sliver_dih_p05_deg"] = float(np.percentile(dih, 5))
        rec["sliver_vol_min"] = float(vol.min())
        rec["sliver_vol_p05"] = float(np.percentile(vol, 5))
    else:
        rec["sliver_dih_min_deg"] = np.nan
        rec["sliver_dih_p05_deg"] = np.nan
        rec["sliver_vol_min"] = np.nan
        rec["sliver_vol_p05"] = np.nan

    # (8) strip-jump ANCHOR watch (pre-reg Appendix A3, Claude audit): the
    # strip's wake-LS convection characteristics run UPSTREAM onto the
    # fuselage surface -- neither the TE (Kutta anchor) nor the far field
    # (pin_gamma anchor) -- so the closure has no inflow data there, and
    # Picard can absorb the near-singular aux rows and still "converge"
    # (B16 precedent: coarse |jump| 53.4 vs Gamma ~ 0.06). jump at an aux
    # node = phi_u - phi_l = -node_side * (phi_ext[aux] - phi_ext[main]).
    if len(strip):
        aux_nodes = np.unique(el[strip].ravel())
        aux_nodes = aux_nodes[mvop.cm.ext_dof_of_node[aux_nodes] >= 0]
    else:
        aux_nodes = np.empty(0, dtype=np.int64)
    if len(aux_nodes):
        aux_dof = mvop.cm.ext_dof_of_node[aux_nodes]
        jump = (-mvop.cm.node_side[aux_nodes].astype(np.float64)
                * (phi_ext[aux_dof] - phi_ext[aux_nodes]))
        gamma_ref = abs(float(np.mean(mvop.te_jump(phi_ext))))
        rec["strip_aux_jump_max"] = float(np.max(np.abs(jump)))
        rec["strip_aux_jump_p95"] = float(np.percentile(np.abs(jump), 95))
        rec["strip_aux_jump_over_gamma"] = (
            rec["strip_aux_jump_max"] / gamma_ref
            if gamma_ref > 1e-12 else np.nan)
    else:
        rec["strip_aux_jump_max"] = np.nan
        rec["strip_aux_jump_p95"] = np.nan
        rec["strip_aux_jump_over_gamma"] = np.nan
    return rec


def te_profile(mesh, mvop, phi_ext):
    """(z, jump) at the TE nodes sorted by z -- the root-section profile for
    the pairwise distortion metric (pre-reg 2.3.3)."""
    z = mesh.nodes[mvop.cm.te_nodes, 2]
    j = np.asarray(mvop.te_jump(phi_ext), dtype=np.float64)
    order = np.argsort(z)
    return z[order], j[order]


def root_distortion(z_j_a, jump_a, z_j_c, jump_c):
    """max |dJump| / max |jump_A| over the inboard quarter of the span."""
    assert np.array_equal(z_j_a, z_j_c), "TE node sets differ (clip must "
    "not touch the TE set)"
    root = z_j_a <= Z_JUNC + ROOT_FRAC * (B_SEMI - Z_JUNC)
    denom = float(np.max(np.abs(jump_a[root]))) if root.any() else 0.0
    if denom <= 1e-12:
        return np.nan
    return float(np.max(np.abs(jump_c[root] - jump_a[root])) / denom)
