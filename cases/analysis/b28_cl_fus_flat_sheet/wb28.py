"""B28 shared harness: the FLAT-sheet fragment LS solves (F side).

F side = the C side (inboard fragment clip, same topology) with ONE variable
moved: the sheet GEOMETRY is ruled flat along +x (the conforming position,
y = 0 through the TE) instead of along the freestream, while the jump
convection u_hat keeps the freestream aim (WakeLevelSet sheet_direction
knob, B28; default None bit-identical, locked by
test_b1_cut_elements.py::TestSheetDirection).

    wls = WakeLevelSet(te_polyline(FUS),
                       direction=(cos a, sin a, 0),      # physical (u_hat)
                       sheet_direction=(1, 0, 0))        # geometric ruling

The C-side comparators are the committed b25 caches (read-only, same code
measures them) -- the F-vs-C root-profile distortion needs the C te_profile.
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
B23 = REPO_ROOT / "cases/analysis/b23_junction_discriminator"
B24 = REPO_ROOT / "cases/analysis/b24_wake_inboard_end"
B25 = REPO_ROOT / "cases/analysis/b25_inboard_fragment_clip"
for p in (REPO_ROOT, B23, B24, B25):
    sys.path.insert(0, str(p))

from wb_common import FUS, SOLVE_KW, load_mesh  # noqa: E402
from wb24 import M_INF  # noqa: E402
from wb25 import measure_f1  # noqa: E402
from pyfp3d.meshgen.fuselage import make_inboard_clip  # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

B25_RESULTS = B25 / "results"


def build_ls_flat(mesh, alpha_deg: float):
    """Fragment clip (C topology) + FLAT sheet geometry (sheet_direction
    = +x); the convection direction stays the freestream."""
    a = np.radians(alpha_deg)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0),
                       sheet_direction=(1.0, 0.0, 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]),
                       inboard_clip=make_inboard_clip(FUS))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return wls, cm, mvop


N_OUTER_BUDGET = 220   # R1 remedy: the FLAT sheet converges geometrically
                       # at ratio ~0.94/outer (C/tilted: ~0.77), so the B25
                       # 80-outer cap is structurally insufficient. Same
                       # physics recipe otherwise (SOLVE_KW verbatim).


def solve_flat(mesh, level: str, alpha_deg: float, m_inf: float = M_INF):
    """Cached F-side solve; a non-converged cache is CONTINUED (phi_init
    warm start, full extended state incl. the wake jump) up to the total
    N_OUTER_BUDGET. Returns (rec, mvop, wls)."""
    wls, cm, mvop = build_ls_flat(mesh, alpha_deg)
    cache = OUT / f"f2_f_{level}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
    phi_init, res_prev, n_done = None, [], 0
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        if bool(d["conv"]):
            rec = dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                       wall_s=float(d["wall"]), n=int(d["n"]),
                       conv=True, gamma=float(d["gamma"]),
                       nlim=int(d["nlim"]), nflr=int(d["nflr"]), cached=True)
            print(f"  [solve F {level} a={alpha_deg}] CACHED", flush=True)
            return rec, mvop, wls
        phi_init, res_prev, n_done = d["phi_ext"], list(d["res"]), int(d["n"])
        print(f"  [solve F {level} a={alpha_deg}] warm start from outer "
              f"{n_done} (res {res_prev[-1]:.2e})", flush=True)
    kw = dict(SOLVE_KW)
    kw["n_outer_max"] = N_OUTER_BUDGET - n_done
    t0 = time.perf_counter()
    r = solve_multivalued_lifting(mvop, mesh, m_inf, alpha_deg=alpha_deg,
                                  phi_init=phi_init, **kw)
    wall_s = time.perf_counter() - t0
    res_full = res_prev + list(r["residual_history"])
    rec = dict(phi_ext=r["phi_ext"], res=res_full,
               wall_s=wall_s, n=n_done + r["n_outer"],
               conv=bool(r["converged"]),
               gamma=float(r["gamma"]), nlim=int(r["n_limited"]),
               nflr=int(r["n_floored"]), cached=False)
    np.savez(cache, phi_ext=r["phi_ext"],
             res=np.asarray(res_full), wall=wall_s,
             n=rec["n"], conv=r["converged"], gamma=r["gamma"],
             nlim=r["n_limited"], nflr=r["n_floored"])
    print(f"  [solve F {level} a={alpha_deg}] conv={r['converged']} "
          f"n_outer={rec['n']} res={res_full[-1]:.2e} "
          f"gamma={r['gamma']:.4f} nlim={r['n_limited']} nflr={r['n_floored']} "
          f"({wall_s:.0f}s)", flush=True)
    return rec, mvop, wls


def load_c_te_profile(mesh, level: str, alpha_deg: float,
                      m_inf: float = M_INF):
    """The committed b25 C-side te profile (read-only cache; the C mvop is
    rebuilt with the unchanged b25 code)."""
    from wb25 import build_ls_clip, te_profile
    cache = B25_RESULTS / f"f1_c_{level}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
    d = np.load(cache, allow_pickle=True)
    _, _, mvop = build_ls_clip(mesh, alpha_deg)
    return te_profile(mesh, mvop, d["phi_ext"])
