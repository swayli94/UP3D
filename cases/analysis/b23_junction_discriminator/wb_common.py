"""Shared solve + measurement harness for the wing-body junction discriminator.

Protocol (PRE_REGISTRATION.md): LS Picard, M0.5, farfield="freestream" +
farfield_aux="pin_gamma" (the B9 measured wing-body recipe with the B17
arbitrated aux pin). Solves cache to results/*.npz (gitignored); committed
evidence = the CSV/PNG each leg writes.

Localization metric: for the top-K M^2 elements we record the centroid, its
distance to the fuselage revolution surface |sqrt(y^2+z^2) - R(x)|, its
signed z offset from the junction station z_junc = r_f, and its distance to
the junction TE point -- enough to say "on the crease / on the wing / on the
fuselage / in the volume" without a curve-distance query.
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import read_mesh                      # noqa: E402
from pyfp3d.meshgen.fuselage import FuselageParams, radius_at  # noqa: E402
from pyfp3d.meshgen.wing3d import x_te                        # noqa: E402
from pyfp3d.meshgen.wingbody import junction_z, te_polyline   # noqa: E402
from pyfp3d.post.surface import planform_area                 # noqa: E402
from pyfp3d.post.unified import wall_forces                   # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_lifting  # noqa: E402
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

FUS = FuselageParams()
Z_JUNC = junction_z(FUS)
TE_JUNC = np.array([x_te(Z_JUNC), 0.0, Z_JUNC])
M_INF = 0.5
ALPHA_REF = 3.06
LS_MESH_DIR = REPO_ROOT / "cases/meshes/onera_m6_wingbody"

SOLVE_KW = dict(farfield="freestream", farfield_aux="pin_gamma",
                n_outer_max=80, tol_residual=1e-7)


def load_mesh(path: Path):
    return read_mesh(str(path))


def build_ls(mesh, alpha_deg: float):
    """WakeLevelSet + CutElementMap + MultivaluedOperator for the given alpha."""
    a = np.radians(alpha_deg)
    wls = WakeLevelSet(te_polyline(FUS), direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return wls, cm, mvop


def solve_ls(mesh, cache_tag: str, alpha_deg: float, m_inf: float = M_INF):
    """Cached LS Picard solve. Returns (rec, mvop) with rec holding phi_ext,
    convergence meta and wall time. Cache key = cache_tag + alpha + m_inf."""
    cache = OUT / f"{cache_tag}_a{alpha_deg:.2f}_m{m_inf:.2f}.npz"
    _, _, mvop = build_ls(mesh, alpha_deg)
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        rec = dict(phi_ext=d["phi_ext"], res=list(d["res"]),
                   wall_s=float(d["wall"]), n=int(d["n"]), conv=bool(d["conv"]),
                   cached=True)
        print(f"  [solve {cache_tag} a={alpha_deg} m={m_inf}] CACHED", flush=True)
        return rec, mvop
    t0 = time.perf_counter()
    r = solve_multivalued_lifting(mvop, mesh, m_inf, alpha_deg=alpha_deg,
                                  **SOLVE_KW)
    wall_s = time.perf_counter() - t0
    rec = dict(phi_ext=r["phi_ext"], res=list(r["residual_history"]),
               wall_s=wall_s, n=r["n_outer"], conv=bool(r["converged"]),
               cached=False)
    np.savez(cache, phi_ext=r["phi_ext"], res=np.asarray(r["residual_history"]),
             wall=wall_s, n=r["n_outer"], conv=r["converged"])
    print(f"  [solve {cache_tag} a={alpha_deg} m={m_inf}] conv={r['converged']} "
          f"n_outer={r['n_outer']} res={r['residual_history'][-1]:.2e} "
          f"({wall_s:.0f}s)", flush=True)
    return rec, mvop


def measure(mesh, mvop, phi_ext, alpha_deg, m_inf=M_INF, topk=8):
    """Discriminator metrics: Mmax + pocket localization + fuselage lift."""
    m2 = np.asarray(mvop.element_mach2(phi_ext, m_inf, 1.4, 1.0))
    cents = mesh.nodes[mesh.elements].mean(axis=1)
    order = np.argsort(m2)[::-1]
    i0 = int(order[0])

    sup = m2 > 1.0
    sup_bbox = {}
    if sup.any():
        sc = cents[sup]
        sup_bbox = dict(n_sup=int(sup.sum()),
                        sup_x_min=float(sc[:, 0].min()),
                        sup_x_max=float(sc[:, 0].max()),
                        sup_y_min=float(sc[:, 1].min()),
                        sup_y_max=float(sc[:, 1].max()),
                        sup_z_min=float(sc[:, 2].min()),
                        sup_z_max=float(sc[:, 2].max()))
    else:
        sup_bbox = dict(n_sup=0, sup_x_min=np.nan, sup_x_max=np.nan,
                        sup_y_min=np.nan, sup_y_max=np.nan,
                        sup_z_min=np.nan, sup_z_max=np.nan)

    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    cl_p = float(wall_forces(mesh, mvop=mvop, phi_ext=phi_ext,
                             alpha_deg=alpha_deg, s_ref=s_ref, m_inf=m_inf,
                             wall_tag="wall")["cl"])
    cl_fus = float(wall_forces(mesh, mvop=mvop, phi_ext=phi_ext,
                               alpha_deg=alpha_deg, s_ref=s_ref, m_inf=m_inf,
                               wall_tag="fuselage")["cl"])

    top = []
    for i in order[:topk]:
        c = cents[int(i)]
        rr = float(np.hypot(c[1], c[2]))
        top.append(dict(elem=int(i), m2=float(m2[int(i)]),
                        mach=float(np.sqrt(m2[int(i)])),
                        x=float(c[0]), y=float(c[1]), z=float(c[2]),
                        dist_fus_surface=abs(rr - radius_at(FUS, float(c[0]))),
                        z_minus_zjunc=float(c[2] - Z_JUNC),
                        dist_te_junc=float(np.linalg.norm(c - TE_JUNC))))

    c0 = cents[i0]
    rec = dict(mmax=float(np.sqrt(m2[i0])), m2_max=float(m2[i0]),
               argmax_elem=i0,
               argmax_x=float(c0[0]), argmax_y=float(c0[1]),
               argmax_z=float(c0[2]),
               argmax_dist_fus=float(abs(np.hypot(c0[1], c0[2])
                                         - radius_at(FUS, float(c0[0])))),
               argmax_z_minus_zjunc=float(c0[2] - Z_JUNC),
               argmax_dist_te_junc=float(np.linalg.norm(c0 - TE_JUNC)),
               cl_p=cl_p, cl_fus=cl_fus,
               cl_fus_over_wing=abs(cl_fus / cl_p) if cl_p else np.nan,
               **sup_bbox)
    return rec, top
