"""
A3 / GA3.5 — C1: does the LS Newton Jacobian match its residual on a 3D mesh?

kimi's 2026-07-17 code review (docs/inspection/20260717-2348-code-review.md,
finding C1) argues that `newton_terms23_side_coo` mis-maps its columns on
MIXED-SIDE PLAIN elements -- uncut elements whose four nodes do not all sit on
one side of the wake level set. They exist only in 3D: elements beyond the
spanwise tip clip, and elements the wake-surface extension passes through
upstream of the TE. `mass_conservation_coo` assembles those with the SIDE-field
density (cut_assembly.py:152-155), so the residual depends on the aux DOFs of
cut nodes through `side_potentials` -- while Terms 2/3 scatter every sensitivity
to MAIN columns (cut_assembly.py:441-456). If the argument holds, J != dR/dphi
there and the LS Newton silently degrades to a quasi-Newton.

B6's FD gate (`tests/test_b6_newton.py::test_newton_jacobian_fd`) cannot see
this: it runs on the quasi-2D NACA mesh, which structurally has no such
elements. This script is that gate's 3D twin. It reuses B6's harness verbatim
(same residual/Jacobian assembly, same central-difference probe) and only
changes the mesh and the probe DIRECTIONS:

  probe "targeted"  -- aux DOFs of cut nodes that touch a mixed-side plain
                       element (where C1 predicts the mismatch)
  probe "control"   -- aux DOFs far from any such element (must stay clean;
                       this is what separates "C1 is real" from "the whole
                       3D Jacobian is off")
  probe "global"    -- all free DOFs, i.e. the B6 gate's own direction

VERDICT LOGIC (pre-registered):
  targeted >> 1e-6 AND control ~ 1e-8  ->  C1 CONFIRMED, localized as claimed
  all three ~ 1e-8                     ->  C1 REFUTED on this mesh
  all three large                      ->  something else is wrong; not C1

This script only MEASURES. Whether to change the column mapping is a separate
decision -- the fix would leave R (hence every converged solution) untouched
but could move the committed step-count trajectories.

Run:  python cases/analysis/c1_ls_jacobian_fd/run_check.py
Cost: ~2-4 min (one coarse M6 Picard seed + 6 residual assemblies).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.constraints.dirichlet import freestream_phi           # noqa: E402
from pyfp3d.kernels.cut_assembly import (                          # noqa: E402
    newton_terms23_side_coo,
    te_kutta_jacobian_coo,
    te_kutta_residual,
)
from pyfp3d.mesh.reader import read_mesh                           # noqa: E402
from pyfp3d.meshgen.wing3d import x_te                             # noqa: E402
from pyfp3d.solve.picard_ls import (                               # noqa: E402
    _farfield_split,
    _neumann_outlet_rhs,
    solve_multivalued_transonic,
)
from pyfp3d.wake import (                                          # noqa: E402
    CutElementMap,
    MultivaluedOperator,
    WakeLevelSet,
)

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MESH_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"

# B7's committed M6 conventions
M_INF = 0.70          # transonic enough for an active pocket (Terms 2/3 live)
ALPHA = 3.06
B_SEMI = 1.1963
UI, GA, C, MC, MCAP, RF = 1.0, 1.4, 1.5, 0.95, 3.0, 0.05


def build(level: str = "coarse"):
    path = MESH_DIR / f"{level}.msh"
    if not path.exists():
        raise SystemExit(
            f"{path} not present (M6 .msh are gitignored). Regenerate with:\n"
            f"  python cases/meshes/onera_m6/generate_onera_m6.py")
    mesh = read_mesh(path)
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return mesh, cm, mvop


def assemble_R_J(mvop, mesh, phi, m_inf):
    """Verbatim from tests/test_b6_newton.py::_assemble_R_J (same objects,
    same row map) so this measures the SHIPPED Newton system, not a variant."""
    cm = mvop.cm
    te_aux = cm.ext_dof_of_node[cm.te_nodes].astype(np.int64)
    ff_nodes = _farfield_split(mesh, ALPHA, UI)[3]
    ff_vals = freestream_phi(mesh.nodes[ff_nodes], ALPHA, UI)
    b = _neumann_outlet_rhs(mesh, ALPHA, UI, mvop.n_total)
    n_total = mvop.n_total
    keep_row = np.ones(n_total)
    keep_row[te_aux] = 0.0
    D = sp.diags(keep_row)
    nonte = np.asarray(mvop._nonte_rows, dtype=np.int64)
    row_map = np.arange(n_total, dtype=np.int64)
    row_map[nonte] = -1
    row_map[te_aux] = cm.te_nodes

    def remap(r, c, d):
        if not len(r):
            return r, c, d
        rr = row_map[r]
        k = rr >= 0
        return rr[k], c[k], d[k]

    phi = phi.copy()
    phi[ff_nodes] = ff_vals
    up, lo = mvop.newton_side_data(phi, m_inf, C, MC, GA, UI, MCAP, RF)
    A = mvop.assemble_matrix(rho_tilde=(up["rho_tilde"], lo["rho_tilde"]),
                             closure="wake_ls", te_kutta="pressure",
                             phi_ext=phi)
    R = A @ phi - b
    R[te_aux] = te_kutta_residual(mvop, phi)
    J = (D @ A).tocoo()
    pr, pc, pd = [J.row], [J.col], [J.data]
    for rr, cc, dd in (te_kutta_jacobian_coo(mvop, phi),
                       remap(*newton_terms23_side_coo(mvop.op, up, UI)),
                       remap(*newton_terms23_side_coo(mvop.op, lo, UI))):
        if len(rr):
            pr.append(rr)
            pc.append(cc)
            pd.append(dd)
    J = sp.coo_matrix((np.concatenate(pd),
                       (np.concatenate(pr), np.concatenate(pc))),
                      shape=(n_total, n_total)).tocsr()
    return R, J, ff_nodes


def mixed_side_plain_elements(mesh, cm):
    """Uncut elements whose nodes straddle the level set.

    "Special" (cut / TE-lower-fan) elements are assembled two-sided by
    construction; everything else is PLAIN and gets one side's density on main
    DOFs. A plain element with mixed node_side is exactly the class C1 names.
    """
    el = np.asarray(mesh.elements, dtype=np.int64)
    is_special = np.zeros(len(el), dtype=bool)
    is_special[cm.cut_elems] = True
    is_special[cm.te_lower_elems] = True
    side_e = cm.node_side[el]
    mixed = side_e.min(axis=1) != side_e.max(axis=1)
    return np.flatnonzero(mixed & ~is_special)


def rel_err(mvop, mesh, phi0, m_inf, J, v, free, eps=1e-7):
    Jv = (J @ v)[free]
    Rp = assemble_R_J(mvop, mesh, phi0 + eps * v, m_inf)[0]
    Rm = assemble_R_J(mvop, mesh, phi0 - eps * v, m_inf)[0]
    fd = ((Rp - Rm) / (2 * eps))[free]
    return float(np.linalg.norm(Jv - fd) / (np.linalg.norm(fd) + 1e-30))


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    mesh, cm, mvop = build("coarse")
    el = np.asarray(mesh.elements, dtype=np.int64)

    msp = mixed_side_plain_elements(mesh, cm)
    print(f"mesh: {len(mesh.nodes)} nodes, {len(el)} tets")
    print(f"cut elements        : {len(cm.cut_elems)}")
    print(f"beyond-tip elements : {len(cm.beyond_tip_elems)}")
    print(f"MIXED-SIDE PLAIN    : {len(msp)}   <-- the C1 class")
    if len(msp) == 0:
        print("\nthe C1 element class is EMPTY on this mesh -- the finding "
              "cannot be exercised here; nothing to measure.")
        return 2

    seed = solve_multivalued_transonic(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                       farfield="freestream",
                                       n_outer_level=200)
    phi0 = seed["phi_ext"]
    R0, J, ff = assemble_R_J(mvop, mesh, phi0, M_INF)
    print(f"seed residual |R|inf = {np.max(np.abs(R0)):.3e}, "
          f"n_nu_active = {mvop.n_nu_active}")

    # --- probe direction sets -------------------------------------------
    # touched: aux DOFs of cut nodes that belong to a mixed-side plain element
    msp_nodes = np.unique(el[msp])
    aux_of = cm.ext_dof_of_node
    touched_aux = aux_of[msp_nodes[aux_of[msp_nodes] >= 0]]
    all_aux = aux_of[aux_of >= 0]
    control_aux = np.setdiff1d(all_aux, touched_aux)
    print(f"aux DOFs: {len(all_aux)} total, {len(touched_aux)} touched by a "
          f"mixed-side plain element, {len(control_aux)} control")

    free = np.flatnonzero(~np.isin(np.arange(mvop.n_total), ff))
    rng = np.random.default_rng(0)
    rows = []
    for name, dofs in (("targeted_aux", touched_aux),
                       ("control_aux", control_aux),
                       ("global_free", free)):
        if len(dofs) == 0:
            rows.append({"probe": name, "n_dofs": 0, "rel_err": float("nan")})
            print(f"  {name:14s}: (empty)")
            continue
        v = np.zeros(mvop.n_total)
        v[dofs] = rng.standard_normal(len(dofs))
        r = rel_err(mvop, mesh, phi0, M_INF, J, v, free)
        rows.append({"probe": name, "n_dofs": int(len(dofs)), "rel_err": r})
        print(f"  {name:14s}: n={len(dofs):6d}  rel err = {r:.3e}")

    # --- eps sweep: is the targeted mismatch a MISSING TERM or just FD
    # noise / non-smoothness (a perturbation flipping upwind branches)?
    # A missing Jacobian term gives an eps-INDEPENDENT relative error; FD
    # non-smoothness does not. This is the discriminator that makes the
    # verdict adversarial rather than merely suggestive.
    v_t = np.zeros(mvop.n_total)
    v_t[touched_aux] = rng.standard_normal(len(touched_aux))
    sweep = []
    for eps in (1e-6, 1e-7, 1e-8):
        r = rel_err(mvop, mesh, phi0, M_INF, J, v_t, free, eps=eps)
        sweep.append((eps, r))
        print(f"  eps={eps:.0e}: targeted rel err = {r:.3e}")
    spread = max(s[1] for s in sweep) / (min(s[1] for s in sweep) + 1e-30)
    print(f"  eps-independence: max/min = {spread:.2f} "
          f"({'stable -> missing term' if spread < 3 else 'eps-sensitive -> suspect FD noise'})")

    tgt = next(r["rel_err"] for r in rows if r["probe"] == "targeted_aux")
    ctl = next(r["rel_err"] for r in rows if r["probe"] == "control_aux")
    if tgt > 1e-6 and ctl < 1e-6:
        verdict = "C1 CONFIRMED (localized to mixed-side plain elements)"
    elif tgt < 1e-6 and ctl < 1e-6:
        verdict = "C1 REFUTED on this mesh (Jacobian exact everywhere probed)"
    else:
        verdict = "INCONCLUSIVE (control also fails -- not the C1 mechanism)"
    print(f"\nVERDICT: {verdict}")

    out = RESULTS / "c1_fd_probes.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mesh", "m_inf", "n_nodes", "n_tets", "n_cut_elems",
                    "n_beyond_tip", "n_mixed_side_plain", "probe", "n_dofs",
                    "rel_err"])
        for r in rows:
            w.writerow(["onera_m6/coarse", M_INF, len(mesh.nodes), len(el),
                        len(cm.cut_elems), len(cm.beyond_tip_elems), len(msp),
                        r["probe"], r["n_dofs"], f"{r['rel_err']:.6e}"])
        for eps, r in sweep:
            w.writerow(["onera_m6/coarse", M_INF, len(mesh.nodes), len(el),
                        len(cm.cut_elems), len(cm.beyond_tip_elems), len(msp),
                        f"targeted_aux_eps{eps:.0e}", len(touched_aux),
                        f"{r:.6e}"])
    (RESULTS / "verdict.txt").write_text(verdict + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
