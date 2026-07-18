"""
B19 Leg B (GB19.6) — is the residual we have the residual we want?

Leg A made the Jacobian exact FOR THE CURRENT R. It did not ask whether that R
is right. For a **mixed-side plain** element -- uncut, nodes straddling the
level set -- the current residual takes

    stiffness  from the MAIN field  (mass_conservation_coo scatters on main DOFs)
    density    from the SIDE field  (newton_side_data builds rho from phi_side)

Those are two different velocity fields for one element, and nothing on record
says that asymmetry was a deliberate modelling decision. The alternative is to
read the density from the main field too, which is what `element_mach2` already
does for these elements (B8's ×5 metric artifact was the same class biting in a
metric instead of in the Jacobian).

**This probe changes nothing.** It measures how far apart the two candidate
densities actually are, on a converged state, so the decision to change R (or
not) is made against a number rather than against an aesthetic argument.
Changing R WOULD move converged answers, so it is a discretization change and
adoption is the user's call, not this script's.

Run: python cases/analysis/c1_ls_jacobian_fd/run_legb_probe.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import read_mesh                        # noqa: E402
from pyfp3d.meshgen.wing3d import x_te                          # noqa: E402
from pyfp3d.physics.isentropic import density_field, limit_q2_field  # noqa: E402
from pyfp3d.solve.picard_ls import solve_multivalued_transonic  # noqa: E402
from pyfp3d.wake import (                                       # noqa: E402
    CutElementMap, MultivaluedOperator, WakeLevelSet,
)

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MESH_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6"
ALPHA, B_SEMI, UI, GA, MCAP = 3.06, 1.1963, 1.0, 1.4, 3.0
M_INF = 0.70


def build(level="coarse"):
    path = MESH_DIR / f"{level}.msh"
    if not path.exists():
        raise SystemExit(f"{path} missing; regenerate with "
                         "cases/meshes/onera_m6/generate_onera_m6.py")
    mesh = read_mesh(path)
    a = np.radians(ALPHA)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, cm, MultivaluedOperator(mesh.nodes, mesh.elements, cm,
                                         levelset=wls)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    mesh, cm, mvop = build()
    el = np.asarray(mesh.elements, dtype=np.int64)

    special = np.zeros(len(el), dtype=bool)
    special[cm.cut_elems] = True
    special[cm.te_lower_elems] = True
    side_e = cm.node_side[el]
    msp = np.flatnonzero((side_e.min(1) != side_e.max(1)) & ~special)
    aux = cm.ext_dof_of_node
    reads_aux = np.array([bool((aux[el[e]] >= 0).any()) for e in msp])
    affected = msp[reads_aux]
    V = mvop.op.V
    print(f"mixed-side plain     : {len(msp)}")
    print(f"  ...that read an aux: {len(affected)}  "
          f"({V[affected].sum() / V.sum() * 100:.4f}% of domain volume)")
    if len(affected) == 0:
        print("nothing to measure")
        return 2

    seed = solve_multivalued_transonic(mvop, mesh, M_INF, alpha_deg=ALPHA,
                                       farfield="freestream", n_outer_level=400)
    phi = seed["phi_ext"]

    # the two candidate fields for these elements
    phi_up, phi_lo = mvop.side_potentials(phi)
    phi_main = phi[: cm.n_main]

    def rho_of(field):
        grad, q2 = mvop.op.velocities(field)
        q2n = q2 / UI ** 2
        return density_field(limit_q2_field(q2n, M_INF, MCAP, GA), M_INF, GA), q2n

    rho_up, q2_up = rho_of(phi_up)
    rho_lo, q2_lo = rho_of(phi_lo)
    rho_mn, q2_mn = rho_of(phi_main)

    # which side does mass_conservation_coo actually use for each?
    plus = cm.node_side[el[affected]].max(axis=1) == 1
    rho_side = np.where(plus, rho_up[affected], rho_lo[affected])
    q2_side = np.where(plus, q2_up[affected], q2_lo[affected])
    rho_main = rho_mn[affected]
    q2_main = q2_mn[affected]

    d_rho = np.abs(rho_side - rho_main)
    rel = d_rho / np.maximum(rho_main, 1e-30)
    d_q2 = np.abs(q2_side - q2_main)
    w = V[affected] / V[affected].sum()

    rows = [
        ("n_affected_elems", len(affected)),
        ("volume_fraction_pct", V[affected].sum() / V.sum() * 100),
        ("d_rho_max", d_rho.max()),
        ("d_rho_mean", d_rho.mean()),
        ("d_rho_volweighted", float(w @ d_rho)),
        ("rel_rho_max_pct", rel.max() * 100),
        ("rel_rho_mean_pct", rel.mean() * 100),
        ("d_q2_max", d_q2.max()),
        ("q2_side_max", q2_side.max()),
        ("q2_main_max", q2_main.max()),
        ("n_elems_rel_gt_1pct", int(np.count_nonzero(rel > 0.01))),
        ("n_elems_rel_gt_10pct", int(np.count_nonzero(rel > 0.10))),
        ("seed_residual", float(np.max(np.abs(seed.get("residual_norm", np.nan))))
         if seed.get("residual_norm") is not None else float("nan")),
    ]
    for k, v in rows:
        print(f"  {k:24s} {v}")

    out = RESULTS / "legb_density_gap.csv"
    with out.open("w", newline="") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["metric", "value"])
        for k, v in rows:
            w_.writerow([k, v])
    print(f"wrote {out}")

    print("\nREADING: this is the size of the modelling asymmetry Leg B asks "
          "about. A gap that is tiny AND confined to a 0.2%-of-volume strip is "
          "an argument for leaving R alone (and recording why); a large one is "
          "an argument for re-deriving it. The decision is the user's either "
          "way -- changing R moves converged answers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
