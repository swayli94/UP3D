"""LE-10: is the controlling variable the supersonic region's CONNECTIVITY?

LE-9 refuted both prior explanations and left a hint from four numbers:

  leg                     entropy  upwind_c  conv   nu_max   n_active   M2max
  ref_entropyON_c1.5      ON       1.5       YES    1.2893   106348     6.426
  ref_entropyOFF_c1.5     OFF      1.5       no     1.2586    48690     5.607
  entropyOFF_c2.0         OFF      2.0       no     1.5570    49072     4.075
  entropyOFF_c3.0         OFF      3.0       no     1.9468     2498     2.571
  taperON_entropyOFF      OFF      1.5       YES    1.2053   105409     4.594

Both refuted:
  * the user's dissipation hypothesis -- upwind_c 1.5 -> 3.0 raised nu_max 55 % (1.259 ->
    1.947), so dissipation genuinely increased, and it still did not converge. Worse, it
    ERASED the shock (n_active 48690 -> 2498, cl_p -11 %).
  * my own "the entropy correction weakens the shock" -- with entropy ON, M2max is 6.426
    against OFF's 5.607, i.e. the peak is STRONGER, not weaker.

What survives is that the two CONVERGING legs -- entropy-ON and taper+entropy-OFF, two
mechanically unrelated fixes -- land within 1 % of each other on n_active (106348 vs 105409)
while all three failures sit at 49k or below. That is suggestive, but n_active is only a COUNT.
The hypothesis it suggests is about TOPOLOGY: the upwind donor walk needs a large, connected
supersonic region, and fails when the region is fragmented.

★ This round tests that directly rather than inferring it from a count: label the supersonic
elements and find the CONNECTED COMPONENTS through the same face adjacency the donor walk uses.

  converging legs have one dominant component, failures fragmented  => topology is the variable
  all legs one dominant component                                   => the count, not topology,
      is what differs, and the hypothesis needs replacing rather than refining

Reported per leg: component count, largest-component size, the fraction of supersonic elements
inside the largest component, and the number of singleton components (a supersonic element with
no supersonic face neighbour -- the worst case for a donor walk).

★ phi is CACHED to .npz this time, so any follow-up analysis is free instead of costing another
re-solve. The previous rounds each paid ~1.2 h to regenerate states that had already been
computed.

Outputs (TRACKED):     bench/gate_results/le10_topology.csv
Outputs (gitignored):  bench/gate_results/le10_cache/<leg>.npz
"""

import csv
import os
import sys
import time
from collections import deque

os.environ.setdefault("NUMBA_NUM_THREADS", "16")
os.environ.setdefault("OMP_NUM_THREADS", "16")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")

import numpy as np                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import run_capability_matrix as cap                                 # noqa: E402
from pyfp3d.constraints.wake import tip_taper_factors               # noqa: E402
from pyfp3d.kernels.upwind import build_face_adjacency               # noqa: E402
from pyfp3d.mesh.metrics import precompute_element_geometry          # noqa: E402
from pyfp3d.mesh.reader import read_mesh                            # noqa: E402
from pyfp3d.mesh.wake_cut import cut_wake                           # noqa: E402
from pyfp3d.meshgen.wing3d import B_SEMI                            # noqa: E402
from pyfp3d.physics.isentropic import mach_number_squared            # noqa: E402
from pyfp3d.post.surface import planform_area                       # noqa: E402
from pyfp3d.post.unified import wall_forces                         # noqa: E402
from pyfp3d.solve.newton import (solve_newton_lifting,              # noqa: E402
                                 solve_newton_transonic)

OUT = os.path.join(HERE, "gate_results")
CACHE = os.path.join(OUT, "le10_cache")
os.makedirs(CACHE, exist_ok=True)
CSV = os.path.join(OUT, "le10_topology.csv")
MP = os.path.join(REPO, "cases", "meshes", "onera_m6_wingbody_conforming",
                  "medium.msh")
M_TARGET, ALPHA = 0.84, 3.06

LEGS = [("ref_entropyON_c1.5", None, True, 1.5),
        ("ref_entropyOFF_c1.5", None, False, 1.5),
        ("entropyOFF_c2.0", None, False, 2.0),
        ("entropyOFF_c3.0", None, False, 3.0),
        ("taperON_entropyOFF", 0.05, False, 1.5)]
KEYS = ["leg", "r_c", "entropy", "upwind_c", "converged", "n_supersonic",
        "n_components", "largest_component", "frac_in_largest", "n_singletons",
        "m2_max", "cl_p", "wall_s", "note"]


def components(adj, mask):
    """Connected components of the masked element set through FACE adjacency -- the same
    adjacency the upwind donor walk traverses, so the components are the walk's own islands."""
    n = len(mask)
    seen = np.zeros(n, dtype=bool)
    sizes = []
    idx = np.flatnonzero(mask)
    for s in idx:
        if seen[s]:
            continue
        q = deque([s]); seen[s] = True; sz = 0
        while q:
            e = q.popleft(); sz += 1
            for nb in adj[e]:
                if nb >= 0 and mask[nb] and not seen[nb]:
                    seen[nb] = True; q.append(nb)
        sizes.append(sz)
    return np.array(sorted(sizes, reverse=True)) if sizes else np.array([0])


def main():
    print(f"M{M_TARGET} / alpha {ALPHA} / wing-body medium -- supersonic topology")
    print("hypothesis: the donor walk needs ONE large connected supersonic "
          "region\n")
    rows = []
    for label, rc, ent, c in LEGS:
        npz = os.path.join(CACHE, f"{label}.npz")
        mc, wc = cut_wake(read_mesh(MP))
        t0 = time.perf_counter()
        if os.path.exists(npz):
            d = np.load(npz)
            phi, conv, m2max, cl = (d["phi"], bool(d["conv"]),
                                    float(d["m2max"]), float(d["cl"]))
            wall = 0.0
            print(f"  {label:22s} (cached)", flush=True)
        else:
            t = (None if rc is None
                 else tip_taper_factors(wc.station_z, B_SEMI, "vanish_smooth",
                                        rc * B_SEMI))
            try:
                seed = solve_newton_lifting(
                    mc, wc, m_inf=cap.WB_MSTART, alpha_deg=ALPHA,
                    **dict(cap.CONF_SEED_KW, entropy_correction=ent,
                           upwind_c=c))
                nk = dict(cap.CONF_RAMP_NK, kutta_estimator="pressure",
                          phi_init=seed["phi"], gamma_init=seed["gamma"],
                          n_picard_seed=0, entropy_correction=ent, upwind_c=c)
                if t is not None:
                    nk["tip_taper"] = t
                r = solve_newton_transonic(mc, wc, m_inf=M_TARGET,
                                           alpha_deg=ALPHA,
                                           m_start=cap.WB_MSTART, dm=cap.DM,
                                           dm_min=0.01, freeze_tol=1e-5,
                                           intermediate_tol=1e-4, newton_kw=nk)
            except Exception as exc:                               # noqa: BLE001
                print(f"  {label:22s} DIED {type(exc).__name__}", flush=True)
                rows.append(dict(leg=label, r_c=rc, entropy=ent, upwind_c=c,
                                 converged=False,
                                 note=f"{type(exc).__name__}: {exc}"))
                continue
            phi = np.asarray(r["phi"])
            m_att = float(r.get("m_last_converged", r.get("m_final", M_TARGET)))
            conv = bool(r["converged"]) and abs(m_att - M_TARGET) < 1e-9
            sref = planform_area(mc.nodes, mc.boundary_faces["wall"])
            cl = wall_forces(mc, phi=phi, alpha_deg=ALPHA, s_ref=sref,
                             m_inf=M_TARGET)["cl"]
            m2max = float(r.get("mach2_max", float("nan")))
            wall = time.perf_counter() - t0
            np.savez_compressed(npz, phi=phi, conv=conv, m2max=m2max, cl=cl)

        #: the supersonic set on the CUT mesh, from the same element gradients the
        #: solver's own switch uses
        B, _V = precompute_element_geometry(mc.nodes, mc.elements)
        g = np.einsum("eaj,ea->ej", B, phi[mc.elements])
        m2 = mach_number_squared(np.einsum("ej,ej->e", g, g), M_TARGET)
        mask = m2 > 1.0
        adj, _ = build_face_adjacency(np.ascontiguousarray(mc.elements))
        sizes = components(adj, mask)
        n_ss = int(mask.sum())
        largest = int(sizes[0]) if n_ss else 0
        singles = int(np.sum(sizes == 1))
        frac = (largest / n_ss) if n_ss else float("nan")
        print(f"  {label:22s} conv={conv!s:5s} n_ss={n_ss:7d} "
              f"comps={len(sizes):5d} largest={largest:7d} "
              f"frac={frac:.4f} singletons={singles:5d} "
              f"M2max={m2max:.3f} cl {cl:.6f}"
              f"{'' if wall == 0 else f' ({wall:.0f}s)'}", flush=True)
        rows.append(dict(leg=label, r_c=rc, entropy=ent, upwind_c=c,
                         converged=conv, n_supersonic=n_ss,
                         n_components=len(sizes), largest_component=largest,
                         frac_in_largest=round(frac, 6) if frac == frac else None,
                         n_singletons=singles, m2_max=m2max, cl_p=cl,
                         wall_s=round(wall, 1), note=""))
    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEYS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {CSV}")
    print("\n=== reading ===")
    ok = [r for r in rows if r.get("converged")]
    bad = [r for r in rows if r.get("n_supersonic") and not r.get("converged")]
    if ok and bad:
        fo = [r["frac_in_largest"] for r in ok if r.get("frac_in_largest")]
        fb = [r["frac_in_largest"] for r in bad if r.get("frac_in_largest")]
        print(f"  frac in largest component: converging {fo}  failing {fb}")
        if fo and fb and min(fo) > max(fb) + 0.05:
            print("  => TOPOLOGY separates them: the converging legs keep the")
            print("  supersonic region connected, the failures fragment it.")
        else:
            print("  => topology does NOT separate them. The hypothesis needs")
            print("  replacing, not refining -- n_active's agreement between the")
            print("  two converging legs is then either coincidence or a proxy")
            print("  for something else still unidentified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
