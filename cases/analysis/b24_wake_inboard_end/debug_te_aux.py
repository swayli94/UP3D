"""B24 debug: locate the TE node(s) that lose their aux DOF at alpha > 0
(the CutElementMap A3 assertion crash in run_e1, B1 medium a=2.0).

Replicates CutElementMap.__init__'s classification verbatim but WITHOUT the
assert, then reports every TE node that ends up outside all cut elements,
with the per-tet reason (no sign change / d_cross / q_cross / fan).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

from wb24 import FUS, LS_MESH_DIR, load_mesh, x_far_of  # noqa: E402
from pyfp3d.meshgen.wingbody import te_polyline, junction_z  # noqa: E402
from pyfp3d.wake.levelset import WakeLevelSet  # noqa: E402
from pyfp3d.wake.cut_elements import _TET_EDGES, _node_min_edge_length  # noqa: E402

level = sys.argv[1] if len(sys.argv) > 1 else "coarse"
alphas = [float(a) for a in sys.argv[2:]] or [2.0]

import os
DELTA = float(os.environ.get("DELTA", "0.0"))

mesh = load_mesh(LS_MESH_DIR / f"{level}.msh")
nodes, el = mesh.nodes, mesh.elements
wall = np.unique(mesh.boundary_faces["wall"])
h_node = _node_min_edge_length(nodes, el)
z0 = junction_z(FUS)

for alpha in alphas:
    a = np.radians(alpha)
    te = te_polyline(FUS, extend="waterline", delta=DELTA,
                     x_far=x_far_of(mesh))
    wls = WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))
    s_raw, d, q = wls.evaluate(nodes)

    te_tol = 1e-3 * h_node
    te_mask = (np.abs(s_raw) < te_tol) & (np.abs(d) < te_tol)
    on_wall = np.zeros(len(nodes), bool)
    on_wall[wall] = True
    te_mask &= on_wall
    te_nodes = np.flatnonzero(te_mask)

    eps = 1e-6 * h_node
    shift = np.abs(s_raw) < eps
    shift |= te_mask
    s = s_raw.copy()
    s[shift] = eps[shift]
    side = np.where(s > 0.0, 1, -1).astype(np.int8)

    side_e = side[el]
    cand = np.flatnonzero(side_e.min(axis=1) != side_e.max(axis=1))
    s_e, d_e, q_e = s[el[cand]], d[el[cand]], q[el[cand]]
    si, sj = s_e[:, _TET_EDGES[:, 0]], s_e[:, _TET_EDGES[:, 1]]
    di, dj = d_e[:, _TET_EDGES[:, 0]], d_e[:, _TET_EDGES[:, 1]]
    qi, qj = q_e[:, _TET_EDGES[:, 0]], q_e[:, _TET_EDGES[:, 1]]
    crossing = (si * sj) < 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(crossing, si / (si - sj), 0.0)
    d_cross = di + t * (dj - di)
    q_cross = qi + t * (qj - qi)
    on_sheet = (crossing & (d_cross > 0.0) & (q_cross >= 0.0)
                & (q_cross <= wls.span_length))
    cut = cand[on_sheet.any(axis=1)]

    is_te = np.zeros(len(nodes), bool)
    is_te[te_nodes] = True
    te_e = is_te[el]
    fan = te_e.any(axis=1) & np.all(te_e | (side[el] == -1), axis=1)
    cut = cut[~fan[cut]]

    cut_nodes = np.unique(el[cut])
    missing = te_nodes[~np.isin(te_nodes, cut_nodes)]
    import os
    if os.environ.get("DUMP_CORNER"):
        corner = te_nodes[nodes[te_nodes, 2] < 0.17]
        missing = np.unique(np.concatenate([missing, corner]))
    print(f"== alpha={alpha}: n_te={len(te_nodes)} n_cand={len(cand)} "
          f"n_cut={len(cut)} missing={len(missing)}")
    for m in missing:
        x = nodes[m]
        print(f"  TE node {m} xyz=({x[0]:.4f},{x[1]:.4f},{x[2]:.4f}) "
              f"s_raw={s_raw[m]:+.3e} d={d[m]:+.3e} q={q[m]:.3f} "
              f"(junction z0={z0:.4f})")
        tets = np.flatnonzero((el == m).any(axis=1))
        print(f"    incident={len(tets)} cand={np.isin(tets, cand).sum()} "
              f"fan={fan[tets].sum()} cut={np.isin(tets, cut).sum()}")
        # dump every node of the incident tets: where do "+" nodes sit?
        nb = np.unique(el[tets])
        k = np.argsort(s_raw[nb])
        print("    nbr nodes (xyz, s_raw, d, q, side) sorted by s_raw:")
        for j in nb[k]:
            x = nodes[j]
            print(f"      {j:7d} ({x[0]:+.4f},{x[1]:+.4f},{x[2]:+.4f}) "
                  f"s={s_raw[j]:+.3e} d={d[j]:+.4f} q={q[j]:.3f} "
                  f"side={side[j]:+d}{' TE' if is_te[j] else ''}")
        for e in tets:
            if side[el[e]].min() == side[el[e]].max():
                continue
            j = int(np.flatnonzero(cand == e)[0])
            nd = el[e]
            print(f"    tet {e} sides={side[nd].tolist()} "
                  f"s_raw={['%+.2e' % v for v in s_raw[nd]]}")
            print(f"      crossing={crossing[j].astype(int).tolist()}")
            print(f"      d_cross={np.round(d_cross[j], 4).tolist()}")
            print(f"      q_cross={np.round(q_cross[j], 4).tolist()} "
                  f"span_length={wls.span_length:.3f}")
            print(f"      on_sheet={on_sheet[j].astype(int).tolist()} "
                  f"fan={bool(fan[e])}")
