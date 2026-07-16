"""
Probe-free conforming Kutta target: wall-adjacent control-volume
recovered-velocity / pressure-equality estimator (P14, routed from A2).

The classic conforming target (constraints/wake.py::kutta_targets) reads the
potential jump at ONE wall probe node per side, one edge off the TE. A2
(cases/analysis/a2_te_kutta_fidelity, closed 2026-07-17) proved that probe
estimator (a) MANUFACTURES the spanwise Gamma(z) jitter from a smooth field
(fixed-Gamma discriminator D = 7.33/25.70 coarse/medium) and (b) enforces
equal *potential* jump where the physical Kutta condition is equal
*pressure*, leaving a 34x/133x larger TE Cp gap than the level-set path's
pressure-equality Kutta.

This module ports the level-set B4 objects (wake/multivalued.py::
_build_te_control_volumes + te_velocities; kernels/cut_assembly.py::
te_kutta_coo) to the conforming CUT mesh:

  * per TE node, an UPPER and a LOWER control volume = the WALL-ADJACENT
    cut-mesh tets referencing that TE node through its slave ("+" = upper)
    resp. master ("-" = lower) copy. Wall-adjacent = the tet OWNS a wall
    boundary face (exact face ownership, not the >= 3-wall-node proxy --
    topologically airtight and measured identical on the committed
    meshes). B4 measured the wall-adjacency restriction is the accuracy:
    full-fan recovery is +11-15% wrong, wall-adjacent < 1%.
  * per side, the volume-weighted P1 velocity over that control volume.
    q_u and q_l come from DIFFERENT element sets, so q_u - q_l is NOT the
    gradient of a jump field and is non-degenerate in Gamma (the same
    reason the LS pressure row can pin Gamma where the wake LS cannot).
  * the pressure-equality residual |q_u|^2 - |q_l|^2 = 0 per station
    (station mean over the station's TE nodes, mirroring kutta_targets'
    mean semantics), which factorizes EXACTLY as

        (q_u + q_l) . (q_u - q_l) = 0.

Two Gamma-derivative flavours, mirroring the LS te_kutta_coo (frozen-mean,
Picard) vs te_kutta_jacobian_coo (exact, Newton) distinction -- they differ
by (q_u - q_l).d(q_u - q_l)/dGamma, which does NOT vanish at the solution
(|q_u| = |q_l| but q_u != q_l):

  * `implied_targets` -- the per-station Gamma* zeroing the FROZEN-MEAN
    linearized residual at fixed phi_red. Returned in Gamma units so the
    existing secant/Aitken drivers consume it as a drop-in replacement for
    kutta_targets. NOTE the outer map Gamma*(Gamma) is rational, not
    affine (phi(Gamma) affine => F quadratic in Gamma): Aitken still
    converges (smooth, weak curvature -- s_bar ~ 2 u_inf dominates) but
    the one-shot affine fixed-point jump property of the probe path does
    not carry over.
  * `newton_rows` -- the EXACT state Jacobian rows (Kp_cut, D):
    Kp_cut = dF/dphi_cut (n_st x n_cut sparse, rebuilt every Newton step
    -- unlike the probe K it is state-dependent) and D = dF/dGamma
    (n_st x n_st dense, banded in practice; the probe path's counterpart
    is exactly -I). D carries only the slave-jump chain: construction
    asserts no control-volume dof is a far-field Dirichlet node, so there
    is no V_red far-field term (and none is silently dropped when the
    caller restricts Kp columns to free dofs).

Everything here is geometry + numpy on tiny TE fans (a few elements per
side per TE node); no Numba, no dependency on PicardOperator -- the
estimator is self-contained on (mesh_cut, wc), so the Laplace driver and
read-only analysis scripts (A2 discriminator rerun) can use it too.
Nothing imports this module unless kutta_estimator="pressure" is requested
(G14.4: the probe default is untouched work- and bit-wise).
"""

from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp

from pyfp3d.mesh.metrics import precompute_element_geometry
from pyfp3d.mesh.wake_cut import (
    WakeCut,
    _all_tet_faces,
    _face_key_view,
    _node_to_elem_csr,
)


class TEControlVolumes:
    """Wall-adjacent upper/lower TE control volumes on a conforming cut mesh.

    Args:
        mesh_cut, wc: output of mesh/wake_cut.cut_wake()
        dirichlet_nodes: optional cut-mesh node ids carrying far-field
            Dirichlet data. If given, construction asserts no control-volume
            dof is one of them (Jacobian-exactness guard, module docstring).
        wall_tag: boundary group of the body surface on the cut mesh.

    Attributes:
        n_te, n_st: TE node / station counts (from wc)
        fan sizes etc. are exposed via `fan_stats()` for diagnostics.
    """

    def __init__(self, mesh_cut, wc: WakeCut, dirichlet_nodes=None,
                 wall_tag: str = "wall"):
        el = np.asarray(mesh_cut.elements, dtype=np.int64)
        nodes = np.asarray(mesh_cut.nodes, dtype=np.float64)
        n_cut = len(nodes)
        self.wc = wc
        self.n_te = len(wc.te_nodes)
        self.n_st = wc.n_stations
        if self.n_te == 0:
            raise ValueError("WakeCut has no TE nodes")

        # -- wall-adjacent elements: EXACT wall-face ownership ---------------
        if wall_tag not in mesh_cut.boundary_faces:
            raise ValueError(f"cut mesh has no '{wall_tag}' face group")
        wall_tris = np.asarray(mesh_cut.boundary_faces[wall_tag],
                               dtype=np.int64)
        all_keys, all_owners = _all_tet_faces(el)
        order = np.argsort(all_keys)
        keys_sorted = all_keys[order]
        owners_sorted = all_owners[order]
        wkeys = _face_key_view(wall_tris)
        lo = np.searchsorted(keys_sorted, wkeys, side="left")
        hi = np.searchsorted(keys_sorted, wkeys, side="right")
        if not np.all(hi - lo == 1):
            bad = int(np.count_nonzero(hi - lo != 1))
            raise AssertionError(
                f"{bad} wall faces not owned by exactly one tet -- "
                "wall tagging inconsistent with the cut mesh"
            )
        wall_adj = np.zeros(len(el), dtype=bool)
        wall_adj[owners_sorted[lo]] = True

        # -- per-TE-node upper (slave copy) / lower (master copy) fans -------
        slave_of = np.full(wc.n_nodes_orig, -1, dtype=np.int64)
        slave_of[wc.master_nodes] = wc.slave_nodes
        te_slaves = slave_of[wc.te_nodes]
        if np.any(te_slaves < 0):
            raise AssertionError(
                "TE node without a slave copy -- P2 duplicated-TE topology "
                "violated (assert_wake_topology should have caught this)"
            )

        offsets, star = _node_to_elem_csr(el, n_cut)

        def _fan(node_id: int) -> np.ndarray:
            e = star[offsets[node_id]: offsets[node_id + 1]]
            return e[wall_adj[e]]

        fans_u = [_fan(int(s)) for s in te_slaves]
        fans_l = [_fan(int(t)) for t in wc.te_nodes]
        empty = [k for k in range(self.n_te)
                 if len(fans_u[k]) == 0 or len(fans_l[k]) == 0]
        if empty:
            raise AssertionError(
                f"TE nodes {empty}: empty wall-adjacent control volume on "
                "at least one side -- check wall tagging / cut topology"
            )

        # side identity: each Kutta probe shares a wall face with its TE
        # node, and that face's owning tet is side-consistent, so the upper
        # probe MUST be a vertex of the upper fan (and lower of lower).
        # This ties the CV sides to the exact convention kutta_targets and
        # Gamma = phi(+) - phi(-) already use (stronger than any geometric
        # centroid check).
        for k in range(self.n_te):
            if wc.kutta_upper[k] not in el[fans_u[k]]:
                raise AssertionError(
                    f"TE node {int(wc.te_nodes[k])}: upper probe not in the "
                    "slave-side fan -- side convention violated"
                )
            if wc.kutta_lower[k] not in el[fans_l[k]]:
                raise AssertionError(
                    f"TE node {int(wc.te_nodes[k])}: lower probe not in the "
                    "master-side fan -- side convention violated"
                )

        # -- flat CSR-style storage + fan-local geometry ---------------------
        def _pack(fans):
            off = np.zeros(self.n_te + 1, dtype=np.int64)
            off[1:] = np.cumsum([len(f) for f in fans])
            elems = (np.concatenate(fans) if off[-1] else
                     np.empty(0, dtype=np.int64))
            seg = np.repeat(np.arange(self.n_te, dtype=np.int64),
                            np.diff(off))
            dofs = el[elems]
            B, V = precompute_element_geometry(nodes, dofs)
            # volume weights normalized per fan
            vsum = np.zeros(self.n_te)
            np.add.at(vsum, seg, V)
            w = V / vsum[seg]
            return {"off": off, "elems": elems, "seg": seg, "dofs": dofs,
                    "B": B, "w": w}

        self._u = _pack(fans_u)
        self._l = _pack(fans_l)

        if dirichlet_nodes is not None and len(dirichlet_nodes) > 0:
            is_dir = np.zeros(n_cut, dtype=bool)
            is_dir[np.asarray(dirichlet_nodes, dtype=np.int64)] = True
            n_hit = int(is_dir[self._u["dofs"]].sum()
                        + is_dir[self._l["dofs"]].sum())
            if n_hit:
                raise AssertionError(
                    f"{n_hit} control-volume dofs are far-field Dirichlet "
                    "nodes -- the dF/dGamma block would be missing the "
                    "V_red far-field term (module docstring); this mesh "
                    "needs that term implemented, refusing to run silently"
                )

        # per-station TE-node counts (station-mean weights)
        self._counts = np.bincount(wc.te_station,
                                   minlength=self.n_st).astype(np.float64)
        # slave-indicator matrix G: phi_cut = T phi_red + G Gamma
        # (same object as WakeConstraint.reduce_operator's G)
        self._G = sp.coo_matrix(
            (np.ones(len(wc.slave_nodes)),
             (wc.slave_nodes, wc.node_station)),
            shape=(n_cut, self.n_st),
        ).tocsr()
        self._n_cut = n_cut

    # -- evaluation -----------------------------------------------------------

    def _side_velocity(self, side: dict, phi_cut: np.ndarray) -> np.ndarray:
        """(n_te, 3) volume-weighted recovered velocity over one side's CVs."""
        x = np.asarray(phi_cut, dtype=np.float64)
        grad = np.einsum("ead,ea->ed", side["B"], x[side["dofs"]])
        q = np.zeros((self.n_te, 3))
        np.add.at(q, side["seg"], side["w"][:, None] * grad)
        return q

    def te_velocities(self, phi_cut) -> Tuple[np.ndarray, np.ndarray]:
        """(q_upper, q_lower), each (n_te, 3) -- the B4 recovery, conforming."""
        return (self._side_velocity(self._u, phi_cut),
                self._side_velocity(self._l, phi_cut))

    def residual_stations(self, phi_cut) -> np.ndarray:
        """(n_st,) pressure-equality residual F_raw = mean_k(|q_u|^2-|q_l|^2).

        Upper = slave ("+") side, so dF/dGamma has one uniform sign per
        mesh, mirroring the probe path's dF/dGamma = -I orientation (the
        sign itself is asserted uniform, never hardcoded)."""
        qu, ql = self.te_velocities(phi_cut)
        f_node = np.einsum("kd,kd->k", qu, qu) - np.einsum("kd,kd->k", ql, ql)
        return np.bincount(self.wc.te_station, weights=f_node,
                           minlength=self.n_st) / self._counts

    def _rows_cut(self, phi_cut, mode: str) -> sp.csr_matrix:
        """(n_st, n_cut) sparse d(station residual)/d(phi_cut).

        mode="exact":  d(|q_u|^2 - |q_l|^2) = 2 q_u.dq_u - 2 q_l.dq_l
        mode="frozen": s_bar.(dq_u - dq_l) at frozen s_bar = q_u + q_l
        (the te_kutta_jacobian_coo vs te_kutta_coo distinction -- exact for
        the Newton elimination + FD guards, frozen for the Picard implied
        target where any nonsingular D preserves the fixed point)."""
        qu, ql = self.te_velocities(phi_cut)
        if mode == "exact":
            vec_u, vec_l = 2.0 * qu, 2.0 * ql
        elif mode == "frozen":
            vec_u = vec_l = qu + ql
        else:
            raise ValueError(f"unknown mode {mode!r}")

        rows, cols, data = [], [], []
        for side, vec, sign in ((self._u, vec_u, 1.0),
                                (self._l, vec_l, -1.0)):
            seg = side["seg"]
            st = self.wc.te_station[seg]
            # d q_side / d phi[dofs[i,a]] = w_i * B[i,a,:]
            contrib = np.einsum("id,iad->ia", vec[seg], side["B"])
            contrib *= (sign * side["w"] / self._counts[st])[:, None]
            rows.append(np.repeat(st, 4))
            cols.append(side["dofs"].ravel())
            data.append(contrib.ravel())
        K = sp.coo_matrix(
            (np.concatenate(data),
             (np.concatenate(rows), np.concatenate(cols))),
            shape=(self.n_st, self._n_cut),
        ).tocsr()
        K.sum_duplicates()
        return K

    def gamma_jacobian(self, phi_cut, mode: str = "exact") -> np.ndarray:
        """(n_st, n_st) dense D = dF/dGamma at this state (via Kp_cut @ G)."""
        return np.asarray((self._rows_cut(phi_cut, mode) @ self._G).todense())

    def newton_rows(self, phi_cut) -> Tuple[sp.csr_matrix, np.ndarray]:
        """(Kp_cut, D): exact dF/dphi_cut (n_st x n_cut sparse) and
        dF/dGamma (n_st x n_st dense). State-dependent -- rebuild every
        Newton step (unlike the probe path's constant K)."""
        Kp = self._rows_cut(phi_cut, "exact")
        D = np.asarray((Kp @ self._G).todense())
        return Kp, D

    def implied_targets(self, phi_cut, gamma) -> np.ndarray:
        """Per-station Gamma* zeroing the frozen-mean residual at fixed
        phi_red: Gamma* = Gamma - M^{-1} F with M = dF_frozen/dGamma.

        At the evaluation state s_bar.(q_u - q_l) == |q_u|^2 - |q_l|^2
        exactly (algebraic identity), so Gamma* = Gamma iff F = 0: the
        fixed point of the target map is exactly the nonlinear pressure
        equality. n_st is <= a few hundred, so the (banded, strictly
        diagonally dominant in practice) M is solved dense."""
        gamma = np.atleast_1d(np.asarray(gamma, dtype=np.float64))
        F = self.residual_stations(phi_cut)
        M = self.gamma_jacobian(phi_cut, mode="frozen")
        try:
            step = np.linalg.solve(M, F)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(
                "frozen-mean dF/dGamma singular -- the recovered two-sided "
                "velocity is degenerate in Gamma at this state (P14 "
                "diagnostic gate G14.1 territory)"
            ) from e
        return gamma - step

    # -- diagnostics -----------------------------------------------------------

    def fan_stats(self) -> dict:
        """Fan-size stats per side (diagnostics; a 1-element fan is legal
        but worth recording -- the CV is then a single-tet gradient)."""
        nu = np.diff(self._u["off"])
        nl = np.diff(self._l["off"])
        return {
            "n_te": self.n_te, "n_st": self.n_st,
            "fan_u_min": int(nu.min()), "fan_u_med": float(np.median(nu)),
            "fan_u_max": int(nu.max()),
            "fan_l_min": int(nl.min()), "fan_l_med": float(np.median(nl)),
            "fan_l_max": int(nl.max()),
        }
