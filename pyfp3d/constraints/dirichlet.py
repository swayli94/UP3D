"""
Far-field Dirichlet data: uniform flow at incidence plus the vortex
far-field correction (design.md Sec 5, incompressible form for P2).

The vortex correction lets the far-field boundary sit at ~15 chords
instead of ~50: a lifting body looks from far away like a point vortex of
the total circulation, so prescribing plain uniform flow there truncates
an O(Gamma/r) velocity and biases cl low. For the quasi-2D M0 cases the
correction is the 2D point vortex, uniform in z:

    phi_v = -(Gamma_total / 2 pi) theta_w,   theta_w in [0, 2 pi)

with theta_w measured around the vortex center from the +x (wake)
direction, so the BRANCH CUT lies on the wake sheet: approaching the cut
from above (y -> 0+, x > x_v) gives theta_w -> 0, from below theta_w ->
2 pi, hence a jump phi(+) - phi(-) = +Gamma_total across the wake --
consistent with the master-slave constraint phi(+) = phi(-) + Gamma, so
duplicated far-field wake nodes need no special-casing beyond evaluating
the MASTER (lower side) with theta_w = 2 pi.

The induced swirl is V_theta = -Gamma/(2 pi r): clockwise for Gamma > 0,
i.e. faster flow above, the correct sense for positive lift with
Gamma = phi_upper(TE) - phi_lower(TE) > 0 and cl = 2 Gamma / (U_inf c).
"""

import numpy as np


def freestream_phi(
    points: np.ndarray, alpha_deg: float = 0.0, u_inf: float = 1.0
) -> np.ndarray:
    """phi_inf = U_inf (x cos(alpha) + y sin(alpha)); quasi-2D, beta = 0."""
    a = np.deg2rad(alpha_deg)
    return u_inf * (points[:, 0] * np.cos(a) + points[:, 1] * np.sin(a))


def vortex_phi_2d(
    points: np.ndarray,
    gamma_total: float,
    center=(0.25, 0.0),
    lower_branch_mask: np.ndarray | None = None,
) -> np.ndarray:
    """2D point-vortex potential with the branch cut along +x from `center`.

    Args:
        points: (n, 3) or (n, 2) evaluation points
        gamma_total: total circulation (span-averaged Gamma for quasi-2D)
        center: vortex location (x_v, y_v), default quarter chord
        lower_branch_mask: boolean mask of points that sit EXACTLY ON the
            cut (y == y_v, x > x_v -- e.g. wake master nodes) and belong to
            the lower (-) side: they get theta_w = 2 pi instead of 0.
            Points off the cut never need it (theta_w is continuous there).

    Returns:
        (n,) phi_v = -(Gamma/2 pi) theta_w
    """
    dx = points[:, 0] - center[0]
    dy = points[:, 1] - center[1]
    theta = np.arctan2(dy, dx)            # (-pi, pi], cut on -x half-axis
    theta_w = np.where(theta < 0.0, theta + 2.0 * np.pi, theta)  # [0, 2 pi)
    if lower_branch_mask is not None:
        on_cut = (dy == 0.0) & (dx > 0.0) & lower_branch_mask
        theta_w = np.where(on_cut, 2.0 * np.pi, theta_w)
    return -(gamma_total / (2.0 * np.pi)) * theta_w


def farfield_dirichlet(
    mesh_cut,
    wc,
    alpha_deg: float,
    gamma_stations: np.ndarray,
    u_inf: float = 1.0,
    vortex_center=(0.25, 0.0),
    farfield_tag: str = "farfield",
):
    """Dirichlet (nodes, values) on the cut mesh's far-field boundary.

    Slave ("+"-side) far-field nodes are included with their upper-side
    branch value for completeness, but the wake constraint eliminates them
    (WakeConstraint.to_reduced_dirichlet drops them); master wake nodes on
    the far field get the lower-side branch (theta_w = 2 pi).

    Args:
        mesh_cut: cut Mesh (boundary_faces[farfield_tag] uses cut-mesh ids)
        wc: WakeCut
        alpha_deg, u_inf: freestream incidence and speed
        gamma_stations: (n_st,) circulation per spanwise station; the 2D
            vortex uses their mean (span-uniform for quasi-2D)
        vortex_center: (x_v, y_v) of the equivalent point vortex

    Returns:
        (dirichlet_nodes, dirichlet_values)
    """
    ff_nodes = np.unique(
        np.asarray(mesh_cut.boundary_faces[farfield_tag], dtype=np.int64)
    )
    pts = mesh_cut.nodes[ff_nodes]
    gamma_total = float(np.mean(np.atleast_1d(gamma_stations)))

    is_master_wake = np.isin(ff_nodes, wc.master_nodes)
    is_slave_wake = np.isin(ff_nodes, wc.slave_nodes)

    values = freestream_phi(pts, alpha_deg, u_inf)
    if gamma_total != 0.0:
        phi_v = vortex_phi_2d(
            pts, gamma_total, vortex_center, lower_branch_mask=is_master_wake
        )
        # Slaves sit at the same coordinates as their masters (on the cut,
        # theta_w = 0 side): arctan2(0, +) = 0 already gives the upper
        # branch, so no extra masking is needed for them.
        del is_slave_wake
        values = values + phi_v
    return ff_nodes, values
