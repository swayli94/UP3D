"""IBL3 solver core: Galerkin P1 finite-element assembly of Drela's six
integral boundary-layer equations on a SurfaceMesh, analytic Jacobian, and a
pseudo-time-globalized Newton solver.

Binding references: Drela AIAA 2013-2437 ("D13") and docs/design_track_v.md
(§2.2 equations, §5 discretization). Equations (steady form, D13 (24)(26)
(28)(29)(31); design doc §2.2):

  R1 x-momentum:   div(Jx) - u div(M) - tau_xw = 0
  R2 z-momentum:   div(Jz) - w div(M) - tau_zw = 0
  R3 kinetic en.:  div(E) - q^2 div(M) - rho Q.grad(q^2) - 2 D = 0
  R4 lateral curv.:div(Ko) + E.grad(psi_i) + (1/2) rho (Q x grad(q^2)).y
                   - rho Qo.grad(q^2) + Dx - 2 Do = 0
  R5/R6 stress:    div(K_tau_c) - S_tau_c = 0   (pinned to CTAU_LAM on laminar
                   nodes, D-TR forced transition)

Defect fluxes M, J, E, Ko, Q, Qo, wall shear and dissipation are assembled
from the closure packet (D13 (62)(63), closures.py). Discretization follows
D13 §III.D faithfully (design doc §5.1): nodal fluxes, strong-form
divergences via constant P1 gradients, 3-point edge-midpoint quadrature
(exact for the linear/quadratic source integrands), Galerkin tent weighting;
only the artificial diffusion is integrated by parts (D13 (74)), acting on
each equation's conserved density (Mx, Mz, e, k_o, k_tau1, k_tau2) with
isotropic h_bar = h I, h = sqrt(2 A) (D-HB). Residual components are taken
in the per-node LOCAL tangent basis (basis_x, basis_zc = basis_x x basis_y)
of D13 SIII.D.1 (the V1 recorded follow-up, implemented for V3): on a
planar +-y-normal surface basis_zc is exactly +z_global, so planar cases
are bit-identical to the pre-V3 global-(x, z) form.

Solution: backward-Euler pseudo-time on the physical conserved densities of
D13 (70)-(72) (Mx - u m, Mz - w m, e - q^2 m, k_o, k_tau1, k_tau2; mass-
lumped, nodal block-diagonal Jacobian, D-PT) with geometric CFL ramp
recovering plain Newton; sparse direct solve (scipy spsolve) on the symbolic
CSR pattern from SurfaceMesh.build_jacobian_pattern (6x6 blocks per node
pair). Assembly kernels are color-locked prange (bit-deterministic),
analytic Jacobian verified against FD (project rule, tests/test_v1_ibl3.py).
"""

import numpy as np
import scipy.sparse.linalg

from pyfp3d.viscous import closures as C
from pyfp3d.viscous.closures import _njit, prange
from pyfp3d.viscous.surface_mesh import QUAD_N

NUNK = 6

# ---------------------------------------------------------------------------
# Packed per-node flux/source table layout (values fv (n,NV), state
# derivatives dv (n,NV,6)); 3-vector fluxes occupy 3 consecutive columns.
# ---------------------------------------------------------------------------
F_M3 = 0    # mass-flux defect vector M            (D13 (62))
F_JX3 = 3   # J . x_hat
F_JZ3 = 6   # J . z_hat
F_E3 = 9    # energy-flux defect vector E
F_KO3 = 12  # curvature-defect flux vector K_o
F_Q3 = 15   # volume-flux defect vector Q
F_QO3 = 18  # curvature volume-flux vector Q_o
F_K13 = 21  # shear-stress flux vector K_tau_1
F_K23 = 24  # shear-stress flux vector K_tau_2
F_TAUX = 27  # wall shear stress components tau_xw, tau_zw (D13 (63))
F_TAUZ = 28
F_DD = 29   # dissipation integrals D, Dx, Do      (D13 (63))
F_DX = 30
F_DO = 31
F_S1 = 32   # stress-equation sources S_tau_1, S_tau_2 (D-CT)
F_S2 = 33
F_DMX = 34  # conserved densities for the diffusion term (D13 (70)-(72)):
F_DMZ = 35  # Mx, Mz, e, k_o, k_tau1, k_tau2
F_DE = 36
F_DKO = 37
F_DKT1 = 38
F_DKT2 = 39
F_DM = 40   # scalar mass defect m = rho * delta_rho (pseudo-time terms)
NV = 41

# divergence-receiving vector fluxes, in equation order M, Jx, Jz, E, Ko, K1, K2
_DIVERGED = (F_M3, F_JX3, F_JZ3, F_E3, F_KO3, F_K13, F_K23)


@_njit(cache=True, fastmath=True)
def _frames(u_e, basis_y, basis_x, basis_zc, q, s1, s2, ucomp, wcomp, q2,
            psi, s1l, s2l):
    """Per-node edge-velocity frames: q, s1 = q_vec/q (fallback basis_x),
    s2 = s1 x n_hat (D13 (38)); velocity components (u, w) = u_e . (basis_x,
    basis_zc) and the flow angle psi = atan2(w, u) in the per-node LOCAL
    tangent basis, where basis_zc = basis_x x basis_y (the right-handed
    local "z"; on a planar +-y-normal surface this is exactly +z_global, so
    planar cases are bit-identical to the pre-V3 global-basis form).

    Also the local components of s1/s2 themselves (s1l/s2l (n,2) =
    s . (basis_x, basis_zc)): the scalar flux contractions (Jx/Jz, tau,
    DMX/DMZ) are per-node local-frame quantities (D13 SIII.D.1, the
    curved-surface local basis projection -- the V1 recorded follow-up,
    implemented for V3). The 3-D s1/s2 vectors are unchanged and still
    carry the vector fluxes.
    """
    n = q.shape[0]
    for i in range(n):
        vx = u_e[i, 0]
        vy = u_e[i, 1]
        vz = u_e[i, 2]
        qi = np.sqrt(vx * vx + vy * vy + vz * vz)
        q[i] = qi
        if qi > 1.0e-12:
            s1x = vx / qi
            s1y = vy / qi
            s1z = vz / qi
        else:
            s1x = basis_x[i, 0]
            s1y = basis_x[i, 1]
            s1z = basis_x[i, 2]
        nx = basis_y[i, 0]
        ny = basis_y[i, 1]
        nz = basis_y[i, 2]
        # s2 = s1 x n (D13 (38); |s1 x n| = 1 up to normalization)
        c1 = s1y * nz - s1z * ny
        c2 = s1z * nx - s1x * nz
        c3 = s1x * ny - s1y * nx
        cm = np.sqrt(c1 * c1 + c2 * c2 + c3 * c3)
        cm = max(cm, 1.0e-30)
        s1[i, 0] = s1x
        s1[i, 1] = s1y
        s1[i, 2] = s1z
        s2[i, 0] = c1 / cm
        s2[i, 1] = c2 / cm
        s2[i, 2] = c3 / cm
        # local-basis velocity components and flow angle
        bxx = basis_x[i, 0]
        bxy = basis_x[i, 1]
        bxz = basis_x[i, 2]
        bzx = basis_zc[i, 0]
        bzy = basis_zc[i, 1]
        bzz = basis_zc[i, 2]
        ui = vx * bxx + vy * bxy + vz * bxz
        wi = vx * bzx + vy * bzy + vz * bzz
        ucomp[i] = ui
        wcomp[i] = wi
        q2[i] = qi * qi
        psi[i] = np.arctan2(wi, ui)
        # local components of s1/s2 (scalar flux contractions)
        s1l[i, 0] = s1x * bxx + s1y * bxy + s1z * bxz
        s1l[i, 1] = s1x * bzx + s1y * bzy + s1z * bzz
        s2l[i, 0] = s2[i, 0] * bxx + s2[i, 1] * bxy + s2[i, 2] * bxz
        s2l[i, 1] = s2[i, 0] * bzx + s2[i, 1] * bzy + s2[i, 2] * bzz


@_njit(cache=True, fastmath=True)
def _nodal_fluxes(states, q, s1, s2, s1l, s2l, rho, outs, douts, c_l, fv, dv):
    """Fill fv (n,NV) and dv (n,NV,6) from the closure packet.

    All formulas are algebraic chains of D13 (62)(63) and D-CT on top of the
    closure outputs; every derivative w.r.t. the six state unknowns follows
    the same chain (douts rows), so the FD-vs-analytic check at residual
    level exercises exactly these tables.

    The vector fluxes are built from the 3-D s1/s2; the scalar contractions
    (Jx/Jz, tau_x/tau_z, DMX/DMZ) use the per-node LOCAL components s1l/s2l
    = s . (basis_x, basis_zc) (D13 SIII.D.1 local basis; "x"/"z" equation
    rows are momentum along basis_x / basis_zc).
    """
    n = states.shape[0]
    for i in prange(n):
        qi = q[i]
        ri = rho[i]
        de = states[i, 0]
        ct1 = states[i, 4]
        ct2 = states[i, 5]
        o = outs[i]
        do = douts[i]
        s1x = s1l[i, 0]
        s1z = s1l[i, 1]
        s2x = s2l[i, 0]
        s2z = s2l[i, 1]
        ds1 = o[C.OUT_DS1]
        ds2 = o[C.OUT_DS2]
        p11 = o[C.OUT_P11]
        p12 = o[C.OUT_P12]
        p22 = o[C.OUT_P22]
        cf1 = o[C.OUT_CF1]
        cf2 = o[C.OUT_CF2]
        ak = o[C.OUT_AK]
        sp1 = o[C.OUT_SP1]
        sp2 = o[C.OUT_SP2]
        sd = o[C.OUT_SD]
        ku1 = o[C.OUT_KU1]
        ku2 = o[C.OUT_KU2]

        rq = ri * qi
        rq2 = rq * qi
        rq3 = rq2 * qi

        # --- 3-vector fluxes, value and derivative rows ---
        for k in range(3):
            a1 = s1[i, k]
            a2 = s2[i, k]
            # M = rho q (ds1 s1 + ds2 s2)
            fv[i, F_M3 + k] = rq * (ds1 * a1 + ds2 * a2)
            # Jx = rho q^2 [ (p11 s1x + p12 s2x) s1 + (p12 s1x + p22 s2x) s2 ]
            ca = p11 * s1x + p12 * s2x
            cb = p12 * s1x + p22 * s2x
            fv[i, F_JX3 + k] = rq2 * (ca * a1 + cb * a2)
            # Jz (z components of the basis vectors)
            cc = p11 * s1z + p12 * s2z
            cd = p12 * s1z + p22 * s2z
            fv[i, F_JZ3 + k] = rq2 * (cc * a1 + cd * a2)
            # E = rho q^3 (ps1 s1 + ps2 s2)
            fv[i, F_E3 + k] = rq3 * (o[C.OUT_PS1] * a1 + o[C.OUT_PS2] * a2)
            # Ko = rho q^3 (tc1 s1 + tc2 s2)
            fv[i, F_KO3 + k] = rq3 * (o[C.OUT_TC1] * a1 + o[C.OUT_TC2] * a2)
            # Q = q (dq1 s1 + dq2 s2)
            fv[i, F_Q3 + k] = qi * (o[C.OUT_DQ1] * a1 + o[C.OUT_DQ2] * a2)
            # Qo = q (dc1 s1 + dc2 s2)
            fv[i, F_QO3 + k] = qi * (o[C.OUT_DC1] * a1 + o[C.OUT_DC2] * a2)
            # K_tau_c = rho q^3 delta ct_c (ku1 s1 + ku2 s2)
            fv[i, F_K13 + k] = rq3 * de * ct1 * (ku1 * a1 + ku2 * a2)
            fv[i, F_K23 + k] = rq3 * de * ct2 * (ku1 * a1 + ku2 * a2)
            for l in range(6):
                dv[i, F_M3 + k, l] = rq * (
                    do[C.OUT_DS1, l] * a1 + do[C.OUT_DS2, l] * a2
                )
                dca = do[C.OUT_P11, l] * s1x + do[C.OUT_P12, l] * s2x
                dcb = do[C.OUT_P12, l] * s1x + do[C.OUT_P22, l] * s2x
                dv[i, F_JX3 + k, l] = rq2 * (dca * a1 + dcb * a2)
                dcc = do[C.OUT_P11, l] * s1z + do[C.OUT_P12, l] * s2z
                dcd = do[C.OUT_P12, l] * s1z + do[C.OUT_P22, l] * s2z
                dv[i, F_JZ3 + k, l] = rq2 * (dcc * a1 + dcd * a2)
                dv[i, F_E3 + k, l] = rq3 * (
                    do[C.OUT_PS1, l] * a1 + do[C.OUT_PS2, l] * a2
                )
                dv[i, F_KO3 + k, l] = rq3 * (
                    do[C.OUT_TC1, l] * a1 + do[C.OUT_TC2, l] * a2
                )
                dv[i, F_Q3 + k, l] = qi * (
                    do[C.OUT_DQ1, l] * a1 + do[C.OUT_DQ2, l] * a2
                )
                dv[i, F_QO3 + k, l] = qi * (
                    do[C.OUT_DC1, l] * a1 + do[C.OUT_DC2, l] * a2
                )
                dk1 = ku1 * a1 + ku2 * a2
                ddk1 = do[C.OUT_KU1, l] * a1 + do[C.OUT_KU2, l] * a2
                dfac = 0.0
                if l == 0:
                    dfac = ct1
                elif l == 4:
                    dfac = de
                dv[i, F_K13 + k, l] = rq3 * (dfac * dk1 + de * ct1 * ddk1)
                dfac2 = 0.0
                if l == 0:
                    dfac2 = ct2
                elif l == 5:
                    dfac2 = de
                dv[i, F_K23 + k, l] = rq3 * (dfac2 * dk1 + de * ct2 * ddk1)

        # --- scalars ---
        fv[i, F_TAUX] = 0.5 * rq2 * (cf1 * s1x + cf2 * s2x)
        fv[i, F_TAUZ] = 0.5 * rq2 * (cf1 * s1z + cf2 * s2z)
        fv[i, F_DD] = rq3 * o[C.OUT_CD]
        fv[i, F_DX] = rq3 * o[C.OUT_CDX]
        fv[i, F_DO] = rq3 * o[C.OUT_CDC]
        # conserved densities (diffusion variables)
        fv[i, F_DMX] = rq * (ds1 * s1x + ds2 * s2x)
        fv[i, F_DMZ] = rq * (ds1 * s1z + ds2 * s2z)
        fv[i, F_DE] = rq2 * o[C.OUT_DQ]
        fv[i, F_DKO] = rq2 * o[C.OUT_DQC]
        fv[i, F_DKT1] = rq2 * de * ct1 * ak
        fv[i, F_DKT2] = rq2 * de * ct2 * ak
        fv[i, F_DM] = ri * o[C.OUT_DRHO]
        for l in range(6):
            dv[i, F_TAUX, l] = 0.5 * rq2 * (
                do[C.OUT_CF1, l] * s1x + do[C.OUT_CF2, l] * s2x
            )
            dv[i, F_TAUZ, l] = 0.5 * rq2 * (
                do[C.OUT_CF1, l] * s1z + do[C.OUT_CF2, l] * s2z
            )
            dv[i, F_DD, l] = rq3 * do[C.OUT_CD, l]
            dv[i, F_DX, l] = rq3 * do[C.OUT_CDX, l]
            dv[i, F_DO, l] = rq3 * do[C.OUT_CDC, l]
            dv[i, F_DMX, l] = rq * (
                do[C.OUT_DS1, l] * s1x + do[C.OUT_DS2, l] * s2x
            )
            dv[i, F_DMZ, l] = rq * (
                do[C.OUT_DS1, l] * s1z + do[C.OUT_DS2, l] * s2z
            )
            dv[i, F_DE, l] = rq2 * do[C.OUT_DQ, l]
            dv[i, F_DKO, l] = rq2 * do[C.OUT_DQC, l]
            dak = do[C.OUT_AK, l]
            d1 = 0.0
            if l == 0:
                d1 = ct1 * ak
            elif l == 4:
                d1 = de * ak
            dv[i, F_DKT1, l] = rq2 * (d1 + de * ct1 * dak)
            d2 = 0.0
            if l == 0:
                d2 = ct2 * ak
            elif l == 5:
                d2 = de * ak
            dv[i, F_DKT2, l] = rq2 * (d2 + de * ct2 * dak)
            dv[i, F_DM, l] = ri * do[C.OUT_DRHO, l]

        # --- stress sources S_c = 2 a1 (P_c - D_c), D-CT; |tau'| is the
        # vector magnitude => sqrt(Ctau_mag), Ctau_mag = sqrt(ct1^2+ct2^2).
        mag2 = ct1 * ct1 + ct2 * ct2
        dm2_4 = 2.0 * ct1
        dm2_5 = 2.0 * ct2
        if mag2 < 1.0e-30:
            mag2 = 1.0e-30
            dm2_4 = 0.0
            dm2_5 = 0.0
        root = np.sqrt(mag2)          # Ctau_mag
        droot_4 = dm2_4 / (2.0 * root)
        droot_5 = dm2_5 / (2.0 * root)
        sroot = np.sqrt(root)         # sqrt(Ctau_mag)
        sroot = max(sroot, 1.0e-15)
        dsr_4 = droot_4 / (2.0 * sroot)
        dsr_5 = droot_5 / (2.0 * sroot)
        a1 = C.A1_BRADSHAW
        P1 = rq3 * ct1 * sp1
        D1 = rq3 / c_l * sroot * ct1 * sd
        fv[i, F_S1] = 2.0 * a1 * (P1 - D1)
        P2 = rq3 * ct2 * sp2
        D2 = rq3 / c_l * sroot * ct2 * sd
        fv[i, F_S2] = 2.0 * a1 * (P2 - D2)
        for l in range(6):
            dsr = 0.0
            if l == 4:
                dsr = dsr_4
            elif l == 5:
                dsr = dsr_5
            dP1 = rq3 * ((1.0 if l == 4 else 0.0) * sp1 + ct1 * do[C.OUT_SP1, l])
            dD1 = rq3 / c_l * (
                dsr * ct1 * sd
                + sroot * ((1.0 if l == 4 else 0.0) * sd + ct1 * do[C.OUT_SD, l])
            )
            dv[i, F_S1, l] = 2.0 * a1 * (dP1 - dD1)
            dP2 = rq3 * ((1.0 if l == 5 else 0.0) * sp2 + ct2 * do[C.OUT_SP2, l])
            dD2 = rq3 / c_l * (
                dsr * ct2 * sd
                + sroot * ((1.0 if l == 5 else 0.0) * sd + ct2 * do[C.OUT_SD, l])
            )
            dv[i, F_S2, l] = 2.0 * a1 * (dP2 - dD2)


@_njit(cache=True)
def _assemble(fv, dv, ucomp, wcomp, q2, psi, rho,
              triangles, gradN, area_tri,
              color_offsets, color_elems,
              veps, veps_s, s1n, basisy, R, Jdata, elem_to_csr, do_jac):
    """Element assembly of the six Galerkin residuals and (optionally) the
    analytic Jacobian into the CSR data array via elem_to_csr.

    Race-free within a color (no shared nodes), so prange over each color
    range is bit-deterministic. R and Jdata are accumulated; the caller
    zeroes them once per Newton step.

    Artificial diffusion (D-HB): isotropic part veps*h_e*(grad G . grad N)
    plus a streamwise tensor part veps_s*h_e*(s1.grad G)*(s1.grad N) with
    s1 the element-averaged edge-direction unit vector. The streamwise part
    damps the centered-Galerkin 2h streamwise mode (GV1.1(e) finding,
    design doc §9 item 4); s1 depends on edge data only, not on the state,
    so the Jacobian stays first-order exact.

    The R4 (Q x grad(q^2)) . n source term is the scalar triple product of
    the interpolated 3-D Q vector, the element q^2 gradient and the
    interpolated LOCAL normal basis_y (D13 SIII.D.1; on a y-normal planar
    surface it reduces operation-for-operation to the pre-V3 global
    formula q3z*gq2x - q3x*gq2z, so planar cases stay bit-identical).
    """
    n_colors = color_offsets.shape[0] - 1
    for c in range(n_colors):
        lo = color_offsets[c]
        hi = color_offsets[c + 1]
        for ei in prange(lo, hi):
            e = color_elems[ei]
            nn = triangles[e]
            area = area_tri[e]
            h_e = np.sqrt(2.0 * area)
            wq = area / 3.0
            gN = gradN[e]  # (3 local nodes, 3 xyz), constant per element

            # ---- strong-form divergences (scalars, constant per element)
            divs = np.zeros(7)
            ddiv = np.zeros((7, 3, 6))  # [flux, local node b, state l]
            for f_i in range(7):
                base = _DIVERGED[f_i]
                acc = 0.0
                for b in range(3):
                    nb = nn[b]
                    for k in range(3):
                        acc += fv[nb, base + k] * gN[b, k]
                        if do_jac:
                            for l in range(6):
                                ddiv[f_i, b, l] += 0.0  # filled below
                divs[f_i] = acc
            if do_jac:
                for f_i in range(7):
                    base = _DIVERGED[f_i]
                    for b in range(3):
                        nb = nn[b]
                        for l in range(6):
                            accd = 0.0
                            for k in range(3):
                                accd += dv[nb, base + k, l] * gN[b, k]
                            ddiv[f_i, b, l] = accd

            # ---- gradients of edge data (constant per element) ----
            gq2 = np.zeros(3)
            gpsi = np.zeros(3)
            for k in range(3):
                gq2[k] = (
                    q2[nn[0]] * gN[0, k]
                    + q2[nn[1]] * gN[1, k]
                    + q2[nn[2]] * gN[2, k]
                )
                gpsi[k] = (
                    psi[nn[0]] * gN[0, k]
                    + psi[nn[1]] * gN[1, k]
                    + psi[nn[2]] * gN[2, k]
                )

            # ---- interpolated 3-D Q vector and local normal at the quad
            # points (for the R4 (Q x grad(q^2)) . n triple product) ----
            qm_tab = np.zeros((3, 3))
            nm_tab = np.zeros((3, 3))
            for m in range(3):
                for b in range(3):
                    nb = nn[b]
                    nmb = QUAD_N[m, b]
                    for k in range(3):
                        qm_tab[m, k] += nmb * fv[nb, F_Q3 + k]
                        nm_tab[m, k] += nmb * basisy[nb, k]

            # ---- diffusion gradients per density column ----
            gdens = np.zeros((6, 3))
            for dcol in range(6):
                base = F_DMX + dcol
                for k in range(3):
                    gdens[dcol, k] = (
                        fv[nn[0], base] * gN[0, k]
                        + fv[nn[1], base] * gN[1, k]
                        + fv[nn[2], base] * gN[2, k]
                    )

            # ---- element stream direction (edge data average, floored) ----
            s1e = np.zeros(3)
            for k in range(3):
                s1e[k] = (s1n[nn[0], k] + s1n[nn[1], k] + s1n[nn[2], k]) / 3.0
            s1m = np.sqrt(s1e[0] * s1e[0] + s1e[1] * s1e[1] + s1e[2] * s1e[2])
            s1m = max(s1m, 1.0e-30)
            s1e[0] /= s1m
            s1e[1] /= s1m
            s1e[2] /= s1m

            # ---- interpolated edge-data velocity weights for the divM
            # source terms: umw[a] = sum_m wq N_a(m) u(m), etc.
            umw = np.zeros(3)
            wmw = np.zeros(3)
            q2mw = np.zeros(3)
            u_m = np.zeros(3)
            w_m = np.zeros(3)
            q2_m = np.zeros(3)
            rho_m = np.zeros(3)
            for m in range(3):
                for b in range(3):
                    nb = nn[b]
                    nmb = QUAD_N[m, b]
                    u_m[m] += nmb * ucomp[nb]
                    w_m[m] += nmb * wcomp[nb]
                    q2_m[m] += nmb * q2[nb]
                    rho_m[m] += nmb * rho[nb]
            for a in range(3):
                for m in range(3):
                    w = wq * QUAD_N[m, a]
                    umw[a] += w * u_m[m]
                    wmw[a] += w * w_m[m]
                    q2mw[a] += w * q2_m[m]

            # ---- residuals ----
            for a in range(3):
                na = nn[a]
                acc = np.zeros(6)
                # divergence terms (Galerkin weight of a constant = wq)
                acc[0] += wq * divs[1] - divs[0] * umw[a]
                acc[1] += wq * divs[2] - divs[0] * wmw[a]
                acc[2] += wq * divs[3] - divs[0] * q2mw[a]
                acc[3] += wq * divs[4]
                acc[4] += wq * divs[5]
                acc[5] += wq * divs[6]
                # diffusion (integrated by parts; constant gradients)
                gab = veps * h_e * area
                gabs = veps_s * h_e * area
                gs1a = gN[a, 0] * s1e[0] + gN[a, 1] * s1e[1] + gN[a, 2] * s1e[2]
                for dcol in range(6):
                    dot = (
                        gdens[dcol, 0] * gN[a, 0]
                        + gdens[dcol, 1] * gN[a, 1]
                        + gdens[dcol, 2] * gN[a, 2]
                    )
                    dots = (
                        gdens[dcol, 0] * s1e[0]
                        + gdens[dcol, 1] * s1e[1]
                        + gdens[dcol, 2] * s1e[2]
                    )
                    acc[dcol] += gab * dot + gabs * dots * gs1a
                # quad-point source integrals
                for m in range(3):
                    w = wq * QUAD_N[m, a]
                    # interpolated flux vectors dotted with constant grads
                    q3g = 0.0
                    qo3g = 0.0
                    e3psi = 0.0
                    taux_m = 0.0
                    tauz_m = 0.0
                    dd_m = 0.0
                    dx_m = 0.0
                    do_m = 0.0
                    s1_m = 0.0
                    s2_m = 0.0
                    for b in range(3):
                        nb = nn[b]
                        nmb = QUAD_N[m, b]
                        for k in range(3):
                            q3g += nmb * fv[nb, F_Q3 + k] * gq2[k]
                            qo3g += nmb * fv[nb, F_QO3 + k] * gq2[k]
                            e3psi += nmb * fv[nb, F_E3 + k] * gpsi[k]
                        taux_m += nmb * fv[nb, F_TAUX]
                        tauz_m += nmb * fv[nb, F_TAUZ]
                        dd_m += nmb * fv[nb, F_DD]
                        dx_m += nmb * fv[nb, F_DX]
                        do_m += nmb * fv[nb, F_DO]
                        s1_m += nmb * fv[nb, F_S1]
                        s2_m += nmb * fv[nb, F_S2]
                    # (Q x grad(q^2)) . n as the scalar triple product of
                    # the interpolated 3-D Q, the element q^2 gradient and
                    # the interpolated local normal; on planar y-normal
                    # meshes nm=(0,1,0) and this reduces operation-for-
                    # operation to the pre-V3 q3z*gq2x - q3x*gq2z.
                    crossy = (
                        (qm_tab[m, 1] * gq2[2] - qm_tab[m, 2] * gq2[1])
                        * nm_tab[m, 0]
                        + (qm_tab[m, 2] * gq2[0] - qm_tab[m, 0] * gq2[2])
                        * nm_tab[m, 1]
                        + (qm_tab[m, 0] * gq2[1] - qm_tab[m, 1] * gq2[0])
                        * nm_tab[m, 2]
                    )
                    acc[0] += w * (-taux_m)
                    acc[1] += w * (-tauz_m)
                    acc[2] += w * (-rho_m[m] * q3g - 2.0 * dd_m)
                    acc[3] += w * (
                        e3psi + 0.5 * rho_m[m] * crossy - rho_m[m] * qo3g
                        + dx_m - 2.0 * do_m
                    )
                    acc[4] += w * (-s1_m)
                    acc[5] += w * (-s2_m)
                R[na, 0] += acc[0]
                R[na, 1] += acc[1]
                R[na, 2] += acc[2]
                R[na, 3] += acc[3]
                R[na, 4] += acc[4]
                R[na, 5] += acc[5]

            # ---- Jacobian ----
            if do_jac:
                for a in range(3):
                    for b in range(3):
                        nb = nn[b]
                        je = np.zeros((6, 6))
                        gab = (
                            gN[a, 0] * gN[b, 0]
                            + gN[a, 1] * gN[b, 1]
                            + gN[a, 2] * gN[b, 2]
                        ) * (veps * h_e * area)
                        gabs = (
                            (gN[a, 0] * s1e[0] + gN[a, 1] * s1e[1]
                             + gN[a, 2] * s1e[2])
                            * (gN[b, 0] * s1e[0] + gN[b, 1] * s1e[1]
                               + gN[b, 2] * s1e[2])
                        ) * (veps_s * h_e * area)
                        for l in range(6):
                            # divergence + diffusion rows
                            je[0, l] += wq * ddiv[1, b, l] - ddiv[0, b, l] * umw[a]
                            je[1, l] += wq * ddiv[2, b, l] - ddiv[0, b, l] * wmw[a]
                            je[2, l] += wq * ddiv[3, b, l] - ddiv[0, b, l] * q2mw[a]
                            je[3, l] += wq * ddiv[4, b, l]
                            je[4, l] += wq * ddiv[5, b, l]
                            je[5, l] += wq * ddiv[6, b, l]
                            for dcol in range(6):
                                je[dcol, l] += (gab + gabs) * dv[nb, F_DMX + dcol, l]
                        # quad-point source derivatives
                        for m in range(3):
                            w = wq * QUAD_N[m, a] * QUAD_N[m, b]
                            if w == 0.0:
                                continue
                            for l in range(6):
                                dq3g = 0.0
                                dqo3g = 0.0
                                de3psi = 0.0
                                for k in range(3):
                                    dq3g += dv[nb, F_Q3 + k, l] * gq2[k]
                                    dqo3g += dv[nb, F_QO3 + k, l] * gq2[k]
                                    de3psi += dv[nb, F_E3 + k, l] * gpsi[k]
                                gxn0 = gq2[1] * nm_tab[m, 2] - gq2[2] * nm_tab[m, 1]
                                gxn1 = gq2[2] * nm_tab[m, 0] - gq2[0] * nm_tab[m, 2]
                                gxn2 = gq2[0] * nm_tab[m, 1] - gq2[1] * nm_tab[m, 0]
                                dcross = (
                                    dv[nb, F_Q3 + 0, l] * gxn0
                                    + dv[nb, F_Q3 + 1, l] * gxn1
                                    + dv[nb, F_Q3 + 2, l] * gxn2
                                )
                                je[0, l] += w * (-dv[nb, F_TAUX, l])
                                je[1, l] += w * (-dv[nb, F_TAUZ, l])
                                je[2, l] += w * (
                                    -rho_m[m] * dq3g - 2.0 * dv[nb, F_DD, l]
                                )
                                je[3, l] += w * (
                                    de3psi + 0.5 * rho_m[m] * dcross
                                    - rho_m[m] * dqo3g
                                    + dv[nb, F_DX, l] - 2.0 * dv[nb, F_DO, l]
                                )
                                je[4, l] += w * (-dv[nb, F_S1, l])
                                je[5, l] += w * (-dv[nb, F_S2, l])
                        # scatter 6x6 block into CSR
                        for ceq in range(6):
                            for l in range(6):
                                Jdata[elem_to_csr[e, a, b, ceq, l]] += je[ceq, l]


# ---------------------------------------------------------------------------
# Solver driver
# ---------------------------------------------------------------------------

class IBL3Solver:
    """Newton/pseudo-time solver for the six IBL3 equations on a SurfaceMesh.

    Parameters
    ----------
    smesh : SurfaceMesh (planar for V1; normals give the wall direction)
    u_e   : (n,3) prescribed edge velocity at surface nodes (tangent)
    rho, mu, mach : scalars or (n,) edge fluid data
    turbulent_flags : (n,) int 0/1 regime per node (D-TR forced transition;
        laminar nodes get their stress rows pinned to CTAU_LAM)
    inflow_mask : (n,) bool Dirichlet nodes; inflow_state : (6,) or (n_in,6)
    c_l : dissipation length constant (D-CT-2); eps_diff : diffusion
        parameter epsilon in V_eps = eps * max(q) (D-HB, recorded in VERDICT)
    eps_diff_s : streamwise-tensor diffusion coefficient (D-HB follow-up,
        implemented after the GV1.1(e) finding): adds
        eps_diff_s * max(q) * h_e * (s1.grad G)(s1.grad N), damping the
        centered-Galerkin 2h streamwise mode. Default 0.02, calibrated on
        the GV1.1(a)/(e) family: the smallest value giving strictly
        decreasing refinement errors is 0.01 (marginal, ds order ~0.5);
        0.02 keeps the strict decrease with H order ~1.0 and real damping
        margin at ~1.3x the absolute error (see VERDICT addendum);
        0.0 recovers the isotropic-only scheme.
    """

    def __init__(self, smesh, u_e, rho, mu, mach, turbulent_flags,
                 inflow_mask, inflow_state, c_l=C.C_L_DEFAULT,
                 eps_diff=0.005, eps_diff_s=0.02):
        sm = smesh
        n = sm.xyz.shape[0]
        self.smesh = sm
        self.n = n
        u_e = np.ascontiguousarray(u_e, dtype=np.float64)
        if u_e.shape != (n, 3):
            raise ValueError(f"u_e shape {u_e.shape}, expected ({n}, 3)")
        if not np.all(np.isfinite(u_e)):
            raise ValueError("u_e contains non-finite values")

        def _arr(v, name):
            a = np.asarray(v, dtype=np.float64)
            if a.ndim == 0:
                a = np.full(n, float(a))
            if a.shape != (n,):
                raise ValueError(f"{name} shape {a.shape}, expected ({n},)")
            return np.ascontiguousarray(a)

        self._rho = _arr(rho, "rho")
        self._mu = _arr(mu, "mu")
        self._mach = _arr(mach, "mach")
        self._flags = np.ascontiguousarray(turbulent_flags, dtype=np.int64)
        if self._flags.shape != (n,):
            raise ValueError("turbulent_flags shape mismatch")
        inflow_mask = np.asarray(inflow_mask, dtype=bool)
        if inflow_mask.shape != (n,):
            raise ValueError("inflow_mask shape mismatch")
        self._inflow_idx = np.where(inflow_mask)[0].astype(np.int64)
        st = np.asarray(inflow_state, dtype=np.float64)
        if st.shape == (NUNK,):
            st = np.tile(st, (self._inflow_idx.size, 1))
        if st.shape != (self._inflow_idx.size, NUNK):
            raise ValueError("inflow_state shape mismatch")
        self._inflow_state = np.ascontiguousarray(st)
        # laminar pin rows (skip nodes already fully Dirichlet)
        lam = np.where(self._flags == 0)[0]
        lam_set = np.setdiff1d(lam, self._inflow_idx)
        self._lam_idx = lam_set.astype(np.int64)
        self.c_l = float(c_l)
        self.eps_diff = float(eps_diff)
        self.eps_diff_s = float(eps_diff_s)

        # frames
        self._q = np.empty(n)
        self._s1 = np.empty((n, 3))
        self._s2 = np.empty((n, 3))
        self._u = np.empty(n)
        self._w = np.empty(n)
        self._q2 = np.empty(n)
        self._psi = np.empty(n)
        self._s1l = np.empty((n, 2))
        self._s2l = np.empty((n, 2))
        # right-handed local "z" (D13 SIII.D.1): on a planar +-y-normal
        # surface this is exactly +z_global, so planar cases are
        # bit-identical to the pre-V3 global-basis form.
        self._basis_zc = np.ascontiguousarray(
            np.cross(sm.basis_x, sm.basis_y))
        _frames(u_e, sm.basis_y, sm.basis_x, self._basis_zc,
                self._q, self._s1, self._s2, self._u, self._w, self._q2,
                self._psi, self._s1l, self._s2l)
        self._veps = self.eps_diff * max(float(np.max(self._q)), 1.0e-30)
        self._veps_s = self.eps_diff_s * max(float(np.max(self._q)), 1.0e-30)

        # Jacobian pattern and diagonal indices
        self.pattern, self.elem_to_csr = sm.build_jacobian_pattern()
        self.Jdata = self.pattern.data
        indptr = self.pattern.indptr
        indices = self.pattern.indices
        nrow = NUNK * n
        diag = np.empty(nrow, dtype=np.int64)
        for r in range(nrow):
            lo = indptr[r]
            hi = indptr[r + 1]
            found = -1
            for p in range(lo, hi):
                if indices[p] == r:
                    found = p
                    break
            if found < 0:
                raise RuntimeError(f"CSR pattern missing diagonal in row {r}")
            diag[r] = found
        self._diag = diag
        # positions of the full diagonal 6x6 block of each node (for the
        # pseudo-time Jacobian, which is nodal block-diagonal)
        diagblk = np.empty((n, NUNK, NUNK), dtype=np.int64)
        for i in range(n):
            for ceq in range(NUNK):
                r = NUNK * i + ceq
                for p in range(indptr[r], indptr[r + 1]):
                    col = indices[p]
                    if NUNK * i <= col < NUNK * i + NUNK:
                        diagblk[i, ceq, col - NUNK * i] = p
        self._diagblk = diagblk

        # work arrays (allocated once; hot path is zero-alloc)
        self.outs = np.empty((n, C.N_OUT))
        self.douts = np.empty((n, C.N_OUT, NUNK))
        self.fv = np.empty((n, NV))
        self.dv = np.empty((n, NV, NUNK))

    # -- internal -----------------------------------------------------------
    def _nodal(self, U):
        C.closure_all(U, self._q, self._rho, self._mu, self._mach,
                      self._flags, self.c_l, self.outs, self.douts)
        _nodal_fluxes(U, self._q, self._s1, self._s2, self._s1l, self._s2l,
                      self._rho, self.outs, self.douts, self.c_l,
                      self.fv, self.dv)

    def _G_current(self):
        """Physical conserved densities G_c at the current fv tables
        (D13 (70)-(72) time terms with steady edge data):
        G = (Mx - u m, Mz - w m, e - q^2 m, k_o, k_tau1, k_tau2)."""
        fv = self.fv
        G = np.empty((self.n, NUNK))
        G[:, 0] = fv[:, F_DMX] - self._u * fv[:, F_DM]
        G[:, 1] = fv[:, F_DMZ] - self._w * fv[:, F_DM]
        G[:, 2] = fv[:, F_DE] - self._q2 * fv[:, F_DM]
        G[:, 3] = fv[:, F_DKO]
        G[:, 4] = fv[:, F_DKT1]
        G[:, 5] = fv[:, F_DKT2]
        return G

    def _dG_current(self):
        """dG_c/d(state) (n,6,6) at the current dv tables."""
        dv = self.dv
        dG = np.empty((self.n, NUNK, NUNK))
        dG[:, 0, :] = dv[:, F_DMX, :] - self._u[:, None] * dv[:, F_DM, :]
        dG[:, 1, :] = dv[:, F_DMZ, :] - self._w[:, None] * dv[:, F_DM, :]
        dG[:, 2, :] = dv[:, F_DE, :] - self._q2[:, None] * dv[:, F_DM, :]
        dG[:, 3, :] = dv[:, F_DKO, :]
        dG[:, 4, :] = dv[:, F_DKT1, :]
        dG[:, 5, :] = dv[:, F_DKT2, :]
        return dG

    def _assemble_into(self, U, R, do_jac):
        self._nodal(U)
        _assemble(self.fv, self.dv, self._u, self._w, self._q2, self._psi,
                  self._rho, self.smesh.triangles, self.smesh.gradN,
                  self.smesh.area_tri, self.smesh.color_offsets,
                  self.smesh.color_elems, self._veps, self._veps_s, self._s1,
                  self.smesh.basis_y,
                  R, self.Jdata, self.elem_to_csr, do_jac)

    def _apply_rows(self, R, U, do_jac):
        """Dirichlet inflow rows and laminar stress-pin rows (D-TR/D-BC)."""
        indptr = self.pattern.indptr
        for ii in range(self._inflow_idx.size):
            i = self._inflow_idx[ii]
            for ceq in range(NUNK):
                R[i, ceq] = U[i, ceq] - self._inflow_state[ii, ceq]
            if do_jac:
                for r in range(NUNK * i, NUNK * i + NUNK):
                    self.Jdata[indptr[r]:indptr[r + 1]] = 0.0
                    self.Jdata[self._diag[r]] = 1.0
        for ii in range(self._lam_idx.size):
            i = self._lam_idx[ii]
            R[i, 4] = U[i, 4] - C.CTAU_LAM
            R[i, 5] = U[i, 5] - C.CTAU_LAM
            if do_jac:
                for r in (NUNK * i + 4, NUNK * i + 5):
                    self.Jdata[indptr[r]:indptr[r + 1]] = 0.0
                    self.Jdata[self._diag[r]] = 1.0

    # -- public -------------------------------------------------------------
    def residual(self, U, dt_inv=None, U_old=None):
        """Steady residual (n,6); optionally with backward-Euler pseudo-time
        on the physical conserved densities (D13 (70)-(72), design D-PT)."""
        U = np.ascontiguousarray(U, dtype=np.float64)
        R = np.zeros((self.n, NUNK))
        self._assemble_into(U, R, False)
        if dt_inv is not None:
            G = self._G_current()
            self._nodal(U_old)
            G_old = self._G_current()
            R += (self.smesh.node_area * dt_inv)[:, None] * (G - G_old)
        self._apply_rows(R, U, False)
        return R

    def residual_jacobian(self, U, dt_inv=None, U_old=None):
        """Residual (n,6) and Jacobian (scipy csr, 6n x 6n)."""
        U = np.ascontiguousarray(U, dtype=np.float64)
        R = np.zeros((self.n, NUNK))
        self.Jdata[:] = 0.0
        self._assemble_into(U, R, True)
        if dt_inv is not None:
            G = self._G_current()
            dG = self._dG_current()
            self._nodal(U_old)
            G_old = self._G_current()
            w = self.smesh.node_area * dt_inv
            R += w[:, None] * (G - G_old)
            # nodal block-diagonal pseudo-time Jacobian
            self.Jdata[self._diagblk] += w[:, None, None] * dG
        self._apply_rows(R, U, True)
        return R, self.pattern

    def solve(self, U0, cfl0=1.0, growth=2.0, cfl_max=1.0e8, tol=1.0e-9,
              max_iter=100, verbose=False):
        """Pseudo-time-globalized Newton. Returns (U, info dict).

        Convergence is judged on the pure steady residual (pseudo-time term
        excluded); the recorded history carries it per accepted iterate.
        Backtracking uses the pseudo-time residual F_pt = R + w(G - G_old)
        as merit function: the step is the Newton step of F_pt (its zero
        level set moves with U_old), so F_pt is the only norm the linearized
        model is guaranteed to descend. Judging the same step on the pure
        steady residual stalls whenever the pseudo-time weight is non-
        negligible (observed on the FS decelerating branch: every step
        rejected from a near-solution seed, cfl collapsing to cfl_min).
        """
        U = np.array(U0, dtype=np.float64, copy=True)
        U_old = U.copy()
        h = np.sqrt(np.maximum(self.smesh.node_area, 1.0e-30))
        cfl = cfl0
        cfl_min = cfl0 * 1.0e-3
        pnorm0 = None
        hist = []
        converged = False
        n_iter = 0
        n_fail = 0
        for it in range(max_iter):
            n_iter = it + 1
            dt_inv = self._q / (cfl * h)
            R, J = self.residual_jacobian(U, dt_inv=dt_inv, U_old=U_old)
            rpure = self.residual(U)
            pnorm = float(np.max(np.abs(rpure)))
            fpnorm = float(np.max(np.abs(R)))
            if pnorm0 is None:
                pnorm0 = max(pnorm, 1.0e-30)
            hist.append(pnorm)
            if verbose:
                print(f"[ibl3] it={it} cfl={cfl:.3e} |R|inf={pnorm:.3e}")
            if pnorm < tol * max(pnorm0, 1.0):
                converged = True
                break
            delta = scipy.sparse.linalg.spsolve(J, -R.ravel())
            step_ok = np.all(np.isfinite(delta))
            omega = 1.0
            U_new = U
            if step_ok:
                dU = delta.reshape(self.n, NUNK)
                # halving backtracking on the pseudo-time residual (the
                # merit function the step actually linearizes)
                step_ok = False
                for _ls in range(12):
                    U_new = U + omega * dU
                    fnew = float(np.max(np.abs(self.residual(
                        U_new, dt_inv=dt_inv, U_old=U_old))))
                    if fnew <= (1.0 - 1.0e-4 * omega) * fpnorm:
                        step_ok = True
                        break
                    omega *= 0.5
            if not step_ok:
                # failed linear solve or no decrease: shrink the pseudo
                # time step and retry from the same iterate
                n_fail += 1
                if n_fail > 10:
                    break
                cfl = max(cfl / 8.0, cfl_min)
                if verbose:
                    print(f"[ibl3]   step rejected, cfl <- {cfl:.3e}")
                continue
            n_fail = 0
            U_old = U
            U = U_new
            cfl = min(cfl * growth, cfl_max)
        rpure = self.residual(U)
        info = {
            "converged": converged,
            "n_iter": n_iter,
            "residual_history": hist,
            "final_residual": float(np.max(np.abs(rpure))),
            "cfl_final": cfl,
            "eps_diff": self.eps_diff,
            "eps_diff_s": self.eps_diff_s,
            "veps": self._veps,
        }
        return U, info
