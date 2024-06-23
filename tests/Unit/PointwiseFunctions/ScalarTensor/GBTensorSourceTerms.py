# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def spacetime_derivative_scalar(lapse, shift, pi_scalar, phi_scalar):
    # n_{0} Pi + gamma^{i}_{0} Phi_{i}; where n_{0} = - lapse
    deriv_0 = -lapse * pi_scalar + np.dot(shift, phi_scalar)

    spacetime_derivative_scalar = np.concatenate(([deriv_0], phi_scalar))

    return spacetime_derivative_scalar


def DDKG_normal_normal_projection(
    lapse,
    shift,
    inverse_spatial_metric,
    phi_scalar,
    d_pi_scalar,
    dt_pi_scalar,
    d_lapse,
):
    lie_pi = dt_pi_scalar - np.dot(shift, d_pi_scalar)
    lie_pi *= 1.0 / lapse

    acc_dot_phi = np.transpose(phi_scalar) @ inverse_spatial_metric @ d_lapse
    acc_dot_phi *= 1.0 / lapse

    result = -lie_pi - acc_dot_phi

    return result


def DDKG_normal_spatial_projection(
    inverse_spatial_metric, extrinsic_curvature, phi_scalar, d_pi_scalar
):
    result = (
        extrinsic_curvature @ inverse_spatial_metric @ phi_scalar
    ) - d_pi_scalar

    return result


def DDKG_spatial_spatial_projection(
    extrinsic_curvature,
    spatial_christoffel_second_kind,
    pi_scalar,
    phi_scalar,
    d_phi_scalar,
):
    chris_times_phi = np.einsum(
        "kij, k->ij", spatial_christoffel_second_kind, phi_scalar
    )

    # Random number test does not impose symmetry of d_phi
    # D_phi = d_phi_scalar - chris_times_phi

    # Symmetrizing
    D_phi = 0.5 * (d_phi_scalar + np.transpose(d_phi_scalar)) - chris_times_phi

    result = D_phi - pi_scalar * extrinsic_curvature

    return result


def DDKG_tensor_from_projections(lapse, shift, nnDDKG, nsDDKG, ssDDKG):
    DDKG = np.zeros((4, 4))

    # 00-component
    DDKG[0, 0] = (
        lapse * lapse * nnDDKG
        + 2.0 * lapse * np.dot(shift, nsDDKG)
        + np.transpose(shift) @ ssDDKG @ shift
    )

    # 0i-component
    DDKG[0, 1:] = lapse * nsDDKG + np.transpose(shift) @ ssDDKG

    # Impose symmetry
    DDKG[1:, 0] = np.transpose(DDKG[0, 1:])

    # ij-component
    DDKG[1:, 1:] = ssDDKG

    return DDKG


def fPsi(psi, eta, zeta):
    return (eta / 8.0) * np.power(psi, 2) + (zeta / 16.0) * np.power(psi, 4)


def first_derivative_fPsi(psi, eta, zeta):
    return (1.0 / 4.0) * (eta * psi + zeta * np.power(psi, 3))


def second_derivative_fPsi(psi, eta, zeta):
    return (1.0 / 4.0) * (eta + 3.0 * zeta * np.power(psi, 2))


def DDFPsi_tensor_from_DDKG_tensor(
    DDKG,
    spacetime_derivative_scalar,
    psi,
    eta,  # first_coupling_psi,
    zeta,  # second_coupling_psi,
):
    result = second_derivative_fPsi(psi, eta, zeta) * np.outer(
        spacetime_derivative_scalar, spacetime_derivative_scalar
    )

    result += first_derivative_fPsi(psi, eta, zeta) * DDKG

    return result


def DDFPsi_normal_normal_projection(
    psi,
    pi_scalar,
    phi_scalar,
    DDKG_normal_normal_projection,
    eta,  #     first_coupling_psi,
    zeta,  #    second_coupling_psi,
):
    result = (
        second_derivative_fPsi(psi, eta, zeta) * np.power(pi_scalar, 2)
        + first_derivative_fPsi(psi, eta, zeta) * DDKG_normal_normal_projection
    )

    return result


def DDFPsi_spatial_normal_projection(
    psi,
    pi_scalar,
    phi_scalar,
    DDKG_spatial_normal_projection,
    eta,  # first_coupling_psi,
    zeta,  # second_coupling_psi,
):
    result = (
        second_derivative_fPsi(psi, eta, zeta) * (-1.0 * pi_scalar * phi_scalar)
        + first_derivative_fPsi(psi, eta, zeta) * DDKG_spatial_normal_projection
    )

    return result


def DDFPsi_spatial_spatial_projection(
    psi,
    pi_scalar,
    phi_scalar,
    DDKG_spatial_spatial_projection,
    eta,  # first_coupling_psi,
    zeta,  # second_coupling_psi,
):
    result = (
        second_derivative_fPsi(psi, eta, zeta)
        * np.outer(phi_scalar, phi_scalar)
        + first_derivative_fPsi(psi, eta, zeta)
        * DDKG_spatial_spatial_projection
    )

    return result


def order_reduced_gb_H_normal_normal_projection(
    inverse_spatial_metric, weyl_electric, ssDDKG
):
    result = np.trace(
        weyl_electric @ inverse_spatial_metric @ ssDDKG @ inverse_spatial_metric
    )

    return result


def compute_S_cross_B(inverse_spatial_metric, weyl_magnetic, ssDDKG):
    # Without the sqrt(gamma) factor

    B_down_up = weyl_magnetic @ inverse_spatial_metric

    S_up_up = inverse_spatial_metric @ ssDDKG @ inverse_spatial_metric

    result = np.zeros(3)
    for l in range(0, 3):
        result += np.cross(S_up_up[l, :], B_down_up[l, :])

    return result


def compute_j_cross_B(inverse_spatial_metric, weyl_magnetic, nsDDKG):
    # Without the sqrt(gamma) factor

    j_up = inverse_spatial_metric @ nsDDKG

    B_down_up = weyl_magnetic @ inverse_spatial_metric

    result = np.zeros((3, 3))
    for j in range(0, 3):
        result[:, j] = np.transpose(np.cross(j_up, B_down_up[j, :]))

    return result


def order_reduced_gb_H_normal_spatial_projection(
    inverse_spatial_metric,
    sqrt_det_spatial_metric,
    weyl_electric,
    nsDDKG,
    S_cross_B,
):
    # Raise index
    j_up = inverse_spatial_metric @ nsDDKG

    # Weyl electric term
    result = weyl_electric @ j_up

    # Weyl magnetic term
    result += sqrt_det_spatial_metric * S_cross_B

    return result


def order_reduced_gb_H_spatial_spatial_projection(
    spatial_metric,
    inverse_spatial_metric,
    sqrt_det_spatial_metric,
    weyl_electric,
    nnDDKG,
    ssDDKG,
    j_cross_B,
    nnH,
):
    trace_S = np.trace(ssDDKG @ inverse_spatial_metric)
    rho = nnDDKG

    S_down_up = ssDDKG @ inverse_spatial_metric
    S_up_up = inverse_spatial_metric @ S_down_up
    S_times_E = S_down_up @ weyl_electric
    trES = nnH  # To avoid recomputing the trace in the cpp function

    # First term
    result = (trace_S + rho) * weyl_electric

    # Second term
    result += -(S_times_E + np.transpose(S_times_E))

    # Third term
    result += sqrt_det_spatial_metric * (j_cross_B + np.transpose(j_cross_B))

    # Fourth term
    result += spatial_metric * trES

    return result


def order_reduced_gb_H_tensor_weyl_part(lapse, shift, nnH, nsH, ssH):
    rho = nnH
    j_vec = nsH
    S = ssH

    result = np.zeros((4, 4))

    # 00-component
    result[0, 0] = (
        np.power(lapse, 2) * rho
        + 2.0 * lapse * np.dot(shift, j_vec)
        + shift @ S @ shift
    )
    # 0i-component
    result[0, 1:] = lapse * j_vec + S @ shift

    # ij-component
    result[1:, 1:] = S

    return result


def order_reduced_Q_tensor(
    spacetime_metric, weyl_electric_scalar, weyl_magnetic_scalar
):
    result = (
        2.0 * (weyl_electric_scalar - weyl_magnetic_scalar) * spacetime_metric
    )

    return result


def order_reduced_gb_H_tensor_ricci_part(
    g, inv_g, T, trace_T, DDKGUpUp, Tdriv, trace_Tdriv
):
    Ricci = T + Tdriv
    trRicci = trace_T + trace_Tdriv

    # First term
    first_term = -0.5 * (
        np.einsum("ac,db,bd->ac", g, Ricci, DDKGUpUp)
        - np.einsum("ad,cb,bd->ac", g, Ricci, DDKGUpUp)
    )

    # Second term
    second_term = 0.5 * (
        np.einsum("bc,da,bd->ac", g, Ricci, DDKGUpUp)
        - np.einsum("bd,ca,bd->ac", g, Ricci, DDKGUpUp)
    )

    # Third term
    third_term = 0.5 * (
        np.einsum("ac,db,bd->ac", g, g, DDKGUpUp)
        - np.einsum("ad,cb,bd->ac", g, g, DDKGUpUp)
    )
    third_term *= (2.0 / 3.0) * trRicci

    result = first_term + second_term + third_term

    return result


def order_reduced_trace_reversed_stress_energy(
    spacetime_metric,
    inverse_spacetime_metric,
    gb_H_tensor_ricci_part,
    gb_H_tensor_weyl_part,
    trace_reversed_canonical_stress_energy,
):
    kappa = 1.0
    H_prefactor = -8.0

    H_SET = (
        H_prefactor * kappa * (gb_H_tensor_weyl_part + gb_H_tensor_ricci_part)
    )
    trace_H_SET = np.trace(inverse_spacetime_metric @ H_SET)

    # Trace-reverse the H tensor
    result = H_SET - 0.5 * trace_H_SET * spacetime_metric

    return result
