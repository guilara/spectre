// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/GBTensorSourceTerms.hpp"

#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"

namespace ScalarTensor {

void spacetime_derivative_scalar(const gsl::not_null<tnsr::a<DataVector, 3>*>
                                     spacetime_derivative_scalar_result,
                                 const Scalar<DataVector>& lapse,
                                 const Scalar<DataVector>& pi_scalar,
                                 const tnsr::i<DataVector, 3>& phi_scalar) {
  // Assemble in symmetric rank-2 4-tensor with lower indices
  // partial_a = Phi_a + n_a Pi
  // with n_0 = - lapse and n_i = 0
  get<0>(*spacetime_derivative_scalar_result) = -get(lapse) * get(pi_scalar);
  for (size_t i = 0; i < 3; ++i) {
    spacetime_derivative_scalar_result->get(i + 1) = phi_scalar.get(i);
  }
}

void DDKG_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> DDKG_normal_normal_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    // Scalar quantities
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_pi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const Scalar<DataVector>& dt_pi_scalar,

    const tnsr::i<DataVector, 3>& d_lapse) {
  // = - L_n Pi from the equations of motion
  tenex::evaluate(
      DDKG_normal_normal_result,
      // - L_n Pi - (1/lapse) Phi^{i} partial_i lapse
      -(1.0 / lapse()) * (dt_pi_scalar() - shift(ti::I) * d_pi_scalar(ti::i)

                          + inverse_spatial_metric(ti::I, ti::J) *
                                phi_scalar(ti::i) * d_lapse(ti::j)

                              )

  );
}

void DDKG_normal_spatial_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> DDKG_normal_spatial_result,

    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,

    // Scalar quantities
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_pi_scalar) {
  tenex::evaluate<ti::i>(
      DDKG_normal_spatial_result, extrinsic_curvature(ti::i, ti::j) *
                                          inverse_spatial_metric(ti::J, ti::K) *
                                          phi_scalar(ti::k)

                                      - d_pi_scalar(ti::i)

  );
}

void DDKG_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> DDKG_spatial_spatial_result,

    // Metric quantities

    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const tnsr::Ijj<DataVector, 3>& spatial_christoffel_second_kind,

    // Scalar quantities
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::ij<DataVector, 3>& d_phi_scalar

    // Provide them with RHS compute tags or from dt<> prefixes
) {
  // Note that D_phi is the covariant derivative and has Christoffel symbols
  tenex::evaluate<ti::i, ti::j>(
      DDKG_spatial_spatial_result,
      -pi_scalar() * extrinsic_curvature(ti::i, ti::j)
          // Note covariant derivative
          // and symmetrize partial derivative of scalar
          + 0.5 * (d_phi_scalar(ti::i, ti::j) + d_phi_scalar(ti::j, ti::i)) -
          spatial_christoffel_second_kind(ti::K, ti::i, ti::j) *
              phi_scalar(ti::k));
}

void DDKG_tensor_from_projections(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDKG_tensor_result,
    // Metric quantities
    const Scalar<DataVector>& lapse,
    // Projections
    const Scalar<DataVector>& nnDDKG, const tnsr::i<DataVector, 3>& nsDDKG,
    const tnsr::ii<DataVector, 3>& ssDDKG) {
  // Assemble in symmetric rank-2 4-tensor with lower indices
  get<0, 0>(*DDKG_tensor_result) = square(get(lapse)) * get(nnDDKG);
  for (size_t i = 0; i < 3; ++i) {
    // nsH with lower indices
    // Check sign
    DDKG_tensor_result->get(0, i + 1) = get(lapse) * nsDDKG.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      // ssH with lower indices
      DDKG_tensor_result->get(i + 1, j + 1) = ssDDKG.get(i, j);
    }
  }
}

void DDKG_trace_minus_eom(
    const gsl::not_null<Scalar<DataVector>*> diagnostic,
    const tnsr::aa<DataVector, 3>& DDKG,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
    const Scalar<DataVector>& dt_pi_scalar,
    const Scalar<DataVector>& scalar_driver) {
  // Check that the trace of the DDKG tensor is consistent with the equation of
  // motion
  // Check that Box Psi - source = 0, where the first term is computed from the
  // DDKG tensor
  tenex::evaluate(diagnostic,
                  // Trace
                  DDKG(ti::a, ti::b) * inverse_spacetime_metric(ti::B, ti::A) -
                      // Equation of motion source term
                      scalar_driver());
}

// template <typename Frame>
void DDKG_tensor_from_projections(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDKG_tensor_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const tnsr::Ijj<DataVector, 3>& spatial_christoffel_second_kind,

    // Scalar quantities
    const Scalar<DataVector>& psi_scalar, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_psi_scalar,
    const tnsr::i<DataVector, 3>& d_pi_scalar,
    const tnsr::ij<DataVector, 3>& d_phi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const Scalar<DataVector>& dt_psi_scalar,
    const Scalar<DataVector>& dt_pi_scalar,
    const tnsr::i<DataVector, 3>& dt_phi_scalar,

    const tnsr::i<DataVector, 3>& d_lapse) {
  // Compute projections of the second derivative tensor
  // nn
  Scalar<DataVector> nnDDKG;
  DDKG_normal_normal_projection(make_not_null(&nnDDKG), lapse, shift,
                                inverse_spatial_metric, phi_scalar, d_pi_scalar,
                                dt_pi_scalar, d_lapse);
  // ns
  tnsr::i<DataVector, 3> nsDDKG;
  DDKG_normal_spatial_projection(make_not_null(&nsDDKG), inverse_spatial_metric,
                                 extrinsic_curvature, phi_scalar, d_pi_scalar);
  // ss
  tnsr::ii<DataVector, 3> ssDDKG;
  DDKG_spatial_spatial_projection(make_not_null(&ssDDKG), extrinsic_curvature,
                                  spatial_christoffel_second_kind, pi_scalar,
                                  phi_scalar, d_phi_scalar);

  // Assemble in symmetric rank-2 4-tensor with lower indices
  get<0, 0>(*DDKG_tensor_result) = square(get(lapse)) * get(nnDDKG);
  for (size_t i = 0; i < 3; ++i) {
    // nsH with lower indices
    // Check sign
    DDKG_tensor_result->get(0, i + 1) = get(lapse) * nsDDKG.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      // ssH with lower indices
      DDKG_tensor_result->get(i + 1, j + 1) = ssDDKG.get(i, j);
    }
  }
}

void DDFPsi_tensor_from_DDKG_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDFPsi_tensor_result,
    const tnsr::aa<DataVector, 3>& DDKG,
    const tnsr::a<DataVector, 3>& spacetime_derivative_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi) {
  const double first_coupling_psi_over_four = first_coupling_psi / 4.0;
  const double second_coupling_psi_over_four = second_coupling_psi / 4.0;

  tenex::evaluate<ti::a, ti::b>(
      DDFPsi_tensor_result,
      // Double prime term
      -first_coupling_psi_over_four * spacetime_derivative_scalar(ti::a) *
              spacetime_derivative_scalar(ti::b) -
          3.0 * second_coupling_psi_over_four * psi() * psi() *
              spacetime_derivative_scalar(ti::a) *
              spacetime_derivative_scalar(ti::b) +
          // Prime term
          (-first_coupling_psi_over_four * psi() -
           second_coupling_psi_over_four * psi() * psi() * psi()) *
              DDKG(ti::a, ti::b)

  );
}

void raise_indices_DDFPsi(
    const gsl::not_null<tnsr::AA<DataVector, 3>*> DDFPsiUpUp_tensor_result,
    const tnsr::aa<DataVector, 3>& DDFPsi,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric) {
  tenex::evaluate<ti::A, ti::B>(DDFPsiUpUp_tensor_result,
                                inverse_spacetime_metric(ti::A, ti::C) *
                                    DDFPsi(ti::c, ti::d) *
                                    inverse_spacetime_metric(ti::D, ti::B));
}

void DDFPsi_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> DDFPsi_normal_normal_result,
    // Scalar quantities
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    // DDKG projections
    const Scalar<DataVector>& DDKG_normal_normal_projection,
    // Coupling function parameters
    const double first_coupling_psi, const double second_coupling_psi) {
  const double first_coupling_psi_over_four = first_coupling_psi / 4.0;
  const double second_coupling_psi_over_four = second_coupling_psi / 4.0;

  tenex::evaluate(DDFPsi_normal_normal_result,
                  // Double prime term
                  -first_coupling_psi_over_four * pi_scalar() * pi_scalar() -
                      3.0 * second_coupling_psi_over_four * psi() * psi() *
                          pi_scalar() * pi_scalar() +
                      // Prime term
                      (-first_coupling_psi_over_four * psi() -
                       second_coupling_psi_over_four * psi() * psi() * psi()) *
                          DDKG_normal_normal_projection()

  );
}

void DDFPsi_spatial_normal_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> DDFPsi_spatial_normal_result,
    // Scalar quantities
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    // DDKG projections
    const tnsr::i<DataVector, 3>& DDKG_spatial_normal_projection,
    // Coupling function parameters
    const double first_coupling_psi, const double second_coupling_psi) {
  const double first_coupling_psi_over_four = first_coupling_psi / 4.0;
  const double second_coupling_psi_over_four = second_coupling_psi / 4.0;

  tenex::evaluate<ti::i>(
      DDFPsi_spatial_normal_result,
      // Double prime term
      -(-first_coupling_psi_over_four * pi_scalar() * phi_scalar(ti::i) -
        3.0 * second_coupling_psi_over_four * psi() * psi() * pi_scalar() *
            phi_scalar(ti::i)) +
          // Prime term
          (-first_coupling_psi_over_four * psi() -
           second_coupling_psi_over_four * psi() * psi() * psi()) *
              DDKG_spatial_normal_projection(ti::i)

  );
}

void DDFPsi_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> DDFPsi_spatial_spatial_result,
    // Scalar quantities
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    // DDKG projections
    const tnsr::ii<DataVector, 3>& DDKG_spatial_spatial_projection,
    // Coupling function parameters
    const double first_coupling_psi, const double second_coupling_psi) {
  const double first_coupling_psi_over_four = first_coupling_psi / 4.0;
  const double second_coupling_psi_over_four = second_coupling_psi / 4.0;
  tenex::evaluate<ti::i, ti::j>(
      DDFPsi_spatial_spatial_result,
      // Double prime term
      -first_coupling_psi_over_four * phi_scalar(ti::i) * phi_scalar(ti::j) -
          3.0 * second_coupling_psi_over_four * psi() * psi() *
              phi_scalar(ti::i) * phi_scalar(ti::j) +
          // Prime term
          (-first_coupling_psi_over_four * psi() -
           second_coupling_psi_over_four * psi() * psi() * psi()) *
              DDKG_spatial_spatial_projection(ti::i, ti::j)

  );
}

void order_reduced_gb_H_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> nnH_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const tnsr::ii<DataVector, 3>& ssDDKG) {
  // Raise indices of the spatial part of the second derivative of the scalar
  tenex::evaluate(nnH_result, weyl_electric(ti::i, ti::j) *
                                  inverse_spatial_metric(ti::J, ti::K) *
                                  ssDDKG(ti::k, ti::l) *
                                  inverse_spatial_metric(ti::L, ti::I));
}

void compute_S_cross_B(
    const gsl::not_null<tnsr::i<DataVector, 3>*> S_cross_B_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& ssDDKG) {
  *S_cross_B_result = make_with_value<tnsr::i<DataVector, 3>>(
      get<0, 0>(inverse_spatial_metric), 0.0);

  for (LeviCivitaIterator<3> levi_civita_it; levi_civita_it; ++levi_civita_it) {
    const auto [i, j, k] = levi_civita_it();
    // S cross B
    // epsilon_{ijk} B_{l}^{k} S^{jl}
    for (size_t l = 0; l < 3; ++l) {
      for (size_t m = 0; m < 3; ++m) {
        for (size_t n = 0; n < 3; ++n) {
          for (size_t p = 0; p < 3; ++p) {
            S_cross_B_result->get(i) +=
                levi_civita_it.sign() *
                (
                    // Raise indices
                    inverse_spatial_metric.get(j, m) * ssDDKG.get(m, n) *
                    inverse_spatial_metric.get(l, n)

                        ) *
                (
                    // Raise index
                    weyl_magnetic.get(l, p) * inverse_spatial_metric.get(p, k));
          }
        }
      }
    }
  }
}

void compute_j_cross_B(
    const gsl::not_null<tnsr::ij<DataVector, 3>*> j_cross_B_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_magnetic,
    const tnsr::i<DataVector, 3>& nsDDKG) {
  *j_cross_B_result = make_with_value<tnsr::ij<DataVector, 3>>(
      get<0, 0>(inverse_spatial_metric), 0.0);

  for (LeviCivitaIterator<3> levi_civita_it; levi_civita_it; ++levi_civita_it) {
    const auto [i, j, k] = levi_civita_it();
    for (size_t l = 0; l < 3; ++l) {
      for (size_t m = 0; m < 3; ++m) {
        for (size_t n = 0; n < 3; ++n) {
          // j cross B
          // Note: For now we don't impose symmetry of this quantity
          j_cross_B_result->get(i, l) +=
              levi_civita_it.sign() *
              (nsDDKG.get(m) * inverse_spatial_metric.get(m, j)) *
              (weyl_magnetic.get(l, n) * inverse_spatial_metric.get(n, k));
        }
      }
    }
  }
}

void order_reduced_gb_H_normal_spatial_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> nsH_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const tnsr::i<DataVector, 3>& nsDDKG,
    const tnsr::i<DataVector, 3>& S_cross_B) {
  tenex::evaluate<ti::i>(
      nsH_result, weyl_electric(ti::i, ti::j) *
                          inverse_spatial_metric(ti::J, ti::K) * nsDDKG(ti::k) +
                      // sqrt(gamma) * epsilon_{ijk} B_{l}^{k} S^{jl}
                      sqrt_det_spatial_metric() * S_cross_B(ti::i));
}

void order_reduced_gb_H_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> ssH_result,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const Scalar<DataVector>& nnDDKG, const tnsr::ii<DataVector, 3>& ssDDKG,
    const tnsr::ij<DataVector, 3>& j_cross_B, const Scalar<DataVector>& nnH) {
  tenex::evaluate<ti::i, ti::j>(
      ssH_result,
      (
          // Trace of ssDDKG
          ssDDKG(ti::k, ti::l) * inverse_spatial_metric(ti::L, ti::K) +

          nnDDKG()) *
              weyl_electric(ti::i, ti::j)
          // -2 * 2 symmetric part E ssDDKG
          - (ssDDKG(ti::j, ti::k) * inverse_spatial_metric(ti::K, ti::L) *
                 weyl_electric(ti::l, ti::i) +
             ssDDKG(ti::i, ti::k) * inverse_spatial_metric(ti::K, ti::L) *
                 weyl_electric(ti::l, ti::j))
          // +2 sqrt(gamma) symmetric part cross prod
          + sqrt_det_spatial_metric() *
                (j_cross_B(ti::i, ti::j) + j_cross_B(ti::j, ti::i))
          // Sum the trace part
          + spatial_metric(ti::i, ti::j) * nnH());
}

void order_reduced_gb_H_tensor_weyl_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const Scalar<DataVector>& lapse, const Scalar<DataVector>& nnH,
    const tnsr::i<DataVector, 3>& nsH, const tnsr::ii<DataVector, 3>& ssH) {
  // Assemble in symmetric rank-2 4-tensor with lower indices
  get<0, 0>(*gb_H_tensor_result) = square(get(lapse)) * get(nnH);
  for (size_t i = 0; i < 3; ++i) {
    // nsH with lower indices
    // Check sign
    gb_H_tensor_result->get(0, i + 1) = get(lapse) * nsH.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      // ssH with lower indices
      gb_H_tensor_result->get(i + 1, j + 1) = ssH.get(i, j);
    }
  }
}

void order_reduced_Q_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> Q_tensor_result,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar) {
  tenex::evaluate<ti::a, ti::b>(
      Q_tensor_result, 2.0 * (weyl_electric_scalar() - weyl_magnetic_scalar()) *
                           spacetime_metric(ti::a, ti::b));
}

// void order_reduced_gb_H_tensor_ricci_part(
//     const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
//     const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
//     const tnsr::aa<DataVector, 3>& T, const tnsr::aa<DataVector, 3>& DDKG) {
//   // g: spacetime metric
//   // inv_g: inverse spacetime metric
//   // T: trace reversed stress energy
//   // 8.0 * M_PI factor in T definition
//   const double one_over_four = 1.0 / 4.0;
//   const double two_over_three = 2.0 / 3.0;

//   tenex::evaluate<ti::a, ti::c>(
//       gb_H_tensor_result, one_over_four *
//                               (-2.0 * g(ti::a, ti::c) * T(ti::b, ti::d) +
//                                g(ti::a, ti::d) * T(ti::b, ti::c) +
//                                g(ti::a, ti::b) * T(ti::c, ti::d) +
//                                g(ti::b, ti::c) * T(ti::a, ti::d) -
//                                2.0 * g(ti::b, ti::d) * T(ti::a, ti::c) +
//                                g(ti::c, ti::d) * T(ti::a, ti::b) +
//                                two_over_three *
//                                    (2.0 * g(ti::a, ti::c) * g(ti::b, ti::d) -
//                                     g(ti::a, ti::d) * g(ti::b, ti::c) -
//                                     g(ti::a, ti::b) * g(ti::c, ti::d)) *
//                                    (
//                                        // Trace T
//                                        T(ti::g, ti::h) * inv_g(ti::H, ti::G))

//                                    ) *
//                               inv_g(ti::B, ti::E) * inv_g(ti::D, ti::F) *
//                               DDKG(ti::e, ti::f));
// }

void order_reduced_gb_H_tensor_ricci_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
    const tnsr::aa<DataVector, 3>& T, const Scalar<DataVector>& trace_T,
    const tnsr::AA<DataVector, 3>& DDKGUpUp) {
  // g: spacetime metric
  // inv_g: inverse spacetime metric
  // T: trace reversed stress energy
  // 8.0 * M_PI factor in T definition

  const double one_over_three = 1.0 / 3.0;

  tenex::evaluate<ti::a, ti::b>(gb_H_tensor_result,
                                (

                                    0.5 * (-g(ti::c, ti::d) * T(ti::a, ti::b) +
                                           g(ti::b, ti::c) * T(ti::a, ti::d) +
                                           g(ti::a, ti::d) * T(ti::b, ti::c) -
                                           g(ti::a, ti::b) * T(ti::c, ti::d)) +

                                    one_over_three *
                                        (-g(ti::a, ti::d) * g(ti::b, ti::c) +
                                         g(ti::a, ti::b) * g(ti::c, ti::d)) *
                                        trace_T()

                                        ) *
                                    DDKGUpUp(ti::C, ti::D));
}

void order_reduced_gb_H_tensor_ricci_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
    const tnsr::aa<DataVector, 3>& T, const Scalar<DataVector>& trace_T,
    const tnsr::AA<DataVector, 3>& DDKGUpUp,
    const tnsr::aa<DataVector, 3>& Tdriv,
    const Scalar<DataVector>& trace_Tdriv) {
  // g: spacetime metric
  // inv_g: inverse spacetime metric
  // T: trace reversed stress energy
  // Tdriv: tensor_driver
  // 8.0 * M_PI factor in T definition

  const double one_over_three = 1.0 / 3.0;

  tenex::evaluate<ti::a, ti::b>(
      gb_H_tensor_result, (

                              0.5 * (-g(ti::c, ti::d) * T(ti::a, ti::b) +
                                     g(ti::b, ti::c) * T(ti::a, ti::d) +
                                     g(ti::a, ti::d) * T(ti::b, ti::c) -
                                     g(ti::a, ti::b) * T(ti::c, ti::d)) +

                              one_over_three *
                                  (-g(ti::a, ti::d) * g(ti::b, ti::c) +
                                   g(ti::a, ti::b) * g(ti::c, ti::d)) *
                                  trace_T() +

                              // Tensor driver part
                              0.5 * (-g(ti::c, ti::d) * Tdriv(ti::a, ti::b) +
                                     g(ti::b, ti::c) * Tdriv(ti::a, ti::d) +
                                     g(ti::a, ti::d) * Tdriv(ti::b, ti::c) -
                                     g(ti::a, ti::b) * Tdriv(ti::c, ti::d)) +

                              one_over_three *
                                  (-g(ti::a, ti::d) * g(ti::b, ti::c) +
                                   g(ti::a, ti::b) * g(ti::c, ti::d)) *
                                  trace_Tdriv()

                                  ) *
                              DDKGUpUp(ti::C, ti::D));
}

void order_reduced_trace_reversed_stress_energy(
    const gsl::not_null<tnsr::aa<DataVector, 3>*>
        order_reduced_trace_reversed_stress_energy_result,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
    const tnsr::aa<DataVector, 3>& gb_H_tensor_ricci_part,
    const tnsr::aa<DataVector, 3>& gb_H_tensor_weyl_part,
    const tnsr::aa<DataVector, 3>& trace_reversed_canonical_stress_energy) {
  // Sum and take the trace-reverse
  // For the Weyl part the trace should be zero, and taking the trace-reverse
  // should not change it
  // 8 pi factor needs to be added here since we absorbed it in the
  // trace reversed part of the canonical stress energy tensor

  // We work in units where set G = 1 / (8 M_PI)
  // const double kappa = 8.0 * M_PI;
  const double kappa = 1.0;

  const double H_tensor_prefactor = -8.0 * kappa;

  // Trace reversed
  // We avoid using tenex::update for now
  tenex::evaluate<ti::a, ti::b>(
      order_reduced_trace_reversed_stress_energy_result,
      H_tensor_prefactor * (gb_H_tensor_ricci_part(ti::a, ti::b) +
                            gb_H_tensor_weyl_part(ti::a, ti::b)) -
          0.5 * H_tensor_prefactor *
              (
                  // Trace of H
                  inverse_spacetime_metric(ti::C, ti::D) *
                  (gb_H_tensor_ricci_part(ti::c, ti::d) +
                   gb_H_tensor_weyl_part(ti::c, ti::d))

                      ) *
              spacetime_metric(ti::a, ti::b));
}

}  // namespace ScalarTensor

namespace fe::ScalarTensorDriver {

void tensor_driver_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> tensor_driver_nn,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::A<DataVector, 3>& spacetime_normal_vector) {
  // n^{a} tensor_driver_{a b} n^{b}
  tenex::evaluate(tensor_driver_nn, spacetime_normal_vector(ti::A) *
                                        tensor_driver(ti::a, ti::b) *
                                        spacetime_normal_vector(ti::B));
}

void tensor_driver_spatial_normal_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> tensor_driver_sn,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::A<DataVector, 3>& spacetime_normal_vector,
    const tnsr::a<DataVector, 3>& spacetime_normal_covector) {
  for (size_t i = 0; i < 3; ++i) {
    tensor_driver_sn->get(i) = 0.0;
    for (size_t a = 0; a < 4; ++a) {
      // n^{a} tensor_driver_{a b} (delta^{b}_{i} - n^{b} n_{i})
      // = n^{a} tensor_driver_{a i}
      tensor_driver_sn->get(i) +=
          spacetime_normal_vector.get(a) * tensor_driver.get(a, i + 1);
    }
  }
}

void tensor_driver_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> tensor_driver_ss,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::A<DataVector, 3>& spacetime_normal_vector,
    const tnsr::a<DataVector, 3>& spacetime_normal_covector) {
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      //  (delta^{a}_{i} - n^{a} n_{i}) * tensor_driver_{a b} *
      //  (delta^{b}_{j} - n^{b} n_{j})
      // = tensor_driver_{ij}
      tensor_driver_ss->get(i, j) = tensor_driver.get(i + 1, j + 1);
    }
  }
}

}  // namespace fe::ScalarTensorDriver
