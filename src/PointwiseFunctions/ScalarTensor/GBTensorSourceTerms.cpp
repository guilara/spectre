// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/GBTensorSourceTerms.hpp"

#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

namespace ScalarTensor {

void DDKG_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> DDKG_normal_normal_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,

    // Scalar quantities

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_pi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const Scalar<DataVector>& dt_pi_scalar) {
  // = - L_n Pi from the equations of motion
  tenex::evaluate(
      DDKG_normal_normal_result,
      -(1.0 / lapse()) * (dt_pi_scalar() - shift(ti::I) * d_pi_scalar(ti::i)));
}

void DDKG_normal_spatial_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> DDKG_normal_spatial_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,

    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,

    // Scalar quantities
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::ij<DataVector, 3>& d_phi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const tnsr::i<DataVector, 3>& dt_phi_scalar) {
  tenex::evaluate<ti::i>(
      DDKG_normal_spatial_result,
      extrinsic_curvature(ti::i, ti::j) * inverse_spatial_metric(ti::J, ti::K) *
              phi_scalar(ti::k)
          // + L_n Phi from the equations of motion
          + (1.0 / lapse()) * (dt_phi_scalar(ti::i) -
                               shift(ti::J) * d_phi_scalar(ti::j, ti::i)));
}

void DDKG_spatial_spatial_projection(
    const gsl::not_null<tnsr::ij<DataVector, 3>*> DDKG_spatial_spatial_result,

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
          + d_phi_scalar(ti::i, ti::j) -
          spatial_christoffel_second_kind(ti::K, ti::i, ti::j) *
              phi_scalar(ti::k));
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
    const tnsr::i<DataVector, 3>& dt_phi_scalar) {
  // Compute projections of the second derivative tensor
  // nn
  Scalar<DataVector> nnDDKG;
  DDKG_normal_normal_projection(make_not_null(&nnDDKG), lapse, shift,
                                d_pi_scalar, dt_pi_scalar);
  // ns
  tnsr::i<DataVector, 3> nsDDKG;
  DDKG_normal_spatial_projection(make_not_null(&nsDDKG), lapse, shift,
                                 inverse_spatial_metric, extrinsic_curvature,
                                 phi_scalar, d_phi_scalar, dt_phi_scalar);
  // ss
  tnsr::ij<DataVector, 3> ssDDKG;
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

/*
void order_reduced_gb_H_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> gb_H_tensor_result,
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> DDKG_tensor_result,
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const tnsr::Ijj<DataVector, 3>& spatial_christoffel_second_kind,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::aa<DataVector, Dim>& spacetime_metric,
    const tnsr::A<DataVector, 3>& normal_spacetime_vector,
    const tnsr::a<DataVector, 3>& normal_spacetime_one_form,
    const tnsr::aa<DataVector>& dd_coupling_function) {
  DDKG_tensor_from_projections(make_not_null(&DDKG_tensor_result),

                               // Metric quantities
                               lapse, shift, spacetime_metric,
                               inverse_spatial_metric, extrinsic_curvature,
                               spatial_christoffel_second_kind,

                               // Scalar quantities
                               psi_scalar, pi_scalar, phi_scalar,

                               // Scalar gradients
                               d_psi_scalar, d_pi_scalar, d_phi_scalar,

                               // Provide them with RHS compute tags
                               dt_psi_scalar, dt_pi_scalar, dt_phi_scalar);

  // Compute projections of the H tensor. Then assemble in 4-tensor
  // Note: Needs computation of the cross products with the levi-civita iterator

  // Preliminaries: Raise indices of the double derivative tensor projections
  tnsr::iJ<DataVector, 3> weyl_magnetic_down_up =
      make_with_value<tnsr::Ij<DataVector, 3>>(get<0, 0>(spacetime_metric),
                                               0.0);
  tnsr::I<DataVector, 3> nsDDKGu =
      make_with_value<tnsr::I<DataVector, 3>>(get<0, 0>(spacetime_metric), 0.0);
  tnsr::iJ<DataVector, 3> ssDDKGdu = make_with_value<tnsr::iJ<DataVector, 3>>(
      get<0, 0>(spacetime_metric), 0.0);
  tnsr::IJ<DataVector, 3> ssDDKGuu = make_with_value<tnsr::IJ<DataVector, 3>>(
      get<0, 0>(spacetime_metric), 0.0);

  tenex::evaluate<ti::i, ti::J>(
      weyl_magnetic_down_up,
      weyl_magnetic(ti::i, ti::l) * inverse_spatial_metric(ti::L, ti::J));
  tenex::evaluate<ti::I>(nsDDKGu,
                         inverse_spatial_metric(ti::I, ti::J) * nsDDKG(ti::j));
  tenex::evaluate<ti::i, ti::J>(
      ssDDKGdu, ssDDKG(ti::i, ti::l) * inverse_spatial_metric(ti::L, ti::J));
  tenex::evaluate<ti::I, ti::J>(
      ssDDKGuu, inverse_spatial_metric(ti::I, ti::L) * ssDDKGdu(ti::l, ti::J));

  // nn
  Scalar<DataVector> nnH =
      make_with_value<Scalar<DataVector>>(get<0, 0>(spacetime_metric), 0.0);
  // Raise indices of the spatial part of the second derivative of the scalar
  tenex::evaluate(nnH, weyl_electric(ti::i, ti::j) ssDDKGuu(ti::I, ti::J));

  // Cross products
  tensor::i<DataVector, 3> ssDDKGuu_cross_Bdu =
      make_with_value<tnsr::i<DataVector, 3>>(get<0, 0>(spacetime_metric), 0.0);
  tensor::ij<DataVector, 3> nsDDKGu_cross_Bdu =
      make_with_value<tnsr::ij<DataVector, 3>>(get<0, 0>(spacetime_metric),
                                               0.0);

  for (LeviCivitaIterator<3> levi_civita_it; levi_civita_it; ++levi_civita_it) {
    const auto [i, j, k] = levi_civita_it();
    // S cross B
    // epsilon_{ijk} B_{l}^{k} S^{jl}
    for (size_t l = 0; l < 3; ++l) {
      ssDDKGuu_cross_Bdu.get(i) += levi_civita_it.sign() * ssDDKGuu.get(j, l) *
                                   weyl_magnetic_down_up.get(l, k);
      // j cross B
      // Note: For now we don't impose symmetry of this quantity
      nsDDKGu_cross_Bdu.get(i, l) += levi_civita_it.sign() * nsDDKGu.get(j) *
                                     weyl_magnetic_down_up.get(l, k);
    }
  }

  // ns
  tnsr::i<DataVector, 3> nsH =
      make_with_value<tnsr::i<DataVector, 3>>(get<0, 0>(spacetime_metric), 0.0);

  tenex::evaluate<ti::i>(
      nsH, weyl_electric(ti::i, ti::j) nsDDKGu(ti::J) +
               // sqrt(gamma) * epsilon_{ijk} B_{l}^{k} S^{jl}
               sqrt_det_spatial_metric() * ssDDKGuu_cross_Bdu(ti::i));

  // ss
  tnsr::ii<DataVector, 3> ssH = make_with_value<tnsr::ii<DataVector, 3>>(
      get<0, 0>(spacetime_metric), 0.0);

  tenex::evaluate<ti::i, ti::j>(
      ssH, (trace_ssDDKG() + nnDDKG()) * weyl_electric(ti::i, ti::j)
               // -2 * 2 symmetric part E ssDDKG
               - (weyl_electric(ti::k, ti::i) * ssDDKGdu(ti::j, ti::K) +
                  weyl_electric(ti::k, ti::j) * ssDDKGdu(ti::i, ti::K))
               // +2 sqrt(gamma) symmetric part cross prod
               + sqrt_det_spatial_metric() * (nsDDKGu_cross_Bdu(ti::i, ti::j) +
                                              nsDDKGu_cross_Bdu(ti::j, ti::i))
               // Sum the trace part
               + spatial_metric(ti::i, ti::j) * nnH());

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
*/

}  // namespace ScalarTensor
