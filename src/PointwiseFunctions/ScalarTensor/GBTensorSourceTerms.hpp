// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

void DDKG_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> DDKG_normal_normal_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,

    // Scalar quantities

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_pi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const Scalar<DataVector>& dt_pi_scalar);

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
    const tnsr::i<DataVector, 3>& dt_phi_scalar);

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
);

void DDKG_tensor_from_projections(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDKG_tensor_result,
    const Scalar<DataVector>& lapse, Scalar<DataVector> nnDDKG,
    tnsr::i<DataVector, 3> nsDDKG, tnsr::ii<DataVector, 3> ssDDKG);

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
    const tnsr::i<DataVector, 3>& dt_phi_scalar);

void order_reduced_gb_H_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> nnH_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const tnsr::ii<DataVector, 3>& ssDDKG);

void compute_S_cross_B(
    const gsl::not_null<tnsr::i<DataVector, 3>*> S_cross_B_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& ssDDKG);

void compute_j_cross_B(
    const gsl::not_null<tnsr::ij<DataVector, 3>*> j_cross_B_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_magnetic,
    const tnsr::i<DataVector, 3>& nsDDKG);

void order_reduced_gb_H_normal_spatial_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> nsH_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const tnsr::i<DataVector, 3>& nsDDKG,
    const tnsr::i<DataVector, 3>& S_cross_B);

void order_reduced_gb_H_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> ssH_result,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const Scalar<DataVector>& nnDDKG, const tnsr::ii<DataVector, 3>& ssDDKG,
    const tnsr::ij<DataVector, 3>& j_cross_B, const Scalar<DataVector>& nnH);

void order_reduced_gb_H_tensor_weyl_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const Scalar<DataVector>& lapse, const Scalar<DataVector>& nnH,
    const tnsr::i<DataVector, 3>& nsH, const tnsr::ii<DataVector, 3>& ssH);

void order_reduced_Q_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> Q_tensor_result,
    const tnsr::aa<DataVector, 3> spacetime_metric,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar);

/*
void order_reduced_gb_H_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> result,
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric const
        tnsr::aa<DataVector, Dim>& spacetime_metric,
    const tnsr::A<DataVector, 3>& normal_spacetime_vector,
    const tnsr::a<DataVector, 3>& normal_spacetime_one_form,
    const tnsr::aa<DataVector, 3>& dd_coupling_function);
*/

namespace Tags {

/*!
 * \brief Compute tag for normal-normal projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct nnDDKGCompute : nnDDKG, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3, Frame>,
      ::Tags::deriv<CurvedScalarWave::Tags::Pi, tmpl::size_t<3>, Frame>,
      ScalarTensor::Tags::RhsPi>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> result,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&,
      const Scalar<DataVector>&) = &DDKG_normal_normal_projection;
  using base = nnDDKG;
};

/*!
 * \brief Compute tag for normal-spatial projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct nsDDKGCompute : nsDDKG, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3, Frame>,
      gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
      gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame>,
      CurvedScalarWave::Tags::Phi<3>,
      ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>, Frame>,
      ScalarTensor::Tags::RhsPhi>;
  using return_type = tnsr::i<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> result,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>&,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::ij<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&) = &DDKG_normal_spatial_projection;
  using base = nsDDKG;
};

/*!
 * \brief Compute tag for normal-spatial projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct ssDDKGCompute : ssDDKG, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame>,
      gr::Tags::SpatialChristoffelSecondKind<DataVector, 3, Frame>,
      CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3>,
      ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>, Frame>>;
  using return_type = tnsr::ii<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ii<DataVector, 3>*> result,
      const tnsr::ii<DataVector, 3>&, const tnsr::Ijj<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const tnsr::ij<DataVector, 3>&) = &DDKG_spatial_spatial_projection;
  using base = ssDDKG;
};

/*!
 * \brief Compute tag for normal-spatial projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct DDKGTensorCompute : DDKGTensor, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>, ScalarTensor::Tags::nnDDKG,
                 ScalarTensor::Tags::nsDDKG, ScalarTensor::Tags::ssDDKG>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> result,
      const Scalar<DataVector>&, Scalar<DataVector>, tnsr::i<DataVector, 3>,
      tnsr::ii<DataVector, 3>) = &DDKG_tensor_from_projections;
  using base = DDKGTensor;
};

/*!
 * \brief Compute tag for normal-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducednnHCompute : OrderReducednnH, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 ScalarTensor::Tags::WeylElectricFull<DataVector, 3, Frame>,
                 ScalarTensor::Tags::ssDDKG>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> result,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) =
      &order_reduced_gb_H_normal_normal_projection;
  using base = OrderReducednnH;
};

/*!
 * \brief Compute S cross B.
 */
template <typename Frame>
struct SCrossBCompute : SCrossB, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::WeylMagnetic<DataVector, 3, Frame>,
                 ScalarTensor::Tags::ssDDKG>;
  using return_type = tnsr::i<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> result,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) = &compute_S_cross_B;
  using base = SCrossB;
};

/*!
 * \brief Compute S cross B.
 */
template <typename Frame>
struct JCrossBCompute : JCrossB, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::WeylMagnetic<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nsDDKG>;
  using return_type = tnsr::ij<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ij<DataVector, 3>*> result,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&) = &compute_j_cross_B;
  using base = JCrossB;
};

/*!
 * \brief Compute tag for normal-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducednsHCompute : OrderReducednsH, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::SqrtDetSpatialMetric<DataVector>,
                 ScalarTensor::Tags::WeylElectricFull<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nsDDKG, ScalarTensor::Tags::SCrossB>;
  using return_type = tnsr::i<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> result,
      const tnsr::II<DataVector, 3>&, const Scalar<DataVector>&,
      const tnsr::ii<DataVector, 3>&, const tnsr::i<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&) =
      &order_reduced_gb_H_normal_spatial_projection;
  using base = OrderReducednsH;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedssHCompute : OrderReducedssH, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::SqrtDetSpatialMetric<DataVector>,
                 ScalarTensor::Tags::WeylElectricFull<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nnDDKG, ScalarTensor::Tags::ssDDKG,
                 ScalarTensor::Tags::JCrossB,
                 ScalarTensor::Tags::OrderReducednnH>;
  using return_type = tnsr::ii<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ii<DataVector, 3>*> result,
      const tnsr::ii<DataVector, 3>& spatial_metric,
      const tnsr::II<DataVector, 3>& inverse_spatial_metric,
      const Scalar<DataVector>& sqrt_det_spatial_metric,
      const tnsr::ii<DataVector, 3>& weyl_electric,
      const Scalar<DataVector>& nnDDKG, const tnsr::ii<DataVector, 3>& ssDDKG,
      const tnsr::ij<DataVector, 3>& j_cross_B, const Scalar<DataVector>& nnH) =
      &order_reduced_gb_H_spatial_spatial_projection;
  using base = OrderReducedssH;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedHTensorCompute : OrderReducedHTensor, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, ScalarTensor::Tags::OrderReducednnH,
      ScalarTensor::Tags::OrderReducednsH, ScalarTensor::Tags::OrderReducedssH>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> result,
      const Scalar<DataVector>&, const Scalar<DataVector>&,
      const tnsr::i<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) = &order_reduced_gb_H_tensor_weyl_part;
  using base = OrderReducedHTensor;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedQTensorCompute : OrderReducedQTensor, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3, Frame>,
                 ScalarTensor::Tags::WeylElectricFullScalar<DataVector>,
                 gr::Tags::WeylMagneticScalar<DataVector>>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> Q_tensor_result,
      const tnsr::aa<DataVector, 3>, const Scalar<DataVector>&,
      const Scalar<DataVector>&) = &order_reduced_Q_tensor;
  using base = OrderReducedQTensor;
};

}  // namespace Tags

}  // namespace ScalarTensor
