// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"

/// \cond
namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl
/// \endcond

namespace ScalarTensor {

template <typename DataType, size_t SpatialDim, typename Frame>
void weyl_electric_ricci(
    const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*>
        weyl_electric_ricci_part,
    const Scalar<DataType>& pi_scalar,
    const tnsr::i<DataType, SpatialDim, Frame>& phi_scalar,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric);

template <typename DataType, size_t SpatialDim, typename Frame>
void contract_electric_parts(
    const gsl::not_null<Scalar<DataType>*> result,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric_vacuum_part,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric_ricci_part,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric);

template <typename DataType, size_t SpatialDim, typename Frame>
void weyl_electric_full(
    const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*>
        weyl_electric_full,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<DataType, SpatialDim, Frame>& extrinsic_curvature,
    const Scalar<DataType>& pi_scalar,
    const tnsr::i<DataType, SpatialDim, Frame>& phi_scalar,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric);

template <typename DataType, size_t SpatialDim, typename Frame>
void weyl_electric_full(
    const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*>
        weyl_electric_full,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<DataType, SpatialDim, Frame>& extrinsic_curvature,
    const Scalar<DataType>& pi_scalar,
    const tnsr::i<DataType, SpatialDim, Frame>& phi_scalar,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_tensor_driver,
    const Scalar<DataType>& trace_tensor_driver);

template <typename DataType, size_t SpatialDim, typename Frame>
void my_trace(
    const gsl::not_null<Scalar<DataType>*> result,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric);

namespace Tags {
/// Compute item for the 4-Ricci part of the electric part of the weyl tensor
///
/// Can be retrieved using gr::Tags::WeylElectricRicci
template <typename DataType, size_t SpatialDim, typename Frame>
struct WeylElectricRicciCompute
    : WeylElectricRicci<DataType, SpatialDim, Frame>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<SpatialDim>,
                 gr::Tags::SpatialMetric<DataType, SpatialDim, Frame>,
                 gr::Tags::InverseSpatialMetric<DataType, SpatialDim, Frame>>;

  using return_type = tnsr::ii<DataType, SpatialDim, Frame>;

  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*> result,
      const Scalar<DataType>&, const tnsr::i<DataType, SpatialDim, Frame>&,
      const tnsr::ii<DataType, SpatialDim, Frame>&,
      const tnsr::II<DataType, SpatialDim, Frame>&)>(
      &weyl_electric_ricci<DataType, SpatialDim, Frame>);

  using base = WeylElectricRicci<DataType, SpatialDim, Frame>;
};

/// Can be retrieved using ScalarTensor::Tags::WeylElectricRicciScalar
template <typename DataType, size_t SpatialDim, typename Frame>
struct WeylElectricRicciScalarComplementCompute
    : WeylElectricRicciScalarComplement<DataType>,
      db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::WeylElectric<DataType, SpatialDim, Frame>,
      ScalarTensor::Tags::WeylElectricRicci<DataType, SpatialDim, Frame>,
      gr::Tags::InverseSpatialMetric<DataType, SpatialDim, Frame>>;

  using return_type = Scalar<DataType>;

  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<DataType>*>,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::II<DataType, SpatialDim, Frame>&)>(
          &ScalarTensor::contract_electric_parts<DataType, SpatialDim, Frame>);

  using base = WeylElectricRicciScalarComplement<DataType>;
};

/// Compute item for the 4-Ricci part of the electric part of the weyl tensor
///
/// Can be retrieved using gr::Tags::WeylElectricRicci
template <typename DataType, size_t SpatialDim, typename Frame>
struct WeylElectricFullCompute : WeylElectricFull<DataType, SpatialDim, Frame>,
                                 db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpatialRicci<DataType, SpatialDim, Frame>,
                 gr::Tags::ExtrinsicCurvature<DataType, SpatialDim, Frame>,
                 CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<SpatialDim>,
                 gr::Tags::SpatialMetric<DataType, SpatialDim, Frame>,
                 gr::Tags::InverseSpatialMetric<DataType, SpatialDim, Frame>,
                 fe::ScalarTensorDriver::Tags::TensorDriverSpatialProjection,
                 fe::ScalarTensorDriver::Tags::TensorDriverTrace>;

  using return_type = tnsr::ii<DataType, SpatialDim, Frame>;

  static constexpr auto function = static_cast<void (*)(
      const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*>
          weyl_electric_full,
      const tnsr::ii<DataType, SpatialDim, Frame>& spatial_ricci,
      const tnsr::ii<DataType, SpatialDim, Frame>& extrinsic_curvature,
      const Scalar<DataType>& pi_scalar,
      const tnsr::i<DataType, SpatialDim, Frame>& phi_scalar,
      const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
      const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
      const tnsr::ii<DataType, SpatialDim, Frame>& spatial_tensor_driver,
      const Scalar<DataType>& trace_tensor_driver)>(
      &weyl_electric_full<DataType, SpatialDim, Frame>);

  using base = WeylElectricFull<DataType, SpatialDim, Frame>;
};

/// Can be retrieved using ScalarTensor::Tags::WeylElectricFullScalar
template <typename DataType, size_t SpatialDim, typename Frame>
struct WeylElectricFullScalarCompute : WeylElectricFullScalar<DataType>,
                                       db::ComputeTag {
  using argument_tags = tmpl::list<
      ScalarTensor::Tags::WeylElectricFull<DataType, SpatialDim, Frame>,
      gr::Tags::InverseSpatialMetric<DataType, SpatialDim, Frame>>;

  using return_type = Scalar<DataType>;

  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<DataType>*>,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::II<DataType, SpatialDim, Frame>&)>(
          &gr::weyl_electric_scalar<DataType, SpatialDim, Frame>);

  using base = WeylElectricFullScalar<DataType>;
};

/// Can be retrieved using ScalarTensor::Tags::WeylElectricFullScalar
template <typename DataType, size_t SpatialDim, typename Frame>
struct WeylElectricFullTraceCompute : WeylElectricFullTrace<DataType>,
                                      db::ComputeTag {
  using argument_tags = tmpl::list<
      ScalarTensor::Tags::WeylElectricFull<DataType, SpatialDim, Frame>,
      gr::Tags::InverseSpatialMetric<DataType, SpatialDim, Frame>>;

  using return_type = Scalar<DataType>;

  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<DataType>*>,
                           const tnsr::ii<DataType, SpatialDim, Frame>&,
                           const tnsr::II<DataType, SpatialDim, Frame>&)>(
          &my_trace<DataType, SpatialDim, Frame>);

  using base = WeylElectricFullTrace<DataType>;
};

}  // namespace Tags
}  // namespace ScalarTensor
