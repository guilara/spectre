// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/WeylElectricRicciPart.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/VectorImpl.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace ScalarTensor {

template <typename DataType, size_t SpatialDim, typename Frame>
void weyl_electric_ricci(
    const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*>
        weyl_electric_ricci_part,
    const Scalar<DataType>& pi_scalar,
    const tnsr::i<DataType, SpatialDim, Frame>& phi_scalar,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric) {
  const double kappa = 8 * M_PI;
  const double one_over_six = 1.0 / 6.0;
  const auto scalar_factor =
      tenex::evaluate(2.0 * pi_scalar() * pi_scalar() +
                      phi_scalar(ti::k) * inverse_spatial_metric(ti::K, ti::L) *
                          phi_scalar(ti::l));
  tenex::evaluate<ti::i, ti::j>(
      weyl_electric_ricci_part,
      -0.5 * kappa * phi_scalar(ti::i) * phi_scalar(ti::j) -
          one_over_six * kappa * scalar_factor() *
              spatial_metric(ti::i, ti::j));
}

template <typename DataType, size_t SpatialDim, typename Frame>
void contract_electric_parts(
    const gsl::not_null<Scalar<DataType>*> result,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric_vacuum_part,
    const tnsr::ii<DataType, SpatialDim, Frame>& weyl_electric_ricci_part,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric) {
  tenex::evaluate(result, (2.0 * weyl_electric_vacuum_part(ti::i, ti::j) +
                           weyl_electric_ricci_part(ti::i, ti::j)) *
                              inverse_spatial_metric(ti::J, ti::K) *
                              weyl_electric_ricci_part(ti::k, ti::l) *
                              inverse_spatial_metric(ti::L, ti::I));
}

template <typename DataType, size_t SpatialDim, typename Frame>
void weyl_electric_full(
    const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*>
        weyl_electric_full,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_ricci,
    const tnsr::ii<DataType, SpatialDim, Frame>& extrinsic_curvature,
    const Scalar<DataType>& pi_scalar,
    const tnsr::i<DataType, SpatialDim, Frame>& phi_scalar,
    const tnsr::ii<DataType, SpatialDim, Frame>& spatial_metric,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric) {
  // Prefactor definitions

  // We work in units where set G = 1 / (8 M_PI)
  // const double kappa = 8.0 * M_PI;
  const double kappa = 1.0;

  const double one_over_six = 1.0 / 6.0;
  const double one_over_three = 1.0 / 3.0;

  // Compute the Weyl Electric scalar
  tenex::evaluate<ti::i, ti::j>(
      weyl_electric_full,
      // Weyl electric vacuum part
      spatial_ricci(ti::i, ti::j) +
          inverse_spatial_metric(ti::K, ti::L) *
              (extrinsic_curvature(ti::k, ti::l) *
                   extrinsic_curvature(ti::i, ti::j) -
               extrinsic_curvature(ti::i, ti::l) *
                   extrinsic_curvature(ti::k, ti::j))
          // Weyl electric Ricci part
          - 0.5 * kappa * phi_scalar(ti::i) * phi_scalar(ti::j) -
          one_over_six * kappa *
              (
                  // Scalar factor
                  2.0 * pi_scalar() * pi_scalar() +
                  phi_scalar(ti::k) * inverse_spatial_metric(ti::K, ti::L) *
                      phi_scalar(ti::l)) *
              spatial_metric(ti::i, ti::j)

  );

  // Remove the trace part
  tenex::update<ti::i, ti::j>(
      weyl_electric_full, (*weyl_electric_full)(ti::i, ti::j)
                              // Remove trace
                              - one_over_three *
                                    (
                                        // Trace
                                        inverse_spatial_metric(ti::K, ti::L) *
                                        (*weyl_electric_full)(ti::k, ti::l)

                                            ) *
                                    spatial_metric(ti::i, ti::j));
}

}  // namespace ScalarTensor

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(2, data)

#define INSTANTIATE(_, data)                                               \
  template void ScalarTensor::weyl_electric_ricci(                         \
      const gsl::not_null<tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>*>  \
          weyl_electric_ricci_part,                                        \
      const Scalar<DTYPE(data)>& pi_scalar,                                \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>& phi_scalar,      \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_metric, \
      const tnsr::II<DTYPE(data), DIM(data), FRAME(data)>&                 \
          inverse_spatial_metric);                                         \
  template void ScalarTensor::contract_electric_parts(                     \
      const gsl::not_null<Scalar<DTYPE(data)>*> result,                    \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>&                 \
          weyl_electric_vacuum_part,                                       \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>&                 \
          weyl_electric_ricci_part,                                        \
      const tnsr::II<DTYPE(data), DIM(data), FRAME(data)>&                 \
          inverse_spatial_metric);                                         \
  template void ScalarTensor::weyl_electric_full(                          \
      const gsl::not_null<tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>*>  \
          weyl_electric_full,                                              \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_ricci,  \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>&                 \
          extrinsic_curvature,                                             \
      const Scalar<DTYPE(data)>& pi_scalar,                                \
      const tnsr::i<DTYPE(data), DIM(data), FRAME(data)>& phi_scalar,      \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& spatial_metric, \
      const tnsr::II<DTYPE(data), DIM(data), FRAME(data)>&                 \
          inverse_spatial_metric);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (double, DataVector),
                        (Frame::Grid, Frame::Inertial))

#undef DIM
#undef DTYPE
#undef FRAME
#undef INSTANTIATE
