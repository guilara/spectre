// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::sgb {

template <typename DataType, size_t SpatialDim, typename Frame>
void f_constraint(
    const gsl::not_null<tnsr::a<DataType, SpatialDim, Frame>*> constraint,
    const tnsr::a<DataType, SpatialDim, Frame>& gauge_function,
    const tnsr::ab<DataType, SpatialDim, Frame>& spacetime_d_gauge_function,
    const tnsr::a<DataType, SpatialDim, Frame>& spacetime_normal_one_form,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::AA<DataType, SpatialDim, Frame>& inverse_spacetime_metric,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& d_pi,
    const tnsr::ijaa<DataType, SpatialDim, Frame>& d_phi,
    const Scalar<DataType>& gamma2,
    const tnsr::iaa<DataType, SpatialDim, Frame>& three_index_constraint,
    const tnsr::aa<DataType, SpatialDim, Frame>& trace_reversed_stress_energy,
    const tnsr::aa<DataType, SpatialDim, Frame>& tensor_driver);

template <typename DataType, size_t SpatialDim, typename Frame>
void f_constraint_add_stress_energy_term(
    const gsl::not_null<tnsr::a<DataType, SpatialDim, Frame>*> constraint,
    const tnsr::AA<DataType, SpatialDim, Frame>& inverse_spacetime_metric,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::a<DataType, SpatialDim, Frame>& spacetime_normal_one_form,
    const tnsr::aa<DataType, SpatialDim, Frame>& trace_reversed_stress_energy,
    const tnsr::aa<DataType, SpatialDim, Frame>& tensor_driver);

namespace Tags {
/*!
 * \brief Compute item to get the F-constraint for the generalized harmonic
 * evolution system with a fixed scalar tensor source.
 *
 * \details See `gh::f_constraint()`. Can be retrieved using
 * `gh::Tags::FConstraint`.
 */
template <size_t SpatialDim, typename Frame>
struct FConstraintCompute
    : gh::Tags::FConstraint<DataVector, SpatialDim, Frame>,
      db::ComputeTag {
  using argument_tags = tmpl::list<
      gh::Tags::GaugeH<DataVector, SpatialDim, Frame>,
      gh::Tags::SpacetimeDerivGaugeH<DataVector, SpatialDim, Frame>,
      gr::Tags::SpacetimeNormalOneForm<DataVector, SpatialDim, Frame>,
      gr::Tags::SpacetimeNormalVector<DataVector, SpatialDim, Frame>,
      gr::Tags::InverseSpatialMetric<DataVector, SpatialDim, Frame>,
      gr::Tags::InverseSpacetimeMetric<DataVector, SpatialDim, Frame>,
      gh::Tags::Pi<DataVector, SpatialDim, Frame>,
      gh::Tags::Phi<DataVector, SpatialDim, Frame>,
      ::Tags::deriv<gh::Tags::Pi<DataVector, SpatialDim, Frame>,
                    tmpl::size_t<SpatialDim>, Frame>,
      ::Tags::deriv<gh::Tags::Phi<DataVector, SpatialDim, Frame>,
                    tmpl::size_t<SpatialDim>, Frame>,
      ::gh::ConstraintDamping::Tags::ConstraintGamma2,
      gh::Tags::ThreeIndexConstraint<DataVector, SpatialDim, Frame>,
      ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, SpatialDim,
                                                    Frame>,
      fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, SpatialDim,
                                                 Frame>>;

  using return_type = tnsr::a<DataVector, SpatialDim, Frame>;

  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<tnsr::a<DataVector, SpatialDim, Frame>*>,
      const tnsr::a<DataVector, SpatialDim, Frame>&,
      const tnsr::ab<DataVector, SpatialDim, Frame>&,
      const tnsr::a<DataVector, SpatialDim, Frame>&,
      const tnsr::A<DataVector, SpatialDim, Frame>&,
      const tnsr::II<DataVector, SpatialDim, Frame>&,
      const tnsr::AA<DataVector, SpatialDim, Frame>&,
      const tnsr::aa<DataVector, SpatialDim, Frame>&,
      const tnsr::iaa<DataVector, SpatialDim, Frame>&,
      const tnsr::iaa<DataVector, SpatialDim, Frame>&,
      const tnsr::ijaa<DataVector, SpatialDim, Frame>&,
      const Scalar<DataVector>&,
      const tnsr::iaa<DataVector, SpatialDim, Frame>&,
      const tnsr::aa<DataVector, SpatialDim, Frame>&,
      const tnsr::aa<DataVector, SpatialDim, Frame>&)>(
      &fe::sgb::f_constraint<DataVector, SpatialDim, Frame>);

  using base = gh::Tags::FConstraint<DataVector, SpatialDim, Frame>;
};
}  // namespace Tags
}  // namespace fe::sgb
