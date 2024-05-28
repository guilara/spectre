// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Characteristics.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup EvolutionSystemsGroup
 * \brief Items related to the scalar tensor driver on a curved background
 */
namespace fe::ScalarTensorDriver {
/// \cond
struct TimeDerivative;
/// \endcond
struct System {
  static constexpr bool is_in_flux_conservative_form = false;
  static constexpr bool has_primitive_and_conservative_vars = false;
  static constexpr size_t volume_dim = 3_st;
  static constexpr bool is_euclidean = false;

  using boundary_conditions_base = BoundaryConditions::BoundaryCondition;
  using boundary_correction_base = BoundaryCorrections::BoundaryCorrection;

  using variables_tag = ::Tags::Variables<
      tmpl::list<Tags::TensorDriver<DataVector, volume_dim>,
                 Tags::Pi<DataVector, volume_dim>, Tags::Psi, Tags::PiScalar>>;
  using flux_variables = tmpl::list<>;
  using gradient_variables =
      tmpl::list<Tags::TensorDriver<DataVector, volume_dim>,
                 Tags::Pi<DataVector, volume_dim>, Tags::Psi, Tags::PiScalar>;

  // Relic alias: needs to be removed once all evolution systems
  // convert to using dg::ComputeTimeDerivative
  using gradients_tags = gradient_variables;

  using spacetime_tag_list = tmpl::list<>;

  using compute_volume_time_derivative_terms = TimeDerivative;

  using compute_largest_characteristic_speed =
      Tags::ComputeLargestCharacteristicSpeed<3, Frame::Inertial>;

  using inverse_spatial_metric_tag =
      gr::Tags::InverseSpatialMetric<DataVector, 3_st>;
};
}  // namespace fe::ScalarTensorDriver
