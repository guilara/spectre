// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/Characteristics.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/TimeDerivativeTerms.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/System.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::DecoupledScalar {
struct System {
  using boundary_conditions_base = BoundaryConditions::BoundaryCondition;
  using boundary_correction_base = BoundaryCorrections::BoundaryCorrection;
  static constexpr bool has_primitive_and_conservative_vars = false;
  static constexpr size_t volume_dim = 3;

  using gh_system = ScalarTensor::System;
  using scalar_system = fe::ScalarDriver::System;

  using variables_tag = ::Tags::Variables<
      tmpl::append<typename gh_system::variables_tag::tags_list,
                   typename scalar_system::variables_tag::tags_list>>;

  using flux_variables = tmpl::append<typename gh_system::flux_variables,
                                      typename scalar_system::flux_variables>;

  using gradient_variables =
      tmpl::append<typename gh_system::gradient_variables,
                   typename scalar_system::gradient_variables>;
  using gradients_tags = gradient_variables;

  static constexpr bool is_in_flux_conservative_form = false;

  using compute_volume_time_derivative_terms = TimeDerivativeTerms;

  using compute_largest_characteristic_speed =
      Tags::ComputeLargestCharacteristicSpeed<>;

  using inverse_spatial_metric_tag =
      typename gh_system::inverse_spatial_metric_tag;
};

}  // namespace fe::DecoupledScalar
