// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/ScalarTensor/Characteristics.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/TimeDerivativeTerms.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup EvolutionSystemsGroup
 * \brief Items related to evolving the first-order scalar tensor system.
 */
namespace ScalarTensor {
struct System {
  using boundary_conditions_base = BoundaryConditions::BoundaryCondition;
  using boundary_correction_base = BoundaryCorrections::BoundaryCorrection;
  static constexpr bool has_primitive_and_conservative_vars = false;
  static constexpr size_t volume_dim = 3;

  using gh_system = gh::System<3_st>;
  using scalar_system = CurvedScalarWave::System<3_st>;

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

}  // namespace ScalarTensor
