// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/BoundaryCorrection.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup ScalarTensorGroup
 * \brief Items related to scalar fields with backreaction in the metric.
 */
namespace ScalarTensor {
/// \cond
struct TimeDerivativeTerms;
/// \endcond

struct System {
  using boundary_conditions_base = BoundaryConditions::BoundaryCondition;
  using boundary_correction_base = BoundaryCorrections::BoundaryCorrection;
  static constexpr bool has_primitive_and_conservative_vars = true;
  static constexpr size_t volume_dim = 3;

  using gh_system = GeneralizedHarmonic::System<3_st>;
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

  // Note: Allocate memory for these using GrTagsForHydro directly
  // in the executable.
  // Need to compute the updates of the missing variables either in
  // TimeDerivativeTerms or ComputeTag
  using spacetime_variables_tag =
      ::Tags::Variables<scalar_system::spacetime_tag_list>;

  using compute_volume_time_derivative_terms = TimeDerivativeTerms;
};

} // namespace ScalarTensor
