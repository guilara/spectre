// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/ScalarTensor/TwoScalars/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup EvolutionSystemsGroup
 * \brief Items related to evolving the first-order scalar tensor system.
 */
namespace ScalarTensor::TwoScalars {
/*!
 * \brief Scalar Tensor system obtained from combining two copies of the
 * CurvedScalarWave and one copy of the gh system.
 */
struct System {
  static constexpr bool has_primitive_and_conservative_vars = false;
  static constexpr size_t volume_dim = 3;

  using gh_system = gh::System<3_st>;
  using scalar_system = CurvedScalarWave::System<3_st>;

  template <size_t ScalarLabel>
  using wrapped_scalar_variables = ::Tags::Variables<
      db::wrap_tags_in<ScalarTensor::TwoScalar::Tags::Csw,
                       typename scalar_system::variables_tag::tags_list,
                       tmpl::size_t<ScalarLabel>>>;

  using variables_tag = ::Tags::Variables<
      tmpl::append<typename gh_system::variables_tag::tags_list,
                   typename wrapped_scalar_variables<1>::tags_list,
                   typename wrapped_scalar_variables<2>::tags_list>>;

  using flux_variables = tmpl::append<typename gh_system::flux_variables,
                                      typename scalar_system::flux_variables,
                                      typename scalar_system::flux_variables>;

  using gradient_variables =
      tmpl::append<typename gh_system::gradient_variables,
                   typename wrapped_scalar_variables<1>::tags_list,
                   typename wrapped_scalar_variables<2>::tags_list>;
  using gradients_tags = gradient_variables;

  static constexpr bool is_in_flux_conservative_form = false;

  using inverse_spatial_metric_tag =
      typename gh_system::inverse_spatial_metric_tag;
};

}  // namespace ScalarTensor::TwoScalars
