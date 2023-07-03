// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/Solutions.hpp"

namespace fe::DecoupledScalar {
/// Base struct for properties common to all FixedDecoupledScalar analytic
/// solutions
struct AnalyticSolution {
  static constexpr size_t volume_dim = 3_st;

  template <typename DataType>
  using tags = tmpl::push_back<
      typename ScalarTensor::AnalyticSolution::template tags<DataType>,
      // Add scalar driver variables here
      fe::ScalarDriver::Tags::Psi, fe::ScalarDriver::Tags::Pi,
      fe::ScalarDriver::Tags::Phi<3_st> >;
};

/*!
 * \ingroup AnalyticSolutionsGroup
 * \brief Holds classes implementing a solution to the ScalarTensor system.
 */
namespace Solutions {}
}  // namespace fe::DecoupledScalar
