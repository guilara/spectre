// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Solutions.hpp"

namespace ScalarTensor {
/// Base struct for properties common to all Scalar Tensor analytic solutions
struct AnalyticSolution {
  static constexpr size_t volume_dim = 3_st;

  template <typename DataType>
  using tags = tmpl::push_back<
      typename gr::AnalyticSolution<3>::template tags<DataType>>;
};

/*!
 * \ingroup AnalyticSolutionsGroup
 * \brief Holds classes implementing a solution to the ScalarTensor system.
 */
namespace Solutions {}
}  // namespace ScalarTensor
