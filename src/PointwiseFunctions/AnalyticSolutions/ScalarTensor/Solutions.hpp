// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Solutions.hpp"

namespace ScalarTensor {
/// Base struct for properties common to all Scalar Tensor analytic solutions
struct AnalyticSolution {
  static constexpr size_t volume_dim = 3_st;

  template <typename DataType>
  using tags = tmpl::push_back<
      typename gr::AnalyticSolution<3>::template tags<DataType>,
      // Add scalar variables here
      CurvedScalarWave::Tags::Psi,
      CurvedScalarWave::Tags::Pi,
      CurvedScalarWave::Tags::Phi<3_st>
      // We can add the gr tags required by CurvedScalarWave here
      // or in the tags for each individual solution
      // gr::Tags::DerivDetSpatialMetric<3_st, Frame, DataType>,
      // gr::Tags::TraceExtrinsicCurvature<DataType>,
      // gr::Tags::SpatialChristoffelFirstKind<3_st, Frame, DataType>,
      // gr::Tags::SpatialChristoffelSecondKind<3_st, Frame, DataType>,
      // gr::Tags::TraceSpatialChristoffelSecondKind<3_st, Frame, DataType>
      >;
};

/*!
 * \ingroup AnalyticSolutionsGroup
 * \brief Holds classes implementing a solution to the ScalarTensor system.
 */
namespace Solutions {}
}  // namespace ScalarTensor
