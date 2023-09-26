// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "PointwiseFunctions/AnalyticData/GeneralRelativity/AnalyticData.hpp"

namespace fe::DecoupledScalar {
/// Base struct for properties common to all Scalar Tensor analytic initial data
struct AnalyticDataBase {
  static constexpr size_t volume_dim = 3_st;

  template <typename DataType>
  using tags = tmpl::push_back<
      typename gr::AnalyticDataBase<volume_dim>::template tags<DataType>,
      CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
      CurvedScalarWave::Tags::Phi<volume_dim>, fe::ScalarDriver::Tags::Psi,
      fe::ScalarDriver::Tags::Pi, fe::ScalarDriver::Tags::Phi<volume_dim>>;
};

/*!
 * \ingroup AnalyticDataGroup
 * \brief Holds classes implementing analytic data for the ScalarTensor
 * system.
 */
namespace AnalyticData {}
}  // namespace fe::DecoupledScalar
