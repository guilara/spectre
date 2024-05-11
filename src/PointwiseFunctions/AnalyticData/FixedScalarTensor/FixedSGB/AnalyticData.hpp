// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "PointwiseFunctions/AnalyticData/GeneralRelativity/AnalyticData.hpp"

namespace fe::sgb {
/// Base struct for properties common to all Scalar Tensor analytic initial data
struct AnalyticDataBase {
  static constexpr size_t volume_dim = 3_st;

  template <typename DataType>
  using tags = tmpl::push_back<
      typename gr::AnalyticDataBase<volume_dim>::template tags<DataType>,
      CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
      CurvedScalarWave::Tags::Phi<volume_dim>,
      fe::ScalarTensorDriver::Tags::Psi, fe::ScalarTensorDriver::Tags::PiScalar,
      fe::ScalarTensorDriver::Tags::TensorDriver<DataType, volume_dim,
                                                 Frame::Inertial>,
      fe::ScalarTensorDriver::Tags::Pi<DataType, volume_dim, Frame::Inertial>>;
};

/*!
 * \ingroup AnalyticDataGroup
 * \brief Holds classes implementing analytic data for the ScalarTensor
 * system.
 */
namespace AnalyticData {}
}  // namespace fe::sgb
