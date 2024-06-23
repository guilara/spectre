// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Helpers/Domain/DomainTestHelpers.hpp"
#include "Helpers/PointwiseFunctions/GeneralRelativity/TestHelpers.hpp"
#include "PointwiseFunctions/ScalarTensor/GBTensorSourceTerms.hpp"

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.ScalarTensor.GBTensorTerms",
                  "[Unit][PointwiseFunctions]") {
  TestHelpers::db::test_simple_tag<ScalarTensor::Tags::SpacetimeDerivScalar>(
      "SpacetimeDerivScalar");
  TestHelpers::db::test_compute_tag<
      ScalarTensor::Tags::SpacetimeDerivScalarCompute<Frame::Inertial>>(
      "SpacetimeDerivScalar");
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/ScalarTensor"};

  pypp::check_with_random_values<1>(
      &ScalarTensor::spacetime_derivative_scalar, "GBTensorSourceTerms",
      {"spacetime_derivative_scalar"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDKG_normal_normal_projection, "GBTensorSourceTerms",
      {"DDKG_normal_normal_projection"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDKG_normal_spatial_projection, "GBTensorSourceTerms",
      {"DDKG_normal_spatial_projection"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDKG_spatial_spatial_projection, "GBTensorSourceTerms",
      {"DDKG_spatial_spatial_projection"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDKG_tensor_from_projections, "GBTensorSourceTerms",
      {"DDKG_tensor_from_projections"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDFPsi_tensor_from_DDKG_tensor, "GBTensorSourceTerms",
      {"DDFPsi_tensor_from_DDKG_tensor"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDFPsi_normal_normal_projection, "GBTensorSourceTerms",
      {"DDFPsi_normal_normal_projection"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDFPsi_spatial_normal_projection, "GBTensorSourceTerms",
      {"DDFPsi_spatial_normal_projection"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::DDFPsi_spatial_spatial_projection, "GBTensorSourceTerms",
      {"DDFPsi_spatial_spatial_projection"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_gb_H_normal_normal_projection,
      "GBTensorSourceTerms", {"order_reduced_gb_H_normal_normal_projection"},
      {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::compute_S_cross_B, "GBTensorSourceTerms",
      {"compute_S_cross_B"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::compute_j_cross_B, "GBTensorSourceTerms",
      {"compute_j_cross_B"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_gb_H_normal_spatial_projection,
      "GBTensorSourceTerms", {"order_reduced_gb_H_normal_spatial_projection"},
      {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_gb_H_spatial_spatial_projection,
      "GBTensorSourceTerms", {"order_reduced_gb_H_spatial_spatial_projection"},
      {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_gb_H_tensor_weyl_part, "GBTensorSourceTerms",
      {"order_reduced_gb_H_tensor_weyl_part"}, {{{1.0e-2, 0.5}}},
      DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_Q_tensor, "GBTensorSourceTerms",
      {"order_reduced_Q_tensor"}, {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_gb_H_tensor_ricci_part,
      "GBTensorSourceTerms", {"order_reduced_gb_H_tensor_ricci_part"},
      {{{1.0e-2, 0.5}}}, DataVector{5});
  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_trace_reversed_stress_energy,
      "GBTensorSourceTerms", {"order_reduced_trace_reversed_stress_energy"},
      {{{1.0e-2, 0.5}}}, DataVector{5});
}
