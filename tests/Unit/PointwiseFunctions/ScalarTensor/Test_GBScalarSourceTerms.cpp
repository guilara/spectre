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
#include "PointwiseFunctions/ScalarTensor/GBScalarSourceTerms.hpp"

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.ScalarTensor.GBScalarTerms",
                  "[Unit][PointwiseFunctions]") {
  TestHelpers::db::test_simple_tag<ScalarTensor::Tags::OrderReducedGBScalar>(
      "OrderReducedGBScalar");
  TestHelpers::db::test_compute_tag<
      ScalarTensor::Tags::OrderReducedGBScalarCompute<Frame::Inertial>>(
      "OrderReducedGBScalar");
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/ScalarTensor"};

  pypp::check_with_random_values<1>(
      &ScalarTensor::order_reduced_gb_scalar_with_tenex<Frame::Inertial>,
      "GBScalarSourceTerms", {"order_reduced_gb_scalar"}, {{{1.0e-2, 0.5}}},
      DataVector{5});
}
