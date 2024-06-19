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
}
