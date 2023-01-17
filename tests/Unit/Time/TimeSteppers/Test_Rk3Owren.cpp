// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Time/TimeSteppers/TimeStepperTestUtils.hpp"
#include "Time/TimeSteppers/Rk3Owren.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"

SPECTRE_TEST_CASE("Unit.Time.TimeSteppers.Rk3Owren", "[Unit][Time]") {
  const TimeSteppers::Rk3Owren stepper{};
  TimeStepperTestUtils::check_substep_properties(stepper);
  TimeStepperTestUtils::integrate_test(stepper, 3, 0, 1.0, 1.0e-9);
  TimeStepperTestUtils::integrate_test(stepper, 3, 0, -1.0, 1.0e-9);
  TimeStepperTestUtils::integrate_test_explicit_time_dependence(stepper, 3, 0,
                                                                -1.0, 1.0e-9);
  TimeStepperTestUtils::integrate_error_test(stepper, 3, 0, 1.0, 1.0e-6, 100,
                                             1.0e-3);
  TimeStepperTestUtils::integrate_error_test(stepper, 3, 0, -1.0, 1.0e-6, 100,
                                             1.0e-3);
  TimeStepperTestUtils::integrate_variable_test(stepper, 3, 0, 1.0e-9);
  TimeStepperTestUtils::stability_test(stepper);
  TimeStepperTestUtils::check_convergence_order(stepper);
  TimeStepperTestUtils::check_dense_output(stepper, 3_st);

  TestHelpers::test_factory_creation<TimeStepper, TimeSteppers::Rk3Owren>(
      "Rk3Owren");
  test_serialization(stepper);
  test_serialization_via_base<TimeStepper, TimeSteppers::Rk3Owren>();
  // test operator !=
  CHECK_FALSE(stepper != stepper);
}
