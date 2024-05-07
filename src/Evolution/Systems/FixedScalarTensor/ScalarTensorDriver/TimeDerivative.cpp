// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedDecoupledScalar/ScalarTensorDriver/TimeDerivative.hpp"

namespace fe::ScalarTensorDriver {
void TimeDerivative::apply() {
  // Use the definition from the CurvedScalarWave system
  CurvedScalarWave::TimeDerivative<3_st>::apply(
      dt_psi, dt_pi, dt_phi,

      result_lapse, result_shift, result_inverse_spatial_metric, result_gamma1,
      result_gamma2,

      d_psi, d_pi, d_phi, pi, phi, lapse, shift, deriv_lapse, deriv_shift,
      upper_spatial_metric, trace_spatial_christoffel,
      trace_extrinsic_curvature, gamma1, gamma2);

  // Add extra terms to the Klein-Gordon equation
  // Make sure all variables called here are in the arguments of apply
  // and in the DataBox
  fe::ScalarDriver::Sources::add_scalar_driver_friction_term_to_dt_pi_scalar(
      dt_pi, pi, lapse, shift, scalar_tau_parameter, scalar_sigma_parameter);

  fe::ScalarDriver::Sources::add_scalar_driver_source_to_dt_pi_scalar(
      dt_pi, scalar_driver_source, lapse);

  // For the tensor driver we want to write a wave equation as well
}
}  // namespace fe::ScalarTensorDriver
