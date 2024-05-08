// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedDecoupledScalar/ScalarTensorDriver/TimeDerivative.hpp"

namespace fe::ScalarTensorDriver {
void TimeDerivative::apply(
    // GH dt variables
    gsl::not_null<tnsr::aa<DataVector, dim>*> dt_spacetime_metric,
    gsl::not_null<tnsr::aa<DataVector, dim>*> dt_pi,
    gsl::not_null<tnsr::iaa<DataVector, dim>*> dt_phi,
    // Scalar dt variables
    gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
    gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
    gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*> dt_phi_scalar,

    // GH temporal variables
    gsl::not_null<Scalar<DataVector>*> lapse,
    gsl::not_null<tnsr::I<DataVector, dim>*> shift,
    gsl::not_null<tnsr::II<DataVector, dim>*> inverse_spatial_metric,
    gsl::not_null<Scalar<DataVector>*> det_spatial_metric,

    // Scalar temporal variables

    // Extra temporal tags
    gsl::not_null<tnsr::aa<DataVector, dim>*> stress_energy,

    // GH spatial derivatives
    const tnsr::iaa<DataVector, dim>& d_spacetime_metric,
    const tnsr::iaa<DataVector, dim>& d_pi,
    const tnsr::ijaa<DataVector, dim>& d_phi,

    // scalar spatial derivatives
    const tnsr::i<DataVector, dim>& d_psi_scalar,
    const tnsr::i<DataVector, dim>& d_pi_scalar,
    const tnsr::ij<DataVector, dim>& d_phi_scalar,

    // GH argument variables
    const tnsr::aa<DataVector, dim>& spacetime_metric,
    const tnsr::aa<DataVector, dim>& pi, const tnsr::iaa<DataVector, dim>& phi,
    const Scalar<DataVector>& gamma0, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2,
    const gh::gauges::GaugeCondition& gauge_condition, const Mesh<dim>& mesh,
    double time,
    const tnsr::I<DataVector, dim, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const std::optional<tnsr::I<DataVector, dim, Frame::Inertial>>&
        mesh_velocity,

    // Scalar argument variables
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, dim>& phi_scalar,
    const Scalar<DataVector>& lapse_scalar,
    const tnsr::I<DataVector, dim>& shift_scalar,
    const tnsr::i<DataVector, dim>& deriv_lapse,
    const tnsr::iJ<DataVector, dim>& deriv_shift,
    const tnsr::II<DataVector, dim>& upper_spatial_metric,
    const tnsr::I<DataVector, dim>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1_scalar,
    const Scalar<DataVector>& gamma2_scalar,

    const Scalar<DataVector>& scalar_source,

    const Scalar<DataVector>& tau_parameter,
    const Scalar<DataVector>& sigma_parameter) {
  // Compute the sourceless part of the RHS of the tensor equation
  const size_t number_of_points = get<0, 0>(*dt_spacetime_metric).size();
  // Need constraint damping on interfaces in DG schemes
  // *temp_gamma1 = gamma1;
  // *temp_gamma2 = gamma2;

  // Compute the spatial metric, determinant, inverse, lapse and shift
  const tnsr::ii<DataVector, Dim> spatial_metric{};
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      make_const_view(make_not_null(&spatial_metric.get(i, j)),
                      spacetime_metric.get(i + 1, j + 1), 0, number_of_points);
    }
  }
  determinant_and_inverse(det_spatial_metric, inverse_spatial_metric,
                          spatial_metric);
  gr::shift(shift, spacetime_metric, *inverse_spatial_metric);
  gr::lapse(lapse, *shift, spacetime_metric);

  // Tensor advection driver
  tenex::evaluate<ti::a, ti::b>(
      dt_spacetime_metric, -(*lapse)() * pi(ti::a, ti::b) +
                               (*shift)(ti::J)*d_spacetime_metric(ti::j, a, b));

  tenex::evaluate<ti::a, ti::b>(dt_pi, (*shift)(ti::J)*d_pi(ti::j, a, b));

  tenex::evaluate<ti::i, ti::a, ti::b>(dt_phi, 0.0 * phi(ti::i, ti::a, ti::b));

  // TODO: Add friction term
  fe::ScalarTensorDriver::Sources::add_tensor_driver_friction_term_to_dt_pi(
      dt_pi, pi, lapse, shift, scalar_tau_parameter, scalar_sigma_parameter);
  // TODO: Add source term
  fe::ScalarTensorDriver::Sources::add_tensor_driver_source_to_dt_pi(
      dt_pi, tensor_driver_source, lapse);

  // *result_gamma1_scalar = gamma1_scalar;
  // *result_gamma2_scalar = gamma2_scalar;

  // Scalar advection driver
  tenex::evaluate(dt_psi_scalar,
                  -lapse() * pi_scalar() + shift(ti::I) * d_psi_scalar(ti::i));

  tenex::evaluate(dt_pi_scalar, shift(ti::I) * d_pi_scalar(ti::i));

  for (size_t index = 0; index < 3_st; ++index) {
    dt_phi_scalar->get(index) = 0.0 * get(lapse) * phi_scalar.get(index);
  }

  // Add friction term
  fe::ScalarDriver::Sources::add_scalar_driver_friction_term_to_dt_pi_scalar(
      dt_pi_scalar, pi_scalar, lapse, shift, scalar_tau_parameter,
      scalar_sigma_parameter);

  // Add source
  fe::ScalarDriver::Sources::add_scalar_driver_source_to_dt_pi_scalar(
      dt_pi_scalar, scalar_driver_source, lapse);
}
}  // namespace fe::ScalarTensorDriver
