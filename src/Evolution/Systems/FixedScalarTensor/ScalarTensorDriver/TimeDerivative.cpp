// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedDecoupledScalar/ScalarTensorDriver/TimeDerivative.hpp"

namespace fe::ScalarTensorDriver {
void TimeDerivative::apply(
    // GH dt variables
    const gsl::not_null<tnsr::aa<DataVector, dim>*> dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, dim>*> dt_pi,
    const gsl::not_null<tnsr::iaa<DataVector, dim>*> dt_phi,
    // Scalar dt variables
    const gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
    const gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
    const gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*>
        dt_phi_scalar,

    // GH temporal variables
    const gsl::not_null<Scalar<DataVector>*> temp_gamma1,
    const gsl::not_null<Scalar<DataVector>*> temp_gamma2,
    const gsl::not_null<tnsr::a<DataVector, dim>*> temp_gauge_function,
    const gsl::not_null<tnsr::ab<DataVector, dim>*>
        temp_spacetime_deriv_gauge_function,
    const gsl::not_null<Scalar<DataVector>*> gamma1gamma2,
    const gsl::not_null<Scalar<DataVector>*> half_half_pi_two_normals,
    const gsl::not_null<Scalar<DataVector>*> normal_dot_gauge_constraint,
    const gsl::not_null<Scalar<DataVector>*> gamma1_plus_1,
    const gsl::not_null<tnsr::a<DataVector, dim>*> pi_one_normal,
    const gsl::not_null<tnsr::a<DataVector, dim>*> gauge_constraint,
    const gsl::not_null<tnsr::i<DataVector, dim>*> half_phi_two_normals,
    const gsl::not_null<tnsr::aa<DataVector, dim>*>
        shift_dot_three_index_constraint,
    const gsl::not_null<tnsr::aa<DataVector, dim>*>
        mesh_velocity_dot_three_index_constraint,
    const gsl::not_null<tnsr::ia<DataVector, dim>*> phi_one_normal,
    const gsl::not_null<tnsr::aB<DataVector, dim>*> pi_2_up,
    const gsl::not_null<tnsr::iaa<DataVector, dim>*> three_index_constraint,
    const gsl::not_null<tnsr::Iaa<DataVector, dim>*> phi_1_up,
    const gsl::not_null<tnsr::iaB<DataVector, dim>*> phi_3_up,
    const gsl::not_null<tnsr::abC<DataVector, dim>*>
        christoffel_first_kind_3_up,
    const gsl::not_null<Scalar<DataVector>*> lapse,
    const gsl::not_null<tnsr::I<DataVector, dim>*> shift,
    const gsl::not_null<tnsr::II<DataVector, dim>*> inverse_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
    const gsl::not_null<tnsr::AA<DataVector, dim>*> inverse_spacetime_metric,
    const gsl::not_null<tnsr::abb<DataVector, dim>*> christoffel_first_kind,
    const gsl::not_null<tnsr::Abb<DataVector, dim>*> christoffel_second_kind,
    const gsl::not_null<tnsr::a<DataVector, dim>*> trace_christoffel,
    const gsl::not_null<tnsr::A<DataVector, dim>*> normal_spacetime_vector,

    // Scalar temporal variables
    const gsl::not_null<Scalar<DataVector>*> result_gamma1_scalar,
    const gsl::not_null<Scalar<DataVector>*> result_gamma2_scalar,

    // Extra temporal tags
    const gsl::not_null<tnsr::aa<DataVector, dim>*> stress_energy,

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

    const Scalar<DataVector>& scalar_source) {
  // Compute the sourceless part of the RHS of the tensor equation
  const size_t number_of_points = get<0, 0>(*dt_spacetime_metric).size();
  // Need constraint damping on interfaces in DG schemes
  *temp_gamma1 = gamma1;
  *temp_gamma2 = gamma2;

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

  *result_gamma1_scalar = gamma1_scalar;
  *result_gamma2_scalar = gamma2_scalar;

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
