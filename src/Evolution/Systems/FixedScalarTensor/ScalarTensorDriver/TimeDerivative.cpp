// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/TimeDerivative.hpp"

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Sources.hpp"

namespace fe::ScalarTensorDriver {
void TimeDerivative::apply(
    // Tensor Driver dt variables
    gsl::not_null<tnsr::aa<DataVector, dim>*> dt_tensor_driver,
    gsl::not_null<tnsr::aa<DataVector, dim>*> dt_pi,

    // Scalar Driver dt variables
    gsl::not_null<Scalar<DataVector>*> dt_scalar_driver,
    gsl::not_null<Scalar<DataVector>*> dt_pi_scalar_driver,

    // Tensor Driver spatial derivatives
    const tnsr::iaa<DataVector, dim>& d_tensor_driver,
    const tnsr::iaa<DataVector, dim>& d_pi,

    // Scalar Driver spatial derivatives
    const tnsr::i<DataVector, dim>& d_scalar_driver,
    const tnsr::i<DataVector, dim>& d_pi_scalar_driver,

    // Tensor Driver argument variables
    const tnsr::aa<DataVector, dim>& tensor_driver,
    const tnsr::aa<DataVector, dim>& pi,

    const tnsr::aa<DataVector, dim, Frame::Inertial>& spacetime_metric,
    const Mesh<dim>& mesh, double time,
    const tnsr::I<DataVector, dim, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const std::optional<tnsr::I<DataVector, dim, Frame::Inertial>>&
        mesh_velocity,

    // Scalar Driver argument variables
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar_driver,

    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, dim>& shift,

    const tnsr::aa<DataVector, dim>& tensor_driver_source,
    const Scalar<DataVector>& scalar_driver_source,

    const Scalar<DataVector>& tau_parameter,
    const Scalar<DataVector>& sigma_parameter) {
  // Compute the sourceless part of the RHS of the tensor equation
  const size_t number_of_points = get<0, 0>(*dt_tensor_driver).size();

  // Tensor advection driver
  tenex::evaluate<ti::a, ti::b>(
      dt_tensor_driver,
      -lapse() * pi(ti::a, ti::b) +
          shift(ti::J) * d_tensor_driver(ti::j, ti::a, ti::b));

  tenex::evaluate<ti::a, ti::b>(dt_pi,
                                shift(ti::J) * d_pi(ti::j, ti::a, ti::b));

  // Add friction term
  fe::ScalarTensorDriver::Sources::add_tensor_driver_friction_term_to_dt_pi(
      dt_pi, pi, lapse, shift, tau_parameter, sigma_parameter);
  // Add source term
  fe::ScalarTensorDriver::Sources::add_tensor_driver_source_to_dt_pi(
      dt_pi, tensor_driver_source, lapse);

  // Scalar advection driver
  tenex::evaluate(dt_scalar_driver, -lapse() * pi_scalar_driver() +
                                        shift(ti::I) * d_scalar_driver(ti::i));

  tenex::evaluate(dt_pi_scalar_driver,
                  shift(ti::I) * d_pi_scalar_driver(ti::i));

  // Add friction term
  fe::ScalarDriver::Sources::add_scalar_driver_friction_term_to_dt_pi_scalar(
      dt_pi_scalar_driver, pi_scalar_driver, lapse, shift, tau_parameter,
      sigma_parameter);

  // Add source
  fe::ScalarDriver::Sources::add_scalar_driver_source_to_dt_pi_scalar(
      dt_pi_scalar_driver, scalar_driver_source, lapse);
}
}  // namespace fe::ScalarTensorDriver
