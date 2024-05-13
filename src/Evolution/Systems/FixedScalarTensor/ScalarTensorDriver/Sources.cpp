// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Sources.hpp"

#include "Utilities/Math.hpp"

namespace fe::ScalarTensorDriver::Sources {

void add_tensor_driver_source_to_dt_pi(
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_pi,
    const tnsr::aa<DataVector, 3>& tensor_driver_source,
    const Scalar<DataVector>& lapse) {
  tenex::update<ti::a, ti::b>(
      dt_pi,
      (*dt_pi)(ti::a, ti::b) + lapse() * tensor_driver_source(ti::a, ti::b));
}

void add_tensor_driver_friction_term_to_dt_pi(
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_pi,
    const tnsr::aa<DataVector, 3>& pi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st>& shift,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  const auto tau_over_sigma =
      tenex::evaluate(scalar_tau_parameter() / scalar_sigma_parameter());
  tenex::update<ti::a, ti::b>(
      dt_pi, (*dt_pi)(ti::a, ti::b) -
                 1.0 * tau_over_sigma() * lapse() * lapse() * pi(ti::a, ti::b));
}

void compute_tensor_driver_source(
    gsl::not_null<tnsr::aa<DataVector, 3>*> tensor_driver_source,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::aa<DataVector, 3>& target_tensor,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  tenex::evaluate<ti::a, ti::b>(
      tensor_driver_source,
      (1.0 / scalar_sigma_parameter()) *
          (tensor_driver(ti::a, ti::b) - target_tensor(ti::a, ti::b)));
}

void compute_target_tensor(
    gsl::not_null<tnsr::aa<DataVector, 3>*> target_tensor,
    const tnsr::aa<DataVector, 3>& tensor_driver) {
  tenex::evaluate<ti::a, ti::b>(target_tensor, tensor_driver(ti::a, ti::b));
}

}  // namespace fe::ScalarTensorDriver::Sources
