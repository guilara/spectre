// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"

namespace fe::ScalarDriver::Sources {

void add_scalar_driver_source_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_source,
    const Scalar<DataVector>& lapse) {
  dt_pi->get() += get(lapse) * scalar_driver_source.get();
}

void add_scalar_driver_friction_term_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_pi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st>& shift, const double scalar_tau_parameter,
    const double scalar_sigma_parameter) {
  dt_pi->get() += -1.0 * (scalar_tau_parameter / scalar_sigma_parameter) *
                  square(lapse.get()) * scalar_driver_pi.get();
}

void compute_scalar_driver_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const double scalar_tau_parameter, const double scalar_sigma_parameter) {
  // Make sure it has the same size
  *scalar_driver_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  scalar_driver_source->get() =
      (1.0 / scalar_sigma_parameter) * (psi.get() - target_psi.get());
}

void compute_target_psi(const gsl::not_null<Scalar<DataVector>*> target_psi,
                        const Scalar<DataVector>& psi) {
  // Make sure it has the same size
  *target_psi = make_with_value<Scalar<DataVector>>(psi, 0.);
  target_psi->get() = psi.get();
}

}  // namespace fe::ScalarDriver::Sources
