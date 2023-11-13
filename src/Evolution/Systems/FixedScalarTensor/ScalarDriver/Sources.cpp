// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"

#include "Utilities/Math.hpp"

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
    const tnsr::I<DataVector, 3_st>& shift,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  dt_pi->get() += -1.0 *
                  (scalar_tau_parameter.get() / scalar_sigma_parameter.get()) *
                  square(lapse.get()) * scalar_driver_pi.get();
}

void compute_scalar_driver_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  // Make sure it has the same size
//   *scalar_driver_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  scalar_driver_source->get() =
      (1.0 / scalar_sigma_parameter.get()) * (psi.get() - target_psi.get());
}

void compute_scalar_driver_source_with_limiter(
    const gsl::not_null<Scalar<DataVector>*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter,
    const double limiter_parameter) {
  // Make sure it has the same size
//   *scalar_driver_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  // Compute scalar driver source
  scalar_driver_source->get() =
      (1.0 / scalar_sigma_parameter.get()) * (psi.get() - target_psi.get());
  // Apply limiter on the magnitude of the source
  double source_value = 0.0;
  double psi_magnitude = 0.0;
  double target_magnitude = 0.0;
  for (size_t index = 0; index < get(psi).size(); ++index) {
    source_value = gsl::at(get(*scalar_driver_source), index);
    psi_magnitude = std::abs(gsl::at(get(psi), index));
    target_magnitude = std::abs(gsl::at(get(target_psi), index));
    // Keep the sign, but limit the magnitude
    get(*scalar_driver_source)[index] =
        sgn(source_value) *
        std::min(std::abs(source_value),
                 limiter_parameter * std::max(psi_magnitude, target_magnitude));
  }
}

void compute_target_psi(const gsl::not_null<Scalar<DataVector>*> target_psi,
                        const Scalar<DataVector>& psi) {
  // Make sure it has the same size
//   *target_psi = make_with_value<Scalar<DataVector>>(psi, 0.);
  target_psi->get() = psi.get();
}

// For the exponential driver
void add_scalar_driver_source_to_dt_psi_scalar_for_exponential_driver(
    gsl::not_null<Scalar<DataVector>*> dt_psi,
    const Scalar<DataVector>& scalar_driver_source) {
  dt_psi->get() += scalar_driver_source.get();
}

void compute_scalar_driver_source_for_exponential_driver(
    const gsl::not_null<Scalar<DataVector>*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  // Make sure it has the same size
  // *scalar_driver_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  scalar_driver_source->get() =
      (-1.0 / scalar_tau_parameter.get()) * (psi.get() - target_psi.get());
}

}  // namespace fe::ScalarDriver::Sources
