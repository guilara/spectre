// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/StressEnergy.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

void add_stress_energy_term_to_dt_pi(
    const gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_pi,
    const tnsr::aa<DataVector, 3_st>& trace_reversed_stress_energy,
    const Scalar<DataVector>& lapse) {
  // We move the M_PI factor to the trace reversed computation
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      dt_pi->get(a, b) -=
          2.0 * get(lapse) * trace_reversed_stress_energy.get(a, b);
    }
  }
}

void trace_reversed_stress_energy(
    const gsl::not_null<tnsr::aa<DataVector, 3_st>*> stress_energy,
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3_st>& phi_scalar,
    const Scalar<DataVector>& lapse) {
  // We work in units where set G = 1 / (8 M_PI)
  // const double kappa = 8.0 * M_PI;
  const double kappa = 1.0;
  get<0, 0>(*stress_energy) = kappa * square(get(lapse) * get(pi_scalar));
  for (size_t i = 0; i < 3; ++i) {
    stress_energy->get(0, i + 1) =
        -kappa * get(lapse) * get(pi_scalar) * phi_scalar.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      stress_energy->get(i + 1, j + 1) =
          kappa * phi_scalar.get(i) * phi_scalar.get(j);
    }
  }
}

void trace_of_trace_reversed_stress_energy(
    const gsl::not_null<Scalar<DataVector>*> trace_of_stress_energy,
    const tnsr::aa<DataVector, 3>& stress_energy,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric) {
  tenex::evaluate(
      trace_of_stress_energy,
      stress_energy(ti::a, ti::b) * inverse_spacetime_metric(ti::B, ti::A));
}

}  // namespace ScalarTensor
