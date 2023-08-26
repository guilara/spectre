// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/TwoScalars/StressEnergy.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor::TwoScalars {

void trace_reversed_stress_energy(
    const gsl::not_null<tnsr::aa<DataVector, 3_st>*> stress_energy,
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3_st>& phi_scalar,
    const Scalar<DataVector>& pi_scalar_2,
    const tnsr::i<DataVector, 3_st>& phi_scalar_2,
    const Scalar<DataVector>& lapse) {
  get<0, 0>(*stress_energy) =
      square(get(lapse) * get(pi_scalar)) + square(get(lapse) * get(pi_scalar));
  for (size_t i = 0; i < 3; ++i) {
    stress_energy->get(0, i + 1) =
        -get(lapse) * get(pi_scalar) * phi_scalar.get(i) -
        get(lapse) * get(pi_scalar_2) * phi_scalar_2.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      stress_energy->get(i + 1, j + 1) =
          phi_scalar.get(i) * phi_scalar.get(j) +
          phi_scalar_2.get(i) * phi_scalar_2.get(j);
    }
  }
}

}  // namespace ScalarTensor::TwoScalars
