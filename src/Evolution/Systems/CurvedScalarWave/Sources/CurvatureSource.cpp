// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Sources/CurvatureSource.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Sources {

Scalar<DataVector> coupling_function_prime(const Scalar<DataVector>& psi,
                                           const double first_coupling_psi,
                                           const double second_coupling_psi) {
  Scalar<DataVector> result{psi.size(), 0.0};
  result.get() =
      first_coupling_psi * psi.get() + second_coupling_psi * square(psi.get());
  return result;
}

void compute_scalar_curvature_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& psi,
    const double first_coupling_psi, const double second_coupling_psi) {
  // We want to control the sign with the 'mass' parameter
  scalar_source->get() = weyl_electric_scalar.get();
  scalar_source->get() *= get(
      coupling_function_prime(psi, first_coupling_psi, second_coupling_psi));
}

}  // namespace CurvedScalarWave::Sources
