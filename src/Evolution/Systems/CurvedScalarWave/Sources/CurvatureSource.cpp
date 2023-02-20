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
  auto result = make_with_value<Scalar<DataVector>>(psi, 1.);
  result.get() *= first_coupling_psi;
  result.get() += second_coupling_psi * psi.get();
  return result;
}

void compute_scalar_curvature_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi) {
  // Make sure it has the same size
  *scalar_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  // Compute the Riemann squared scalar in vacuum
  scalar_source->get() =
      8.0 * (weyl_electric_scalar.get() - weyl_magnetic_scalar.get());
  // Multiply by the source coupling function
  scalar_source->get() *= get(
      coupling_function_prime(psi, first_coupling_psi, second_coupling_psi));
}

}  // namespace CurvedScalarWave::Sources
