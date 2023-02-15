// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Sources/CurvatureSource.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Sources {

void compute_scalar_curvature_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar, const double mass_psi) {
  // We want to control the sign with the 'mass' parameter
  scalar_source->get() = mass_psi * weyl_electric_scalar.get();
}

}  // namespace CurvedScalarWave::Sources
