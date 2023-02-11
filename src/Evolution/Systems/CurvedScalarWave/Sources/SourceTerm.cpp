// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Sources/SourceTerm.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Sources {

void compute_scalar_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi, const double mass_psi) {
  scalar_source->get() = square(mass_psi) * psi.get();
}

void compute_scalar_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi) {
  // Mass
  const double mass_psi = 1.0;
  scalar_source->get() = square(mass_psi) * psi.get();
}

void add_scalar_source_to_dt_pi(const gsl::not_null<Scalar<DataVector>*> dt_pi,
                                const Scalar<DataVector>& scalar_source,
                                const Scalar<DataVector>& lapse) {
  dt_pi->get() += get(lapse) * scalar_source.get();
}

}  // namespace CurvedScalarWave::Sources
