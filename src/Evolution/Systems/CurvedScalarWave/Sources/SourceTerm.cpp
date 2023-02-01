// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Sources/SourceTerm.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Sources {

void compute_scalar_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi) {
  // Mass term
  scalar_source->get() = /* factor times */ psi.get();
}

void add_scalar_source_to_dt_pi(const gsl::not_null<Scalar<DataVector>*> dt_pi,
                                const Scalar<DataVector>& scalar_source,
                                const Scalar<DataVector>& lapse) {
  // Source
  dt_pi->get() += /* factor times */ get(lapse) * scalar_source.get();
}

}  // namespace CurvedScalarWave::Sources
