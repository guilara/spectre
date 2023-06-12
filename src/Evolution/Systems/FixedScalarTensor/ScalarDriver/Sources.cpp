// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"

namespace fe::ScalarDriver::Sources {

void compute_scalar_driver_source(const gsl::not_null<return_type*> result,
                                  const Scalar<DataVector>& psi,
                                  const Scalar<DataVector>& target_psi) {
  // Essentially psi - target_psi
}

}  // namespace fe::ScalarDriver::Sources
