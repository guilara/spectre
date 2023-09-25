// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Diagnostics.hpp"

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarDriver {

void driver_tracking_diagnostic(
    const gsl::not_null<Scalar<DataVector>*> diagnostic,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const double scalar_tau_parameter, const double scalar_sigma_parameter) {
  // Make sure it has the same size
//   *diagnostic = make_with_value<Scalar<DataVector>>(psi, 0.);
  diagnostic->get() = psi.get() - target_psi.get();
}

}  // namespace fe::ScalarDriver
