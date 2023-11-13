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
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  // Make sure it has the same size
//   *diagnostic = make_with_value<Scalar<DataVector>>(psi, 0.);
  diagnostic->get() = psi.get() - target_psi.get();
}

void shift_minus_mesh_velocity(
    const gsl::not_null<tnsr::I<DataVector, 3>*> result,
    const tnsr::I<DataVector, 3>& shift,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity) {
  if (mesh_velocity.has_value()) {
    tenex::evaluate<ti::I>(result, shift(ti::I) - mesh_velocity.value()(ti::I));
  } else {
    tenex::evaluate<ti::I>(result, shift(ti::I));
  }
}

}  // namespace fe::ScalarDriver
