// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Diagnostics.hpp"

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarTensorDriver {

void tensor_driver_tracking_diagnostic(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> diagnostic,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::aa<DataVector, 3>& target_tensor,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter) {
  tenex::evaluate<ti::a, ti::b>(
      diagnostic, tensor_driver(ti::a, ti::b) - target_tensor(ti::a, ti::b));
}

}  // namespace fe::ScalarTensorDriver
