// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor::Sources {

void add_scalar_source_to_dt_pi_scalar(gsl::not_null<Scalar<DataVector>*> dt_pi,
                                const Scalar<DataVector>& scalar_source,
                                const Scalar<DataVector>& lapse) {
  dt_pi->get() += get(lapse) * scalar_source.get();
}

}  // namespace ScalarTensor::Sources
