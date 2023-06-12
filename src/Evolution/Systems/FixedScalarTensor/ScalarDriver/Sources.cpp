// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"

namespace fe::ScalarDriver::Sources {

void add_scalar_driver_source_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_source,
    const Scalar<DataVector>& lapse) {
  dt_pi->get() += get(lapse) * scalar_source.get();
}

void compute_scalar_driver_source(
    const gsl::not_null<return_type*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi) {
  // Make sure it has the same size
  *scalar_driver_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  scalar_driver_source->get() = psi.get() - target_psi.get();
}

}  // namespace fe::ScalarDriver::Sources
