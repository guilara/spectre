// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/ScalarCharge.hpp"

template <typename Frame>
void ScalarTensor::scalar_charge_integrand(
    const gsl::not_null<Scalar<DataVector>*> result,
    const tnsr::i<DataVector, 3, Frame>& phi,
    const tnsr::I<DataVector, 3, Frame>& unit_normal_vector,
    const Scalar<DataVector>& area_element) {
  destructive_resize_components(result, get(area_element).size());
  for (auto& component : *result) {
    component = 0.0;
  }
  // Project the scalar gradient on the normal vector
  get(*result) = phi.get(i) * get<0>(unit_normal_vector);
  for (size_t j = 1; j < 3; ++j) {
    get(*result) += phi.get(j) * unit_normal_vector.get(j);
  }
  // Multiply by area element (Need this?)
  get(*result) *= area_element.get();
  // Multiply by integral prefactor
  get(*result) /= -4.0 * M_PI;
}
