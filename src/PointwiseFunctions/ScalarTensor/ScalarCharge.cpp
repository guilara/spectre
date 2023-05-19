// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/ScalarCharge.hpp"

#include <cmath>

#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

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
  get(*result) = phi.get(0) * get<0>(unit_normal_vector);
  for (size_t j = 1; j < 3; ++j) {
    get(*result) += phi.get(j) * unit_normal_vector.get(j);
  }
  // Multiply by integral prefactor
  get(*result) /= -4.0 * M_PI;
}

#define FRAME(data) BOOST_PP_TUPLE_ELEM(0, data)
#define INSTANTIATE(_, data)                                         \
  template void ScalarTensor::scalar_charge_integrand<FRAME(data)>(  \
      const gsl::not_null<Scalar<DataVector>*> result,               \
      const tnsr::i<DataVector, 3, FRAME(data)>& phi,                \
      const tnsr::I<DataVector, 3, FRAME(data)>& unit_normal_vector, \
      const Scalar<DataVector>& area_element);
GENERATE_INSTANTIATIONS(INSTANTIATE, (Frame::Inertial))
#undef INSTANTIATE
#undef FRAME
