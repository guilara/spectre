// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/Spectral/PowerMonitors.hpp"

#include <cmath>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace PowerMonitors {

template <size_t Dim>
void compute_power_monitor(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,
    const Mesh<Dim>& mesh) {
  Parallel::printf("Inside compute monitor function \n");
  Parallel::printf("get_size(pi) = %u \n", get_size(get(pi)));
  destructive_resize_components(result, get_size(get(pi)));
  get(*result) += square(get(pi));
}

} // namespace PowerMonitors

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                   \
  template void PowerMonitors::compute_power_monitor(                          \
      gsl::not_null<Scalar<DataVector>*> result, const Scalar<DataVector>& pi, \
      const tnsr::i<DataVector, DIM(data), Frame::Inertial>& phi,              \
      const Mesh< DIM(data) >& mesh);                                          \

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
