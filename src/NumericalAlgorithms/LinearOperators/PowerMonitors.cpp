// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"

#include <cmath>

#include "DataStructures/ApplyMatrices.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/StripeIterator.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace PowerMonitors {

template <size_t Dim>
void compute_power_monitor(const gsl::not_null<Scalar<DataVector>*> result,
                           const Scalar<DataVector>& pi,
                           const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,
                           const Mesh<Dim>& mesh) {
  // Set result size
  destructive_resize_components(result, get_size(get(pi)));

  // Extract nodal coefficients
  DataVector test_data_vector = pi.get();

  std::array<DataVector, Dim> check_power_monitor_array_function =
      power_monitor_array(test_data_vector, mesh);

  size_t sliced_dim = 0;
  size_t n_stripe = mesh.extents(sliced_dim);
  for (size_t index = 0; index < n_stripe; ++index) {
    for (SliceIterator si(mesh.extents(), sliced_dim, index); si;
         ++si) {
      // Fill the slice with the value of the power_monitor_array
      get(*result)[si.volume_offset()] =
          check_power_monitor_array_function[sliced_dim].data()[index];
    }
  }

}  // compute_power_monitor

template <size_t Dim>
std::array<DataVector, Dim> power_monitor_array(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  // Result Data vectors
  std::array<DataVector, Dim> result;

  // Get modal coefficients
  const ModalVector modal_coefficients = to_modal_coefficients(
    input_data_vector, mesh);

  double slice_sum = 0.0;
  size_t n_slice = 0;
  size_t n_stripe = 0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {

    n_slice = mesh.extents().slice_away(sliced_dim).product();
    n_stripe = mesh.extents(sliced_dim);

    result[sliced_dim].destructive_resize(n_stripe);

    for (size_t index = 0; index < n_stripe; ++index) {
      slice_sum = 0.0;
      // We might want to move SliceIterator to the middle loop
      for (SliceIterator si(mesh.extents(), sliced_dim, index); si;
           ++si) {
        slice_sum += square(modal_coefficients[si.volume_offset()]);
      }
      slice_sum /= n_slice;
      slice_sum = sqrt(slice_sum);

      result[sliced_dim].data()[index] = slice_sum;

      if (index >= result[sliced_dim].size()) {
        Parallel::printf("Error in index of DataVector in array \n");
      }
    }
  }

  return result;

}

}  // namespace PowerMonitors

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                   \
  template void PowerMonitors::compute_power_monitor(                          \
      gsl::not_null<Scalar<DataVector>*> result, const Scalar<DataVector>& pi, \
      const tnsr::i<DataVector, DIM(data), Frame::Inertial>& phi,              \
      const Mesh<DIM(data)>& mesh);                                            \
  template std::array<DataVector, DIM(data)>                                   \
  PowerMonitors::power_monitor_array(const DataVector& input_data_vector,      \
                                     const Mesh<DIM(data)>& mesh);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
