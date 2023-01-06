// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/Spectral/PowerMonitors.hpp"

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
void compute_power_monitor(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,
    const Mesh<Dim>& mesh) {

  // Set result size
  destructive_resize_components(result, get_size(get(pi)));

  std::array<DataVector, Dim> check_power_monitor_array_function
      = power_monitor_array(pi.get(), mesh);

  for (size_t check_dim = 0; check_dim < Dim; ++check_dim) {
    Parallel::printf("Check entry of power_monitor_array_function: %lf \n",
                     check_power_monitor_array_function[check_dim].data()[0]);
  }

  } // compute_power_monitor

// New function
template <size_t Dim>
std::array<DataVector, Dim> power_monitor_array(
    const DataVector& input_data_vector,
    const Mesh<Dim>& mesh) {

  // Result Data vectors
  std::array<DataVector, Dim> array_of_power_monitors;
  Parallel::printf("array first entry size = %u \n",
                   array_of_power_monitors[0].size());

  Parallel::printf("Inside power_monitor_array function \n");

  size_t my_number_of_grid_points = mesh.number_of_grid_points();
  Parallel::printf("number of (total) gridpoints: %u \n",
                   my_number_of_grid_points);

  // Get modal coefficients
  const ModalVector mod_coeffs =
      to_modal_coefficients(input_data_vector, mesh);
  Parallel::printf("Size of mod_coeffs Modal Vector: %u \n", mod_coeffs.size());

  // <<<<<<<<<<<<<<<<<<<
  // Here we compute the marginalized coefficients for each dimension
  // and repeat for all dims
  // <<<<<<<<<<<<<<<<<<<

  double my_slice_sum{0.0};
  size_t num_elems_slice;
  size_t num_elems_stripe;
  size_t counter = 0;
  size_t data_per_dim_counter = 0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {
    data_per_dim_counter = 0;
    Parallel::printf("Computing in dim =  %u\n", sliced_dim);
    num_elems_slice = mesh.extents().slice_away(sliced_dim).product();
    num_elems_stripe = mesh.extents(sliced_dim);
    Parallel::printf("num_elems_slice =  %u\n", num_elems_slice);
    Parallel::printf("num_elems_stripe =  %u\n", num_elems_stripe);

    array_of_power_monitors[sliced_dim].destructive_resize(num_elems_stripe);
    Parallel::printf("array first entry size = %u \n",
                     array_of_power_monitors[sliced_dim].size());

    for (size_t index = 0; index < num_elems_stripe; ++index) {
      Parallel::printf("Summing in slice with index %u\n", index);
      my_slice_sum = 0.0;
      counter = 0;
      // We might want to move SliceIterator to the middle loop
      // Check Nils D. code
      for (SliceIterator si(mesh.extents(), sliced_dim, index); si;
           ++si, ++counter) {
        // Check index structure of mod_coeffs
        my_slice_sum += square(mod_coeffs[si.volume_offset()]);
      }
      my_slice_sum /= num_elems_slice;
      my_slice_sum = sqrt(my_slice_sum);
      ++data_per_dim_counter;
      Parallel::printf("Final my_slice_sum (avg) = %lf \n", my_slice_sum);

      array_of_power_monitors[sliced_dim].data()[index] = my_slice_sum;
      Parallel::printf("Check Data Vector array entry filled: %lf \n",
                       array_of_power_monitors[sliced_dim].data()[index]);
      if (index >= array_of_power_monitors[sliced_dim].size()) {
        Parallel::printf("Error in index of DataVector in array \n");
      }
    }
    Parallel::printf("data in this dimension = %u \n", data_per_dim_counter);
  }

  return array_of_power_monitors;

} // power_monitor_array

} // namespace PowerMonitors

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                   \
  template void PowerMonitors::compute_power_monitor(                          \
      gsl::not_null<Scalar<DataVector>*> result, const Scalar<DataVector>& pi, \
      const tnsr::i<DataVector, DIM(data), Frame::Inertial>& phi,              \
      const Mesh<DIM(data)>& mesh);                                            \
  template std::array<DataVector, DIM(data)>                                   \
    PowerMonitors::power_monitor_array(                                        \
      const DataVector& input_data_vector,                                     \
      const Mesh<DIM(data)>& mesh);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
