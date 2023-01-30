// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"

#include <array>
#include <cmath>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/SliceIterator.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace PowerMonitors {

template <size_t Dim>
void power_monitors(const gsl::not_null<std::array<DataVector, Dim>*> result,
                const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  // Get modal coefficients
  const ModalVector modal_coefficients =
      to_modal_coefficients(input_data_vector, mesh);

  double slice_sum = 0.0;
  size_t n_slice = 0;
  size_t n_stripe = 0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {
    n_slice = mesh.extents().slice_away(sliced_dim).product();
    n_stripe = mesh.extents(sliced_dim);

    gsl::at(*result, sliced_dim).destructive_resize(n_stripe);

    for (size_t index = 0; index < n_stripe; ++index) {
      slice_sum = 0.0;
      for (SliceIterator si(mesh.extents(), sliced_dim, index); si; ++si) {
        slice_sum += square(modal_coefficients[si.volume_offset()]);
      }
      slice_sum /= n_slice;
      slice_sum = sqrt(slice_sum);

      gsl::at(*result, sliced_dim)[index] = slice_sum;
    }
  }
}

template <size_t Dim>
std::array<DataVector, Dim> power_monitors(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  std::array<DataVector, Dim> result{};
  power_monitors(make_not_null(&result), input_data_vector, mesh);
  return result;
}

template <size_t Dim>
void truncation_error(gsl::not_null<std::array<double, Dim>*> result,
                      const DataVector& input_data_vector,
                      const Mesh<Dim>& mesh) {
  // Temporary: compute power monitors
  auto pm_array = power_monitors<Dim>(input_data_vector, mesh);

  // Compute relative truncation error in each dimension
  size_t number_of_points_per_dim = 0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {

    // Number of power monitors in the current dimension
    number_of_points_per_dim = gsl::at(pm_array, sliced_dim).size();

    // Compute weighted average and total sum in the current dimension
    double weighted_average = 0.0;
    double weight_sum = 0.0;
    double weight_value = 0.0;
    for (size_t index = 0; index < number_of_points_per_dim; ++index) {
      // Compute current weight
      weight_value =
          exp(-square(index - number_of_points_per_dim + 0.5));
      // Add weighted power monitor
      // (Need to check log argument or add floor)
      weighted_average += weight_value *
                          log10(gsl::at(pm_array, sliced_dim)[index]);
      // Add term to weighted sum
      weight_sum += weight_value;
    }
    weighted_average /= weight_sum;

    // Maximum between the first two power monitors
    // (Need to check log argument or add floor)
    double first_term = log10(std::max(gsl::at(pm_array, sliced_dim)[0],
                                  gsl::at(pm_array, sliced_dim)[1]));

    // Compute relative truncation error
    gsl::at(*result, sliced_dim) = first_term - weighted_average;
  }
}

template <size_t Dim>
std::array<DataVector, Dim> truncation_error(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  std::array<double, Dim> result{};
  truncation_error(make_not_null(&result), input_data_vector, mesh);
  return result;
}

void maximum_of_variable(gsl::not_null<double*> result,
                         const DataVector& input_data_vector) {
  //
  double max_value = 0.0;
  for (auto value : input_data_vector) {
    max_value += square(value);
  }
  max_value /= input_data_vector.size();
  *result = sqrt(max_value);
}

}  // namespace PowerMonitors

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                   \
  template std::array<DataVector, DIM(data)>                                   \
  PowerMonitors::power_monitors(const DataVector& input_data_vector,           \
                                     const Mesh<DIM(data)>& mesh);             \
  template void PowerMonitors::power_monitors(                                 \
    const gsl::not_null<std::array<DataVector, DIM(data)>*> result,            \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template void PowerMonitors::truncation_error(                               \
    const gsl::not_null<std::array<double, DIM(data)>*> result,                \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
