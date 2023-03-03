// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
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
  for (size_t d = 0; d < Dim; ++d) {
    n_slice = mesh.extents().slice_away(d).product();
    n_stripe = mesh.extents(d);

    gsl::at(*result, d).destructive_resize(n_stripe);

    for (size_t i = 0; i < n_stripe; ++i) {
      slice_sum = 0.0;
      for (SliceIterator si(mesh.extents(), d, i); si; ++si) {
        slice_sum += square(modal_coefficients[si.volume_offset()]);
      }
      slice_sum /= n_slice;
      slice_sum = sqrt(slice_sum);

      gsl::at(*result, d)[i] = slice_sum;
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

double relative_truncation_error(const DataVector& input_power_monitors,
                                 const size_t num_modes) {
  ASSERT(num_modes <= input_power_monitors.size(),
         "Number of modes needs less or equal than the number of input power "
         "monitors");
  ASSERT(2_st <= num_modes,
         "Number of modes needs to be larger or equal than 2.");
  // Compute weighted average and total sum in the current dimension
  double weighted_average = 0.0;
  double weight_sum = 0.0;
  double weight_value = 0.0;
  for (size_t i = 0; i < num_modes; ++i) {
    // Compute current weight
    weight_value = exp(-square(i - static_cast<double>(num_modes) + 0.5));
    // Add weighted power monitor
    const auto mode = gsl::at(input_power_monitors, i);
    if (not equal_within_roundoff(mode, 0.0)) {
      weighted_average += weight_value * log10(mode);
    }
    else {
      // If zero within roundoff set a maximum resolution
      weighted_average += - 34.0 * weight_value;
    }
    // Add term to weighted sum
    weight_sum += weight_value;
  }
  weighted_average /= weight_sum;

  // Maximum between the first two power monitors
  double leading_term = std::max(gsl::at(input_power_monitors, 0_st),
                                 gsl::at(input_power_monitors, 1_st));
  if (not equal_within_roundoff(leading_term, 0.0)) {
    leading_term = log10(leading_term);
  }
  else {
    // If zero within roundoff set a maximum resolution
    leading_term = - 34.0;
  }

  // // Compute relative truncation error
  // return leading_term - weighted_average;
  // Compute relative truncation error
  // Pileup modes may change the expected sign from positive to negative
  // Prevent this by taking the maximum with zero
  return std::max(leading_term - weighted_average, 0.0);
}

template <size_t Dim>
std::array<double, Dim> relative_truncation_error(
    const std::array<DataVector, Dim>& input_power_monitors) {
  std::array<double, Dim> result{};
  // Compute relative truncation error in each dimension
  for (size_t d = 0; d < Dim; ++d) {
    const auto& modes = gsl::at(input_power_monitors, d);
    // Compute relative truncation error
    gsl::at(result, d) =
        relative_truncation_error(modes, modes.size());
  }
  return result;
}

template <size_t Dim>
std::array<double, Dim> relative_truncation_error(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  ASSERT(2_st <= input_data_vector.size(),
         "Number of data points needs to be larger or equal than 2.");
  const auto pm_array = power_monitors<Dim>(input_data_vector, mesh);
  return relative_truncation_error(pm_array);
}

template <size_t Dim>
std::array<double, Dim> truncation_error(const DataVector& input_data_vector,
                                       const Mesh<Dim>& mesh) {
   return truncation_error<Dim>(
      relative_truncation_error<Dim>(input_data_vector, mesh),
      input_data_vector);
}

template <size_t Dim>
std::array<double, Dim> truncation_error(
    const std::array<double, Dim>& relative_truncation_error,
    const DataVector& input_data_vector){
  std::array<double, Dim> result{};
  // Use infinity norm
  const double umax = max(abs(input_data_vector));
  for (size_t d = 0; d < Dim; ++d) {
    gsl::at(result, d) =
        umax * pow(10.0, -1.0 * gsl::at(relative_truncation_error, d));
  }
  return result;
}

template <size_t Dim>
double truncation_error_max(const Scalar<DataVector>& input_scalar,
                            const Mesh<Dim>& mesh) {
  const std::array<double, Dim> error_in_all_dimensions =
      truncation_error(input_scalar.get(), mesh);
  const auto result_it = std::max_element(error_in_all_dimensions.begin(),
                              error_in_all_dimensions.end());
  const double result = *result_it;
  return result;
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
  template std::array<double, DIM(data)>                                       \
    PowerMonitors::relative_truncation_error(                                  \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template std::array<double, DIM(data)> PowerMonitors::truncation_error(      \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template std::array<double, DIM(data)> PowerMonitors::truncation_error(      \
    const std::array<double, DIM(data)>& relative_truncation_error,            \
    const DataVector& input_data_vector);                                      \
  template double PowerMonitors::truncation_error_max(                         \
    const Scalar<DataVector>& input_scalar, const Mesh<DIM(data)>& mesh);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
