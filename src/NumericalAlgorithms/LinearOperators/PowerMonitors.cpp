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

#include "Parallel/Printf.hpp"

namespace PowerMonitors {

double detail::relative_truncation_error_impl(
    const DataVector& input_power_monitors, const size_t upperBound) {
  // upperBound needs to be smaller than the number of power monitors in the
  // current dimension, i.e.
  // upperBound < gsl::at(input_power_monitors, sliced_dim).size()

  // Compute weighted average and total sum in the current dimension
  double weighted_average = 0.0;
  double weight_sum = 0.0;
  double weight_value = 0.0;
  for (size_t index = 0; index < upperBound; ++index) {
    // Compute current weight
    weight_value = exp(-square(
        index - static_cast<double>(upperBound) + 0.5)
        );
    Parallel::printf("weight = %lf \n ", weight_value);
    // Add weighted power monitor
    // (Need to check log argument or add floor)
    weighted_average +=
        weight_value * gsl::at(input_power_monitors, index);
    // Add term to weighted sum
    weight_sum += weight_value;
  }
  weighted_average /= weight_sum;

  // Maximum between the first two power monitors
  // (Need to check log argument or add floor)
  double leading_term = std::max(gsl::at(input_power_monitors, 0_st),
                                     gsl::at(input_power_monitors, 1_st));

  // Compute relative truncation error
  return leading_term - weighted_average;
}

template <size_t Dim>
void power_monitors(const gsl::not_null<std::array<DataVector, Dim>*> result,
                const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  // Get modal coefficients
  const ModalVector modal_coefficients =
      to_modal_coefficients(input_data_vector, mesh);

  double slice_sum = 0.0;
  size_t n_slice = 0;
  size_t n_stripe = 0;
  const double log_floor = 1.0e-16;
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
      slice_sum = log10(std::max(slice_sum, log_floor));

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
void relative_truncation_error(
    const gsl::not_null<std::array<double, Dim>*> result,
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  // Temporary: compute power monitors
  const auto pm_array = power_monitors<Dim>(input_data_vector, mesh);

  // Compute relative truncation error in each dimension
  size_t number_of_points_per_dim = 0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {
    // Number of power monitors per dimension
    number_of_points_per_dim = gsl::at(pm_array, sliced_dim).size();
    // Compute relative truncation error
    gsl::at(*result, sliced_dim) = detail::relative_truncation_error_impl(
        gsl::at(pm_array, sliced_dim), number_of_points_per_dim);
  }
}

template <size_t Dim>
std::array<double, Dim> relative_truncation_error(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  std::array<double, Dim> result{};
  relative_truncation_error(make_not_null(&result), input_data_vector, mesh);
  return result;
}

double maximum_of_variable(const DataVector& input_data_vector) {
  // double max_value = 0.0;
  // for (auto value : input_data_vector) {
  //   // Use L2 norm
  //   max_value += square(value);
  // }
  // max_value /= input_data_vector.size();
  // double result = sqrt(max_value);
  // Use infinity norm
  return max(abs(input_data_vector));
}

template <size_t Dim>
void error_estimate(const gsl::not_null<std::array<double, Dim>*> result,
                    const DataVector& input_data_vector,
                    const Mesh<Dim>& mesh) {
  // Define tolerances
  const double atol = 0.0;  // Check where to get this
  const double rtol = 1.0;  // Check where to get this

  // Get relative truncation error
  auto trunc_error = relative_truncation_error<Dim>(input_data_vector, mesh);
  auto umax = maximum_of_variable(input_data_vector);

  double error = 0.0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {
    double exponent = - 1.0 * gsl::at(trunc_error, sliced_dim);
    error = umax * pow(10.0, exponent);
    error /= (atol + umax * rtol);
    gsl::at(*result, sliced_dim) = error;
  }
}

template <size_t Dim>
std::array<double, Dim> error_estimate(const DataVector& input_data_vector,
                                       const Mesh<Dim>& mesh) {
  std::array<double, Dim> result{};
  error_estimate(make_not_null(&result), input_data_vector, mesh);
  return result;
}

template <size_t Dim>
void absolute_truncation_error_estimate(
    const gsl::not_null<std::array<double, Dim>*> result,
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  // Temporary: compute power monitors
  const auto pm_array = power_monitors<Dim>(input_data_vector, mesh);

  // Compute the relative truncation error with all power monitors and with
  // one less
  size_t number_of_points_per_dim = 0;
  for (size_t sliced_dim = 0; sliced_dim < Dim; ++sliced_dim) {
    // Number of power monitors
    number_of_points_per_dim = gsl::at(pm_array, sliced_dim).size();
    // Can also compare with
    double relative_truncation_error =
        detail::relative_truncation_error_impl(gsl::at(pm_array, sliced_dim),
                                               number_of_points_per_dim);
    double relative_truncation_error_with_one_fewer_mode =
        detail::relative_truncation_error_impl(gsl::at(pm_array, sliced_dim),
                                               number_of_points_per_dim - 1);
    // Return the minimum of these two
    gsl::at(*result, sliced_dim) =
        std::min(relative_truncation_error,
                 relative_truncation_error_with_one_fewer_mode);
  }
}

template <size_t Dim>
std::array<double, Dim> absolute_truncation_error_estimate(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh) {
  std::array<double, Dim> result{};
  absolute_truncation_error_estimate(make_not_null(&result), input_data_vector,
                                     mesh);
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
  template void PowerMonitors::relative_truncation_error(                      \
    const gsl::not_null<std::array<double, DIM(data)>*> result,                \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template std::array<double, DIM(data)> PowerMonitors::error_estimate(        \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template void PowerMonitors::error_estimate(                                 \
    const gsl::not_null<std::array<double, DIM(data)>*> result,                \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template std::array<double, DIM(data)>                                       \
  PowerMonitors::absolute_truncation_error_estimate(                           \
    const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);         \
  template void PowerMonitors::absolute_truncation_error_estimate(             \
      const gsl::not_null<std::array<double, DIM(data)>*> result,              \
      const DataVector& input_data_vector, const Mesh<DIM(data)>& mesh);

  GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
