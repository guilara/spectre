// Distributed under the MIT License.
// See LICENSE.txt for details.

// \file
// Tests of power monitors.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestCreation.hpp"
#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace {

void test_power_monitors_impl() {
  size_t number_of_points_per_dimension = 4;
  size_t number_of_points = pow<2>(number_of_points_per_dimension);

  // Test a constant function
  const DataVector test_data_vector = DataVector{number_of_points, 1.0};

  const Mesh<2_st> mesh{number_of_points_per_dimension,
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};

  auto test_power_monitors =
      PowerMonitors::power_monitors<2_st>(test_data_vector, mesh);

  // The only non-zero modal coefficient of a constant is the one corresponding
  // to the first Legendre polynomial
  DataVector check_data_vector =
      DataVector{number_of_points_per_dimension, 0.0};
  check_data_vector[0] = 1.0 / sqrt(number_of_points_per_dimension);

  const std::array<DataVector, 2> expected_power_monitors{check_data_vector,
                                                          check_data_vector};

  CHECK_ITERABLE_APPROX(test_power_monitors, expected_power_monitors);
}

void test_power_monitors_second_impl() {
  size_t number_of_points_per_dimension = 4;

  const Mesh<2_st> mesh{number_of_points_per_dimension,
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};

  const auto logical_coords = logical_coordinates(mesh);

  // Build a test function containing only one Legendre basis function
  // per dimension
  size_t x_mode = 0;
  size_t y_mode = 1;
  std::array<size_t, 2> coeff = {x_mode, y_mode};

  DataVector u_nodal(mesh.number_of_grid_points(), 1.0);
  for (size_t dim = 0; dim < 2; ++dim) {
    u_nodal *=
        Spectral::compute_basis_function_value<Spectral::Basis::Legendre>(
            gsl::at(coeff, dim), logical_coords.get(dim));
  }

  auto test_power_monitors =
      PowerMonitors::power_monitors<2_st>(u_nodal, mesh);

  // The only non-zero modal coefficient of a constant is the one corresponding
  // to the specified Legendre polynomial

  // In the x direction
  DataVector check_data_vector_x =
      DataVector{number_of_points_per_dimension, 0.0};
  check_data_vector_x[x_mode] = 1.0 / sqrt(number_of_points_per_dimension);

  // In the y direction
  DataVector check_data_vector_y =
      DataVector{number_of_points_per_dimension, 0.0};
  check_data_vector_y[y_mode] = 1.0 / sqrt(number_of_points_per_dimension);

  // We compare against the expected array
  const std::array<DataVector, 2> expected_power_monitors{
      check_data_vector_x, check_data_vector_y};

  CHECK_ITERABLE_APPROX(test_power_monitors, expected_power_monitors);
}

void test_relative_truncation_error_impl() {
  // We check that the relative truncation error is consistent for a basis
  // function in a mesh with a given number of collocation points and another
  // with one less point
  size_t number_of_points_per_dimension = 10;

  const Mesh<2_st> mesh{number_of_points_per_dimension,
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  const auto logical_coords = logical_coordinates(mesh);

  // We will compare with a mesh with one less point
  const Mesh<2_st> mesh_one_less{number_of_points_per_dimension - 1,
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  const auto logical_coords_one_less = logical_coordinates(mesh_one_less);

  // Build a test function containing only one Legendre basis function
  // per dimension
  size_t x_mode = 0;
  size_t y_mode = 1;
  std::array<size_t, 2> coeff = {x_mode, y_mode};

  DataVector u_nodal(mesh.number_of_grid_points(), 1.0);
  DataVector u_nodal_one_less(mesh_one_less.number_of_grid_points(), 1.0);
  for (size_t dim = 0; dim < 2; ++dim) {
    u_nodal *=
        Spectral::compute_basis_function_value<Spectral::Basis::Legendre>(
            gsl::at(coeff, dim), logical_coords.get(dim));
    u_nodal_one_less *=
        Spectral::compute_basis_function_value<Spectral::Basis::Legendre>(
            gsl::at(coeff, dim), logical_coords_one_less.get(dim));
  }

  auto test_power_monitors = PowerMonitors::power_monitors<2_st>(u_nodal, mesh);
  auto test_power_monitors_one_less =
      PowerMonitors::power_monitors<2_st>(u_nodal_one_less, mesh_one_less);

  // Compute the truncation error for u_nodal with one less mode
  auto test_relative_truncation_error_one_mode_less_x = pow(10.0,
      -1.0 * PowerMonitors::relative_truncation_error(
          gsl::at(test_power_monitors, 0_st), mesh.extents(0) - 1));
  auto test_relative_truncation_error_one_mode_less_y = pow(10.0,
      -1.0 * PowerMonitors::relative_truncation_error(
          gsl::at(test_power_monitors, 1_st), mesh.extents(1) - 1));

  // Compute the relative truncation error for u_nodal_one_less with all its
  // modes
  auto test_relative_truncation_error_log10 =
      PowerMonitors::relative_truncation_error<2_st>(u_nodal_one_less,
                                                     mesh_one_less);
  auto test_relative_truncation_error_x =
      pow(10.0, -1.0 * test_relative_truncation_error_log10[0_st]);
  auto test_relative_truncation_error_y =
      pow(10.0, -1.0 * test_relative_truncation_error_log10[1_st]);

  // Check consistency between the two truncation errors
  CHECK_ITERABLE_APPROX(test_relative_truncation_error_x,
        test_relative_truncation_error_one_mode_less_x);
  CHECK_ITERABLE_APPROX(test_relative_truncation_error_y,
        test_relative_truncation_error_one_mode_less_y);

  // Test truncation error
  // Get error estimates in both dimensions
  auto test_truncation_error =
      PowerMonitors::truncation_error<2_st>(u_nodal, mesh);

  // Compare with the result from the relative truncation error
  auto relative_truncation_error =
      PowerMonitors::relative_truncation_error<2_st>(u_nodal, mesh);
  auto relative_truncation_error_x =
      pow(10.0, -1.0 * relative_truncation_error[0]);
  auto relative_truncation_error_y =
      pow(10.0, -1.0 * relative_truncation_error[1]);
  std::array<double, 2> check_vector = {relative_truncation_error_x,
                                        relative_truncation_error_y};

  CHECK_ITERABLE_APPROX(test_truncation_error, check_vector);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Numerical.LinearOperators.PowerMonitors",
                  "[NumericalAlgorithms][LinearOperators][Unit]") {
  test_power_monitors_impl();
  test_power_monitors_second_impl();
  test_relative_truncation_error_impl();
}
