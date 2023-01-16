// Distributed under the MIT License.
// See LICENSE.txt for details.

// \file
// Tests of power monitors.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "Framework/TestCreation.hpp"
#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Math.hpp"

namespace {

void test_power_monitor_array_impl() {
  size_t number_of_points_per_dimension = 4;
  size_t number_of_points = pow<2>(number_of_points_per_dimension);

  // Test a constant function
  const DataVector test_data_vector = DataVector{number_of_points, 1.0};

  const Mesh<2_st> mesh{number_of_points_per_dimension,
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};

  auto test_power_monitor_array =
        PowerMonitors::power_monitor_array<2_st>(test_data_vector, mesh);

  // The only non-zero modal coefficient of a constant is the one corresponding
  // to the first Legendre polynomial
  DataVector check_data_vector =
      DataVector{number_of_points_per_dimension, 0.0};
  check_data_vector[0] = 1.0 / sqrt(number_of_points_per_dimension);

  const std::array<DataVector, 2> check_power_monitor_array{
      check_data_vector,
      check_data_vector};

  CHECK_ITERABLE_APPROX(test_power_monitor_array, check_power_monitor_array);
}

void test_power_monitor_array_second_impl() {
  size_t number_of_points_per_dimension = 4;
  size_t number_of_points = pow<2>(number_of_points_per_dimension);

  const Mesh<2_st> mesh{number_of_points_per_dimension,
                        Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};

  const auto logical_coords = logical_coordinates(mesh);
  DataVector u_nodal_expected(mesh.number_of_grid_points(), 0.0);

  // Build a test function containing only one Legendre basis function
  double basis_factor = 1.0;
  size_t x_mode = 0;
  size_t y_mode = 1;
  std::vector<std::array<size_t, 2>> coeffs_to_include = {{x_mode, y_mode}};
  auto coeff = coeffs_to_include[0];

  DataVector basis_function =
          basis_factor *
          Spectral::compute_basis_function_value<Spectral::Basis::Legendre>(
              coeff[0], get<0>(logical_coords));
  for (size_t dim = 1; dim < 2; ++dim) {
    basis_function *=
        1.0 *
        Spectral::compute_basis_function_value<Spectral::Basis::Legendre>(
            gsl::at(coeff, dim), logical_coords.get(dim));
  }
  u_nodal_expected += basis_function;

  auto test_power_monitor_array =
      PowerMonitors::power_monitor_array<2_st>(u_nodal_expected, mesh);

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

  // We test against this array
  const std::array<DataVector, 2> check_power_monitor_array{check_data_vector_x,
                                                        check_data_vector_y};

  CHECK_ITERABLE_APPROX(test_power_monitor_array, check_power_monitor_array);
}

void test_power_monitor_array() {
  test_power_monitor_array_impl();
  test_power_monitor_array_second_impl();
}

} // namespace

SPECTRE_TEST_CASE("Unit.Numerical.LinearOperators.PowerMonitors",
                  "[NumericalAlgorithms][LinearOperators][Unit]") {
    test_power_monitor_array();
}
