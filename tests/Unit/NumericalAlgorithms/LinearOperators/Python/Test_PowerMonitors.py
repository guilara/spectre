#Distributed under the MIT License.
#See LICENSE.txt for details.

from spectre.NumericalAlgorithms.LinearOperators import (
    power_monitors, relative_truncation_error, error_estimate)
from spectre.Spectral import (Mesh2D, Basis, Quadrature)
from spectre.DataStructures import DataVector

import numpy as np
import unittest


class TestPowerMonitors(unittest.TestCase):
    # Check the case for a constant function where the power monitors
    # should be given by the first basis function
    def test_power_monitors(self):
        num_points_per_dimension = 4

        extent = num_points_per_dimension
        basis = Basis.Legendre
        quadrature = Quadrature.GaussLobatto
        mesh = Mesh2D(extent, basis, quadrature)

        np_vec = np.ones(mesh.number_of_grid_points())
        test_vec = DataVector(np_vec)

        test_array = power_monitors(test_vec, mesh)
        np_test_array = np.asarray(test_array)

        # Define a floor before taking the logarithm
        log_floor = 1.0e-16

        check_vec_0 = log_floor * np.ones(num_points_per_dimension)
        check_vec_0[0] = 1.0 / np.sqrt(num_points_per_dimension)
        check_vec_0 = np.log10(check_vec_0)

        check_vec_1 = log_floor * np.ones(num_points_per_dimension)
        check_vec_1[0] = 1.0 / np.sqrt(num_points_per_dimension)
        check_vec_1 = np.log10(check_vec_1)

        np_check_array = np.array([check_vec_0, check_vec_1])

        np.testing.assert_allclose(np_test_array, np_check_array, 1e-12, 1e-12)

    # Check that the truncation error for a constant unit function is
    # consistent with numerical floor when we use a high order basis
    def test_relative_truncation_error(self):
        num_points_per_dimension = 3

        extent = num_points_per_dimension
        basis = Basis.Legendre
        quadrature = Quadrature.GaussLobatto
        mesh = Mesh2D(extent, basis, quadrature)

        np_vec = np.ones(mesh.number_of_grid_points())
        test_vec = DataVector(np_vec)

        test_relative_truncation_error_exponent = np.asarray(
            relative_truncation_error(test_vec, mesh))

        test_relative_truncation_error = np.power(
            10.0 * np.ones(2), -1.0 * test_relative_truncation_error_exponent)

        expected_truncation_error = 1.0e-16 * np.ones(2)

        np.testing.assert_allclose(test_relative_truncation_error,
                                   expected_truncation_error, 1e-12, 1e-12)

    # Check that the truncation error for a constant unit function is
    # consistent with numerical floor when we use a high order basis
    def test_error_estimate(self):
        num_points_per_dimension = 3

        extent = num_points_per_dimension
        basis = Basis.Legendre
        quadrature = Quadrature.GaussLobatto
        mesh = Mesh2D(extent, basis, quadrature)

        np_vec = np.ones(mesh.number_of_grid_points())
        test_vec = DataVector(np_vec)

        test_error_estimate = np.asarray(error_estimate(test_vec, mesh))

        expected_truncation_error = 1.0e-16 * np.ones(2)

        np.testing.assert_allclose(test_error_estimate,
                                   expected_truncation_error, 1e-12, 1e-12)


if __name__ == '__main__':
    unittest.main()
