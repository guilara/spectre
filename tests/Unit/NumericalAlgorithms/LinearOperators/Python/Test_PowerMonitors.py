#Distributed under the MIT License.
#See LICENSE.txt for details.

from spectre.LinearOperators import power_monitor_array
from spectre.Spectral import (Mesh2D, Basis, Quadrature)
from spectre.DataStructures import DataVector

import numpy as np
import unittest


class TestPowerMonitors(unittest.TestCase):
    # Check the case for a constant function where the power monitors
    # should be given by the first basis function
    def test_power_monitor_array(self):
        elements_per_dimension = 4

        extent = elements_per_dimension
        basis = Basis.Legendre
        quadrature = Quadrature.GaussLobatto
        mesh = Mesh2D(extent, basis, quadrature)

        np_vec = np.ones(mesh.number_of_grid_points())
        test_vec = DataVector(np_vec)

        test_array = power_monitor_array(test_vec, mesh)
        np_test_array = np.asarray(test_array)

        check_vec_0 = np.zeros(elements_per_dimension)
        check_vec_0[0] = 1.0 / np.sqrt(elements_per_dimension)

        check_vec_1 = np.zeros(elements_per_dimension)
        check_vec_1[0] = 1.0 / np.sqrt(elements_per_dimension)

        np_check_array = np.array([check_vec_0, check_vec_1])

        np.testing.assert_allclose(np_test_array, np_check_array, 1e-12, 1e-12)


if __name__ == '__main__':
    unittest.main()
