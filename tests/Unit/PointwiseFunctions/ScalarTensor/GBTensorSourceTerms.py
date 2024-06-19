# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def spacetime_derivative_scalar(lapse, pi_scalar, phi_scalar):
    # n_{0} Pi; where n_{0} = - lapse
    deriv_0 = -lapse * pi_scalar

    spacetime_derivative_scalar = np.concatenate(([deriv_0], phi_scalar))

    return spacetime_derivative_scalar
