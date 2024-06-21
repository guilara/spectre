# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def spacetime_derivative_scalar(lapse, shift, pi_scalar, phi_scalar):
    # n_{0} Pi + gamma^{i}_{0} Phi_{i}; where n_{0} = - lapse
    deriv_0 = -lapse * pi_scalar + np.dot(shift, phi_scalar)

    spacetime_derivative_scalar = np.concatenate(([deriv_0], phi_scalar))

    return spacetime_derivative_scalar
