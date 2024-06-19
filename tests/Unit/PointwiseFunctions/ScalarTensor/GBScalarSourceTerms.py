# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def order_reduced_gb_scalar(
    weyl_electric_scalar,
    weyl_magnetic_scalar,
    trace_reversed_stress_energy,
    trace_of_trace_reversed_stress_energy,
    inverse_spacetime_metric,
    tensor_driver,
    trace_of_tensor_driver,
):
    result = 8.0 * (weyl_electric_scalar - weyl_magnetic_scalar)

    # Sum the rhs tensors
    SET = trace_reversed_stress_energy + tensor_driver
    # Sum the traces
    trace_of_SET = (
        trace_of_trace_reversed_stress_energy + trace_of_tensor_driver
    )

    # Raise the second index
    SET_down_up = np.einsum("ij,jk->ik", SET, inverse_spacetime_metric)

    result += -2.0 * np.einsum("ij,ji", SET_down_up, SET_down_up)

    result += (2.0 / 3.0) * trace_of_SET * trace_of_SET

    return result
