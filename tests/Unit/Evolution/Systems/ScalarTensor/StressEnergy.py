# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def trace_reversed_stress_energy(spacetime_metric, shift, lapse):

    return np.zeros(spacetime_metric.shape)


def add_stress_energy_term_to_dt_pi(trace_reversed_stress_energy, lapse):
    return 0.1234 - 16.0 * np.pi * lapse * trace_reversed_stress_energy
