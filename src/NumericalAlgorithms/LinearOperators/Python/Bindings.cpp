// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>

#include "NumericalAlgorithms/LinearOperators/Python/PowerMonitors.hpp"

PYBIND11_MODULE(_PyLinearOperators, m) {  // NOLINT
  PowerMonitors::py_bindings::bind_power_monitor_array(m);
}
