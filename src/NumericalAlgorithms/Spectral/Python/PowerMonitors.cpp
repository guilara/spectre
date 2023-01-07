// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/Spectral/PowerMonitors.hpp"

#include <pybind11/pybind11.h>

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace py = pybind11;

namespace PowerMonitors::py_bindings {

void bind_power_monitor_array(py::module& m) {
  m.def(
      "power_monitor_array",
      [](const DataVector& input_data_vector, const Mesh<1>& mesh) {
        return PowerMonitors::power_monitor_array(input_data_vector, mesh);
      },
      py::arg("DataVector"), py::arg("mesh"),
      py::return_value_policy::reference, "PowerMonitors array ");
  m.def(
      "power_monitor_array",
      [](const DataVector& input_data_vector, const Mesh<2>& mesh) {
        return PowerMonitors::power_monitor_array(input_data_vector, mesh);
      },
      py::arg("DataVector"), py::arg("mesh"),
      py::return_value_policy::reference, "PowerMonitors array ");
  m.def(
      "power_monitor_array",
      [](const DataVector& input_data_vector, const Mesh<3>& mesh) {
        return PowerMonitors::power_monitor_array(input_data_vector, mesh);
      },
      py::arg("DataVector"), py::arg("mesh"),
      py::return_value_policy::reference, "PowerMonitors array ");
}

}  // namespace PowerMonitors::py_bindings
