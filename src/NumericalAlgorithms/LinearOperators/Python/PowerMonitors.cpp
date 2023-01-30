// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/Python/PowerMonitors.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace py = pybind11;

namespace PowerMonitors::py_bindings {

namespace {
template <size_t Dim>
void bind_power_monitors_impl(py::module& m) {  // NOLINT
  m.def("power_monitors",
        py::overload_cast<const DataVector&, const Mesh<Dim>&> (
          &power_monitors<Dim>),
        py::arg("input_data_vector"), py::arg("mesh"));
}

template <size_t Dim>
void bind_relative_truncation_error_impl(py::module& m) {  // NOLINT
  m.def("relative_truncation_error",
        py::overload_cast<const DataVector&, const Mesh<Dim>&>(
            &relative_truncation_error<Dim>),
        py::arg("input_data_vector"), py::arg("mesh"));
}

template <size_t Dim>
void bind_error_estimate_impl(py::module& m) {  // NOLINT
  m.def("error_estimate",
        py::overload_cast<const DataVector&, const Mesh<Dim>&>(
            &error_estimate<Dim>),
        py::arg("input_data_vector"), py::arg("mesh"));
}

}  // namespace

void bind_power_monitors(py::module& m) {
  bind_power_monitors_impl<1>(m);
  bind_power_monitors_impl<2>(m);
  bind_power_monitors_impl<3>(m);
  bind_relative_truncation_error_impl<1>(m);
  bind_relative_truncation_error_impl<2>(m);
  bind_relative_truncation_error_impl<3>(m);
  bind_error_estimate_impl<1>(m);
  bind_error_estimate_impl<2>(m);
  bind_error_estimate_impl<3>(m);
}

}  // namespace PowerMonitors::py_bindings
