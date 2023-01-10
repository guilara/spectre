# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "PyLinearOperators")

spectre_python_add_module(
  LinearOperators
  LIBRARY_NAME ${LIBRARY}
  SOURCES
  Bindings.cpp
  PowerMonitors.cpp
)

spectre_python_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  PowerMonitors.hpp
  )

spectre_python_link_libraries(
  ${LIBRARY}
  PRIVATE
  LinearOperators
  Spectral
  pybind11::module
)
