# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarTensorDriver)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Characteristics.cpp
  Sources.cpp
  TimeDerivative.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Characteristics.hpp
  Sources.hpp
  System.hpp
  Tags.hpp
  TimeDerivative.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Domain
  ErrorHandling
  Serialization
  Utilities
  CurvedScalarWave
  )

  add_subdirectory(BoundaryConditions)
  add_subdirectory(BoundaryCorrections)
