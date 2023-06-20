# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarDriver)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Characteristics.cpp
  Sources.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Characteristics.hpp
  Initialize.hpp
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
  GeneralRelativity
  Serialization
  Utilities
  CurvedScalarWave
  )

add_subdirectory(BoundaryConditions)
add_subdirectory(BoundaryCorrections)
