# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarDriver)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Diagnostics.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Diagnostics.hpp
  System.hpp
  Tags.hpp
  TimeDerivativeTerms.hpp
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
  ScalarTensor
  ScalarTensorSources
  )

add_subdirectory(BoundaryConditions)
add_subdirectory(BoundaryCorrections)
