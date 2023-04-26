# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarTensorSolutions)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  # KerrSchildScalar.cpp
  MinkowskiZeroScalar.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  # KerrSchildScalar.hpp
  MinkowskiZeroScalar.hpp
  Solutions.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  ErrorHandling
  GeneralRelativitySolutions
  CurvedScalarWave
  InitialDataUtilities
  Options
  )
