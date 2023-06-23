# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY FixedDecoupledScalarSolutions)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  MinkowskiZeroScalar.cpp
  KerrSchildScalar.cpp
  KerrSphericalHarmonic.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  MinkowskiZeroScalar.hpp
  KerrSchildScalar.hpp
  KerrSphericalHarmonic.hpp
  Solutions.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  ErrorHandling
  GeneralRelativitySolutions
  CurvedScalarWave
  ScalarDriver
  InitialDataUtilities
  Options
  Serialization
  )
