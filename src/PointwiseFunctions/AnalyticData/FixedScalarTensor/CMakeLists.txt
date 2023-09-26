# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY FixedDecoupledScalarAnalyticData)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  MinkowskiZeroScalar.cpp
  KerrSphericalHarmonic.cpp
  SphericalKerrSchildSH.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  AnalyticData.hpp
  MinkowskiZeroScalar.hpp
  KerrSphericalHarmonic.hpp
  SphericalKerrSchildSH.hpp
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
