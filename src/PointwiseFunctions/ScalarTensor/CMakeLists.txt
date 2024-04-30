# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarTensorPointwise)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  GBScalarSourceTerms.cpp
  GBTensorSourceTerms.cpp
  ScalarCharge.cpp
  WeylElectricRicciPart.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  GBScalarSourceTerms.hpp
  GBTensorSourceTerms.hpp
  ScalarCharge.hpp
  WeylElectricRicciPart.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  INTERFACE
  ErrorHandling
  )
