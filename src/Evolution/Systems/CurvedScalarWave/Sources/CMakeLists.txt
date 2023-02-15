# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY CurvedScalarWaveSources)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  CurvatureSource.cpp
  SourceTerm.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  CurvatureSource.hpp
  SourceTerm.hpp
  Tags.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  CurvedScalarWave
  )
