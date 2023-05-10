# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarTensorSources)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ScalarSource.cpp
  )
ourceTerm
spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ScalarSource.hpp
  Tags.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  ScalarTensor
  )
