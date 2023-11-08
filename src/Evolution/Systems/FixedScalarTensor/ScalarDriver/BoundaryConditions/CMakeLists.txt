# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  AnalyticConstant.cpp
  BoundaryCondition.cpp
  ConstraintPreservingSphericalRadiation.cpp
  DemandOutgoingCharSpeeds.cpp
)

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  AnalyticConstant.hpp
  BoundaryCondition.hpp
  ConstraintPreservingSphericalRadiation.hpp
  DemandOutgoingCharSpeeds.hpp
  Factory.hpp
)
