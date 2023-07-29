# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_GhScalarTensorSolutions")

set(LIBRARY_SOURCES
  Test_InstantiateWrappedGr.cpp
  )

add_test_library(
  ${LIBRARY}
  "PointwiseFunctions/AnalyticSolutions/GhScalarTensor"
  "${LIBRARY_SOURCES}"
  ""
  )

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  GeneralRelativitySolutions
  GhScalarTensorSolutions
  ScalarTensorSolutions
  Options
  Utilities
  )
