# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_ScalarTensorSolutions")

set(LIBRARY_SOURCES
  Test_InstantiateWrappedGr.cpp
  )

add_test_library(
  ${LIBRARY}
  "PointwiseFunctions/AnalyticSolutions/GhScalarTensor"
  "${LIBRARY_SOURCES}"
  "GeneralRelativitySolutions;GhScalarTensorSolutions;ScalarTensorSolutions;Options;Utilities"
  )
