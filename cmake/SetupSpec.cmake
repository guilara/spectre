# Distributed under the MIT License.
# See LICENSE.txt for details.

find_package(SpEC)

if (NOT SpEC_FOUND)
  # Make SpEC scripts available in Python independently of whether the SpEC
  # exporter has been found. These can be used until we have ported them to
  # SpECTRE.
  if (SPEC_ROOT)
    set(PYTHONPATH "${SPEC_ROOT}/Support/Python:${PYTHONPATH}")
  endif()
  return()
endif()

# Make SpEC scripts available in Python. These can be used until we have ported
# them to SpECTRE.
set(PYTHONPATH "${SPEC_ROOT}/Support/Python:${PYTHONPATH}")

file(APPEND
  "${CMAKE_BINARY_DIR}/BuildInfo.txt"
  "SpEC exporter: ${SPEC_EXPORTER_ROOT}\n"
  )
