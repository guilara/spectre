# Distributed under the MIT License.
# See LICENSE.txt for details.

function(add_scalartensor_executable INITIAL_DATA_NAME INITIAL_DATA LIBS_TO_LINK)
  add_spectre_parallel_executable(
    "EvolveScalarTensor${INITIAL_DATA_NAME}"
    EvolveScalarTensor
    Evolution/Executables/ScalarTensor
    "EvolutionMetavars<${INITIAL_DATA}>"
    "${LIBS_TO_LINK}"
    )
endfunction(add_scalartensor_executable)

# function(add_scalartensor_with_horizon_executable INITIAL_DATA_NAME INITIAL_DATA LIBS_TO_LINK)
#   add_spectre_parallel_executable(
#     "EvolveScalarTensor${INITIAL_DATA_NAME}"
#     EvolveScalarTensorWithHorizon
#     Evolution/Executables/ScalarTensor
#     "EvolutionMetavars<${INITIAL_DATA}>"
#     "${LIBS_TO_LINK}"
#     )
# endfunction(add_scalartensor_with_horizon_executable)

set(LIBS_TO_LINK
  CoordinateMaps
  CurvedScalarWave
  DiscontinuousGalerkin
  DomainCreators
  EventsAndDenseTriggers
  EventsAndTriggers
  Evolution
  GeneralRelativitySolutions
  GeneralizedHarmonic
  ScalarTensor
  Informer
  Limiters
  LinearOperators
  MathFunctions
  Observer
  Options
  Parallel
  Time
  Utilities
  )

add_scalartensor_executable(
  Minksowski
  GeneralizedHarmonic::Solutions::WrappedGr<gr::Solutions::Minkowski>
  "${LIBS_TO_LINK}"
  )

# add_scalartensor_with_horizon_executable(
#   BondiHoyleAccretion
#   GeneralizedHarmonic::Solutions::WrappedGr<grmhd::AnalyticData::BondiHoyleAccretion>
#   "${LIBS_TO_LINK}"
#   )
