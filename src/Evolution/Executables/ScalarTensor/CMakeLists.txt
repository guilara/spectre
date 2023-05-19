# Distributed under the MIT License.
# See LICENSE.txt for details.

function(add_scalar_tensor_executable_without_horizon
    DIM LIBS_TO_LINK)
  add_spectre_parallel_executable(
    "EvolveScalarTensor${DIM}D"
    EvolveScalarTensor
    Evolution/Executables/ScalarTensor
    "EvolutionMetavars<${DIM}, false>"
    "${LIBS_TO_LINK}"
    )
endfunction(add_scalar_tensor_executable_without_horizon)

set(LIBS_TO_LINK
  CoordinateMaps
  DiscontinuousGalerkin
  Domain
  DomainCreators
  EventsAndDenseTriggers
  EventsAndTriggers
  Evolution
  GeneralRelativity
  GeneralizedHarmonic
  Informer
  Interpolation
  LinearOperators
  MathFunctions
  Observer
  Options
  Parallel
  ParallelInterpolation
  PhaseControl
  Serialization
  Time
  Utilities
  #
  CurvedScalarWave
  CurvedWaveEquationAnalyticData
  ScalarTensor
  ScalarTensorPointwise
  ScalarTensorSolutions
  ScalarTensorSources
  GhScalarTensorSolutions
  )

add_scalar_tensor_executable_without_horizon(
  3
  "${LIBS_TO_LINK};GeneralRelativitySolutions"
  )

add_spectre_parallel_executable(
  EvolveScalarTensorSingleBlackHole
  EvolveScalarTensorSingleBlackHole
  Evolution/Executables/ScalarTensor
  "EvolutionMetavars<3, false>"
  "${LIBS_TO_LINK};ApparentHorizons;ControlSystem;GeneralRelativitySolutions"
  )

# add_scalartensor_with_horizon_executable(
#   BondiHoyleAccretion
#   GeneralizedHarmonic::Solutions::WrappedGr<grmhd::AnalyticData::BondiHoyleAccretion>
#   "${LIBS_TO_LINK}"
#   )
