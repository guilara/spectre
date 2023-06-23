# Distributed under the MIT License.
# See LICENSE.txt for details.

function(add_scalar_tensor_executable_without_horizon
    DIM LIBS_TO_LINK)
  add_spectre_parallel_executable(
    "EvolveFixedDecoupledScalar${DIM}D"
    EvolveFixedDecoupledScalar
    Evolution/Executables/FixedScalarTensor/FixedDecoupledScalar
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
  ScalarDriver
  FixedDecoupledScalar
  FixedDecoupledScalarSolutions
  GhFixedDecoupledScalarSolutions
  )

add_scalar_tensor_executable_without_horizon(
  3
  "${LIBS_TO_LINK};GeneralRelativitySolutions"
  )

add_spectre_parallel_executable(
  EvolveFixedDecoupledScalarSingleBlackHole
  EvolveFixedDecoupledScalarSingleBlackHole
  Evolution/Executables/FixedScalarTensor/FixedDecoupledScalar
  "EvolutionMetavars<3, false>"
  "${LIBS_TO_LINK};ApparentHorizons;ControlSystem;GeneralRelativitySolutions"
  )
