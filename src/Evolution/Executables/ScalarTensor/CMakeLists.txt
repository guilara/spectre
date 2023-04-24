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
  Time
  Utilities
  #
  CurvedScalarWave
  CurvedWaveEquationAnalyticData
  ScalarTensor
  ScalarTensorSolutions
  GhScalarTensorSolutions
  )

add_scalar_tensor_executable_without_horizon(
  3
  "${LIBS_TO_LINK};GeneralRelativitySolutions"
)
