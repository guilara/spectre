# Distributed under the MIT License.
# See LICENSE.txt for details.

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
  ScalarTensorAnalyticData
  ScalarTensorSources
  GhScalarTensorAnalyticData
  ScalarDriver
  FixedDecoupledScalar
  FixedDecoupledScalarAnalyticData
  GhFixedScalarTensorAnalyticData
  )

add_spectre_parallel_executable(
  EvolveFixedDecoupledScalarSingleBlackHole
  EvolveFixedDecoupledScalarSingleBlackHole
  Evolution/Executables/FixedScalarTensor/FixedDecoupledScalar
  "EvolutionMetavars<3, false>"
  "${LIBS_TO_LINK};ApparentHorizons;ControlSystem;GeneralRelativitySolutions"
  )
