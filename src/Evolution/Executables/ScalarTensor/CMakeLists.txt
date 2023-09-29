# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  CoordinateMaps
  CurvedScalarWave
  DiscontinuousGalerkin
  Domain
  DomainCreators
  EventsAndDenseTriggers
  EventsAndTriggers
  Evolution
  GeneralRelativity
  GeneralRelativityAnalyticData
  GeneralRelativitySolutions
  GeneralizedHarmonic
  GhScalarTensorAnalyticData
  Informer
  Interpolation
  LinearOperators
  MathFunctions
  Observer
  Options
  Parallel
  ParallelInterpolation
  PhaseControl
  ScalarTensor
  ScalarTensorAnalyticData
  ScalarTensorPointwise
  Serialization
  Time
  Utilities
  )

add_spectre_parallel_executable(
  EvolveScalarTensorSingleBlackHole
  EvolveScalarTensorSingleBlackHole
  Evolution/Executables/ScalarTensor
  "EvolutionMetavars"
  "${LIBS_TO_LINK};ApparentHorizonFinder;ControlSystem"
  )

add_spectre_parallel_executable(
  "EvolveScalarTensorBinaryBlackHole"
  EvolveScalarTensorBinaryBlackHole
  Evolution/Executables/ScalarTensor
  "EvolutionMetavars"
  "${LIBS_TO_LINK};ApparentHorizonFinder;Cce;ControlSystem;Deadlock;EvolutionTriggers"
)


if (TARGET SpEC::Exporter)
  target_compile_definitions(
    EvolveScalarTensorBinaryBlackHole PRIVATE HAS_SPEC_EXPORTER)
endif()
