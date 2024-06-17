// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <vector>

#include "ControlSystem/Actions/LimitTimeStep.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "Domain/Creators/Factory1D.hpp"
#include "Domain/Creators/Factory2D.hpp"
#include "Domain/Creators/Factory3D.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/FunctionsOfTimeAreReady.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsCharacteristicSpeeds.hpp"
#include "Evolution/Actions/RunEventsAndDenseTriggers.hpp"
#include "Evolution/ComputeTags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Factory.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Initialization/NonconservativeSystem.hpp"
#include "Evolution/Initialization/SetVariables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/SetInitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Equations.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Gauges.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/SetPiAndPhiFromConstraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Initialize.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/InterpolationTarget.hpp"
//
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
//
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/CurvedScalarWave/Initialize.hpp"
#include "Evolution/Systems/CurvedScalarWave/PsiSquared.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
//
#include "Evolution/Systems/ScalarTensor/Actions/InitializeConstraintGammas.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/ProductOfConditions.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/ProductOfCorrections.hpp"
#include "Evolution/Systems/ScalarTensor/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Initialize.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/StressEnergy.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
//
//
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/BoundaryConditions/ProductOfConditions.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/BoundaryCorrections/ProductOfCorrections.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/Initialize.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/Tags.hpp"
//
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Diagnostics.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Sources.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
//
#include "Evolution/Tags/Filter.hpp"
#include "Evolution/TypeTraits.hpp"
#include "IO/Importers/Actions/RegisterWithElementDataReader.hpp"
#include "IO/Importers/ElementDataReader.hpp"
#include "IO/Observer/Actions/ObserverRegistration.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/ExponentialFilter.hpp"
#include "NumericalAlgorithms/LinearOperators/FilterAction.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseControl/CheckpointAndExitAfterWallclock.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "Parallel/PhaseControl/Factory.hpp"
#include "Parallel/PhaseControl/VisitAndReturn.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Actions/AddSimpleTags.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
//
#include "ParallelAlgorithms/Actions/RandomizeVariables.hpp"
//
#include "Evolution/Actions/RunEventsAndDenseTriggers.hpp"
#include "ParallelAlgorithms/Actions/MemoryMonitor/ContributeMemoryData.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ErrorOnFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/IgnoreFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/Events/Factory.hpp"
#include "ParallelAlgorithms/Events/MonitorMemory.hpp"
#include "ParallelAlgorithms/Events/ObserveTimeStep.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/CleanUpInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolationTargetReceiveVars.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceivePoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceiveVolumeData.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorRegisterElement.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/TryToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Events/Interpolate.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Interpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/SphericalKerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
//
#include "PointwiseFunctions/AnalyticData/GhScalarTensor/Factory.hpp"
//
#include "PointwiseFunctions/AnalyticData/GhFixedScalarTensor/Factory.hpp"
//
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintGammas.hpp"
#include "PointwiseFunctions/GeneralRelativity/Psi4Real.hpp"
#include "PointwiseFunctions/GeneralRelativity/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
//
#include "PointwiseFunctions/ScalarTensor/ScalarCharge.hpp"
//
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Actions/RecordTimeStepperData.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/Actions/UpdateU.hpp"
#include "Time/StepChoosers/Cfl.hpp"
#include "Time/StepChoosers/Constant.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/StepChoosers/Increase.hpp"
#include "Time/StepChoosers/PreventRapidIncrease.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/StepChoosers/StepToTimes.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {

struct Inertial;
}  // namespace Frame
namespace PUP {
class er;
}  // namespace PUP
namespace Parallel {
template <typename Metavariables>
class CProxy_GlobalCache;
}  // namespace Parallel
/// \endcond

template <typename EvolutionMetavarsDerived>
struct FixedScalarTensorTemplateBase;

namespace detail {
template <bool UseNumericalInitialData>
constexpr auto make_default_phase_order() {
  if constexpr (UseNumericalInitialData) {
    //   if constexpr (false) {
    // Register needs to be before InitializeTimeStepperHistory so that CCE is
    // properly registered when the self-start happens
    return std::array{Parallel::Phase::Initialization,
                      Parallel::Phase::RegisterWithElementDataReader,
                      Parallel::Phase::ImportInitialData,
                      Parallel::Phase::InitializeInitialDataDependentQuantities,
                      Parallel::Phase::Register,
                      Parallel::Phase::InitializeTimeStepperHistory,
                      Parallel::Phase::Evolve,
                      Parallel::Phase::Exit};
  } else {
    return std::array{Parallel::Phase::Initialization,
                      Parallel::Phase::InitializeInitialDataDependentQuantities,
                      Parallel::Phase::Register,
                      Parallel::Phase::InitializeTimeStepperHistory,
                      Parallel::Phase::Evolve,
                      Parallel::Phase::Exit};
  }
}

template <size_t VolumeDim>
struct ObserverTags {
  //   static constexpr size_t volume_dim = VolumeDim;
  static constexpr size_t volume_dim = 3_st;
  //   using system = gh::System<volume_dim>;
  //   using system_scalar_tensor = ScalarTensor::System;
  using system = fe::sgb::System;

  using variables_tag = typename system::variables_tag;
  //   using variables_tag_scalar_tensor =
  //       typename system_scalar_tensor::variables_tag;

  using analytic_solution_fields = typename variables_tag::tags_list;
  //   using analytic_solution_fields_scalar_tensor =
  //       typename variables_tag_scalar_tensor::tags_list;

  //   using initial_data_list =
  //       gh::Solutions::all_solutions<volume_dim>;
  //   using initial_data_list =
  //       tmpl::list<gh::Solutions::WrappedGr<
  //           gr::Solutions::Minkowski<volume_dim>>>;
  using initial_data_list = gh::fe::sgb::AnalyticData::all_analytic_data;
  using analytic_compute = evolution::Tags::AnalyticSolutionsCompute<
      volume_dim, analytic_solution_fields, false, initial_data_list>;
  //   using analytic_compute_scalar_tensor =
  //       evolution::Tags::AnalyticSolutionsCompute<
  //           volume_dim, analytic_solution_fields_scalar_tensor, false,
  //           initial_data_list>;
  using deriv_compute = ::Tags::DerivCompute<
      variables_tag, domain::Tags::Mesh<volume_dim>,
      domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                    Frame::Inertial>,
      typename system::gradient_variables>;
  //   using deriv_compute_scalar_tensor = ::Tags::DerivCompute<
  //       variables_tag_scalar_tensor,
  //       domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
  //                                     Frame::Inertial>,
  //       typename system_scalar_tensor::gradient_variables>;

  //   using error_compute = Tags::ErrorsCompute<analytic_solution_fields>;
  //   using error_compute_scalar_tensor =
  //       Tags::ErrorsCompute<analytic_solution_fields_scalar_tensor>;

  //   using error_tags = db::wrap_tags_in<Tags::Error,
  //   analytic_solution_fields>; using error_tags_scalar_tensor =
  //       db::wrap_tags_in<Tags::Error,
  //       analytic_solution_fields_scalar_tensor>;

  using observe_fields = tmpl::append<
      tmpl::push_back<
          tmpl::append<
              // GH field tags
              system::gh_system::gh_system::variables_tag::tags_list,
              // Driver tags
              system::scalar_system::variables_tag::tags_list>,
          // Csw field tags
          ScalarTensor::Tags::CswCompute<CurvedScalarWave::Tags::Psi>,
          ScalarTensor::Tags::CswCompute<CurvedScalarWave::Tags::Pi>,
          ScalarTensor::Tags::CswCompute<
              CurvedScalarWave::Tags::Phi<volume_dim>>,
          //   system::scalar_system::variables_tag::tags_list,
          // (These gauge tags need subsequent tags to compile [why?])
          gh::Tags::GaugeH<DataVector, volume_dim, Frame::Inertial>,
          gh::Tags::SpacetimeDerivGaugeH<DataVector, volume_dim,
                                         Frame::Inertial>,
          // 3 plus 1 Tags and derivatives
          gr::Tags::SpatialMetric<DataVector, volume_dim, Frame::Inertial>,
          gr::Tags::DetSpatialMetric<DataVector>,
          gr::Tags::InverseSpatialMetric<DataVector, volume_dim,
                                         Frame::Inertial>,
          gr::Tags::Shift<DataVector, volume_dim, Frame::Inertial>,
          gr::Tags::Lapse<DataVector>,
          //   gr::Tags::SqrtDetSpatialMetricCompute<volume_dim,
          //   Frame::Inertial,
          //                                         DataVector>,
          gr::Tags::SqrtDetSpatialMetric<DataVector>,
          gr::Tags::SpacetimeNormalOneForm<DataVector, volume_dim,
                                           Frame::Inertial>,
          gr::Tags::SpacetimeNormalVector<DataVector, volume_dim,
                                          Frame::Inertial>,
          gr::Tags::InverseSpacetimeMetric<DataVector, volume_dim,
                                           Frame::Inertial>,
          ::Tags::deriv<
              gr::Tags::SpatialMetric<DataVector, volume_dim, Frame::Inertial>,
              tmpl::size_t<volume_dim>, Frame::Inertial>,
          gr::Tags::SpatialChristoffelFirstKind<DataVector, volume_dim,
                                                Frame::Inertial>,
          gr::Tags::SpatialChristoffelSecondKind<DataVector, volume_dim,
                                                 Frame::Inertial>,
          // 3 plus 1 variables used by CSW
          gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, volume_dim,
                                                      Frame::Inertial>,
          gr::Tags::ExtrinsicCurvature<DataVector, volume_dim, Frame::Inertial>,
          gr::Tags::TraceExtrinsicCurvature<DataVector>,
          // More 3 plus 1 variables
          //   ::Tags::DerivTensorCompute<
          //       gr::Tags::SpatialChristoffelSecondKind<
          //           volume_dim, ::Frame::Inertial, DataVector>,
          //    ::domain::Tags::InverseJacobian<volume_dim,
          //    Frame::ElementLogical,
          //                                       Frame::Inertial>>,
          ::Tags::deriv<gr::Tags::SpatialChristoffelSecondKind<
                            DataVector, volume_dim, Frame::Inertial>,
                        tmpl::size_t<volume_dim>, Frame::Inertial>,
          //   gr::Tags::SpatialRicciCompute<volume_dim, ::Frame::Inertial,
          //                                 DataVector>,
          gr::Tags::SpatialRicci<DataVector, volume_dim, Frame::Inertial>,
          //   gr::Tags::SpatialRicciScalarCompute<volume_dim,
          //   ::Frame::Inertial,
          //                                       DataVector>,
          gr::Tags::SpatialRicciScalar<DataVector>,
          // Compute the constraints of GH
          gh::Tags::GaugeConstraintCompute<volume_dim, Frame::Inertial>,
          gh::Tags::TwoIndexConstraintCompute<volume_dim, Frame::Inertial>,
          gh::Tags::ThreeIndexConstraintCompute<volume_dim, Frame::Inertial>,
          // Compute the constraints of CSW
          ScalarTensor::Tags::CswOneIndexConstraintCompute<volume_dim>,
          ScalarTensor::Tags::CswTwoIndexConstraintCompute<volume_dim>,
          // GH constraint norms
          ::Tags::PointwiseL2NormCompute<gh::Tags::GaugeConstraint<
              DataVector, volume_dim, Frame::Inertial>>,
          ::Tags::PointwiseL2NormCompute<gh::Tags::TwoIndexConstraint<
              DataVector, volume_dim, Frame::Inertial>>,
          ::Tags::PointwiseL2NormCompute<gh::Tags::ThreeIndexConstraint<
              DataVector, volume_dim, Frame::Inertial>>,
          // CSW constraint norms
          ::Tags::PointwiseL2NormCompute<ScalarTensor::Tags::Csw<
              CurvedScalarWave::Tags::OneIndexConstraint<volume_dim>>>,
          ::Tags::PointwiseL2NormCompute<ScalarTensor::Tags::Csw<
              CurvedScalarWave::Tags::TwoIndexConstraint<volume_dim>>>,
          // Damping parameters
          gh::ConstraintDamping::Tags::ConstraintGamma0,
          gh::ConstraintDamping::Tags::ConstraintGamma1,
          gh::ConstraintDamping::Tags::ConstraintGamma2,
          ScalarTensor::Tags::CswCompute<
              CurvedScalarWave::Tags::ConstraintGamma1>,
          ScalarTensor::Tags::CswCompute<
              CurvedScalarWave::Tags::ConstraintGamma2>,
          //   fe::ScalarDriver::Tags::ConstraintGamma1,
          //   fe::ScalarDriver::Tags::ConstraintGamma2,
          // Scalar Driver parameters
          fe::ScalarTensorDriver::Tags::SigmaParameter,
          fe::ScalarTensorDriver::Tags::TauParameter,
          // Sources
          ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, volume_dim,
                                                        Frame::Inertial>,
          //   ScalarTensor::Sources::Tags::ScalarSourceCompute,
          ScalarTensor::Tags::ScalarSource,
          // More diagnostic compute tags
          ScalarTensor::Tags::GBScalarCompute<DataVector>,
          ScalarTensor::Tags::CouplingFunctionDerivativeCompute<DataVector>,
          // Driver quantities
          fe::ScalarTensorDriver::Tags::TargetScalar,
          fe::ScalarTensorDriver::Tags::TargetTensor<DataVector, volume_dim,
                                                     Frame::Inertial>,
          ::Tags::PointwiseL2NormCompute<
              fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, volume_dim,
                                                         Frame::Inertial>>,
          ::Tags::PointwiseL2NormCompute<
              fe::ScalarTensorDriver::Tags::TargetTensor<DataVector, volume_dim,
                                                         Frame::Inertial>>,
          fe::ScalarTensorDriver::Tags::TensorDriverSource<
              DataVector, volume_dim, Frame::Inertial>,
          ::Tags::PointwiseL2NormCompute<
              fe::ScalarTensorDriver::Tags::TensorDriverSource<
                  DataVector, volume_dim, Frame::Inertial>>,
          fe::ScalarTensorDriver::Tags::ScalarDriverSource,
          fe::ScalarTensorDriver::Tags::TensorTrackingDiagnosticCompute<
              Frame::Inertial, DataVector>,
          ::Tags::PointwiseL2NormCompute<
              fe::ScalarTensorDriver::Tags::TensorTrackingDiagnostic<
                  DataVector, volume_dim, Frame::Inertial>>,
          fe::ScalarTensorDriver::Tags::ScalarTrackingDiagnosticCompute<
              Frame::Inertial, DataVector>,
          ::Tags::PointwiseL2NormCompute<
              ScalarTensor::Tags::TraceReversedStressEnergy<
                  DataVector, volume_dim, Frame::Inertial>>,
          // ::Tags::PointwiseL2NormCompute<
          //     ScalarTensor::Tags::OrderReducedTraceReversedStressEnergy>,
          //   ::Tags::PointwiseL2NormCompute<
          //       fe::ScalarTensorDriver::Tags::ScalarTrackingDiagnostic>,
          // Compute the constraints of CSW
          //   fe::ScalarDriver::Tags::OneIndexConstraintCompute,
          //   fe::ScalarDriver::Tags::TwoIndexConstraintCompute,
          // Driver constraint norms
          //   ::Tags::PointwiseL2NormCompute<
          //       fe::ScalarDriver::Tags::OneIndexConstraint>,
          //   ::Tags::PointwiseL2NormCompute<
          //       fe::ScalarDriver::Tags::TwoIndexConstraint>,
          // Coordinates

          // BRT
          ScalarTensor::Tags::WeylElectricFullTraceCompute<
              DataVector, volume_dim, Frame::Inertial>,
          ScalarTensor::Tags::OrderReducedGBScalar,
          ScalarTensor::Tags::DDFPsiTensor,
          ScalarTensor::Tags::OrderReducedHTensor,
          ScalarTensor::Tags::OrderReducedHTensorRicciPart,
          ScalarTensor::Tags::OrderReducedTraceReversedStressEnergy,
          ScalarTensor::Tags::EomFromDDKGTensorDiagnosticCompute<
              Frame::Inertial>,

          ::domain::Tags::Coordinates<volume_dim, Frame::Grid>,
          ::domain::Tags::Coordinates<volume_dim, Frame::Inertial>>,
      //   error_tags,
      // The 4-index constraint is only implemented in 3d
      tmpl::conditional_t<
          volume_dim == 3,
          tmpl::list<
              gh::Tags::FourIndexConstraintCompute<3, Frame::Inertial>,
              gh::Tags::FConstraintCompute<3, Frame::Inertial>,
              ::Tags::PointwiseL2NormCompute<gh::Tags::FConstraint<
                  DataVector, volume_dim, Frame::Inertial>>,
              ::Tags::PointwiseL2NormCompute<gh::Tags::FourIndexConstraint<
                  DataVector, volume_dim, Frame::Inertial>>,
              gh::Tags::ConstraintEnergyCompute<3, Frame::Inertial>,
              //   ::Tags::DerivTensorCompute<
              //       gr::Tags::ExtrinsicCurvature<3, Frame::Inertial>,
              //       ::domain::Tags::InverseJacobian<
              //           volume_dim, Frame::ElementLogical,
              //           Frame::Inertial>>,
              ::Tags::deriv<gr::Tags::ExtrinsicCurvature<DataVector, volume_dim,
                                                         Frame::Inertial>,
                            tmpl::size_t<volume_dim>, Frame::Inertial>,
              gr::Tags::GradExtrinsicCurvature<DataVector, 3, Frame::Inertial>,
              //   gr::Tags::WeylElectricCompute<3, Frame::Inertial,
              //   DataVector>,
              gr::Tags::WeylElectric<DataVector, 3, Frame::Inertial>,
              gr::Tags::WeylElectricScalar<DataVector>,
              gr::Tags::WeylMagneticScalar<DataVector>,
              gr::Tags::Psi4RealCompute<Frame::Inertial>>,
          tmpl::list<>>>;
  using non_tensor_compute_tags = tmpl::list<
      ::Events::Tags::ObserverMeshCompute<volume_dim>,
      ::Events::Tags::ObserverCoordinatesCompute<volume_dim, Frame::Inertial>,
      ::Events::Tags::ObserverInverseJacobianCompute<
          volume_dim, Frame::ElementLogical, Frame::Inertial>,
      ::Events::Tags::ObserverJacobianCompute<volume_dim, Frame::ElementLogical,
                                              Frame::Inertial>,
      ::Events::Tags::ObserverDetInvJacobianCompute<Frame::ElementLogical,
                                                    Frame::Inertial>,
      ::Events::Tags::ObserverMeshVelocityCompute<volume_dim, Frame::Inertial>,
      //   analytic_compute, error_compute,
      gh::gauges::Tags::GaugeAndDerivativeCompute<volume_dim>>;

  using field_observations =
      dg::Events::field_observations<volume_dim, observe_fields,
                                     non_tensor_compute_tags>;

  // We collect here all the tags needed for interpolation in all surfaces
  using scalar_charge_vars_to_interpolate_to_target = tmpl::list<
      gr::Tags::SpatialMetric<DataVector, volume_dim, Frame::Inertial>,
      gr::Tags::InverseSpatialMetric<DataVector, volume_dim, Frame::Inertial>,
      CurvedScalarWave::Tags::Phi<volume_dim>, CurvedScalarWave::Tags::Psi,
      fe::ScalarTensorDriver::Tags::Psi>;

  using scalar_charge_compute_items_on_target = tmpl::list<
      ylm::Tags::ThetaPhiCompute<::Frame::Inertial>,
      ylm::Tags::RadiusCompute<::Frame::Inertial>,
      ylm::Tags::RhatCompute<::Frame::Inertial>,
      ylm::Tags::InvJacobianCompute<::Frame::Inertial>,
      ylm::Tags::JacobianCompute<::Frame::Inertial>,
      ylm::Tags::DxRadiusCompute<::Frame::Inertial>,
      ylm::Tags::NormalOneFormCompute<::Frame::Inertial>,
      ylm::Tags::OneOverOneFormMagnitudeCompute<DataVector, volume_dim,
                                                ::Frame::Inertial>,
      ylm::Tags::UnitNormalOneFormCompute<::Frame::Inertial>,
      ylm::Tags::UnitNormalVectorCompute<::Frame::Inertial>,
      gr::surfaces::Tags::AreaElementCompute<::Frame::Inertial>,
      ScalarTensor::StrahlkorperScalar::Tags::ScalarChargeIntegrandCompute,
      gr::surfaces::Tags::SurfaceIntegralCompute<
          ScalarTensor::StrahlkorperScalar::Tags::ScalarChargeIntegrand,
          ::Frame::Inertial>,
      gr::surfaces::Tags::SurfaceIntegralCompute<CurvedScalarWave::Tags::Psi,
                                                 ::Frame::Inertial>,
      gr::surfaces::Tags::SurfaceIntegralCompute<
          fe::ScalarTensorDriver::Tags::Psi, ::Frame::Inertial>,
      CurvedScalarWave::Tags::PsiSquaredCompute,
      //   fe::ScalarDriver::Tags::PsiSquaredCompute,
      gr::surfaces::Tags::SurfaceIntegralCompute<
          CurvedScalarWave::Tags::PsiSquared, ::Frame::Inertial>
      //       ,
      //   gr::surfaces::Tags::SurfaceIntegralCompute<
      //       fe::ScalarDriver::Tags::PsiSquared, ::Frame::Inertial>
      >;

  using scalar_charge_surface_obs_tags = tmpl::list<
      gr::surfaces::Tags::SurfaceIntegralCompute<
          ScalarTensor::StrahlkorperScalar::Tags::ScalarChargeIntegrand,
          ::Frame::Inertial>,
      gr::surfaces::Tags::SurfaceIntegralCompute<CurvedScalarWave::Tags::Psi,
                                                 ::Frame::Inertial>,
      gr::surfaces::Tags::SurfaceIntegralCompute<
          fe::ScalarTensorDriver::Tags::Psi, ::Frame::Inertial>,
      gr::surfaces::Tags::SurfaceIntegralCompute<
          CurvedScalarWave::Tags::PsiSquared, ::Frame::Inertial>
      //       ,
      //   gr::surfaces::Tags::SurfaceIntegralCompute<
      //       fe::ScalarDriver::Tags::PsiSquared, ::Frame::Inertial>
      >;
};

// template <size_t VolumeDim, bool LocalTimeStepping>
template <size_t VolumeDim, bool LocalTimeStepping,
          bool UseNumericalInitialData>
struct FactoryCreation : tt::ConformsTo<Options::protocols::FactoryCreation> {
  //   static constexpr size_t volume_dim = VolumeDim;
  static constexpr size_t volume_dim = 3_st;
  //   using system = gh::System<volume_dim>;
  //   using system_scalar_tensor = ScalarTensor::System;
  using system = fe::sgb::System;

  //   using initial_data_list =
  //       gh::Solutions::all_solutions<volume_dim>;
  //   using initial_data_list =
  //       tmpl::list<gh::Solutions::WrappedGr<
  //           gr::Solutions::Minkowski<volume_dim>>>;
  using initial_data_list = gh::fe::sgb::AnalyticData::all_analytic_data;
  using factory_classes = tmpl::map<
      tmpl::pair<DenseTrigger, DenseTriggers::standard_dense_triggers>,
      tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
      tmpl::pair<
          Event,
          tmpl::flatten<tmpl::list<
              Events::Completion, Events::MonitorMemory<volume_dim>,
              typename detail::ObserverTags<volume_dim>::field_observations,
              Events::time_events<system>>>>,
      //   tmpl::pair<gh::BoundaryConditions::BoundaryCondition<
      //                  volume_dim>,
      //              gh::BoundaryConditions::
      //                  standard_boundary_conditions<volume_dim>>,
      //
      tmpl::pair<fe::sgb::BoundaryConditions::BoundaryCondition,
                 fe::sgb::BoundaryConditions::standard_boundary_conditions>,
      //
      tmpl::pair<gh::gauges::GaugeCondition, gh::gauges::all_gauges>,
      tmpl::pair<evolution::initial_data::InitialData,
                 //  gh::Solutions::all_solutions<volume_dim>
                 //  initial_data_list
                 tmpl::conditional_t<UseNumericalInitialData,
                                     tmpl::list<gh::NumericInitialData>,
                                     initial_data_list>>,
      tmpl::pair<LtsTimeStepper, TimeSteppers::lts_time_steppers>,
      //   tmpl::pair<PhaseChange,
      //              tmpl::list<PhaseControl::VisitAndReturn<
      //                             Parallel::Phase::LoadBalancing>,
      //                        PhaseControl::CheckpointAndExitAfterWallclock>>,
      tmpl::pair<PhaseChange, PhaseControl::factory_creatable_classes>,
      tmpl::pair<StepChooser<StepChooserUse::LtsStep>,
                 StepChoosers::standard_step_choosers<system>>,
      tmpl::pair<
          StepChooser<StepChooserUse::Slab>,
          StepChoosers::standard_slab_choosers<system, LocalTimeStepping>>,
      tmpl::pair<TimeSequence<double>,
                 TimeSequences::all_time_sequences<double>>,
      tmpl::pair<TimeSequence<std::uint64_t>,
                 TimeSequences::all_time_sequences<std::uint64_t>>,
      tmpl::pair<TimeStepper, TimeSteppers::time_steppers>,
      tmpl::pair<Trigger, tmpl::append<Triggers::logical_triggers,
                                       Triggers::time_triggers>>>;
};
}  // namespace detail

template <template <size_t, bool> class EvolutionMetavarsDerived,
          size_t VolumeDim, bool UseNumericalInitialData>
struct FixedScalarTensorTemplateBase<
    EvolutionMetavarsDerived<VolumeDim, UseNumericalInitialData>> {
  //   using derived_metavars =
  //       EvolutionMetavarsDerived<VolumeDim, UseNumericalInitialData>;
  //   static constexpr size_t volume_dim = VolumeDim;
  using derived_metavars =
      EvolutionMetavarsDerived<3_st, UseNumericalInitialData>;
  static constexpr size_t volume_dim = 3_st;
  //   using system = gh::System<volume_dim>;
  //   using system_scalar = CurvedScalarWave::System<volume_dim>;
  //   using system_combined = ScalarTensor::System;
  using system = fe::sgb::System;
  using TimeStepperBase = LtsTimeStepper;

  static constexpr bool local_time_stepping =
      TimeStepperBase::local_time_stepping;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}

  //   using factory_creation =
  //       detail::FactoryCreation<volume_dim, local_time_stepping>;
  using factory_creation =
      detail::FactoryCreation<volume_dim, local_time_stepping,
                              UseNumericalInitialData>;

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::push_back<
          tmpl::at<typename factory_creation::factory_classes, Event>>>;

  using initialize_initial_data_dependent_quantities_actions = tmpl::list<
      // I think these are variables than can be retrieved from ID
      // gh::Actions::InitializeGhAnd3Plus1Variables<volume_dim>,
      Initialization::Actions::AddComputeTags<
          fe::sgb::Initialization::scalar_tensor_3plus1_compute_tags<
              volume_dim>>,
      tmpl::conditional_t<
          UseNumericalInitialData,
          // Until we read numerical data for the scalar
          // we set them to some analytical profile given some numerical data
          // for the metric quantities
          //   Initialization::Actions::AddSimpleTags<
          //       ScalarTensor::Actions::InitializeEvolvedScalarVariables>,
          tmpl::list<>, tmpl::list<>>,
      Actions::MutateApply<gh::gauges::SetPiAndPhiFromConstraints<volume_dim>>,
      // Initialization::Actions::GrTagsForHydro<system>,
      //   Initialization::Actions::AddSimpleTags<
      //       ScalarTensor::Initialization::
      //           InitializeConstraintDampingGammasGaussian>,
      //   Initialization::Actions::AddSimpleTags<
      //       fe::ScalarDriver::Initialization::
      //           InitializeConstraintDampingGammasGaussian>,
      Parallel::Actions::TerminatePhase>;

  // A tmpl::list of tags to be added to the GlobalCache by the
  // metavariables
  using const_global_cache_tags = tmpl::list<
      gh::gauges::Tags::GaugeCondition,
      evolution::initial_data::Tags::InitialData,
      gh::ConstraintDamping::Tags::DampingFunctionGamma0<volume_dim,
                                                         Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma1<volume_dim,
                                                         Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma2<volume_dim,
                                                         Frame::Grid>,
      ScalarTensor::ConstraintDamping::Tags::DampingFunctionGamma1<volume_dim,
                                                                   Frame::Grid>,
      ScalarTensor::ConstraintDamping::Tags::DampingFunctionGamma2<volume_dim,
                                                                   Frame::Grid>,
      //   fe::DecoupledScalar::ConstraintDamping::Tags::DampingFunctionGamma1<
      //       volume_dim, Frame::Grid>,
      //   fe::DecoupledScalar::ConstraintDamping::Tags::DampingFunctionGamma2<
      //       volume_dim, Frame::Grid>,
      // Scalar driver parameters
      fe::sgb::ConstraintDamping::Tags::DampingFunctionSigmaParameter<
          volume_dim, Frame::Grid>,
      fe::sgb::ConstraintDamping::Tags::DampingFunctionTauParameter<
          volume_dim, Frame::Grid>,
      // Source parameters
      ScalarTensor::Tags::ScalarMass,
      ScalarTensor::Tags::ScalarFirstCouplingParameter,
      ScalarTensor::Tags::ScalarSecondCouplingParameter>;
  //   ,
  // Scalar driver parameters
  //   fe::ScalarDriver::Tags::ScalarSigmaParameter,
  //       fe::ScalarDriver::Tags::ScalarTauParameter
  //   ,
  //  fe::ScalarDriver::Tags::DriverLimiterParameter,
  // Constraint damping
  //   fe::ScalarDriver::Tags::AmplitudeConstraintGamma2,
  //   fe::ScalarDriver::Tags::SigmaConstraintGamma2,
  //   fe::ScalarDriver::Tags::OffsetConstraintGamma2,
  //   ScalarTensor::Tags::AmplitudeConstraintGamma2,
  //   ScalarTensor::Tags::SigmaConstraintGamma2,
  //   ScalarTensor::Tags::OffsetConstraintGamma2
  //   >;

  using dg_registration_list =
      tmpl::list<observers::Actions::RegisterEventsWithObservers>;

  static constexpr auto default_phase_order =
      detail::make_default_phase_order<UseNumericalInitialData>();

  template <typename ControlSystems>
  using step_actions = tmpl::list<
      evolution::dg::Actions::ComputeTimeDerivative<
          volume_dim, system, AllStepChoosers, local_time_stepping>,
      tmpl::conditional_t<
          local_time_stepping,
          tmpl::list<evolution::Actions::RunEventsAndDenseTriggers<tmpl::list<
                         ::domain::CheckFunctionsOfTimeAreReadyPostprocessor,
                         evolution::dg::ApplyBoundaryCorrections<
                             local_time_stepping, system, volume_dim, true>>>,
                     evolution::dg::Actions::ApplyLtsBoundaryCorrections<
                         system, volume_dim, false>,
                     // We allow for separate filtering of the system variables
                     dg::Actions::Filter<Filters::Exponential<0>,
                                         system::gh_system::gh_system::
                                             variables_tag::tags_list>,
                     dg::Actions::Filter<Filters::Exponential<1>,
                                         system::gh_system::scalar_system::
                                             variables_tag::tags_list>,
                     dg::Actions::Filter<
                         Filters::Exponential<2>,
                         system::scalar_system::variables_tag::tags_list>>
          // tmpl::list<>
          ,
          tmpl::list<
              evolution::dg::Actions::ApplyBoundaryCorrectionsToTimeDerivative<
                  system, volume_dim, false>,
              Actions::RecordTimeStepperData<system>,
              evolution::Actions::RunEventsAndDenseTriggers<tmpl::list<>>,
              control_system::Actions::LimitTimeStep<ControlSystems>,
              Actions::UpdateU<system>,
              // We allow for separate filtering of the system variables
              dg::Actions::Filter<
                  Filters::Exponential<0>,
                  system::gh_system::gh_system::variables_tag::tags_list>,
              dg::Actions::Filter<
                  Filters::Exponential<1>,
                  system::gh_system::scalar_system::variables_tag::tags_list>,
              dg::Actions::Filter<
                  Filters::Exponential<2>,
                  system::scalar_system::variables_tag::tags_list>>
          // tmpl::list<>
          >>;

  // For labeling the yaml option for RandomizeVariables
  struct RandomizeInitialGuess {};

  template <bool UseControlSystems>
  using initialization_actions = tmpl::list<
      Initialization::Actions::InitializeItems<
          Initialization::TimeStepping<derived_metavars, TimeStepperBase>,
          evolution::dg::Initialization::Domain<volume_dim, UseControlSystems>,
          Initialization::TimeStepperHistory<derived_metavars>>,
      Initialization::Actions::NonconservativeSystem<system>,
      //
      //   Initialization::Actions::NonconservativeSystem<system_scalar>,
      //
      std::conditional_t<
          UseNumericalInitialData, tmpl::list<>,
          evolution::Initialization::Actions::SetVariables<
              domain::Tags::Coordinates<volume_dim, Frame::ElementLogical>>>,
      // Random noise system::variables_tag
      //   Actions::RandomizeVariables<typename system::variables_tag,
      //                               RandomizeInitialGuess>,
      Initialization::Actions::AddComputeTags<::Tags::DerivCompute<
          typename system::variables_tag, domain::Tags::Mesh<volume_dim>,
          domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                        Frame::Inertial>,
          typename system::gradient_variables>>,
      // gh::Actions::InitializeGhAnd3Plus1Variables<volume_dim>,
      Initialization::Actions::AddComputeTags<
          tmpl::push_back<StepChoosers::step_chooser_compute_tags<
              FixedScalarTensorTemplateBase, local_time_stepping>>>,
      ::evolution::dg::Initialization::Mortars<volume_dim, system>,
      evolution::Actions::InitializeRunEventsAndDenseTriggers,
      Parallel::Actions::TerminatePhase>;
};
