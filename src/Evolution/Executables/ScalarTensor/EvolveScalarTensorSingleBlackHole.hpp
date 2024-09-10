// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "ControlSystem/Actions/InitializeMeasurements.hpp"
#include "ControlSystem/Actions/PrintCurrentMeasurement.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Size/Factory.hpp"
#include "ControlSystem/ControlErrors/Size/State.hpp"
#include "ControlSystem/Measurements/SingleHorizon.hpp"
#include "ControlSystem/Metafunctions.hpp"
#include "ControlSystem/Systems/Shape.hpp"
#include "ControlSystem/Systems/Size.hpp"
#include "ControlSystem/Trigger.hpp"
#include "Domain/FunctionsOfTime/OutputTimeBounds.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Evolution/Actions/RunEventsAndTriggers.hpp"
#include "Evolution/Deadlock/PrintDgElementArray.hpp"
#include "Evolution/DiscontinuousGalerkin/InboxTags.hpp"
#include "Evolution/Executables/ScalarTensor/ScalarTensorBase.hpp"
#include "Evolution/Systems/Cce/Callbacks/DumpBondiSachsOnWorldtube.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/SetInitialData.hpp"
#include "Evolution/Systems/ScalarTensor/Actions/SetInitialData.hpp"
#include "Options/FactoryHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "Parallel/Protocols/RegistrationMetavariables.hpp"
#include "ParallelAlgorithms/Actions/FunctionsOfTimeAreReady.hpp"
#include "ParallelAlgorithms/Amr/Projectors/CopyFromCreatorOrLeaveAsIs.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ErrorOnFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/IgnoreFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeExcisionBoundaryVolumeQuantities.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeExcisionBoundaryVolumeQuantities.tpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeHorizonVolumeQuantities.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeHorizonVolumeQuantities.tpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/InterpolationTarget.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ObserveCenters.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/CleanUpInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/ElementInitInterpPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolationTargetReceiveVars.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceivePoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceiveVolumeData.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorRegisterElement.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/TryToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveSurfaceData.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveYlms.hpp"
#include "ParallelAlgorithms/Interpolation/Events/Interpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Interpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/ChangeSlabSize/Action.hpp"
#include "Time/ChangeSlabSize/Tags.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/Tags/StepperErrors.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"

struct EvolutionMetavars : public ScalarTensorTemplateBase<EvolutionMetavars> {
  using st_base = ScalarTensorTemplateBase<EvolutionMetavars>;
  using typename st_base::initialize_initial_data_dependent_quantities_actions;
  using typename st_base::system;

  static constexpr size_t volume_dim = 3_st;

  static constexpr Options::String help{
      "Evolve the Einstein field equations in GH gauge coupled to a scalar "
      "field \n"
      "on a domain with a single horizon and corresponding excised region"};

  template <typename Frame>
  struct Ah : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;
    using tags_to_observe = ::ah::tags_for_observing<Frame>;
    using surface_tags_to_observe = ::ah::surface_tags_for_observing;
    using compute_vars_to_interpolate = ah::ComputeHorizonVolumeQuantities;
    using vars_to_interpolate_to_target =
        ::ah::vars_to_interpolate_to_target<volume_dim, Frame>;
    using compute_items_on_target =
        ::ah::compute_items_on_target<volume_dim, Frame>;
    using compute_target_points =
        intrp::TargetPoints::ApparentHorizon<Ah, Frame>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::FindApparentHorizon<Ah, Frame>>;
    using horizon_find_failure_callback =
        intrp::callbacks::IgnoreFailedApparentHorizon;
    using post_horizon_find_callbacks = tmpl::list<
        intrp::callbacks::ObserveTimeSeriesOnSurface<tags_to_observe, Ah>,
        intrp::callbacks::ObserveSurfaceData<surface_tags_to_observe, Ah,
                                             Frame>,
        // Needs to be Frame::Grid or Frame::Distorted
        ::ah::callbacks::ObserveCenters<Ah, ::Frame::Distorted>>;
  };

  using AhA = Ah<::Frame::Distorted>;

  struct ExcisionBoundaryA
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;
    using tags_to_observe =
        tmpl::list<gr::Tags::Lapse<DataVector>,
                   gr::Tags::Shift<DataVector, 3, Frame::Grid>>;
    using compute_vars_to_interpolate =
        ah::ComputeExcisionBoundaryVolumeQuantities;
    using vars_to_interpolate_to_target = tags_to_observe;
    using compute_items_on_source = tmpl::list<>;
    using compute_items_on_target = tmpl::list<>;
    using compute_target_points =
        intrp::TargetPoints::Sphere<ExcisionBoundaryA, ::Frame::Grid>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::ObserveSurfaceData<
            tags_to_observe, ExcisionBoundaryA, ::Frame::Grid>>;
    // run_callbacks
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  template <size_t SphereNumber>
  struct SphericalSurfaceTmp
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    using vars_to_interpolate_to_target =
        detail::ObserverTags::scalar_charge_vars_to_interpolate_to_target;
    using compute_items_on_target =
        detail::ObserverTags::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurfaceTmp<SphereNumber>,
                                    ::Frame::Inertial>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::ObserveTimeSeriesOnSurface<
                       detail::ObserverTags::scalar_charge_surface_obs_tags,
                       SphericalSurfaceTmp<SphereNumber>>,
                   intrp::callbacks::ObserveYlms<
                       CurvedScalarWave::Tags::Psi,
                       SphericalSurfaceTmp<SphereNumber>, ::Frame::Inertial>>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
    static std::string name() {
      return "SphericalSurface" + std::to_string(SphereNumber);
    }
  };

  using SphericalSurface = SphericalSurfaceTmp<1>;
  using SphericalSurface2 = SphericalSurfaceTmp<2>;
  using SphericalSurface3 = SphericalSurfaceTmp<3>;
  using SphericalSurface4 = SphericalSurfaceTmp<4>;
  using SphericalSurface5 = SphericalSurfaceTmp<5>;
  using SphericalSurface6 = SphericalSurfaceTmp<6>;

  using control_systems =
      tmpl::list<control_system::Systems::Shape<
                     ::domain::ObjectLabel::None, 2,
                     control_system::measurements::SingleHorizon<
                         ::domain::ObjectLabel::None>>,
                 control_system::Systems::Size<::domain::ObjectLabel::None, 2>>;

  static constexpr bool use_control_systems =
      tmpl::size<control_systems>::value > 0;

  struct BondiSachs;

  using interpolation_target_tags = tmpl::push_back<
      control_system::metafunctions::interpolation_target_tags<control_systems>,
      AhA, ExcisionBoundaryA, SphericalSurface, SphericalSurface2,
      SphericalSurface3, SphericalSurface4, SphericalSurface5,
      SphericalSurface6, BondiSachs>;
  using interpolator_source_vars = ::ah::source_vars<volume_dim>;
  using source_vars_no_deriv =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, volume_dim>,
                 gh::Tags::Pi<DataVector, volume_dim>,
                 gh::Tags::Phi<DataVector, volume_dim>>;

  using scalar_charge_interpolator_source_vars =
      detail::ObserverTags::scalar_charge_vars_to_interpolate_to_target;

  struct BondiSachs : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    static std::string name() { return "BondiSachsInterpolation"; }
    using temporal_id = ::Tags::Time;
    using vars_to_interpolate_to_target = source_vars_no_deriv;
    using compute_target_points =
        intrp::TargetPoints::Sphere<BondiSachs, ::Frame::Inertial>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::DumpBondiSachsOnWorldtube<BondiSachs>>;
    using compute_items_on_target = tmpl::list<>;
    template <typename Metavariables>
    using interpolating_component = typename Metavariables::st_dg_element_array;
  };

  // The interpolator_source_vars need to be the same in both the
  // Interpolate event and the InterpolateWithoutInterpComponent event.  The
  // Interpolate event interpolates to the horizon, and the
  // InterpolateWithoutInterpComponent event interpolates to the excision
  // boundary. Every Target gets the same interpolator_source_vars, so they need
  // to be made the same. Otherwise a static assert is triggered.
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = Options::add_factory_classes<
        // Restrict to monotonic time steppers in LTS to avoid control
        // systems deadlocking.
        tmpl::insert<
            tmpl::erase<typename st_base::factory_creation::factory_classes,
                        LtsTimeStepper>,
            tmpl::pair<LtsTimeStepper,
                       TimeSteppers::monotonic_lts_time_steppers>>,
        tmpl::pair<
            Event,
            tmpl::flatten<tmpl::list<
                intrp::Events::Interpolate<volume_dim, AhA,
                                           interpolator_source_vars>,
                control_system::metafunctions::control_system_events<
                    control_systems>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, BondiSachs, source_vars_no_deriv>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, ExcisionBoundaryA, interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, SphericalSurface,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, SphericalSurface2,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, SphericalSurface3,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, SphericalSurface4,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, SphericalSurface5,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, SphericalSurface6,
                    scalar_charge_interpolator_source_vars>>>>,
        tmpl::pair<DenseTrigger,
                   control_system::control_system_triggers<control_systems>>,
        tmpl::pair<control_system::size::State,
                   control_system::size::States::factory_creatable_states>>;
  };

  using typename st_base::const_global_cache_tags;

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::push_back<
          tmpl::at<typename factory_creation::factory_classes, Event>,
          typename AhA::post_horizon_find_callbacks>>;

  using dg_registration_list =
      tmpl::push_back<typename st_base::dg_registration_list,
                      intrp::Actions::RegisterElementWithInterpolator>;

  using step_actions = typename st_base::template step_actions<control_systems>;

  using initialization_actions = tmpl::push_back<
      tmpl::pop_back<typename st_base::template initialization_actions<
          use_control_systems>>,
      control_system::Actions::InitializeMeasurements<control_systems>,
      intrp::Actions::ElementInitInterpPoints<
          intrp::Tags::InterpPointInfo<EvolutionMetavars>>,
      tmpl::back<typename st_base::template initialization_actions<
          use_control_systems>>>;

  using st_dg_element_array = DgElementArray<
      EvolutionMetavars,
      tmpl::flatten<tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::RegisterWithElementDataReader,
              tmpl::list<importers::Actions::RegisterWithElementDataReader,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::ImportInitialData,
              tmpl::list<ScalarTensor::Actions::SetInitialData,
                         ScalarTensor::Actions::ReceiveNumericInitialData,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeInitialDataDependentQuantities,
              initialize_initial_data_dependent_quantities_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              SelfStart::self_start_procedure<step_actions, system>>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<Parallel::Phase::CheckDomain,
                                 tmpl::list<::amr::Actions::SendAmrDiagnostics,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::list<
                  ::domain::Actions::CheckFunctionsOfTimeAreReady<volume_dim>,
                  evolution::Actions::RunEventsAndTriggers,
                  Actions::ChangeSlabSize, step_actions, Actions::AdvanceTime,
                  PhaseControl::Actions::ExecutePhaseChange>>>>>;

  using gh_dg_element_array = st_dg_element_array;

  struct amr : tt::ConformsTo<::amr::protocols::AmrMetavariables> {
    using element_array = st_dg_element_array;

    using projectors = tmpl::list<
        Initialization::ProjectTimeStepping<volume_dim>,
        evolution::dg::Initialization::ProjectDomain<volume_dim>,
        Initialization::ProjectTimeStepperHistory<EvolutionMetavars>,
        ::amr::projectors::ProjectVariables<volume_dim,
                                            typename system::variables_tag>,
        evolution::dg::Initialization::ProjectMortars<EvolutionMetavars>,
        evolution::Actions::ProjectRunEventsAndDenseTriggers,
        ::amr::projectors::DefaultInitialize<
            Initialization::Tags::InitialTimeDelta,
            Initialization::Tags::InitialSlabSize<st_base::local_time_stepping>,
            ::domain::Tags::InitialExtents<volume_dim>,
            ::domain::Tags::InitialRefinementLevels<volume_dim>,
            evolution::dg::Tags::Quadrature,
            Tags::StepperErrors<typename system::variables_tag>,
            SelfStart::Tags::InitialValue<typename system::variables_tag>,
            SelfStart::Tags::InitialValue<Tags::TimeStep>,
            SelfStart::Tags::InitialValue<Tags::Next<Tags::TimeStep>>,
            evolution::dg::Tags::BoundaryData<volume_dim>>,
        ::amr::projectors::CopyFromCreatorOrLeaveAsIs<tmpl::push_back<
            typename control_system::Actions::InitializeMeasurements<
                control_systems>::simple_tags,
            intrp::Tags::InterpPointInfo<EvolutionMetavars>,
            Tags::ChangeSlabSize::NumberOfExpectedMessages,
            Tags::ChangeSlabSize::NewSlabSize>>>;
  };

  //   template <typename ParallelComponent>
  //   struct registration_list {
  //     using type = std::conditional_t<
  //         std::is_same_v<ParallelComponent, st_dg_element_array>,
  //         dg_registration_list, tmpl::list<>>;
  //   };

  struct registration
      : tt::ConformsTo<Parallel::protocols::RegistrationMetavariables> {
    using element_registrars =
        tmpl::map<tmpl::pair<st_dg_element_array, dg_registration_list>>;
  };

  using control_components =
      control_system::control_components<EvolutionMetavars, control_systems>;

  static void run_deadlock_analysis_simple_actions(
      Parallel::GlobalCache<EvolutionMetavars>& cache,
      const std::vector<std::string>& deadlocked_components) {
    const auto& functions_of_time =
        Parallel::get<::domain::Tags::FunctionsOfTime>(cache);

    const std::string time_bounds =
        ::domain::FunctionsOfTime::output_time_bounds(functions_of_time);

    Parallel::printf("%s\n", time_bounds);

    if (alg::count(deadlocked_components,
                   pretty_type::name<st_dg_element_array>()) == 1) {
      tmpl::for_each<control_components>([&cache](auto component_v) {
        using component = tmpl::type_from<decltype(component_v)>;
        Parallel::simple_action<
            control_system::Actions::PrintCurrentMeasurement>(
            Parallel::get_parallel_component<component>(cache));
      });

      if constexpr (Parallel::is_dg_element_collection_v<gh_dg_element_array>) {
        Parallel::threaded_action<Parallel::Actions::SimpleActionOnElement<
            deadlock::PrintElementInfo, true>>(
            Parallel::get_parallel_component<gh_dg_element_array>(cache));
      } else {
        Parallel::simple_action<deadlock::PrintElementInfo>(
            Parallel::get_parallel_component<gh_dg_element_array>(cache));
      }
    }
  }

  using component_list = tmpl::flatten<tmpl::list<
      ::amr::Component<EvolutionMetavars>,
      observers::Observer<EvolutionMetavars>,
      observers::ObserverWriter<EvolutionMetavars>,
      mem_monitor::MemoryMonitor<EvolutionMetavars>,
      importers::ElementDataReader<EvolutionMetavars>, st_dg_element_array,
      intrp::Interpolator<EvolutionMetavars>,
      control_system::control_components<EvolutionMetavars, control_systems>,
      tmpl::transform<interpolation_target_tags,
                      tmpl::bind<intrp::InterpolationTarget,
                                 tmpl::pin<EvolutionMetavars>, tmpl::_1>>>>;
};
