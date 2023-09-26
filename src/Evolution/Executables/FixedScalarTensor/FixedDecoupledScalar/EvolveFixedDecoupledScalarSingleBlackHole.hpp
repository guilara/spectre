// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <vector>

#include "ControlSystem/Actions/InitializeMeasurements.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/Event.hpp"
#include "ControlSystem/Measurements/SingleHorizon.hpp"
#include "ControlSystem/Systems/Shape.hpp"
#include "ControlSystem/Trigger.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/FunctionsOfTimeAreReady.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Evolution/Actions/RunEventsAndTriggers.hpp"
// #include
// "Evolution/Executables/GeneralizedHarmonic/GeneralizedHarmonicBase.hpp"
// #include "Evolution/Executables/ScalarTensor/ScalarTensorBase.hpp"
//
#include "Evolution/Systems/GeneralizedHarmonic/Actions/SetInitialData.hpp"
//
// #include "Evolution/Systems/ScalarTensor/Actions/NumericInitialData.hpp"
#include "Evolution/Systems/ScalarTensor/Actions/SetInitialData.hpp"
//
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/RegisterDerived.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/RegisterDerivedWithCharm.hpp"
//
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/RegisterDerived.hpp"
//
#include "Evolution/Executables/FixedScalarTensor/FixedDecoupledScalar/FixedScalarTensorBase.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/RegisterDerived.hpp"
//
#include "Options/FactoryHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeExcisionBoundaryVolumeQuantities.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeExcisionBoundaryVolumeQuantities.tpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeHorizonVolumeQuantities.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeHorizonVolumeQuantities.tpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/CleanUpInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/ElementInitInterpPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolationTargetReceiveVars.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceivePoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceiveVolumeData.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorRegisterElement.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/TryToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ErrorOnFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/IgnoreFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveSurfaceData.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Events/Interpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Interpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/ErrorHandling/SegfaultHandler.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

template <size_t VolumeDim, bool UseNumericalInitialData>
// struct EvolutionMetavars
//     : public ScalarTensorTemplateBase<
//           EvolutionMetavars<VolumeDim, UseNumericalInitialData>> {
struct EvolutionMetavars
    : public FixedScalarTensorTemplateBase<
          EvolutionMetavars<3_st, UseNumericalInitialData>> {
  using st_base = FixedScalarTensorTemplateBase<EvolutionMetavars>;
  using typename st_base::initialize_initial_data_dependent_quantities_actions;
  using typename st_base::system;
  //   static constexpr size_t volume_dim = VolumeDim;
  static constexpr size_t volume_dim = 3_st;

  static constexpr Options::String help{
      "Evolve the Einstein field equations in GH gauge coupled to a scalar "
      "field with an extra scalar driver \n"
      "on a domain with a single horizon and corresponding excised region"};

  struct AhA : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;
    using tags_to_observe = ::ah::tags_for_observing<Frame::Inertial>;
    using surface_tags_to_observe = ::ah::surface_tags_for_observing;
    using compute_vars_to_interpolate = ah::ComputeHorizonVolumeQuantities;
    using vars_to_interpolate_to_target =
        ::ah::vars_to_interpolate_to_target<volume_dim, ::Frame::Inertial>;
    using compute_items_on_target =
        ::ah::compute_items_on_target<volume_dim, Frame::Inertial>;
    using compute_target_points =
        intrp::TargetPoints::ApparentHorizon<AhA, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::FindApparentHorizon<AhA, ::Frame::Inertial>;
    using horizon_find_failure_callback =
        intrp::callbacks::IgnoreFailedApparentHorizon;
    using post_horizon_find_callbacks = tmpl::list<
        intrp::callbacks::ObserveTimeSeriesOnSurface<tags_to_observe, AhA>,
        intrp::callbacks::ObserveSurfaceData<surface_tags_to_observe, AhA,
                                             ::Frame::Inertial>>;
  };

  struct ExcisionBoundaryA
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;
    using tags_to_observe =
        tmpl::list<gr::Tags::Lapse<DataVector>,
                   gh::ConstraintDamping::Tags::ConstraintGamma1,
                   gh::CharacteristicSpeedsOnStrahlkorper<Frame::Grid>>;
    using compute_vars_to_interpolate =
        ah::ComputeExcisionBoundaryVolumeQuantities;
    using vars_to_interpolate_to_target =
        tmpl::list<gr::Tags::Lapse<DataVector>,
                   gr::Tags::Shift<DataVector, 3, Frame::Grid>,
                   gr::Tags::SpatialMetric<DataVector, 3, Frame::Grid>,
                   gh::ConstraintDamping::Tags::ConstraintGamma1>;
    using compute_items_on_source = tmpl::list<>;
    using compute_items_on_target = tmpl::append<tmpl::list<
        gr::Tags::DetAndInverseSpatialMetricCompute<DataVector, 3, Frame::Grid>,
        ylm::Tags::OneOverOneFormMagnitudeCompute<DataVector, 3, Frame::Grid>,
        ylm::Tags::UnitNormalOneFormCompute<Frame::Grid>,
        gh::CharacteristicSpeedsOnStrahlkorperCompute<3, Frame::Grid>>>;
    using compute_target_points =
        intrp::TargetPoints::Sphere<ExcisionBoundaryA, ::Frame::Grid>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveSurfaceData<tags_to_observe, ExcisionBoundaryA,
                                             ::Frame::Grid>;
    // run_callbacks
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  struct SphericalSurface
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    // Note: These need to be the same as in `interpolator_source_vars`.
    // For now, all interpolator targets in this executable need the same
    // tags here
    using vars_to_interpolate_to_target =
        detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;
    // Most of these tags are required to compute the unit normal
    using compute_items_on_target =
        detail::ObserverTags<3_st>::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurface, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<
            detail::ObserverTags<3_st>::scalar_charge_surface_obs_tags,
            SphericalSurface>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  struct SphericalSurface2
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    // Note: These need to be the same as in `interpolator_source_vars`.
    // For now, all interpolator targets in this executable need the same
    // tags here
    using vars_to_interpolate_to_target =
        detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;
    // Most of these tags are required to compute the unit normal
    using compute_items_on_target =
        detail::ObserverTags<3_st>::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurface2, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<
            detail::ObserverTags<3_st>::scalar_charge_surface_obs_tags,
            SphericalSurface2>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  struct SphericalSurface3
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    // Note: These need to be the same as in `interpolator_source_vars`.
    // For now, all interpolator targets in this executable need the same
    // tags here
    using vars_to_interpolate_to_target =
        detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;
    // Most of these tags are required to compute the unit normal
    using compute_items_on_target =
        detail::ObserverTags<3_st>::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurface3, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<
            detail::ObserverTags<3_st>::scalar_charge_surface_obs_tags,
            SphericalSurface3>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  struct SphericalSurface4
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    // Note: These need to be the same as in `interpolator_source_vars`.
    // For now, all interpolator targets in this executable need the same
    // tags here
    using vars_to_interpolate_to_target =
        detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;
    // Most of these tags are required to compute the unit normal
    using compute_items_on_target =
        detail::ObserverTags<3_st>::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurface4, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<
            detail::ObserverTags<3_st>::scalar_charge_surface_obs_tags,
            SphericalSurface4>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  struct SphericalSurface5
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    // Note: These need to be the same as in `interpolator_source_vars`.
    // For now, all interpolator targets in this executable need the same
    // tags here
    using vars_to_interpolate_to_target =
        detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;
    // Most of these tags are required to compute the unit normal
    using compute_items_on_target =
        detail::ObserverTags<3_st>::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurface5, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<
            detail::ObserverTags<3_st>::scalar_charge_surface_obs_tags,
            SphericalSurface5>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  struct SphericalSurface6
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;

    // Note: These need to be the same as in `interpolator_source_vars`.
    // For now, all interpolator targets in this executable need the same
    // tags here
    using vars_to_interpolate_to_target =
        detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;
    // Most of these tags are required to compute the unit normal
    using compute_items_on_target =
        detail::ObserverTags<3_st>::scalar_charge_compute_items_on_target;
    using compute_target_points =
        intrp::TargetPoints::Sphere<SphericalSurface6, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<
            detail::ObserverTags<3_st>::scalar_charge_surface_obs_tags,
            SphericalSurface6>;
    template <typename metavariables>
    using interpolating_component = typename metavariables::st_dg_element_array;
  };

  using control_systems = tmpl::list<control_system::Systems::Shape<
      ::domain::ObjectLabel::None, 2,
      control_system::measurements::SingleHorizon<
          ::domain::ObjectLabel::None>>>;

  static constexpr bool use_control_systems =
      tmpl::size<control_systems>::value > 0;

  //   using interpolation_target_tags = tmpl::list<AhA, ExcisionBoundaryA>;
  //   using interpolation_target_tags = tmpl::list<AhA>;
  //   using interpolation_target_tags = tmpl::list<AhA, SphericalSurface,
  //                                     SphericalSurface2, SphericalSurface3,
  //                                     SphericalSurface4, SphericalSurface5,
  //                                     SphericalSurface6>;
  //   using interpolator_source_vars_excision_boundary = tmpl::list<
  //       gr::Tags::SpacetimeMetric<volume_dim, Frame::Inertial>,
  //       gh::ConstraintDamping::Tags::ConstraintGamma1>;
  //   using interpolator_source_vars = tmpl::remove_duplicates<
  //       tmpl::append<interpolator_source_vars_excision_boundary,
  //                    ::ah::source_vars<volume_dim>>>;

  using interpolation_target_tags = tmpl::push_back<
      control_system::metafunctions::interpolation_target_tags<control_systems>,
      AhA, ExcisionBoundaryA, SphericalSurface, SphericalSurface2,
      SphericalSurface3, SphericalSurface4, SphericalSurface5,
      SphericalSurface6>;
  using interpolator_source_vars = ::ah::source_vars<volume_dim>;

  using scalar_charge_interpolator_source_vars =
      detail::ObserverTags<3_st>::scalar_charge_vars_to_interpolate_to_target;

  // The interpolator_source_vars need to be the same in both the
  // Interpolate event and the InterpolateWithoutInterpComponent event.  The
  // Interpolate event interpolates to the horizon, and the
  // InterpolateWithoutInterpComponent event interpolates to the excision
  // boundary. Every Target gets the same interpolator_source_vars, so they need
  // to be made the same. Otherwise a static assert is triggered.
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    // using factory_classes = Options::add_factory_classes<
    //     typename st_base::factory_creation::factory_classes,
    //     tmpl::pair<Event,
    //                tmpl::list<intrp::Events::Interpolate<
    //                               3, AhA, interpolator_source_vars>,
    //                         intrp::Events::InterpolateWithoutInterpComponent<
    //                               3, ExcisionBoundaryA, EvolutionMetavars,
    //                               interpolator_source_vars>>>>;
    using factory_classes = Options::add_factory_classes<
        typename st_base::factory_creation::factory_classes,
        tmpl::pair<
            Event,
            tmpl::flatten<tmpl::list<
                intrp::Events::Interpolate<3, AhA, interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, ExcisionBoundaryA, interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, SphericalSurface,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, SphericalSurface2,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, SphericalSurface3,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, SphericalSurface4,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, SphericalSurface5,
                    scalar_charge_interpolator_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, SphericalSurface6,
                    scalar_charge_interpolator_source_vars>>>>,
        tmpl::pair<DenseTrigger,
                   control_system::control_system_triggers<control_systems>>>;
  };

  using typename st_base::const_global_cache_tags;

  //   using observed_reduction_data_tags =
  //       observers::collect_reduction_data_tags<tmpl::push_back<
  //           tmpl::at<typename factory_creation::factory_classes, Event>,
  //           typename AhA::post_horizon_find_callbacks,
  //           typename ExcisionBoundaryA::post_interpolation_callback>>;
  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::push_back<
          tmpl::at<typename factory_creation::factory_classes, Event>,
          typename AhA::post_horizon_find_callbacks>>;

  using dg_registration_list =
      tmpl::push_back<typename st_base::dg_registration_list,
                      intrp::Actions::RegisterElementWithInterpolator>;

  using typename st_base::step_actions;

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
          tmpl::conditional_t<
              UseNumericalInitialData,
              //   tmpl::list<>,
              tmpl::list<
                  Parallel::PhaseActions<
                      Parallel::Phase::RegisterWithElementDataReader,
                      tmpl::list<
                          importers::Actions::RegisterWithElementDataReader,
                          Parallel::Actions::TerminatePhase>>,
                  Parallel::PhaseActions<
                      Parallel::Phase::ImportInitialData,
                      tmpl::list<
                          gh::Actions::SetInitialData,
                          gh::Actions::ReceiveNumericInitialData,
                          Parallel::Actions::TerminatePhase>>>,
              tmpl::list<>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeInitialDataDependentQuantities,
              initialize_initial_data_dependent_quantities_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              SelfStart::self_start_procedure<step_actions,
                                              //  tmpl::list<>,
                                              system>>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::list<::domain::Actions::CheckFunctionsOfTimeAreReady,
                         evolution::Actions::RunEventsAndTriggers,
                         Actions::ChangeSlabSize, step_actions,
                         Actions::AdvanceTime,
                         PhaseControl::Actions::ExecutePhaseChange>>>>>;

  template <typename ParallelComponent>
  struct registration_list {
    using type = std::conditional_t<
        std::is_same_v<ParallelComponent, st_dg_element_array>,
        dg_registration_list, tmpl::list<>>;
  };

  //   using component_list = tmpl::flatten<tmpl::list<
  //       observers::Observer<EvolutionMetavars>,
  //       observers::ObserverWriter<EvolutionMetavars>,
  //       std::conditional_t<UseNumericalInitialData, tmpl::list<>,
  //                          //
  //                          importers::ElementDataReader<EvolutionMetavars>,
  //                          tmpl::list<>>,
  //       st_dg_element_array, intrp::Interpolator<EvolutionMetavars>,
  //       intrp::InterpolationTarget<EvolutionMetavars, AhA>,
  //       intrp::InterpolationTarget<EvolutionMetavars, ExcisionBoundaryA>>>;
  using component_list = tmpl::flatten<tmpl::list<
      observers::Observer<EvolutionMetavars>,
      observers::ObserverWriter<EvolutionMetavars>,
      mem_monitor::MemoryMonitor<EvolutionMetavars>,
      std::conditional_t<UseNumericalInitialData,
                         //  tmpl::list<>,
                         importers::ElementDataReader<EvolutionMetavars>,
                         tmpl::list<>>,
      st_dg_element_array, intrp::Interpolator<EvolutionMetavars>,
      control_system::control_components<EvolutionMetavars, control_systems>,
      tmpl::transform<interpolation_target_tags,
                      tmpl::bind<intrp::InterpolationTarget,
                                 tmpl::pin<EvolutionMetavars>, tmpl::_1>>
      //   intrp::InterpolationTarget<EvolutionMetavars, AhA>,
      //   intrp::InterpolationTarget<EvolutionMetavars, SphericalSurface>,
      //   intrp::InterpolationTarget<EvolutionMetavars, SphericalSurface2>,
      //   intrp::InterpolationTarget<EvolutionMetavars, SphericalSurface3>,
      //   intrp::InterpolationTarget<EvolutionMetavars, SphericalSurface4>,
      //   intrp::InterpolationTarget<EvolutionMetavars, SphericalSurface5>,
      //   intrp::InterpolationTarget<EvolutionMetavars, SphericalSurface6>
      >>;
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling, &disable_openblas_multithreading,
    &domain::creators::time_dependence::register_derived_with_charm,
    &domain::FunctionsOfTime::register_derived_with_charm,
    // &gh::BoundaryCorrections::register_derived_with_charm,
    // &ScalarTensor::BoundaryCorrections::register_derived_with_charm,
    &fe::DecoupledScalar::BoundaryCorrections::register_derived_with_charm,
    &domain::creators::register_derived_with_charm,
    &gh::ConstraintDamping::register_derived_with_charm,
    &register_factory_classes_with_charm<metavariables>};

static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions, &enable_segfault_handler};
