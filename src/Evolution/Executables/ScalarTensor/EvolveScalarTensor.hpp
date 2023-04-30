// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <vector>

#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
// #include
// "Evolution/Executables/GeneralizedHarmonic/GeneralizedHarmonicBase.hpp"
#include "Evolution/Executables/ScalarTensor/ScalarTensorBase.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/RegisterDerived.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/RegisterDerivedWithCharm.hpp"
//
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/RegisterDerived.hpp"
//
#include "Options/Options.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/PhaseControl/PhaseControlTags.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/GaugeWave.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"

template <size_t VolumeDim, bool UseNumericalInitialData>
struct EvolutionMetavars
    : public ScalarTensorTemplateBase<
          EvolutionMetavars<3_st, UseNumericalInitialData>> {
//   using gh_base = GeneralizedHarmonicTemplateBase<
//       EvolutionMetavars<VolumeDim, UseNumericalInitialData>>;
  using st_base = ScalarTensorTemplateBase<
      EvolutionMetavars<3_st, UseNumericalInitialData>>;
  using typename st_base::const_global_cache_tags;
  using typename st_base::dg_registration_list;
  using typename st_base::initialization_actions;
  using typename st_base::initialize_initial_data_dependent_quantities_actions;
  using typename st_base::observed_reduction_data_tags;
  using typename st_base::step_actions;
  using typename st_base::system;

  using st_dg_element_array = DgElementArray<
      EvolutionMetavars,
      tmpl::flatten<tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,
          tmpl::conditional_t<UseNumericalInitialData, tmpl::list<>,
              tmpl::list<>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeInitialDataDependentQuantities,
              initialize_initial_data_dependent_quantities_actions>,
          Parallel::PhaseActions<Parallel::Phase::InitializeTimeStepperHistory,
                                 SelfStart::self_start_procedure<
                                       step_actions,
                                    //  tmpl::list<>,
                                      system>>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              //   tmpl::list<>
              tmpl::list<
                         Actions::RunEventsAndTriggers,
                         Actions::ChangeSlabSize,
                         step_actions,
                         Actions::AdvanceTime,
                         PhaseControl::Actions::ExecutePhaseChange
                         >
                         >
                         >
                         >
                         >;

  template <typename ParallelComponent>
  struct registration_list {
    using type = std::conditional_t<
        std::is_same_v<ParallelComponent, st_dg_element_array>,
        dg_registration_list, tmpl::list<>>;
  };

  using component_list = tmpl::flatten<tmpl::list<
      observers::Observer<EvolutionMetavars>,
      observers::ObserverWriter<EvolutionMetavars>,
      std::conditional_t<UseNumericalInitialData, tmpl::list<>, tmpl::list<>>,
      st_dg_element_array>>;

  static constexpr Options::String help{
      "Evolve the Einstein field equations in GH gauge coupled to a scalar "
      "field \n"};
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling,
    &setup_memory_allocation_failure_reporting,
    &disable_openblas_multithreading,
    &domain::creators::time_dependence::register_derived_with_charm,
    &domain::FunctionsOfTime::register_derived_with_charm,
    // &GeneralizedHarmonic::BoundaryCorrections::register_derived_with_charm,
    &ScalarTensor::BoundaryCorrections::register_derived_with_charm,
    &domain::creators::register_derived_with_charm,
    &GeneralizedHarmonic::ConstraintDamping::register_derived_with_charm,
    &Parallel::register_factory_classes_with_charm<metavariables>};

static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
