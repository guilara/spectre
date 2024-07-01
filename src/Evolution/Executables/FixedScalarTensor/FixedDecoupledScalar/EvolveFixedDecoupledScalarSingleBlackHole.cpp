// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Executables/FixedScalarTensor/FixedDecoupledScalar/EvolveFixedDecoupledScalarSingleBlackHole.hpp"

#include <vector>

#include "ControlSystem/ControlErrors/Size/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/RegisterDerived.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/RegisterDerivedWithCharm.hpp"
#include "Parallel/CharmMain.tpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

// Parameters chosen in CMakeLists.txt
using metavariables = EvolutionMetavars<DIM, USE_NUMERICAL_ID>;

extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<metavariables>();
  Parallel::charmxx::register_init_node_and_proc(
      {&domain::creators::register_derived_with_charm,
       &domain::creators::time_dependence::register_derived_with_charm,
       &domain::FunctionsOfTime::register_derived_with_charm,
       &fe::DecoupledScalar::BoundaryCorrections::register_derived_with_charm,
       &gh::ConstraintDamping::register_derived_with_charm,
       &control_system::size::register_derived_with_charm,
       &register_factory_classes_with_charm<metavariables>},
      {});
}
