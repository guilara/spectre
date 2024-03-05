// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Executables/FixedScalarTensor/FixedDecoupledScalar/EvolveFixedDecoupledScalarBinaryBlackHole.hpp"

#include <vector>

#include "ControlSystem/ControlErrors/Size/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/RegisterDerived.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/ConstraintDamping/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/ScalarTensor/ConstraintDamping/RegisterDerivedWithCharm.hpp"
#include "Parallel/CharmMain.tpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<EvolutionMetavars>();
  Parallel::charmxx::register_init_node_and_proc(
      {&domain::creators::register_derived_with_charm,
       &domain::creators::time_dependence::register_derived_with_charm,
       &domain::FunctionsOfTime::register_derived_with_charm,
       &fe::DecoupledScalar::BoundaryCorrections::register_derived_with_charm,
       &gh::ConstraintDamping::register_derived_with_charm,
       &ScalarTensor::ConstraintDamping::register_derived_with_charm,
       &fe::DecoupledScalar::ConstraintDamping::register_derived_with_charm,
       &control_system::size::register_derived_with_charm,
       &register_factory_classes_with_charm<EvolutionMetavars>},
      {});
}
