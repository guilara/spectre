// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/RegisterDerived.hpp"

#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/Factory.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace fe::DecoupledScalar {
void register_derived_with_charm() {
  register_derived_classes_with_charm<BoundaryCorrection>();
}
}  // namespace fe::DecoupledScalar
