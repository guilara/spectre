// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/FixedSGB/BoundaryCorrections/RegisterDerived.hpp"

#include "Evolution/Systems/FixedScalarTensor/FixedSGB/BoundaryCorrections/Factory.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace fe::sgb::BoundaryCorrections {
void register_derived_with_charm() {
  register_derived_classes_with_charm<BoundaryCorrection>();
}
}  // namespace fe::sgb::BoundaryCorrections
