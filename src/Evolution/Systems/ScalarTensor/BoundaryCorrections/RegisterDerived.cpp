// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/RegisterDerived.hpp"

#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/Factory.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"

namespace ScalarTensor::BoundaryCorrections {
void register_derived_with_charm() {
  Parallel::register_derived_classes_with_charm<BoundaryCorrection>();
}
}  // namespace ScalarTensor::BoundaryCorrections
