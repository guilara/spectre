// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/ConstraintPreservingSphericalRadiation.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeString.hpp"

namespace fe::ScalarDriver::BoundaryConditions {

std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
ConstraintPreservingSphericalRadiation::get_clone() const {
  return std::make_unique<ConstraintPreservingSphericalRadiation>(*this);
}

void ConstraintPreservingSphericalRadiation::pup(PUP::er& p) {
  BoundaryCondition::pup(p);
  p | csw_constraint_preserving_;
}

ConstraintPreservingSphericalRadiation::ConstraintPreservingSphericalRadiation(
    CkMigrateMessage* const msg)
    : BoundaryCondition(msg) {}

ConstraintPreservingSphericalRadiation::ConstraintPreservingSphericalRadiation()
    : {
  csw_constraint_preserving_ = CurvedScalarWave::BoundaryConditions::
      ConstraintPreservingSphericalRadiation<3>();
}

// NOLINTNEXTLINE
PUP::able::PUP_ID ConstraintPreservingSphericalRadiation::my_PUP_ID = 0;

}  // namespace fe::ScalarDriver::BoundaryConditions
