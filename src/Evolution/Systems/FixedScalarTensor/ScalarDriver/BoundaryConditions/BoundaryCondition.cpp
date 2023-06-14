// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensor/BoundaryConditions/BoundaryCondition.hpp"

#include <pup.h>

#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace fe::ScalarDriver::BoundaryConditions {
BoundaryCondition<3_st>::BoundaryCondition(CkMigrateMessage* const msg)
    : domain::BoundaryConditions::BoundaryCondition(msg) {}

void BoundaryCondition<3_st>::pup(PUP::er& p) {
  domain::BoundaryConditions::BoundaryCondition::pup(p);
}

}  // namespace fe::ScalarDriver::BoundaryConditions
