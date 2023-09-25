// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/AnalyticConstant.hpp"

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
AnalyticConstant::get_clone() const {
  return std::make_unique<AnalyticConstant>(*this);
}

void AnalyticConstant::pup(PUP::er& p) {
  BoundaryCondition::pup(p);
}

AnalyticConstant::AnalyticConstant(CkMigrateMessage* const msg)
    : BoundaryCondition(msg) {}

// NOLINTNEXTLINE
PUP::able::PUP_ID AnalyticConstant::my_PUP_ID = 0;

}  // namespace fe::ScalarDriver::BoundaryConditions
