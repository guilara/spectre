// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/ConstraintPreservingSphericalRadiation.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarDriver::BoundaryConditions {
/// Typelist of standard BoundaryConditions
using standard_boundary_conditions =
    tmpl::list<DemandOutgoingCharSpeeds,
               domain::BoundaryConditions::Periodic<BoundaryCondition>>;
}  // namespace fe::ScalarDriver::BoundaryConditions
