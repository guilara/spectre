// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/AnalyticConstant.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarTensorDriver::BoundaryConditions {
/// Typelist of standard BoundaryConditions
using standard_boundary_conditions =
    tmpl::list<AnalyticConstant, DemandOutgoingCharSpeeds,
               //    ConstraintPreservingSphericalRadiation,
               domain::BoundaryConditions::Periodic<BoundaryCondition>>;
}  // namespace fe::ScalarTensorDriver::BoundaryConditions
