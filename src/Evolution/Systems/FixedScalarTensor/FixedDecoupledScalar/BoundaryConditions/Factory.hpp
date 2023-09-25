// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryConditions/Factory.hpp"
//
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryConditions/BoundaryCondition.hpp"
//
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/Factory.hpp"
//
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/AnalyticConstant.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
//
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryConditions/ProductOfConditions.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::DecoupledScalar::BoundaryConditions {

namespace detail {

template <typename DerivedGhCondition, typename DerivedScalarCondition>
using ProductOfConditionsIfConsistent = tmpl::conditional_t<
    (DerivedGhCondition::bc_type ==
     evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds) xor
        (DerivedScalarCondition::bc_type ==
         evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds),
    tmpl::list<>,
    ProductOfConditions<DerivedGhCondition, DerivedScalarCondition>>;

template <typename GhList, typename ScalarList>
struct AllProductConditions;

template <typename GhList, typename... ScalarConditions>
struct AllProductConditions<GhList, tmpl::list<ScalarConditions...>> {
  using type = tmpl::flatten<tmpl::list<tmpl::transform<
      GhList, tmpl::bind<ProductOfConditionsIfConsistent, tmpl::_1,
                         tmpl::pin<ScalarConditions>>>...>>;
};

}  // namespace detail

// remove the periodic BCs from the creatable classes of the
// individual systems; for the remaining conditions, include a
// `ProductOfConditions` for each pair with compatible `bc_type`s.
/// Typelist of standard BoundaryConditions

// For now, we only support a subset of the available boundary conditions
// We copy this list as in the factory of BCs for each system
// but omit periodic boundary conditions
/// Typelist of standard BoundaryConditions

// We write here the product of boundary conditions as in the Scalar Tensor
// system
using subset_standard_boundary_conditions_gh =
    detail::remove_periodic_conditions_t<
        typename ScalarTensor::BoundaryConditions::
            standard_boundary_conditions>;

using subset_standard_boundary_conditions_scalar =
    tmpl::list<fe::ScalarDriver::BoundaryConditions::AnalyticConstant,
               fe::ScalarDriver::BoundaryConditions::DemandOutgoingCharSpeeds>;

using standard_boundary_conditions =
    tmpl::push_back<typename detail::AllProductConditions<
                        subset_standard_boundary_conditions_gh,
                        subset_standard_boundary_conditions_scalar>::type,
                    domain::BoundaryConditions::Periodic<BoundaryCondition>>;

}  // namespace fe::DecoupledScalar::BoundaryConditions
