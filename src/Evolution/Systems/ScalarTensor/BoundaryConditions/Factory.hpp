// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Factory.hpp"
//
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Bjorhus.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletMinkowski.hpp"
//
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Factory.hpp"
//
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/AnalyticConstant.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/ConstraintPreservingSphericalRadiation.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
//
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/ProductOfConditions.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarTensor::BoundaryConditions {

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

template <typename ClassList>
struct remove_periodic_conditions {
using type = tmpl::remove_if<
    ClassList,
    std::is_base_of<domain::BoundaryConditions::MarkAsPeriodic, tmpl::_1>>;
};

template <typename ClassList>
using remove_periodic_conditions_t =
    typename remove_periodic_conditions<ClassList>::type;

}  // namespace detail

// remove the periodic BCs from the creatable classes of the
// individual systems; for the remaining conditions, include a
// `ProductOfConditions` for each pair with compatible `bc_type`s.
/// Typelist of standard BoundaryConditions
// using standard_boundary_conditions = tmpl::push_back<
//     typename detail::AllProductConditions<
//         detail::remove_periodic_conditions_t<
//             typename GeneralizedHarmonic::BoundaryConditions::
//                 standard_boundary_conditions<3_st>>,
//         detail::remove_periodic_conditions_t<
//             typename CurvedScalarWave::BoundaryConditions::
//                 standard_boundary_conditions<3_st>>>::type,
//     domain::BoundaryConditions::Periodic<BoundaryCondition>>;

// For now, we only support a subset of the available boundary conditions
// We copy this list as in the factory of BCs for each system
// but omit periodic boundary conditions
/// Typelist of standard BoundaryConditions
using subset_standard_boundary_conditions_gh = tmpl::list<
  // GeneralizedHarmonic::BoundaryConditions::ConstraintPreservingBjorhus<3_st>,
    GeneralizedHarmonic::BoundaryConditions::DemandOutgoingCharSpeeds<3_st>,
    GeneralizedHarmonic::BoundaryConditions::DirichletAnalytic<3_st>
    // ,
    // GeneralizedHarmonic::BoundaryConditions::DirichletMinkowski<3_st>
    >;

using subset_standard_boundary_conditions_scalar = tmpl::list<
    CurvedScalarWave::BoundaryConditions::AnalyticConstant<3_st>,
    // CurvedScalarWave::BoundaryConditions::
    //     ConstraintPreservingSphericalRadiation<3_st>,
    CurvedScalarWave::BoundaryConditions::DemandOutgoingCharSpeeds<3_st>>;
using standard_boundary_conditions =
    tmpl::push_back<typename detail::AllProductConditions<
                        subset_standard_boundary_conditions_gh,
                        subset_standard_boundary_conditions_scalar>::type,
                    domain::BoundaryConditions::Periodic<BoundaryCondition>>;

}  // namespace ScalarTensor::BoundaryConditions
