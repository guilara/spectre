// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Factory.hpp"
// Add Scalar headers
#include "Utilities/TMPL.hpp"

namespace ScalarTensor::BoundaryConditions {

namespace detail {

template </* Add templated variables */
          typename DerivedGhCondition, typename DerivedScalarCondition>
using ProductOfConditionsIfConsistent = tmpl::conditional_t<
    (DerivedGhCondition::bc_type ==
     evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds) xor
        (DerivedScalarCondition::bc_type ==
         evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds),
    tmpl::list<>,
    ProductOfConditions<DerivedGhCondition, DerivedScalarCondition>>;

template </* Add templated variables */
          typename GhList, typename ScalarList>
struct AllProductConditions {};

template </* Add templated variables */
          typename GhList, typename... ScalarConditions>
struct AllProductConditions</* Add templated variables */
                            GhList, tmpl::list<ScalarConditions...>> {
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
using standard_boundary_conditions = tmpl::push_back<
    typename detail::AllProductConditions<
        detail::remove_periodic_conditions_t<
            typename GeneralizedHarmonic::BoundaryConditions::
                standard_boundary_conditions<3_st>>,
        detail::remove_periodic_conditions_t<
            typename CurvedScalarWave::BoundaryConditions::
                standard_boundary_conditions>>::type,
    domain::BoundaryConditions::Periodic<BoundaryCondition>>;

}  // namespace ScalarTensor::BoundaryConditions
