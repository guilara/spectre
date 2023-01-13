// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <utility>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivativeHelpers.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
// Missing Curved Scalar and Scalar Tensor headers
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace ScalarTensor::BoundaryConditions {
namespace detail {

// This defines evolution::BoundaryConditions::Type bc_type
// for different templated variables
// Have not included Ghost and TimeDerivative BCs
template </* Add templated variables */>
struct UnionOfBcTypes {};

template <evolution::BoundaryConditions::Type GhBcType>
struct UnionOfBcTypes<
    GhBcType, evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds> {
  static_assert(
      GhBcType == evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds,
      "If either boundary condition in `ProductOfConditions` has "
      "`Type::DemandOutgoingCharSpeeds`, both must have "
      "`Type::DemandOutgoingCharSpeeds`");
};

template <evolution::BoundaryConditions::Type ScalarBcType>
struct UnionOfBcTypes<
    evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds,
    ScalarBcType> {
  static_assert(
      ScalarBcType ==
          evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds,
      "If either boundary condition in `ProductOfConditions` has "
      "`Type::DemandOutgoingCharSpeeds`, both must have "
      "`Type::DemandOutgoingCharSpeeds`");
};

template <>
struct UnionOfBcTypes<
    evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds,
    evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds> {
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds;
};

// Need to understand if we need this
// This could combine ghost for one and something else for the other
struct ScalarDoNothingGhostCondition {
  // What does this do?
  std::optional<std::string> dg_ghost();
};

// Need to understand if we need this
struct GhDoNothingGhostCondition {
  // What does this do?
  std::optional<std::string> dg_ghost();
};

template <
    typename DerivedGhCondition, typename DerivedScalarCondition,
    typename GhEvolvedTagList, typename ScalarEvolvedTagList,
    typename GhFluxTagList, typename ScalarFluxTagList,
    typename GhInteriorEvolvedTagList,
    typename ScalarInteriorEvolvedTagList,
    typename DeduplicatedInteriorEvolvedTagList,
    // Comment this ?
    typename GhInteriorPrimitiveTagList,
    typename ScalarInteriorPrimitiveTgs,
    //
    typename GhInteriorTempTagList,
    typename ScalarInteriorTempTagList,
    typename DeduplicatedTempTagList, typename GhInteriorDtTagList,
    typename ScalarInteriorDtTagList, typename GhInteriorDerivTagList,
    typename ScalarInteriorDerivTagList, typename GhGridlessTagList,
    typename ScalarGridlessTagList,
    typename DeduplicatedGridlessTagList>
struct ProductOfConditionsImpl;

template <typename DerivedGhCondition, typename DerivedScalarCondition,
          typename... GhEvolvedTags, typename... ScalarEvolvedTags,
          typename... GhFluxTags, typename... ScalarFluxTags,
          typename... GhInteriorEvolvedTags,
          typename... ScalarInteriorEvolvedTags,
          typename... DeduplicatedInteriorEvolvedTags,
          // Comment this ?
          typename... GhInteriorPrimitiveTags,
          typename... ScalarInteriorPrimitiveTags,
          //
          typename... GhInteriorTempTags, typename... ScalarInteriorTempTags,
          typename... DeduplicatedTempTags, typename... GhInteriorDtTags,
          typename... ScalarInteriorDtTags, typename... GhInteriorDerivTags,
          typename... ScalarInteriorDerivTags, typename... GhGridlessTags,
          typename... ScalarGridlessTags, typename... DeduplicatedGridlessTags>
struct ProductOfConditionsImpl<
    DerivedGhCondition, DerivedScalarCondition, tmpl::list<GhEvolvedTags...>,
    tmpl::list<ScalarEvolvedTags...>, tmpl::list<GhFluxTags...>,
    tmpl::list<ScalarFluxTags...>, tmpl::list<GhInteriorEvolvedTags...>,
    tmpl::list<ScalarInteriorEvolvedTags...>,
    tmpl::list<DeduplicatedInteriorEvolvedTags...>,
    // Comment this?
    tmpl::list<GhInteriorPrimitiveTags...>,
    tmpl::list<ScalarInteriorPrimitiveTags...>,
    //
    tmpl::list<GhInteriorTempTags...>, tmpl::list<ScalarInteriorTempTags...>,
    tmpl::list<DeduplicatedTempTags...>, tmpl::list<GhInteriorDtTags...>,
    tmpl::list<ScalarInteriorDtTags...>, tmpl::list<GhInteriorDerivTags...>,
    tmpl::list<ScalarInteriorDerivTags...>, tmpl::list<GhGridlessTags...>,
    tmpl::list<ScalarGridlessTags...>,
    tmpl::list<DeduplicatedGridlessTags...>> {

// Puts together this function from both systems
// CurvedScalarWave seems not to have dg_ghost() defined
// template </* Add templated variables */>
// static std::optional<std::string> dg_ghost(/* Add variables */);

template </* Add templated variables */
          typename... GridlessVariables>
static std::optional<std::string> dg_demand_outgoing_char_speeds(
    /* Add variables */
    const DerivedGhCondition& gh_condition,
    const DerivedScalarCondition& scalar_condition,
    const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& normal_vector,
    const typename DeduplicatedInteriorEvolvedTags::type&...
        int_evolved_variables,
    // Comment this ?
    const typename GhInteriorPrimitiveTags::type&... gh_int_prim_variables,
    const typename ScalarInteriorPrimitiveTags::type&...
        scalar_int_prim_variables,
    //
    const typename DeduplicatedTempTags::type&... temp_variables,
    const typename GhInteriorDtTags::type&... gh_int_dt_variables,
    const typename ScalarInteriorDtTags::type&... scalar_int_dt_variables,
    const typename GhInteriorDerivTags::type&... gh_int_deriv_variables,
    const typename ScalarInteriorDerivTags::type&...
        scalar_int_deriv_variables,
    const GridlessVariables&... gridless_variables) {
  using gridless_tags_and_types =
      tmpl::map<tmpl::pair<DeduplicatedGridlessTags, GridlessVariables>...>;

  tuples::TaggedTuple<
      Tags::detail::TemporaryReference<
          DeduplicatedGridlessTags,
          tmpl::at<gridless_tags_and_types, DeduplicatedGridlessTags>>...,
      Tags::detail::TemporaryReference<DeduplicatedTempTags>...,
      Tags::detail::TemporaryReference<DeduplicatedInteriorEvolvedTags>...>
      shuffle_refs{gridless_variables..., temp_variables...,
                   int_evolved_variables...};
  // DemandOutgoingCharSpeeds condition is only valid if both boundary
  // conditions are DemandOutgoingCharSpeeds, so we directly apply both. A
  // static_assert elsewhere is triggered if only one boundary condition is
  // DemandOutgoingCharSpeeds.
  auto gh_string = gh_condition.dg_demand_outgoing_char_speeds(
      face_mesh_velocity, normal_covector, normal_vector,
      tuples::get<Tags::detail::TemporaryReference<GhInteriorEvolvedTags>>(
          shuffle_refs)...,
      gh_int_prim_variables...,
      tuples::get<Tags::detail::TemporaryReference<GhInteriorTempTags>>(
          shuffle_refs)...,
      gh_int_dt_variables..., gh_int_deriv_variables...,
      tuples::get<Tags::detail::TemporaryReference<
          GhGridlessTags, tmpl::at<gridless_tags_and_types, GhGridlessTags>>>(
          shuffle_refs)...);
  auto scalar_string = scalar_condition.dg_demand_outgoing_char_speeds(
      face_mesh_velocity, normal_covector, normal_vector,
      tuples::get<
          Tags::detail::TemporaryReference<ScalarInteriorEvolvedTags>>(
          shuffle_refs)...,
      scalar_int_prim_variables...,
      tuples::get<Tags::detail::TemporaryReference<ScalarInteriorTempTags>>(
          shuffle_refs)...,
      scalar_int_dt_variables..., scalar_int_deriv_variables...,
      tuples::get<Tags::detail::TemporaryReference<
          ScalarGridlessTags,
          tmpl::at<gridless_tags_and_types, ScalarGridlessTags>>>(
          shuffle_refs)...);
  if (not gh_string.has_value()) {
    return scalar_string;
  }
  if (not scalar_string.has_value()) {
    return gh_string;
  }
  return gh_string.value() + ";" + scalar_string.value();
}

// template </* Add templated variables */>
// static std::optional<std::string> dg_time_derivative(/* Add variables */);
};
}  // namespace detail

/*!
 * \brief Apply a boundary condition to the combined Generalized Harmonic (GH)
 * and scalar wave system using the boundary conditions defined separately
 * for the GH and Curved Scalar systems.
 */
template <typename DerivedGhCondition, typename DerivedScalarCondition>
class ProductOfConditions final : public BoundaryCondition {
 public:
    // Add missing using statements and ProductOfConditionsImpl function calls

// Removed Ghost and TimeDerivative BCs in the using statements
  using dg_interior_evolved_variables_tags =
      tmpl::remove_duplicates<tmpl::append<
          typename DerivedGhCondition::dg_interior_evolved_variables_tags,
          typename DerivedScalarCondition::dg_interior_evolved_variables_tags>>;

  // Comment this line?
  using dg_interior_primitive_variables_tags = tmpl::append<
      tmpl::list<>,
      typename DerivedScalarCondition::dg_interior_primitive_variables_tags>;
  //

  using dg_interior_temporary_tags = tmpl::remove_duplicates<tmpl::append<
      typename DerivedGhCondition::dg_interior_temporary_tags,
      typename DerivedScalarCondition::dg_interior_temporary_tags>>;

  using dg_gridless_tags = tmpl::remove_duplicates<
      tmpl::append<typename DerivedGhCondition::dg_gridless_tags,
                   typename DerivedScalarCondition::dg_gridless_tags>>;

  using dg_interior_dt_vars_tags = tmpl::append<
      evolution::dg::Actions::detail::get_dt_vars_from_boundary_condition<
          DerivedGhCondition>,
      evolution::dg::Actions::detail::get_dt_vars_from_boundary_condition<
          DerivedScalarCondition>>;

  using dg_interior_deriv_vars_tags = tmpl::append<
      evolution::dg::Actions::detail::get_deriv_vars_from_boundary_condition<
          DerivedGhCondition>,
      evolution::dg::Actions::detail::get_deriv_vars_from_boundary_condition<
          DerivedScalarCondition>>;

  using product_of_conditions_impl = detail::ProductOfConditionsImpl<
      DerivedGhCondition, DerivedScalarCondition,
      typename GeneralizedHarmonic::System<3_st>::variables_tag::tags_list,
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list,
      db::wrap_tags_in<
          ::Tags::Flux,
          typename GeneralizedHarmonic::System<3_st>::flux_variables,
          tmpl::size_t<3_st>, Frame::Inertial>,
      db::wrap_tags_in<::Tags::Flux,
                       typename CurvedScalarWave::System<3_st>::flux_variables,
                       tmpl::size_t<3_st>, Frame::Inertial>,
      typename DerivedGhCondition::dg_interior_evolved_variables_tags,
      typename DerivedScalarCondition::dg_interior_evolved_variables_tags,
      dg_interior_evolved_variables_tags,
      // Comment this line?
      tmpl::list<>,
      typename DerivedScalarCondition::dg_interior_primitive_variables_tags,
      //
      typename DerivedGhCondition::dg_interior_temporary_tags,
      typename DerivedScalarCondition::dg_interior_temporary_tags,
      dg_interior_temporary_tags,
      evolution::dg::Actions::detail::get_dt_vars_from_boundary_condition<
          DerivedGhCondition>,
      evolution::dg::Actions::detail::get_dt_vars_from_boundary_condition<
          DerivedScalarCondition>,
      evolution::dg::Actions::detail::get_deriv_vars_from_boundary_condition<
          DerivedGhCondition>,
      evolution::dg::Actions::detail::get_deriv_vars_from_boundary_condition<
          DerivedScalarCondition>,
      typename DerivedGhCondition::dg_gridless_tags,
      typename DerivedScalarCondition::dg_gridless_tags, dg_gridless_tags>;

  static std::string name() {
    return "Product" + pretty_type::name<DerivedGhCondition>() + "And" +
           pretty_type::name<DerivedScalarCondition>();
  }

  struct GhCondition {
    using type = DerivedGhCondition;
    static std::string name() {
      return "GeneralizedHarmonic" + pretty_type::name<DerivedGhCondition>();
    }
    static constexpr Options::String help{
        "The Generalized Harmonic part of the product boundary condition"};
  };

  struct ScalarCondition {
    using type = DerivedScalarCondition;
    static std::string name() {
      return "Scalar" + pretty_type::name<DerivedScalarCondition>();
    }
    static constexpr Options::String help{
        "The Scalar part of the product boundary condition"};
  };

  using options = tmpl::list<GhCondition, ScalarCondition>;

  static constexpr Options::String help = {
      "Direct product of a GH and Curved Scalar boundary conditions. "
      "See the documentation for the two individual boundary conditions for "
      "further details."};

  ProductOfConditions() = default;
  ProductOfConditions(DerivedGhCondition gh_condition,
                      DerivedScalarCondition scalar_condition)
      : derived_gh_condition_{gh_condition},
        derived_scalar_condition_{scalar_condition} {}
  ProductOfConditions(const ProductOfConditions&) = default;
  ProductOfConditions& operator=(const ProductOfConditions&) = default;
  ProductOfConditions(ProductOfConditions&&) = default;
  ProductOfConditions& operator=(ProductOfConditions&&) = default;
  ~ProductOfConditions() override = default;

  /// \cond
  explicit ProductOfConditions(CkMigrateMessage* msg)
      : BoundaryCondition(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, ProductOfConditions);
  /// \endcond

  void pup(PUP::er& p) override;

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  // template <typename... Args>
  // std::optional<std::string> dg_ghost(Args&&... args);

  template <typename... Args>
  std::optional<std::string> dg_demand_outgoing_char_speeds(
      Args&&... args) const {
    return product_of_conditions_impl::dg_demand_outgoing_char_speeds(
        derived_gh_condition_, derived_scalar_condition_,
        std::forward<Args>(args)...);

  // template <typename... Args>
  // std::optional<std::string> dg_time_derivative(Args && ... args);

  }

 private:
  DerivedGhCondition derived_gh_condition_;
  DerivedScalarCondition derived_scalar_condition_;
};

template <typename DerivedGhCondition, typename DerivedScalarCondition>
void ProductOfConditions<DerivedGhCondition, DerivedScalarCondition>::pup(
    PUP::er& p) {
  BoundaryCondition::pup(p);
  p | derived_gh_condition_;
  p | derived_scalar_condition_;
}

template <typename DerivedGhCondition, typename DerivedScalarCondition>
auto ProductOfConditions<DerivedGhCondition,
                         DerivedScalarCondition>::get_clone() const
    -> std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> {
  return std::make_unique<ProductOfConditions>(*this);
}

/// \cond
template <typename DerivedGhCondition, typename DerivedScalarCondition>
PUP::able::PUP_ID ProductOfConditions<DerivedGhCondition,
                                      DerivedScalarCondition>::my_PUP_ID =
    0;  // NOLINT
/// \endcond

template <typename DerivedGhCondition, typename DerivedScalarCondition>
ProductOfConditions(DerivedGhCondition gh_condition,
                    DerivedScalarCondition scalar_condition)
    -> ProductOfConditions<DerivedGhCondition, DerivedScalarCondition>;

}  // namespace ScalarTensor::BoundaryConditions
