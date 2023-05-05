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
//
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
//
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

template <evolution::BoundaryConditions::Type GhBcType,
          evolution::BoundaryConditions::Type ScalarBcType>
struct UnionOfBcTypes {
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;
};

// This defines evolution::BoundaryConditions::Type bc_type
// for different templated variables
// For now only Ghost type for exterior boundaries and DemandOutgoing for
// interior boundaries
template <>
struct UnionOfBcTypes<evolution::BoundaryConditions::Type::Ghost,
                      evolution::BoundaryConditions::Type::Ghost> {
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;
};

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

}  // namespace detail

/*!
 * \brief Apply a boundary condition to the combined Generalized Harmonic (GH)
 * and scalar wave system using the boundary conditions defined separately
 * for the GH and Curved Scalar systems.
 */
template <typename DerivedGhCondition, typename DerivedScalarCondition>
class ProductOfConditions final : public BoundaryCondition {
 public:
  static constexpr evolution::BoundaryConditions::Type bc_type =
      detail::UnionOfBcTypes<DerivedGhCondition::bc_type,
                             DerivedScalarCondition::bc_type>::bc_type;

  // Removed Ghost and TimeDerivative BCs in the using statements
  using dg_interior_evolved_variables_tags =
      tmpl::remove_duplicates<tmpl::append<
          typename DerivedGhCondition::dg_interior_evolved_variables_tags,
          typename DerivedScalarCondition::dg_interior_evolved_variables_tags>>;

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

  std::optional<std::string> dg_demand_outgoing_char_speeds(
      // GH arguments
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          outward_directed_normal_covector,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>&
          outward_directed_normal_vector,

      // Scalar arguments
    //   const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    //       face_mesh_velocity,
    //   const tnsr::i<DataVector, Dim>& normal_covector,
    //   const tnsr::I<DataVector, Dim>& /*normal_vector*/,

      // GH
      const Scalar<DataVector>& gamma_1,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,

      // Scalar
      const Scalar<DataVector>& gamma1_scalar
    //,const Scalar<DataVector>& lapse, const tnsr::I<DataVector, Dim>& shift
      ) const {
    // DemandOutgoingCharSpeeds condition is only valid if both boundary
    // conditions are DemandOutgoingCharSpeeds, so we directly apply both. A
    // static_assert elsewhere is triggered if only one boundary condition is
    // DemandOutgoingCharSpeeds.
    auto gh_string = derived_gh_condition_.dg_demand_outgoing_char_speeds(
            face_mesh_velocity,
            outward_directed_normal_covector,
            outward_directed_normal_vector,

            gamma_1, lapse, shift);
    auto scalar_string =
        derived_scalar_condition_.dg_demand_outgoing_char_speeds(
            face_mesh_velocity,
            outward_directed_normal_covector,
            outward_directed_normal_vector,

            gamma1_scalar, lapse, shift);
    if (not gh_string.has_value()) {
      return scalar_string;
    }
    if (not scalar_string.has_value()) {
      return gh_string;
    }
    return gh_string.value() + ";" + scalar_string.value();
  }

  // A temporal solution to using pack expansions (variadic arguments) is to
  // use constexpr if statements checking for the different boundary
  // conditions passed to dg_ghost and passing the right arguments to each.
  // Update: Seems that I cannot make conditional the definition of a class
  // member function

  // Boundary conditions for Dirichlet-Minkowski/Constant
  std::optional<std::string> dg_ghost(
      // GH evolved variables
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          spacetime_metric,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*> pi,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*> phi,
      // Scalar evolved variables. Change names
      const gsl::not_null<Scalar<DataVector>*> psi_scalar,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_scalar,
      // GH temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1,
      const gsl::not_null<Scalar<DataVector>*> gamma2,
      const gsl::not_null<Scalar<DataVector>*> lapse,
      const gsl::not_null<tnsr::I<DataVector, 3_st, Frame::Inertial>*> shift,
      // Scalar temporary variables. Change names
      // const gsl::not_null<Scalar<DataVector>*> lapse,
      // const gsl::not_null<tnsr::I<DataVector, Dim>*> shift,
      const gsl::not_null<Scalar<DataVector>*> gamma1_scalar,
      const gsl::not_null<Scalar<DataVector>*> gamma2_scalar,
      // Inverse metric
      const gsl::not_null<tnsr::II<DataVector, 3_st, Frame::Inertial>*>
          inv_spatial_metric,
      // Mesh variables
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& normal_vector,
      // GH interior variables
      const Scalar<DataVector>& interior_gamma1,
      const Scalar<DataVector>& interior_gamma2,
      // Scalar interior variables
      const tnsr::II<DataVector, 3_st, Frame::Inertial>&
          inverse_spatial_metric_interior,
      const Scalar<DataVector>& gamma1_interior_scalar,
      const Scalar<DataVector>& gamma2_interior_scalar,
      const Scalar<DataVector>& lapse_interior,
      const tnsr::I<DataVector, 3_st>& shift_interior) const {
    // Note: Check that CurvedScalarWave does not update GH variables
    // to a different value. If it does, invert the order of application of the
    // corrections first, so that the GH update is applied at last

    // GeneralizedHarmonic::BoundaryConditions::DirichletMinkowski
    auto gh_string = derived_gh_condition_.dg_ghost(
        spacetime_metric,
        pi,
        phi,
        gamma1,
        gamma2,
        lapse,
        shift,
        inv_spatial_metric,
        face_mesh_velocity,
        normal_covector,
        normal_vector,
        interior_gamma1,
        interior_gamma2
        );

    // CurvedScalarWave::BoundaryConditions::AnalyticConstant
    auto scalar_string = derived_scalar_condition_.dg_ghost(
      // Change names
        psi_scalar,
        pi_scalar,
        phi_scalar,
        lapse,
        shift,
        gamma1_scalar,
        gamma2_scalar,
        inv_spatial_metric,
        face_mesh_velocity,
        normal_covector,
        normal_vector,
        inverse_spatial_metric_interior,
        gamma1_interior_scalar,
        gamma2_interior_scalar,
        lapse_interior,
        shift_interior);
    if (not gh_string.has_value()) {
      return scalar_string;
    }
    if (not scalar_string.has_value()) {
      return gh_string;
    }
    return gh_string.value() + ";" + scalar_string.value();
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
