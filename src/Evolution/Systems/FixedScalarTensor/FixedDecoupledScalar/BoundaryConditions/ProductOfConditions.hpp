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
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/System.hpp"
//
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace fe::DecoupledScalar::BoundaryConditions {
namespace detail {

template <evolution::BoundaryConditions::Type GhBcType,
          evolution::BoundaryConditions::Type ScalarBcType>
struct UnionOfBcTypes {
  static_assert(
          (GhBcType ==
           evolution::BoundaryConditions::Type::GhostAndTimeDerivative) or
          (ScalarBcType ==
           evolution::BoundaryConditions::Type::GhostAndTimeDerivative),
      "TimeDerivative boundary conditions are not yet supported.");
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;
};

template <evolution::BoundaryConditions::Type GhBcType>
struct UnionOfBcTypes<GhBcType,
                      evolution::BoundaryConditions::Type::TimeDerivative> {
  static_assert(GhBcType == evolution::BoundaryConditions::Type::TimeDerivative,
                "If either boundary condition in `ProductOfConditions` has "
                "`Type::TimeDerivative`, both must have "
                "`Type::TimeDerivative`");
};

template <evolution::BoundaryConditions::Type ScalarBcType>
struct UnionOfBcTypes<evolution::BoundaryConditions::Type::TimeDerivative,
                      ScalarBcType> {
  static_assert(ScalarBcType ==
                    evolution::BoundaryConditions::Type::TimeDerivative,
                "If either boundary condition in `ProductOfConditions` has "
                "`Type::TimeDerivative`, both must have "
                "`Type::TimeDerivative`");
};

template <>
struct UnionOfBcTypes<evolution::BoundaryConditions::Type::TimeDerivative,
                      evolution::BoundaryConditions::Type::TimeDerivative> {
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::TimeDerivative;
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
      return "ScalarTensor" + pretty_type::name<DerivedGhCondition>();
    }
    static constexpr Options::String help{
        "The ScalarTensor part of the product boundary condition"};
  };

  struct ScalarCondition {
    using type = DerivedScalarCondition;
    static std::string name() {
      return "ScalarDriver" + pretty_type::name<DerivedScalarCondition>();
    }
    static constexpr Options::String help{
        "The ScalarDriver part of the product boundary condition"};
  };

  using options = tmpl::list<GhCondition, ScalarCondition>;

  static constexpr Options::String help = {
      "Direct product of a ScalarTensor and ScalarDriver boundary conditions. "
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
      // GH
      const Scalar<DataVector>& gamma_1, const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
      // Scalar
      const Scalar<DataVector>& gamma1_scalar,
      // Scalar Driver
      const Scalar<DataVector>& gamma1_scalar_driver) const {
    // DemandOutgoingCharSpeeds condition is only valid if both boundary
    // conditions are DemandOutgoingCharSpeeds, so we directly apply both. A
    // static_assert elsewhere is triggered if only one boundary condition is
    // DemandOutgoingCharSpeeds.
    auto gh_string = derived_gh_condition_.dg_demand_outgoing_char_speeds(
        // GH arguments
        face_mesh_velocity, outward_directed_normal_covector,
        outward_directed_normal_vector,
        // GH
        gamma_1, lapse, shift,
        // Scalar
        gamma1_scalar);
    auto scalar_string =
        derived_scalar_condition_.dg_demand_outgoing_char_speeds(
            face_mesh_velocity, outward_directed_normal_covector,
            outward_directed_normal_vector,

            gamma1_scalar_driver, lapse, shift);
    if (not gh_string.has_value()) {
      return scalar_string;
    }
    if (not scalar_string.has_value()) {
      return gh_string;
    }
    return gh_string.value() + ";" + scalar_string.value();
  }

  // We overload dg_ghost for the different analytic boundary conditions

  // Boundary conditions for Dirichlet-Analytic/Constant
  std::optional<std::string> dg_ghost(
      // GH evolved variables
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          spacetime_metric,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*> pi,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*> phi,
      // Scalar evolved variables
      const gsl::not_null<Scalar<DataVector>*> psi_scalar,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_scalar,
      // Scalar Driver evolved variables
      const gsl::not_null<Scalar<DataVector>*> psi_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar_driver,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_scalar_driver,
      // GH temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1,
      const gsl::not_null<Scalar<DataVector>*> gamma2,
      const gsl::not_null<Scalar<DataVector>*> lapse,
      const gsl::not_null<tnsr::I<DataVector, 3_st, Frame::Inertial>*> shift,
      // Scalar temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1_scalar,
      const gsl::not_null<Scalar<DataVector>*> gamma2_scalar,
      // Scalar temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> gamma2_scalar_driver,
      // Inverse metric
      const gsl::not_null<tnsr::II<DataVector, 3_st, Frame::Inertial>*>
          inv_spatial_metric,
      // Mesh variables
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& normal_vector,
      // GH interior variables
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& coords,
      const Scalar<DataVector>& interior_gamma1,
      const Scalar<DataVector>& interior_gamma2,
      // Scalar interior variables
      const tnsr::II<DataVector, 3_st, Frame::Inertial>&
          inverse_spatial_metric_interior,
      const Scalar<DataVector>& gamma1_interior_scalar,
      const Scalar<DataVector>& gamma2_interior_scalar,
      const Scalar<DataVector>& lapse_interior,
      const tnsr::I<DataVector, 3_st>& shift_interior,
      // Scalar interior variables
      const Scalar<DataVector>& gamma1_interior_scalar_driver,
      const Scalar<DataVector>& gamma2_interior_scalar_driver,
      const double time) const {
    // Note: Check that CurvedScalarWave does not update GH variables
    // to a different value. If it does, invert the order of application of the
    // corrections first, so that the GH update is applied at last

    // gh::BoundaryConditions::DirichletAnalytic
    auto gh_string = derived_gh_condition_.dg_ghost(
        // GH evolved variables
        spacetime_metric, pi, phi,
        // Scalar evolved variables. Change names
        psi_scalar, pi_scalar, phi_scalar,
        // GH temporary variables
        gamma1, gamma2, lapse, shift,
        // Scalar temporary variables
        gamma1_scalar, gamma2_scalar,
        // Inverse metric
        inv_spatial_metric,
        // Mesh variables
        face_mesh_velocity, normal_covector, normal_vector,
        // GH interior variables
        coords, interior_gamma1, interior_gamma2,
        // Scalar interior variables
        inverse_spatial_metric_interior, gamma1_interior_scalar,
        gamma2_interior_scalar, lapse_interior, shift_interior, time);

    // CurvedScalarWave::BoundaryConditions::AnalyticConstant
    auto scalar_string = derived_scalar_condition_.dg_ghost(
        // Change names
        psi_scalar_driver, pi_scalar_driver, phi_scalar_driver, lapse, shift,
        gamma1_scalar_driver, gamma2_scalar_driver, inv_spatial_metric,
        face_mesh_velocity, normal_covector, normal_vector,
        inverse_spatial_metric_interior, gamma1_interior_scalar_driver,
        gamma2_interior_scalar_driver, lapse_interior, shift_interior);
    if (not gh_string.has_value()) {
      return scalar_string;
    }
    if (not scalar_string.has_value()) {
      return gh_string;
    }
    return gh_string.value() + ";" + scalar_string.value();
  }

  // Boundary conditions for Dirichlet-Minkowski/Constant
  std::optional<std::string> dg_ghost(
      // GH evolved variables
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          spacetime_metric,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*> pi,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*> phi,
      // Scalar evolved variables
      const gsl::not_null<Scalar<DataVector>*> psi_scalar,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_scalar,
      // Scalar driver evolved variables
      const gsl::not_null<Scalar<DataVector>*> psi_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar_driver,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_scalar_driver,
      // GH temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1,
      const gsl::not_null<Scalar<DataVector>*> gamma2,
      const gsl::not_null<Scalar<DataVector>*> lapse,
      const gsl::not_null<tnsr::I<DataVector, 3_st, Frame::Inertial>*> shift,
      // Scalar temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1_scalar,
      const gsl::not_null<Scalar<DataVector>*> gamma2_scalar,
      // Scalar driver temporary variables
      const gsl::not_null<Scalar<DataVector>*> gamma1_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> gamma2_scalar_driver,
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
      const tnsr::I<DataVector, 3_st>& shift_interior,
      // Scalar driver interior variables
      const Scalar<DataVector>& gamma1_interior_scalar_driver,
      const Scalar<DataVector>& gamma2_interior_scalar_driver) const {
    // Note: Check that CurvedScalarWave does not update GH variables
    // to a different value. If it does, invert the order of application of the
    // corrections first, so that the GH update is applied at last

    // gh::BoundaryConditions::DirichletMinkowski
    auto gh_string = derived_gh_condition_.dg_ghost(
        // GH evolved variables
        spacetime_metric, pi, phi,
        // Scalar evolved variables
        psi_scalar, pi_scalar, phi_scalar,
        // GH temporary variables
        gamma1, gamma2, lapse, shift,
        // Scalar temporary variables
        gamma1_scalar, gamma2_scalar,
        // Inverse metric
        inv_spatial_metric,
        // Mesh variables
        face_mesh_velocity, normal_covector, normal_vector,
        // GH interior variables
        interior_gamma1, interior_gamma2,
        // Scalar interior variables
        inverse_spatial_metric_interior, gamma1_interior_scalar,
        gamma2_interior_scalar, lapse_interior, shift_interior);

    // CurvedScalarWave::BoundaryConditions::AnalyticConstant
    auto scalar_string = derived_scalar_condition_.dg_ghost(
        // Change names
        psi_scalar_driver, pi_scalar_driver, phi_scalar_driver, lapse, shift,
        gamma1_scalar_driver, gamma2_scalar_driver, inv_spatial_metric,
        face_mesh_velocity, normal_covector, normal_vector,
        inverse_spatial_metric_interior, gamma1_interior_scalar_driver,
        gamma2_interior_scalar_driver, lapse_interior, shift_interior);
    if (not gh_string.has_value()) {
      return scalar_string;
    }
    if (not scalar_string.has_value()) {
      return gh_string;
    }
    return gh_string.value() + ";" + scalar_string.value();
  }

  std::optional<std::string> dg_time_derivative(
      // GH
      gsl::not_null<tnsr::aa<DataVector, dim, Frame::Inertial>*>
          dt_spacetime_metric_correction,
      gsl::not_null<tnsr::aa<DataVector, dim, Frame::Inertial>*>
          dt_pi_correction,
      gsl::not_null<tnsr::iaa<DataVector, dim, Frame::Inertial>*>
          dt_phi_correction,
      // Scalar
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar_correction,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar_correction,
      gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*>
          dt_phi_scalar_correction,
      // Scalar Driver
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar_driver_correction,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar_driver_correction,
      gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*>
          dt_phi_scalar_driver_correction,

      const std::optional<tnsr::I<DataVector, dim, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, dim, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, dim, Frame::Inertial>& normal_vector,

      // c.f. GH dg_interior_evolved_variables_tags
      const tnsr::aa<DataVector, dim, Frame::Inertial>& spacetime_metric,
      const tnsr::aa<DataVector, dim, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, dim, Frame::Inertial>& phi,
      // Scalar evolved variables
      const Scalar<DataVector>& psi_scalar,
      const tnsr::i<DataVector, dim, Frame::Inertial>& phi_scalar,
      // Scalar Driver evolved variables
      const Scalar<DataVector>& psi_scalar_driver,
      const tnsr::i<DataVector, dim, Frame::Inertial>& phi_scalar_driver,

      // c.f. GH dg_interior_temporary_tags
      const tnsr::I<DataVector, dim, Frame::Inertial>& coords,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, dim, Frame::Inertial>& shift,
      const tnsr::AA<DataVector, dim, Frame::Inertial>&
          inverse_spacetime_metric,
      const tnsr::A<DataVector, dim, Frame::Inertial>&
          spacetime_unit_normal_vector,
      const tnsr::iaa<DataVector, dim, Frame::Inertial>& three_index_constraint,
      const tnsr::a<DataVector, dim, Frame::Inertial>& gauge_source,
      const tnsr::ab<DataVector, dim, Frame::Inertial>&
          spacetime_deriv_gauge_source,
      // Scalar interior temporary
      const Scalar<DataVector>& gamma1_scalar,
      const Scalar<DataVector>& gamma2_scalar,
      // Scalar Driver interior temporary
      const Scalar<DataVector>& gamma1_scalar_driver,
      const Scalar<DataVector>& gamma2_scalar_driver,

      // c.f. dg_interior_dt_vars_tags
      const tnsr::aa<DataVector, dim, Frame::Inertial>&
          logical_dt_spacetime_metric,
      const tnsr::aa<DataVector, dim, Frame::Inertial>& logical_dt_pi,
      const tnsr::iaa<DataVector, dim, Frame::Inertial>& logical_dt_phi,
      // Scalar interior dt tags
      const Scalar<DataVector>& logical_dt_psi_scalar,
      const Scalar<DataVector>& logical_dt_pi_scalar,
      const tnsr::i<DataVector, dim>& logical_dt_phi_scalar,
      // Scalar Driver interior dt tags
      const Scalar<DataVector>& logical_dt_psi_scalar_driver,
      const Scalar<DataVector>& logical_dt_pi_scalar_driver,
      const tnsr::i<DataVector, dim>& logical_dt_phi_scalar_driver,

      // c.f. GH dg_interior_deriv_vars_tags
      const tnsr::iaa<DataVector, dim, Frame::Inertial>& d_spacetime_metric,
      const tnsr::iaa<DataVector, dim, Frame::Inertial>& d_pi,
      const tnsr::ijaa<DataVector, dim, Frame::Inertial>& d_phi,
      // Scalar deriv vars
      const tnsr::i<DataVector, dim, Frame::Inertial>& d_psi_scalar,
      const tnsr::i<DataVector, dim, Frame::Inertial>& d_pi_scalar,
      const tnsr::ij<DataVector, dim, Frame::Inertial>& d_phi_scalar,
      // Scalar Driver deriv vars
      const tnsr::i<DataVector, dim, Frame::Inertial>& d_psi_scalar_driver,
      const tnsr::i<DataVector, dim, Frame::Inertial>& d_pi_scalar_driver,
      const tnsr::ij<DataVector, dim, Frame::Inertial>& d_phi_scalar_driver)
      const {
    // GH Bjorus boundary condition
    auto gh_string = derived_gh_condition_.dg_time_derivative(
        // GH
        dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
        // Scalar
        dt_psi_scalar_correction, dt_pi_scalar_correction,
        dt_phi_scalar_correction, face_mesh_velocity, normal_covector,
        normal_vector,
        // c.f. GH dg_interior_evolved_variables_tags
        spacetime_metric, pi, phi,
        // Scalar evolved variables
        psi_scalar, phi_scalar,
        // c.f. GH dg_interior_temporary_tags
        coords, gamma1, gamma2, lapse, shift, inverse_spacetime_metric,
        spacetime_unit_normal_vector, three_index_constraint, gauge_source,
        spacetime_deriv_gauge_source,
        // Scalar interior temporary
        gamma1_scalar, gamma2_scalar,
        // c.f. dg_interior_dt_vars_tags
        logical_dt_spacetime_metric, logical_dt_pi, logical_dt_phi,
        // Scalar interior dt tags
        logical_dt_psi_scalar, logical_dt_pi_scalar, logical_dt_phi_scalar,
        // c.f. GH dg_interior_deriv_vars_tags
        d_spacetime_metric, d_pi, d_phi,
        // Scalar deriv vars
        d_psi_scalar, d_pi_scalar, d_phi_scalar);

    // Scalar ConstraintPreservingSphericalRadiation boundary conditions
    auto scalar_string = derived_scalar_condition_.dg_time_derivative(
        dt_psi_scalar_driver_correction, dt_pi_scalar_driver_correction,
        dt_phi_scalar_driver_correction, face_mesh_velocity, normal_covector,
        normal_vector, psi_scalar_driver, phi_scalar_driver, coords,
        gamma1_scalar_driver, gamma2_scalar_driver, lapse, shift,
        logical_dt_psi_scalar_driver, logical_dt_pi_scalar_driver,
        logical_dt_phi_scalar_driver, d_psi_scalar_driver, d_pi_scalar_driver,
        d_phi_scalar_driver);

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

}  // namespace fe::DecoupledScalar::BoundaryConditions
