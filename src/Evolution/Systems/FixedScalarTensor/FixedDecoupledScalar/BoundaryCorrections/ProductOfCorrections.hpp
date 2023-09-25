// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>

#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

namespace fe::DecoupledScalar::BoundaryCorrections {

/*!
 * \brief Apply a boundary condition to the combined Generalized Harmonic (GH)
 * and scalar field system using boundary corrections defined separately.
 */
template <typename DerivedGhCorrection, typename DerivedScalarCorrection>
class ProductOfCorrections final : public BoundaryCorrection {
 public:
  // Need to remove duplicates in the following unions of tags
  using dg_package_field_tags = tmpl::remove_duplicates<
      tmpl::append<typename DerivedGhCorrection::dg_package_field_tags,
                   typename DerivedScalarCorrection::dg_package_field_tags>>;

  using dg_package_data_temporary_tags = tmpl::remove_duplicates<tmpl::append<
      typename DerivedGhCorrection::dg_package_data_temporary_tags,
      typename DerivedScalarCorrection::dg_package_data_temporary_tags>>;

  using dg_package_data_volume_tags = tmpl::remove_duplicates<tmpl::append<
      typename DerivedGhCorrection::dg_package_data_volume_tags,
      typename DerivedScalarCorrection::dg_package_data_volume_tags>>;

  static std::string name() {
    return "Product" + pretty_type::name<DerivedGhCorrection>() +
           "ScalarTensor" + "And" +
           pretty_type::name<DerivedScalarCorrection>() + "ScalarDriver";
  }

  struct GhCorrection {
    using type = DerivedGhCorrection;
    static std::string name() {
      // The boundary conditions for the different systems cannot have the same
      // name, otherwise errors during option parsing will occur
      return pretty_type::name<DerivedGhCorrection>() + "ScalarTensor";
    }
    static constexpr Options::String help{
        "The ScalarTensor part of the product boundary condition"};
  };
  struct ScalarCorrection {
    using type = DerivedScalarCorrection;
    static std::string name() {
      // The boundary conditions for the different systems cannot have the same
      // name, otherwise errors during option parsing will occur
      return pretty_type::name<DerivedScalarCorrection>() + "ScalarDriver";
    }
    static constexpr Options::String help{
        "The Scalar Driver part of the product boundary condition"};
  };

  using options = tmpl::list<GhCorrection, ScalarCorrection>;

  static constexpr Options::String help = {
      "Direct product of a ScalarTensor and ScalarDriver boundary correction. "
      "See the documentation for the two individual boundary corrections for "
      "further details."};

  ProductOfCorrections() = default;
  ProductOfCorrections(DerivedGhCorrection gh_correction,
                       DerivedScalarCorrection scalar_correction)
      : derived_gh_correction_{gh_correction},
        derived_scalar_correction_{scalar_correction} {}
  ProductOfCorrections(const ProductOfCorrections&) = default;
  ProductOfCorrections& operator=(const ProductOfCorrections&) = default;
  ProductOfCorrections(ProductOfCorrections&&) = default;
  ProductOfCorrections& operator=(ProductOfCorrections&&) = default;
  ~ProductOfCorrections() override = default;

  /// \cond
  explicit ProductOfCorrections(CkMigrateMessage* msg)
      : BoundaryCorrection(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ProductOfCorrections);  // NOLINT
  /// \endcond
  void pup(PUP::er& p) override {
    BoundaryCorrection::pup(p);
    p | derived_gh_correction_;
    p | derived_scalar_correction_;
  }

  std::unique_ptr<BoundaryCorrection> get_clone() const override {
    return std::make_unique<ProductOfCorrections>(*this);
  }

  double dg_package_data(
      // GH packaged fields
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_v_spacetime_metric,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_v_zero,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_v_plus,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_v_minus,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_n_times_v_plus,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_n_times_v_minus,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          packaged_char_speed_gamma2_v_spacetime_metric,
      const gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds,
      // Scalar packaged fields
      const gsl::not_null<Scalar<DataVector>*> packaged_v_psi_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          packaged_v_zero_scalar,
      const gsl::not_null<Scalar<DataVector>*> packaged_v_plus_scalar,
      const gsl::not_null<Scalar<DataVector>*> packaged_v_minus_scalar,
      const gsl::not_null<Scalar<DataVector>*> packaged_gamma2_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          packaged_interface_unit_normal_scalar,
      const gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds_scalar,
      // Scalar driver packaged fields
      const gsl::not_null<Scalar<DataVector>*> packaged_v_psi_scalar_driver,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          packaged_v_zero_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> packaged_v_plus_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> packaged_v_minus_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*> packaged_gamma2_scalar_driver,
    //   gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
    //       packaged_interface_unit_normal_scalar_driver,
      const gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds_scalar_driver,
      // GH variables
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>& spacetime_metric,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>& phi,
      // Scalar variables
      const Scalar<DataVector>& psi_scalar, const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi_scalar,
      // Scalar variables
      const Scalar<DataVector>& psi_scalar_driver,
      const Scalar<DataVector>& pi_scalar_driver,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi_scalar_driver,
      // GH fluxes
      // Scalar fluxes
      // GH temporaries
      const Scalar<DataVector>& constraint_gamma1,
      const Scalar<DataVector>& constraint_gamma2,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
      // Scalar temporaries (without repeating tags)
      // const Scalar<DataVector>& lapse,
      // const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,

      const Scalar<DataVector>& constraint_gamma1_scalar,
      const Scalar<DataVector>& constraint_gamma2_scalar,
      // Scalar driver temporaries (without repeating tags)
      const Scalar<DataVector>& constraint_gamma1_scalar_driver,
      const Scalar<DataVector>& constraint_gamma2_scalar_driver,
      // Mesh variables
      const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
      const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
          mesh_velocity,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity
      // GH volume quantities
      // Scalar volume quantities
      // Scalar driver volume quantities
  ) const {
    // Note: Check that CurvedScalarWave does not update GH variables
    // to a different value. If it does, invert the order of application of the
    // corrections first, so that the GH update is applied at last
    const double gh_correction_result = derived_gh_correction_.dg_package_data(
        // GH packaged fields
        packaged_char_speed_v_spacetime_metric, packaged_char_speed_v_zero,
        packaged_char_speed_v_plus, packaged_char_speed_v_minus,
        packaged_char_speed_n_times_v_plus, packaged_char_speed_n_times_v_minus,
        packaged_char_speed_gamma2_v_spacetime_metric, packaged_char_speeds,
        // Scalar packaged fields
        packaged_v_psi_scalar, packaged_v_zero_scalar, packaged_v_plus_scalar,
        packaged_v_minus_scalar, packaged_gamma2_scalar,
        packaged_interface_unit_normal_scalar, packaged_char_speeds_scalar,
        // GH variables
        spacetime_metric, pi, phi,
        // Scalar variables
        psi_scalar, pi_scalar, phi_scalar,
        // GH fluxes
        // Scalar fluxes
        // GH temporaries
        constraint_gamma1, constraint_gamma2, lapse, shift,
        // Scalar temporaries (without repeating tags)
        constraint_gamma1_scalar, constraint_gamma2_scalar,
        // Mesh variables
        normal_covector, normal_vector, mesh_velocity, normal_dot_mesh_velocity
        // GH volume quantities
        // Scalar volume quantities
    );

    const double scalar_correction_result =
        derived_scalar_correction_.dg_package_data(
            // Scalar packaged variables
            packaged_v_psi_scalar_driver, packaged_v_zero_scalar_driver,
            packaged_v_plus_scalar_driver, packaged_v_minus_scalar_driver,
            packaged_gamma2_scalar_driver,
            // packaged_interface_unit_normal_scalar_driver,
            packaged_interface_unit_normal_scalar,

            packaged_char_speeds_scalar_driver,
            // Scalar variables
            psi_scalar_driver, pi_scalar_driver, phi_scalar_driver,
            // Scalar temporaries
            lapse, shift,

            constraint_gamma1_scalar_driver, constraint_gamma2_scalar_driver,
            // Scalar mesh variables
            normal_covector, normal_vector, mesh_velocity,
            normal_dot_mesh_velocity);
    return std::max(gh_correction_result, scalar_correction_result);
  }

  void dg_boundary_terms(
      // GH boundary corrections
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          boundary_correction_spacetime_metric,
      const gsl::not_null<tnsr::aa<DataVector, 3_st, Frame::Inertial>*>
          boundary_correction_pi,
      const gsl::not_null<tnsr::iaa<DataVector, 3_st, Frame::Inertial>*>
          boundary_correction_phi,
      // Scalar boundary corrections
      const gsl::not_null<Scalar<DataVector>*> psi_boundary_correction_scalar,
      const gsl::not_null<Scalar<DataVector>*> pi_boundary_correction_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_boundary_correction_scalar,
      // Scalar boundary corrections
      const gsl::not_null<Scalar<DataVector>*>
          psi_boundary_correction_scalar_driver,
      const gsl::not_null<Scalar<DataVector>*>
          pi_boundary_correction_scalar_driver,
      const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_boundary_correction_scalar_driver,
      // GH internal packages field tags
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>&
          char_speed_v_spacetime_metric_int,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>& char_speed_v_zero_int,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>& char_speed_v_plus_int,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>& char_speed_v_minus_int,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>&
          char_speed_normal_times_v_plus_int,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>&
          char_speed_normal_times_v_minus_int,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>&
          char_speed_constraint_gamma2_v_spacetime_metric_int,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int,
      // Scalar internal packaged field tags
      const Scalar<DataVector>& v_psi_int_scalar,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero_int_scalar,
      const Scalar<DataVector>& v_plus_int_scalar,
      const Scalar<DataVector>& v_minus_int_scalar,
      const Scalar<DataVector>& gamma2_int_scalar,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          interface_unit_normal_int_scalar,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int_scalar,
      // Scalar driver internal packaged field tags
      const Scalar<DataVector>& v_psi_int_scalar_driver,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          v_zero_int_scalar_driver,
      const Scalar<DataVector>& v_plus_int_scalar_driver,
      const Scalar<DataVector>& v_minus_int_scalar_driver,
      const Scalar<DataVector>& gamma2_int_scalar_driver,
      //   const tnsr::i<DataVector, 3_st, Frame::Inertial>&
      //       interface_unit_normal_int_scalar_driver,
      const tnsr::a<DataVector, 3, Frame::Inertial>&
          char_speeds_int_scalar_driver,
      // GH external packaged fields
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>&
          char_speed_v_spacetime_metric_ext,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>& char_speed_v_zero_ext,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>& char_speed_v_plus_ext,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>& char_speed_v_minus_ext,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>&
          char_speed_normal_times_v_plus_ext,
      const tnsr::iaa<DataVector, 3_st, Frame::Inertial>&
          char_speed_normal_times_v_minus_ext,
      const tnsr::aa<DataVector, 3_st, Frame::Inertial>&
          char_speed_constraint_gamma2_v_spacetime_metric_ext,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext,
      // Scalar external packaged fields
      const Scalar<DataVector>& v_psi_ext_scalar,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero_ext_scalar,
      const Scalar<DataVector>& v_plus_ext_scalar,
      const Scalar<DataVector>& v_minus_ext_scalar,
      const Scalar<DataVector>& gamma2_ext_scalar,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          interface_unit_normal_ext_scalar,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext_scalar,
      // Scalar driver external packaged fields
      const Scalar<DataVector>& v_psi_ext_scalar_driver,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          v_zero_ext_scalar_driver,
      const Scalar<DataVector>& v_plus_ext_scalar_driver,
      const Scalar<DataVector>& v_minus_ext_scalar_driver,
      const Scalar<DataVector>& gamma2_ext_scalar_driver,
      //   const tnsr::i<DataVector, 3_st, Frame::Inertial>&
      //       interface_unit_normal_ext_scalar_driver,
      const tnsr::a<DataVector, 3, Frame::Inertial>&
          char_speeds_ext_scalar_driver,
      // DG formulation
      const dg::Formulation dg_formulation) const {
    // Note: Check that CurvedScalarWave does not update GH variables
    // to a different value. If it does, invert the order of application of the
    // corrections first, so that the GH update is applied at last
    derived_gh_correction_.dg_boundary_terms(
        // GH boundary corrections
        boundary_correction_spacetime_metric, boundary_correction_pi,
        boundary_correction_phi,
        // Scalar boundary corrections
        psi_boundary_correction_scalar, pi_boundary_correction_scalar,
        phi_boundary_correction_scalar,
        // GH internal packages field tags
        char_speed_v_spacetime_metric_int, char_speed_v_zero_int,
        char_speed_v_plus_int, char_speed_v_minus_int,
        char_speed_normal_times_v_plus_int, char_speed_normal_times_v_minus_int,
        char_speed_constraint_gamma2_v_spacetime_metric_int, char_speeds_int,
        // Scalar internal packaged field tags
        v_psi_int_scalar, v_zero_int_scalar, v_plus_int_scalar,
        v_minus_int_scalar, gamma2_int_scalar, interface_unit_normal_int_scalar,
        char_speeds_int_scalar,
        // GH external packaged fields
        char_speed_v_spacetime_metric_ext, char_speed_v_zero_ext,
        char_speed_v_plus_ext, char_speed_v_minus_ext,
        char_speed_normal_times_v_plus_ext, char_speed_normal_times_v_minus_ext,
        char_speed_constraint_gamma2_v_spacetime_metric_ext, char_speeds_ext,
        // Scalar external packaged fields
        v_psi_ext_scalar, v_zero_ext_scalar, v_plus_ext_scalar,
        v_minus_ext_scalar, gamma2_ext_scalar, interface_unit_normal_ext_scalar,
        char_speeds_ext_scalar,
        // DG formulation
        dg_formulation);

    derived_scalar_correction_.dg_boundary_terms(
        // scalar_boundary_corrections...,
        psi_boundary_correction_scalar_driver,
        pi_boundary_correction_scalar_driver,
        phi_boundary_correction_scalar_driver,

        // scalar_internal_packaged_fields...,
        v_psi_int_scalar_driver, v_zero_int_scalar_driver,
        v_plus_int_scalar_driver, v_minus_int_scalar_driver,
        gamma2_int_scalar_driver, interface_unit_normal_int_scalar,
        char_speeds_int_scalar_driver,

        // scalar_external_packaged_fields...,
        v_psi_ext_scalar_driver, v_zero_ext_scalar_driver,
        v_plus_ext_scalar_driver, v_minus_ext_scalar_driver,
        gamma2_ext_scalar_driver, interface_unit_normal_ext_scalar,
        char_speeds_ext_scalar_driver,

        dg_formulation);
  }

  const DerivedGhCorrection& gh_correction() const {
    return derived_gh_correction_;
  }

  const DerivedScalarCorrection& scalar_correction() const {
    return derived_scalar_correction_;
  }

 private:
  DerivedGhCorrection derived_gh_correction_;
  DerivedScalarCorrection derived_scalar_correction_;
};

/// \cond
template <typename DerivedGhCorrection, typename DerivedScalarCorrection>
PUP::able::PUP_ID ProductOfCorrections<DerivedGhCorrection,
                                       DerivedScalarCorrection>::my_PUP_ID =
    0;  // NOLINT
/// \endcond
}  // namespace fe::DecoupledScalar::BoundaryCorrections
