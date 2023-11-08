// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
//
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
//
#include "Evolution/Systems/CurvedScalarWave/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Characteristics.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace fe::ScalarDriver::BoundaryCorrections {
class UpwindPenalty final : public BoundaryCorrection {
 private:
  struct CharSpeedsTensor : db::SimpleTag {
    using type = tnsr::a<DataVector, 3, Frame::Inertial>;
  };

  CurvedScalarWave::BoundaryCorrections::UpwindPenalty<3_st>
      boundary_correction_for_scalar_ =
          CurvedScalarWave::BoundaryCorrections::UpwindPenalty<3_st>();

 public:
  using options = tmpl::list<>;
  static constexpr Options::String help = {
      "Computes the UpwindPenalty boundary correction term for the scalar wave "
      "system in curved spacetime."};

  UpwindPenalty() = default;
  UpwindPenalty(const UpwindPenalty&) = default;
  UpwindPenalty& operator=(const UpwindPenalty&) = default;
  UpwindPenalty(UpwindPenalty&&) = default;
  UpwindPenalty& operator=(UpwindPenalty&&) = default;
  ~UpwindPenalty() override = default;

  /// \cond
  explicit UpwindPenalty(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(UpwindPenalty);  // NOLINT
  /// \endcond
  void pup(PUP::er& p) override;  // NOLINT

  std::unique_ptr<BoundaryCorrection> get_clone() const override;

  using dg_package_field_tags =
      tmpl::list<Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus, Tags::VMinus,
                 Tags::ConstraintGamma2,
                 ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<
                     3_st, Frame::Inertial>>,
                 CharSpeedsTensor>;
  using dg_package_data_temporary_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3_st>,
                 Tags::ConstraintGamma1, Tags::ConstraintGamma2>;
  using dg_package_data_volume_tags = tmpl::list<>;

  double dg_package_data(
      gsl::not_null<Scalar<DataVector>*> packaged_v_psi,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          packaged_v_zero,
      gsl::not_null<Scalar<DataVector>*> packaged_v_plus,
      gsl::not_null<Scalar<DataVector>*> packaged_v_minus,
      gsl::not_null<Scalar<DataVector>*> packaged_gamma2,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          packaged_interface_unit_normal,
      gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds,

      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,

      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
      const Scalar<DataVector>& constraint_gamma1,
      const Scalar<DataVector>& constraint_gamma2,

      const tnsr::i<DataVector, 3_st, Frame::Inertial>& interface_unit_normal,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>&
          interface_unit_normal_vector,
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          mesh_velocity,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
    // Use the CurvedScalarWave routines
    return boundary_correction_for_scalar_.dg_package_data(
        packaged_v_psi, packaged_v_zero, packaged_v_plus, packaged_v_minus,
        packaged_gamma2, packaged_interface_unit_normal, packaged_char_speeds,

        psi, pi, phi,

        lapse, shift, constraint_gamma1, constraint_gamma2,

        interface_unit_normal, interface_unit_normal_vector, mesh_velocity,
        normal_dot_mesh_velocity);
  }

  void dg_boundary_terms(
      gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,
      gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
          phi_boundary_correction,

      const Scalar<DataVector>& v_psi_int,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero_int,
      const Scalar<DataVector>& v_plus_int,
      const Scalar<DataVector>& v_minus_int,
      const Scalar<DataVector>& gamma2_int,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          interface_unit_normal_int,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int,

      const Scalar<DataVector>& v_psi_ext,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero_ext,
      const Scalar<DataVector>& v_plus_ext,
      const Scalar<DataVector>& v_minus_ext,
      const Scalar<DataVector>& gamma2_ext,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&
          interface_unit_normal_ext,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext,
      dg::Formulation dg_formulation) const {
    // Use the CurvedScalarWave routines
    CurvedScalarWave::BoundaryCorrections::UpwindPenalty<3_st>
        boundary_correction_for_scalar;
    boundary_correction_for_scalar_.dg_boundary_terms(
        psi_boundary_correction, pi_boundary_correction,
        phi_boundary_correction,

        v_psi_int, v_zero_int, v_plus_int, v_minus_int, gamma2_int,

        interface_unit_normal_int, char_speeds_int,

        v_psi_ext, v_zero_ext, v_plus_ext, v_minus_ext, gamma2_ext,

        interface_unit_normal_ext, char_speeds_ext, dg_formulation);

    // Since the principal part of the psi equation doesn't really mix
    // with that of the gradient equations, we simply turn off the boundary
    // corrections for the gradient fields. In this way we keep the shift in the
    // equation, and the mesh velocity correction.

    // Set to zero the pi and phi boundary corrections
    get(*pi_boundary_correction) = 0.0 * get(*psi_boundary_correction);
    for (size_t i = 0; i < 3_st; ++i) {
      phi_boundary_correction->get(i) = 0.0 * interface_unit_normal_int.get(i);
    }
  }
};
}  // namespace fe::ScalarDriver::BoundaryCorrections
