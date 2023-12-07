// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/String.hpp"
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

namespace ScalarWave::BoundaryCorrections {
/*!
 * \brief Computes a simple penalty applied directly to the evolved fields and
 * without a characteristic decomposition.
 */
template <size_t Dim>
class SimplePenalty final : public BoundaryCorrection<Dim> {
 private:
  struct NormalTimesVPlus : db::SimpleTag {
    using type = tnsr::i<DataVector, Dim, Frame::Inertial>;
  };
  struct NormalTimesVMinus : db::SimpleTag {
    using type = tnsr::i<DataVector, Dim, Frame::Inertial>;
  };
  struct Gamma2VPsi : db::SimpleTag {
    using type = Scalar<DataVector>;
  };
  struct CharSpeedsTensor : db::SimpleTag {
    using type = tnsr::i<DataVector, 3, Frame::Inertial>;
  };
  double penalty_factor_ = std::numeric_limits<double>::signaling_NaN();

 public:
  struct PenaltyFactor {
    using type = double;
    static constexpr Options::String help = {"Amplitude of penalty factor"};
  };

  using options = tmpl::list<PenaltyFactor>;
  static constexpr Options::String help = {
      "Computes the SimplePenalty boundary correction term for the scalar wave "
      "system."};

  SimplePenalty(double penalty_factor);
  SimplePenalty() = default;
  SimplePenalty(const SimplePenalty&) = default;
  SimplePenalty& operator=(const SimplePenalty&) = default;
  SimplePenalty(SimplePenalty&&) = default;
  SimplePenalty& operator=(SimplePenalty&&) = default;
  ~SimplePenalty() override = default;

  /// \cond
  explicit SimplePenalty(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(SimplePenalty);  // NOLINT
  /// \endcond
  void pup(PUP::er& p) override;  // NOLINT

  std::unique_ptr<BoundaryCorrection<Dim>> get_clone() const override;

  using dg_package_field_tags =
      tmpl::list<Tags::VPsi, Tags::VZero<Dim>, Tags::VPlus, Tags::VMinus,
                 NormalTimesVPlus, NormalTimesVMinus, Gamma2VPsi,
                 CharSpeedsTensor>;
  using dg_package_data_temporary_tags = tmpl::list<Tags::ConstraintGamma2>;
  using dg_package_data_volume_tags = tmpl::list<>;

  double dg_package_data(
      gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_psi,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          packaged_char_speed_v_zero,
      gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_plus,
      gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_minus,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          packaged_char_speed_n_times_v_plus,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          packaged_char_speed_n_times_v_minus,
      gsl::not_null<Scalar<DataVector>*> packaged_char_speed_gamma2_v_psi,
      gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds,

      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,

      const Scalar<DataVector>& constraint_gamma2,

      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*mesh_velocity*/,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const;

  void dg_boundary_terms(
      gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,
      gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          phi_boundary_correction,

      const Scalar<DataVector>& char_speed_v_psi_int,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& char_speed_v_zero_int,
      const Scalar<DataVector>& char_speed_v_plus_int,
      const Scalar<DataVector>& char_speed_v_minus_int,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          char_speed_normal_times_v_plus_int,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          char_speed_normal_times_v_minus_int,
      const Scalar<DataVector>& char_speed_constraint_gamma2_v_psi_int,
      const tnsr::i<DataVector, 3, Frame::Inertial>& char_speeds_int,

      const Scalar<DataVector>& char_speed_v_psi_ext,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& char_speed_v_zero_ext,
      const Scalar<DataVector>& char_speed_v_plus_ext,
      const Scalar<DataVector>& char_speed_v_minus_ext,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          char_speed_minus_normal_times_v_plus_ext,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          char_speed_minus_normal_times_v_minus_ext,
      const Scalar<DataVector>& char_speed_constraint_gamma2_v_psi_ext,
      const tnsr::i<DataVector, 3, Frame::Inertial>& char_speeds_ext,
      dg::Formulation /*dg_formulation*/) const;
};
}  // namespace ScalarWave::BoundaryCorrections
