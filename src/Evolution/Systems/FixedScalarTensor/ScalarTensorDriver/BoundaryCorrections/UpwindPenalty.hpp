// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
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

namespace fe::ScalarTensorDriver::BoundaryCorrections {
/*!
 * \brief Computes the generalized harmonic upwind multipenalty boundary
 * correction.
 */
class UpwindPenalty final : public BoundaryCorrection {
 private:
  struct CharSpeedsTensor : db::SimpleTag {
    using type = tnsr::a<DataVector, 3, Frame::Inertial>;
  };

 public:
  using options = tmpl::list<>;
  static constexpr Options::String help = {
      "Computes the UpwindPenalty boundary correction term for the scalar "
      "tensor driver system."};

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
      tmpl::list<Tags::VScalarDriver<DataVector>, VPiScalar<DataVector>,
                 Tags::VTensorDriver<DataVector, 3>, Tags::VPi<DataVector, 3>,
                 //  ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<
                 //      3_st, Frame::Inertial>>,
                 CharSpeedsTensor>;
  using dg_package_data_temporary_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>>;
  using dg_package_data_primitive_tags = tmpl::list<>;
  using dg_package_data_volume_tags = tmpl::list<>;

  double dg_package_data(
      gsl::not_null<Scalar<DataVector>*> packaged_v_scalar_driver,
      gsl::not_null<Scalar<DataVector>*> packaged_v_pi_scalar,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
          packaged_v_tensor_driver,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> packaged_v_pi,
      //   gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
      //       packaged_interface_unit_normal,
      gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds,

      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,

      const tnsr::aa<DataVector, 3, Frame::Inertial>& tensor_driver,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,

      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3, Frame::Inertial>& shift,

      const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
      const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
          mesh_velocity,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const;

  void dg_boundary_terms(
      gsl::not_null<Scalar<DataVector, 3, Frame::Inertial>*>
          boundary_correction_scalar_driver,
      gsl::not_null<Scalar<DataVector, 3, Frame::Inertial>*>
          boundary_correction_pi_scalar,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
          boundary_correction_tensor_driver,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
          boundary_correction_pi,

      const Scalar<DataVector>& char_speed_v_scalar_driver_int,
      const Scalar<DataVector>& char_speed_v_pi_scalar_int,
      const tnsr::aa<DataVector, 3, Frame::Inertial>&
          char_speed_v_tensor_driver_int,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& char_speed_v_pi_int,

      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int,

      const Scalar<DataVector>& char_speed_v_scalar_driver_ext,
      const Scalar<DataVector>& char_speed_v_pi_scalar_ext,
      const tnsr::aa<DataVector, 3, Frame::Inertial>&
          char_speed_v_tensor_driver_ext,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& char_speed_v_pi_ext,

      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext,
      dg::Formulation /*dg_formulation*/) const;
};

bool operator==(const UpwindPenalty& lhs, const UpwindPenalty& rhs);

bool operator!=(const UpwindPenalty& lhs, const UpwindPenalty& rhs);
}  // namespace fe::ScalarTensorDriver::BoundaryCorrections
