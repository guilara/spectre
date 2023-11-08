// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FaceNormal.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/ConstraintPreservingSphericalRadiation.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
/// \endcond

namespace fe::ScalarDriver::BoundaryConditions {
/// A `BoundaryCondition` that imposes the scalar to be the zero at the
/// outer boundary.
class ConstraintPreservingSphericalRadiation final : public BoundaryCondition {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Constraint-preserving boundary conditions with a second order "
      "Bayliss-Turkel radiation boundary condition."};

  ConstraintPreservingSphericalRadiation();
  ConstraintPreservingSphericalRadiation() = default;
  /// \cond
  ConstraintPreservingSphericalRadiation(
      ConstraintPreservingSphericalRadiation&&) = default;
  ConstraintPreservingSphericalRadiation& operator=(
      ConstraintPreservingSphericalRadiation&&) = default;
  ConstraintPreservingSphericalRadiation(
      const ConstraintPreservingSphericalRadiation&) = default;
  ConstraintPreservingSphericalRadiation& operator=(
      const ConstraintPreservingSphericalRadiation&) = default;
  /// \endcond
  ~ConstraintPreservingSphericalRadiation() override = default;

  explicit ConstraintPreservingSphericalRadiation(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition,
      ConstraintPreservingSphericalRadiation);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::TimeDerivative;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags =
      tmpl::list<Tags::Psi, Tags::Phi<3>>;
  using dg_interior_temporary_tags =
      tmpl::list<domain::Tags::Coordinates<3, Frame::Inertial>,
                 Tags::ConstraintGamma1, Tags::ConstraintGamma2,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>>;
  using dg_interior_dt_vars_tags =
      tmpl::list<::Tags::dt<Tags::Psi>, ::Tags::dt<Tags::Pi>,
                 ::Tags::dt<Tags::Phi<3>>>;
  using dg_interior_deriv_vars_tags =
      tmpl::list<::Tags::deriv<Tags::Psi, tmpl::size_t<3>, Frame::Inertial>,
                 ::Tags::deriv<Tags::Pi, tmpl::size_t<3>, Frame::Inertial>,
                 ::Tags::deriv<Tags::Phi<3>, tmpl::size_t<3>, Frame::Inertial>>;
  using dg_gridless_tags = tmpl::list<>;

  std::optional<std::string> dg_time_derivative(
      gsl::not_null<Scalar<DataVector>*> dt_psi_correction,
      gsl::not_null<Scalar<DataVector>*> dt_pi_correction,
      gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*> dt_phi_correction,
      const std::optional<tnsr::I<DataVector, 3>>& face_mesh_velocity,
      const tnsr::i<DataVector, 3>& normal_covector,
      const tnsr::I<DataVector, 3>& normal_vector,
      const Scalar<DataVector>& psi, const tnsr::i<DataVector, 3>& phi,
      const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
      const Scalar<DataVector>& logical_dt_psi,
      const Scalar<DataVector>& logical_dt_pi,
      const tnsr::i<DataVector, 3>& logical_dt_phi,
      const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
      const tnsr::ij<DataVector, 3>& d_phi) const {
    // Use the boundary condition from CurvedScalarWave
    auto fe_string = csw_constraint_preserving_.dg_time_derivative(
        dt_psi_correction, dt_pi_correction, dt_phi_correction,
        face_mesh_velocity, normal_covector, normal_vector, psi, phi, coords,
        gamma1, gamma2, lapse, shift, logical_dt_psi, logical_dt_pi,
        logical_dt_phi, d_psi, d_pi, d_phi);

    if (not fe_string.has_value()) {
      return fe_string;
    }
    return std::nullopt;
  }

 private:
  CurvedScalarWave::BoundaryConditions::ConstraintPreservingSphericalRadiation<
      3_st>
      csw_constraint_preserving_;
};
}  // namespace fe::ScalarDriver::BoundaryConditions
