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
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
/// \endcond

namespace fe::ScalarTensorDriver::BoundaryConditions {
/// A `BoundaryCondition` that imposes the scalar to be the zero at the
/// outer boundary.
class AnalyticConstant final : public BoundaryCondition {
 public:
  struct Amplitude {
    using type = double;
    static constexpr Options::String help = {
        "Amplitude of the scalar at the boundary"};
  };
  using options = tmpl::list<Amplitude>;
  static constexpr Options::String help{
      "Boundary conditions which check that all characteristic "
      "fields are outflowing."};

  AnalyticConstant(double amplitude);
  AnalyticConstant() = default;
  /// \cond
  AnalyticConstant(AnalyticConstant&&) = default;
  AnalyticConstant& operator=(AnalyticConstant&&) = default;
  AnalyticConstant(const AnalyticConstant&) = default;
  AnalyticConstant& operator=(const AnalyticConstant&) = default;
  /// \endcond
  ~AnalyticConstant() override = default;

  explicit AnalyticConstant(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, AnalyticConstant);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags = tmpl::list<
      gr::Tags::InverseSpatialMetric<DataVector, 3_st, Frame::Inertial>,
      //   Tags::ConstraintGamma1, Tags::ConstraintGamma2,
      gr::Tags::Lapse<DataVector>,
      gr::Tags::Shift<DataVector, 3_st, Frame::Inertial>>;
  using dg_interior_dt_vars_tags = tmpl::list<>;
  using dg_interior_deriv_vars_tags = tmpl::list<>;
  using dg_gridless_tags = tmpl::list<>;

  std::optional<std::string> dg_ghost(
      const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
          tensor_driver,
      const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,
      const gsl::not_null<Scalar<DataVector>*> psi,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar,

      const gsl::not_null<Scalar<DataVector>*> lapse,
      const gsl::not_null<tnsr::I<DataVector, 3_st>*> shift,
      //   const gsl::not_null<Scalar<DataVector>*> gamma1,
      //   const gsl::not_null<Scalar<DataVector>*> gamma2,
      const gsl::not_null<tnsr::II<DataVector, 3_st, Frame::Inertial>*>
          inverse_spatial_metric,
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, 3_st>& normal_covector,
      const tnsr::I<DataVector, 3_st>& normal_vector,
      //   const Scalar<DataVector>& psi_interior,
      //   const Scalar<DataVector>& pi_interior,
      //   const tnsr::i<DataVector, 3_st>& phi_interior,
      const tnsr::II<DataVector, 3_st, Frame::Inertial>&
          inverse_spatial_metric_interior,
      //   const Scalar<DataVector>& gamma1_interior,
      //   const Scalar<DataVector>& gamma2_interior,
      const Scalar<DataVector>& lapse_interior,
      const tnsr::I<DataVector, 3_st>& shift_interior) const;

 private:
  double amplitude_ = std::numeric_limits<double>::signaling_NaN();
};
}  // namespace fe::ScalarTensorDriver::BoundaryConditions
