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
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Characteristics.hpp"
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
/// A `BoundaryCondition` that only verifies that all characteristic speeds are
/// directed out of the domain; no boundary data is altered by this boundary
/// condition.
class DemandOutgoingCharSpeeds final : public BoundaryCondition {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Boundary conditions which check that all characteristic "
      "fields are outflowing."};
  DemandOutgoingCharSpeeds() = default;
  /// \cond
  DemandOutgoingCharSpeeds(DemandOutgoingCharSpeeds&&) = default;
  DemandOutgoingCharSpeeds& operator=(DemandOutgoingCharSpeeds&&) = default;
  DemandOutgoingCharSpeeds(const DemandOutgoingCharSpeeds&) = default;
  DemandOutgoingCharSpeeds& operator=(const DemandOutgoingCharSpeeds&) =
      default;
  /// \endcond
  ~DemandOutgoingCharSpeeds() override = default;

  explicit DemandOutgoingCharSpeeds(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, DemandOutgoingCharSpeeds);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags =
      tmpl::list<Tags::ConstraintGamma1, gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, 3_st>>;
  using dg_interior_dt_vars_tags = tmpl::list<>;
  using dg_interior_deriv_vars_tags = tmpl::list<>;
  using dg_gridless_tags = tmpl::list<>;

  std::optional<std::string> dg_demand_outgoing_char_speeds(
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, 3_st>& normal_covector,
      const tnsr::I<DataVector, 3_st>& normal_vector,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st>& shift) const {
    // Use the boundary condition from CurvedScalarWave
    // auto fe_string = csw_boundary_instance_.dg_demand_outgoing_char_speeds(
    //     face_mesh_velocity, normal_covector, normal_vector, gamma1, lapse,
    //     shift);
    // if (not fe_string.has_value()) {
    //   return fe_string;
    // }
    // return std::nullopt;

    tnsr::a<DataVector, 3, Frame::Inertial> char_speeds{lapse.size()};

    // These are the modified char speeds for the advection equation
    characteristic_speeds(make_not_null(&char_speeds), gamma1, lapse, shift,
                          normal_covector);

    if (face_mesh_velocity.has_value()) {
      const auto face_speed = dot_product(normal_covector, *face_mesh_velocity);
      for (auto& char_speed : char_speeds) {
        char_speed -= get(face_speed);
      }
    }
    for (size_t i = 0; i < char_speeds.size(); ++i) {
      if (min(char_speeds[i]) < 0.) {
        return MakeString{}
               << "Detected negative characteristic speed at boundary with "
                  "outgoing char speeds boundary conditions specified. The "
                  "speed is "
               << min(char_speeds[i]) << " for index " << i
               << ". To see which characteristic field this corresponds to, "
                  "check the function `characteristic_speeds` in "
                  "Evolution/Systems/CurvedScalarWave/Characteristics.hpp.";
      }
    }
    return std::nullopt;  // LCOV_EXCL_LINE
  }

 private:
  CurvedScalarWave::BoundaryConditions::DemandOutgoingCharSpeeds<3_st>
      csw_boundary_instance_ =
          CurvedScalarWave::BoundaryConditions::DemandOutgoingCharSpeeds<
              3_st>();
};
}  // namespace fe::ScalarDriver::BoundaryConditions
