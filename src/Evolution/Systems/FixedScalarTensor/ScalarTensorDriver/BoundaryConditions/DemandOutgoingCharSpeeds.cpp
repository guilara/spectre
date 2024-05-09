// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <pup.h>
#include <string>

#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Characteristics.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarTensorDriver::BoundaryConditions {
DemandOutgoingCharSpeeds::DemandOutgoingCharSpeeds(CkMigrateMessage* const msg)
    : BoundaryCondition(msg) {}

std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
DemandOutgoingCharSpeeds::get_clone() const {
  return std::make_unique<DemandOutgoingCharSpeeds>(*this);
}

void DemandOutgoingCharSpeeds::pup(PUP::er& p) { BoundaryCondition::pup(p); }

std::optional<std::string>
DemandOutgoingCharSpeeds::dg_demand_outgoing_char_speeds(
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, 3, Frame::Inertial>&
        outward_directed_normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>&
    /*outward_directed_normal_vector*/,

    const Scalar<DataVector>& gamma_1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift) {
  const auto char_speeds = characteristic_speeds(
      gamma_1, lapse, shift, outward_directed_normal_covector,
      face_mesh_velocity);
  Scalar<DataVector> normal_dot_mesh_velocity;
  if (face_mesh_velocity.has_value()) {
    normal_dot_mesh_velocity = dot_product(outward_directed_normal_covector,
                                           face_mesh_velocity.value());
  }
  double min_speed = std::numeric_limits<double>::signaling_NaN();
  for (size_t i = 0; i < char_speeds.size(); ++i) {
    if (face_mesh_velocity.has_value()) {
      min_speed = min(gsl::at(char_speeds, i) - get(normal_dot_mesh_velocity));
    } else {
      min_speed = min(gsl::at(char_speeds, i));
    }
    if (min_speed < 0.0) {
      return {MakeString{}
              << "DemandOutgoingCharSpeeds boundary condition violated with "
                 "speed index "
              << i << " ingoing: " << min_speed
              << "\n speed: " << gsl::at(char_speeds, i)
              << "\nn_i: " << outward_directed_normal_covector
              << "\n"
                 "See gh::characteristic_speeds for the "
                 "index ordering of characteristic speeds\n"};
    }
  }
  return std::nullopt;
}

// NOLINTNEXTLINE
PUP::able::PUP_ID DemandOutgoingCharSpeeds::my_PUP_ID = 0;

}  // namespace fe::ScalarTensorDriver::BoundaryConditions
