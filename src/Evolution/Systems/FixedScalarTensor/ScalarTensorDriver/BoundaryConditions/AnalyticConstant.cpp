// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryConditions/AnalyticConstant.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeString.hpp"

namespace fe::ScalarTensorDriver::BoundaryConditions {

std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
AnalyticConstant::get_clone() const {
  return std::make_unique<AnalyticConstant>(*this);
}

void AnalyticConstant::pup(PUP::er& p) {
  BoundaryCondition::pup(p);
  p | amplitude_;
}

AnalyticConstant::AnalyticConstant(CkMigrateMessage* const msg)
    : BoundaryCondition(msg) {}

AnalyticConstant::AnalyticConstant(const double amplitude)
    : amplitude_(amplitude) {}

std::optional<std::string> AnalyticConstant::dg_ghost(
    const gsl::not_null<Scalar<DataVector>*> psi,
    const gsl::not_null<Scalar<DataVector>*> pi_scalar,
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
        tensor_driver,
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,

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
    const tnsr::I<DataVector, 3_st>& shift_interior) const {
  // Use the boundary condition from CurvedScalarWave
  *psi = make_with_value<Scalar<DataVector>>(lapse_interior, amplitude_);
  *pi_scalar = make_with_value<Scalar<DataVector>>(lapse_interior, 0.0);
  *lapse = lapse_interior;
  *shift = shift_interior;
  *inverse_spatial_metric = inverse_spatial_metric_interior;
  // *gamma1 = gamma1_interior;
  // *gamma2 = gamma2_interior;

  *tensor_driver =
      make_with_value<tnsr::aa<DataVector, 3>>(lapse_interior, amplitude_);
  *pi = make_with_value<tnsr::aa<DataVector, 3>>(lapse_interior, amplitude_);

  return {};
}

// NOLINTNEXTLINE
PUP::able::PUP_ID AnalyticConstant::my_PUP_ID = 0;

}  // namespace fe::ScalarTensorDriver::BoundaryConditions
