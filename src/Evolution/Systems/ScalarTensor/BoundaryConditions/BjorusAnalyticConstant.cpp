// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/BoundaryConditions/BjorusAnalyticConstant.hpp"

#include <cstddef>
#include <memory>
#include <pup.h>

namespace ScalarTensor::BoundaryConditions {
BjorusAnalyticConstant::BjorusAnalyticConstant(
    gh::BoundaryConditions::detail::ConstraintPreservingBjorhusType type,
    const double amplitude_scalar)
    : constraint_preserving_(type), amplitude_scalar_(amplitude_scalar) {}

// LCOV_EXCL_START
BjorusAnalyticConstant::BjorusAnalyticConstant(CkMigrateMessage* const msg)
    : BoundaryCondition(msg) {}
// LCOV_EXCL_STOP

std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
BjorusAnalyticConstant::get_clone() const {
  return std::make_unique<BjorusAnalyticConstant>(*this);
}

void BjorusAnalyticConstant::pup(PUP::er& p) {
  BoundaryCondition::pup(p);
  p | constraint_preserving_;
  p | amplitude_scalar_;
}

std::optional<std::string> BjorusAnalyticConstant::dg_ghost(
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
        spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*> phi,

    const gsl::not_null<Scalar<DataVector>*> psi_scalar,
    const gsl::not_null<Scalar<DataVector>*> pi_scalar,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*> phi_scalar,

    const gsl::not_null<Scalar<DataVector>*> gamma1,
    const gsl::not_null<Scalar<DataVector>*> gamma2,
    const gsl::not_null<Scalar<DataVector>*> lapse,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> shift,
    const gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*>
        inv_spatial_metric,

    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,

    const tnsr::aa<DataVector, 3, Frame::Inertial>& interior_spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& interior_pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& interior_phi,

    const Scalar<DataVector>& psi_scalar_interior,
    const Scalar<DataVector>& pi_scalar_interior,
    const tnsr::i<DataVector, 3>& phi_scalar_interior,

    const tnsr::I<DataVector, 3, Frame::Inertial>& /*coords*/,
    const Scalar<DataVector>& interior_gamma1,
    const Scalar<DataVector>& interior_gamma2,
    const Scalar<DataVector>& interior_lapse,
    const tnsr::I<DataVector, 3>& interior_shift,
    const tnsr::II<DataVector, 3>& interior_inv_spatial_metric,
    const tnsr::AA<DataVector, 3,
                   Frame::Inertial>& /*inverse_spacetime_metric*/,
    const tnsr::A<DataVector, 3, Frame::Inertial>&
    /*spacetime_unit_normal_vector*/,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& /*three_index_constraint*/,
    const tnsr::a<DataVector, 3, Frame::Inertial>& /*gauge_source*/,
    const tnsr::ab<DataVector, 3, Frame::Inertial>&
    /*spacetime_deriv_gauge_source*/,

    // const Scalar<DataVector>& gamma1_scalar_interior,
    // const Scalar<DataVector>& gamma2_scalar_interior,

    // c.f. dg_interior_dt_vars_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>&
    /*logical_dt_spacetime_metric*/,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& /*logical_dt_pi*/,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& /*logical_dt_phi*/,

    const Scalar<DataVector>& logical_dt_psi_scalar,
    const Scalar<DataVector>& logical_dt_pi_scalar,
    const tnsr::i<DataVector, 3, Frame::Inertial>& logical_dt_phi_scalar,

    // c.f. dg_interior_deriv_vars_tags
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& /*d_spacetime_metric*/,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& /*d_pi*/,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& /*d_phi*/,

    const tnsr::i<DataVector, 3, Frame::Inertial>& d_psi_scalar,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_pi_scalar,
    const tnsr::ij<DataVector, 3, Frame::Inertial>& d_phi_scalar) const {
  *gamma1 = interior_gamma1;
  *gamma2 = interior_gamma2;
  *spacetime_metric = interior_spacetime_metric;
  *pi = interior_pi;
  *phi = interior_phi;
  *lapse = interior_lapse;
  *shift = interior_shift;
  *inv_spatial_metric = interior_inv_spatial_metric;

  // Implement as in CurvedScalarWave::AnalyticConstant
  *psi_scalar =
      make_with_value<Scalar<DataVector>>(interior_gamma1, amplitude_scalar_);
  *pi_scalar = make_with_value<Scalar<DataVector>>(interior_gamma1, 0.0);
  *phi_scalar = make_with_value<tnsr::i<DataVector, 3, Frame::Inertial>>(
      interior_gamma1, 0.0);

  return {};
}

std::optional<std::string> BjorusAnalyticConstant::dg_time_derivative(
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
        dt_spacetime_metric_correction,
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
        dt_pi_correction,
    const gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*>
        dt_phi_correction,

    const gsl::not_null<Scalar<DataVector>*> dt_psi_scalar_correction,
    const gsl::not_null<Scalar<DataVector>*> dt_pi_scalar_correction,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        dt_phi_scalar_correction,

    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
    // c.f. dg_interior_evolved_variables_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,

    const Scalar<DataVector>& psi_scalar, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3, Frame::Inertial>& phi_scalar,
    // c.f. dg_interior_primitive_variables_tags

    // c.f. dg_interior_temporary_tags
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::II<DataVector, 3>& /*interior_inv_spatial_metric*/,
    const tnsr::AA<DataVector, 3, Frame::Inertial>& inverse_spacetime_metric,
    const tnsr::A<DataVector, 3, Frame::Inertial>& spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& three_index_constraint,
    const tnsr::a<DataVector, 3, Frame::Inertial>& gauge_source,
    const tnsr::ab<DataVector, 3, Frame::Inertial>&
        spacetime_deriv_gauge_source,

    // const Scalar<DataVector>& gamma1_scalar_interior,
    // const Scalar<DataVector>& gamma2_scalar_interior,

    // c.f. dg_interior_dt_vars_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>& logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& logical_dt_pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& logical_dt_phi,

    const Scalar<DataVector>& logical_dt_psi_scalar,
    const Scalar<DataVector>& logical_dt_pi_scalar,
    const tnsr::i<DataVector, 3, Frame::Inertial>& logical_dt_phi_scalar,

    // c.f. dg_interior_deriv_vars_tags
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi,

    const tnsr::i<DataVector, 3, Frame::Inertial>& d_psi_scalar,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_pi_scalar,
    const tnsr::ij<DataVector, 3, Frame::Inertial>& d_phi_scalar) const {
  get(*dt_psi_scalar_correction) = 0.0;
  get(*dt_pi_scalar_correction) = 0.0;
  for (size_t i = 0; i < 3; ++i) {
    dt_phi_scalar_correction->get(i) = 0.0;
  }

  return constraint_preserving_.dg_time_derivative(
      dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
      face_mesh_velocity, normal_covector, normal_vector, spacetime_metric, pi,
      phi, coords, gamma1, gamma2, lapse, shift, inverse_spacetime_metric,
      spacetime_unit_normal_vector, three_index_constraint, gauge_source,
      spacetime_deriv_gauge_source, logical_dt_spacetime_metric, logical_dt_pi,
      logical_dt_phi, d_spacetime_metric, d_pi, d_phi);
}

// NOLINTNEXTLINE
PUP::able::PUP_ID BjorusAnalyticConstant::my_PUP_ID = 0;
}  // namespace ScalarTensor::BoundaryConditions
