// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/BoundaryCorrections/UpwindPenalty.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Characteristics.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"  // IWYU pragma: keep

namespace fe::ScalarTensorDriver::BoundaryCorrections {

UpwindPenalty::UpwindPenalty(CkMigrateMessage* msg) : BoundaryCorrection(msg) {}

std::unique_ptr<BoundaryCorrection> UpwindPenalty::get_clone() const {
  return std::make_unique<UpwindPenalty>(*this);
}

void UpwindPenalty::pup(PUP::er& p) { BoundaryCorrection::pup(p); }

double UpwindPenalty::dg_package_data(
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
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>& mesh_velocity,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
  // Compute the char speeds without the mesh movement, then add the mesh
  // movement.

  // TODO: Add this function to Characteristics.hpp
  characteristic_fields(packaged_v_scalar_driver, packaged_v_pi_scalar,
                        packaged_v_tensor_driver, packaged_v_pi,
                        psi, pi_scalar, tensor_driver, pi, normal_covector);

  characteristic_speeds(packaged_char_speeds, lapse, shift, normal_covector,
                        mesh_velocity);
  if (normal_dot_mesh_velocity.has_value()) {
    get<0>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<1>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<2>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<3>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
  }

  return max(max(get<0>(*packaged_char_speeds), get<1>(*packaged_char_speeds),
                 get<2>(*packaged_char_speeds), get<3>(*packaged_char_speeds)));
}

void UpwindPenalty::dg_boundary_terms(
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
    dg::Formulation /*dg_formulation*/) const {
  const size_t num_pts = char_speeds_int[0].size();
  Variables<tmpl::list<::Tags::TempScalar<0>, ::Tags::TempScalar<1>,
                       ::Tags::TempScalar<2>, ::Tags::TempScalar<3>,
                       ::Tags::TempScalar<4>, ::Tags::TempScalar<5>,
                       ::Tags::TempScalar<6>, ::Tags::TempScalar<7>>>
      buffer(num_pts);
  DataVector& weighted_lambda_scalar_driver_int =
      get(get<::Tags::TempScalar<0>>(buffer));
  weighted_lambda_scalar_driver_int = step_function(-char_speeds_int[0]);
  DataVector& weighted_lambda_scalar_driver_ext =
      get(get<::Tags::TempScalar<1>>(buffer));
  weighted_lambda_scalar_driver_ext = -step_function(char_speeds_ext[0]);

  DataVector& weighted_lambda_pi_scalar_int =
      get(get<::Tags::TempScalar<2>>(buffer));
  weighted_lambda_pi_scalar_int = step_function(-char_speeds_int[1]);
  DataVector& weighted_lambda_pi_scalar_ext =
      get(get<::Tags::TempScalar<3>>(buffer));
  weighted_lambda_pi_scalar_ext = -step_function(char_speeds_ext[1]);

  DataVector& weighted_lambda_tensor_driver_int =
      get(get<::Tags::TempScalar<4>>(buffer));
  weighted_lambda_tensor_driver_int = step_function(-char_speeds_int[2]);
  DataVector& weighted_lambda_tensor_driver_ext =
      get(get<::Tags::TempScalar<5>>(buffer));
  weighted_lambda_tensor_driver_ext = -step_function(char_speeds_ext[2]);

  DataVector& weighted_lambda_pi_int = get(get<::Tags::TempScalar<6>>(buffer));
  weighted_lambda_pi_int = step_function(-char_speeds_int[3]);
  DataVector& weighted_lambda_pi_ext = get(get<::Tags::TempScalar<7>>(buffer));
  weighted_lambda_pi_ext = -step_function(char_speeds_ext[3]);

  boundary_correction_scalar_driver->get() =
      weighted_lambda_scalar_driver_ext * char_speed_v_scalar_driver_ext.get() -
      weighted_lambda_scalar_driver_int * char_speed_v_scalar_driver_int.get();

  boundary_correction_pi_scalar->get() =
      weighted_lambda_pi_scalar_ext * char_speed_v_pi_scalar_ext.get() -
      weighted_lambda_pi_scalar_int * char_speed_v_pi_scalar_int.get();

  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      boundary_correction_tensor_driver->get(a, b) =
          weighted_lambda_tensor_driver_ext *
              char_speed_v_tensor_driver_ext.get(a, b) -
          weighted_lambda_tensor_driver_int *
              char_speed_v_tensor_driver_int.get(a, b);
      boundary_correction_pi->get(a, b) =
          weighted_lambda_pi_ext * char_speed_v_pi_ext.get(a, b) -
          weighted_lambda_pi_int * char_speed_v_pi_int.get(a, b);
    }
  }
}

bool operator==(const UpwindPenalty& /*lhs*/, const UpwindPenalty& /*rhs*/) {
  return true;
}

bool operator!=(const UpwindPenalty& lhs, const UpwindPenalty& rhs) {
  return not(lhs == rhs);
}

// NOLINTNEXTLINE
PUP::able::PUP_ID UpwindPenalty::my_PUP_ID = 0;

}  // namespace fe::ScalarTensorDriver::BoundaryCorrections
