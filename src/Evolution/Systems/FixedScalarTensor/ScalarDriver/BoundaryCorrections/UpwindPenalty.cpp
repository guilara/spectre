// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryCorrections/UpwindPenalty.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Characteristics.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"  // IWYU pragma: keep

namespace fe::ScalarDriver::BoundaryCorrections {
UpwindPenalty::UpwindPenalty(CkMigrateMessage* msg) : BoundaryCorrection(msg) {}

std::unique_ptr<BoundaryCorrection> UpwindPenalty::get_clone() const {
  return std::make_unique<UpwindPenalty>(*this);
}

void UpwindPenalty::pup(PUP::er& p) { BoundaryCorrection::pup(p); }

double UpwindPenalty::dg_package_data(
    gsl::not_null<Scalar<DataVector>*> packaged_v_psi,
    gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> packaged_v_zero,
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
  *packaged_gamma2 = constraint_gamma2;
  *packaged_interface_unit_normal = interface_unit_normal;
  //   CurvedScalarWave::characteristic_fields(
  //       packaged_v_psi, packaged_v_zero, packaged_v_plus, packaged_v_minus,
  //       constraint_gamma2, psi, pi, phi, interface_unit_normal,
  //       interface_unit_normal_vector);
  //   CurvedScalarWave::characteristic_speeds(packaged_char_speeds,
  //                                           constraint_gamma1, lapse, shift,
  //                                           interface_unit_normal);
  characteristic_fields(packaged_v_psi, packaged_v_zero, packaged_v_plus,
                        packaged_v_minus, constraint_gamma2, psi, pi, phi,
                        interface_unit_normal, interface_unit_normal_vector);
  characteristic_speeds(packaged_char_speeds, constraint_gamma1, lapse, shift,
                        interface_unit_normal);
  if (normal_dot_mesh_velocity.has_value()) {
    get<0>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<1>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<2>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<3>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
  }

  return max(max(get<0>(*packaged_char_speeds), get<1>(*packaged_char_speeds),
                 get<2>(*packaged_char_speeds)));
}

void UpwindPenalty::dg_boundary_terms(
    gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,
    gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
    gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>
        phi_boundary_correction,

    const Scalar<DataVector>& v_psi_int,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero_int,
    const Scalar<DataVector>& v_plus_int, const Scalar<DataVector>& v_minus_int,
    const Scalar<DataVector>& gamma2_int,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& interface_unit_normal_int,
    const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int,

    const Scalar<DataVector>& v_psi_ext,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero_ext,
    const Scalar<DataVector>& v_plus_ext, const Scalar<DataVector>& v_minus_ext,
    const Scalar<DataVector>& gamma2_ext,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& interface_unit_normal_ext,
    const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext,
    dg::Formulation dg_formulation) const {
  // The implementations assumes that the external unit normal vector is
  // exactly negative the internal unit normal vector and that the internal
  // gamma2 is equal to the external gamma2. For Gauss quadrature this will not
  // be exactly true due to interpolation error but that should not matter.

  // Psi stays the same for the advection driver
  get(*psi_boundary_correction) =
      -step_function(get<0>(char_speeds_ext)) * get<0>(char_speeds_ext) *
          get(v_psi_ext) -
      step_function(-get<0>(char_speeds_int)) * get<0>(char_speeds_int) *
          get(v_psi_int) -
      get(gamma2_int) * (-step_function(get<2>(char_speeds_ext)) *
                             get<2>(char_speeds_ext) * get(v_plus_ext) -
                         step_function(-get<2>(char_speeds_int)) *
                             get<2>(char_speeds_int) * get(v_plus_int));

  //   auto& temp_1 = get<3 - 1>(*phi_boundary_correction);
  //   temp_1 = -step_function(get<2>(char_speeds_ext)) *
  //   get<2>(char_speeds_ext) *
  //                get(v_plus_ext) -
  //            step_function(-get<3>(char_speeds_int)) *
  //            get<3>(char_speeds_int) *
  //                get(v_minus_int);

  // in 2+ dimensions the calculation is done without any memory allocations
  //   DataVector temp_2{};
  //   if constexpr (3 > 1) {
  //     temp_2.set_data_ref(make_not_null(&get<3 -
  //     2>(*phi_boundary_correction)));
  //   }
  //   temp_2 = step_function(-get<2>(char_speeds_int)) *
  //   get<2>(char_speeds_int) *
  //                get(v_plus_int) +
  //            step_function(get<3>(char_speeds_ext)) * get<3>(char_speeds_ext)
  //            *
  //                get(v_minus_ext);
  get(*pi_boundary_correction) = -step_function(get<2>(char_speeds_ext)) *
                                     get<2>(char_speeds_ext) * get(v_plus_ext) -
                                 step_function(-get<2>(char_speeds_int)) *
                                     get<2>(char_speeds_int) * get(v_plus_int);

  //   temp_1 = -0.5 * (temp_1 + temp_2);
  for (size_t i = 0; i < 3; ++i) {
    phi_boundary_correction->get(i) =
        (-step_function(get<3>(char_speeds_ext)) * get<3>(char_speeds_ext) *
             get(v_minus_ext) -
         step_function(-get<3>(char_speeds_int)) * get<3>(char_speeds_int) *
             get(v_minus_int)) *
            interface_unit_normal_int.get(i) -
        step_function(get<1>(char_speeds_ext)) * get<1>(char_speeds_ext) *
            v_zero_ext.get(i) -
        step_function(-get<1>(char_speeds_int)) * get<1>(char_speeds_int) *
            v_zero_int.get(i);
  }
}

// NOLINTNEXTLINE
PUP::able::PUP_ID UpwindPenalty::my_PUP_ID = 0;

}  // namespace fe::ScalarDriver::BoundaryCorrections
