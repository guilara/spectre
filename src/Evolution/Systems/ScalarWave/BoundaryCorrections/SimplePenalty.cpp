// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/BoundaryCorrections/SimplePenalty.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"  // IWYU pragma: keep

namespace ScalarWave::BoundaryCorrections {
template <size_t Dim>
SimplePenalty<Dim>::SimplePenalty(CkMigrateMessage* msg)
    : BoundaryCorrection<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<BoundaryCorrection<Dim>> SimplePenalty<Dim>::get_clone() const {
  return std::make_unique<SimplePenalty>(*this);
}

template <size_t Dim>
void SimplePenalty<Dim>::pup(PUP::er& p) {
  BoundaryCorrection<Dim>::pup(p);
  p | penalty_factor_;
}

template <size_t Dim>
SimplePenalty<Dim>::SimplePenalty(const double penalty_factor)
    : penalty_factor_(penalty_factor) {}

template <size_t Dim>
double SimplePenalty<Dim>::dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_psi,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_char_speed_v_zero,
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_plus,
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_minus,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_char_speed_n_times_v_plus,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_char_speed_n_times_v_minus,
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_gamma2_v_psi,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        packaged_char_speeds,

    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,

    const Scalar<DataVector>& constraint_gamma2,

    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    /*mesh_velocity*/,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
  if (normal_dot_mesh_velocity.has_value()) {
    get<0>(*packaged_char_speeds) = -get(*normal_dot_mesh_velocity);
    get<1>(*packaged_char_speeds) = 1.0 - get(*normal_dot_mesh_velocity);
    get<2>(*packaged_char_speeds) = -1.0 - get(*normal_dot_mesh_velocity);
  } else {
    get<0>(*packaged_char_speeds) = 0.0;
    get<1>(*packaged_char_speeds) = 1.0;
    get<2>(*packaged_char_speeds) = -1.0;
  }

  // Computes the contribution to the boundary correction from one side of the
  // interface.
  //
  // Note: when SimplePenalty::dg_boundary_terms() is called, an Element passes
  // in its own packaged data to fill the interior fields, and its neighbor's
  // packaged data to fill the exterior fields. This introduces a sign flip for
  // each normal used in computing the exterior fields.
  get(*packaged_char_speed_gamma2_v_psi) = get(constraint_gamma2) * get(psi);
  {
    // Use v_psi allocation as n^i Phi_i
    dot_product(packaged_char_speed_v_psi, normal_covector, phi);
    const auto& normal_dot_phi = get(*packaged_char_speed_v_psi);

    for (size_t i = 0; i < Dim; ++i) {
      packaged_char_speed_v_zero->get(i) =
          get<0>(*packaged_char_speeds) *
          (phi.get(i) - normal_covector.get(i) * normal_dot_phi);
    }

    get(*packaged_char_speed_v_plus) =
        get<1>(*packaged_char_speeds) *
        (get(pi) + normal_dot_phi - get(*packaged_char_speed_gamma2_v_psi));
    get(*packaged_char_speed_v_minus) =
        get<2>(*packaged_char_speeds) *
        (get(pi) - normal_dot_phi - get(*packaged_char_speed_gamma2_v_psi));
  }

  for (size_t d = 0; d < Dim; ++d) {
    packaged_char_speed_n_times_v_plus->get(d) =
        get(*packaged_char_speed_v_plus) * normal_covector.get(d);
    packaged_char_speed_n_times_v_minus->get(d) =
        get(*packaged_char_speed_v_minus) * normal_covector.get(d);
  }

  get(*packaged_char_speed_v_psi) = get<0>(*packaged_char_speeds) * get(psi);
  get(*packaged_char_speed_gamma2_v_psi) *= get<0>(*packaged_char_speeds);

  const auto result =
      max(max(get<0>(*packaged_char_speeds), get<1>(*packaged_char_speeds),
              get<2>(*packaged_char_speeds)));

  // To avoid ambiguities in the penalties package the evolved fields themselves
  get(*packaged_char_speed_v_psi) = get(psi);
  get(*packaged_char_speed_v_plus) = get(pi);
  for (size_t d = 0; d < Dim; ++d) {
    packaged_char_speed_v_zero->get(d) = phi.get(d);
  }

  return result;
}

template <size_t Dim>
void SimplePenalty<Dim>::dg_boundary_terms(
    const gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,
    const gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
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
    dg::Formulation /*dg_formulation*/) const {
  // For now we reconstruct write the evolved fields using the characteristic
  // fields. But the latter are not really needed. If we package the evolved
  // fields, we could simplify the expressions.
  //   psi_boundary_correction->get() =
  //       penalty_factor_ * (get(char_speed_v_psi_ext) -
  //       get(char_speed_v_psi_int));
  psi_boundary_correction->get() = 0.0;

  get(*pi_boundary_correction) = penalty_factor_ * (get(char_speed_v_plus_ext) -
                                                    get(char_speed_v_plus_int));

  for (size_t d = 0; d < Dim; ++d) {
    get(*pi_boundary_correction) +=
        penalty_factor_ *
        (char_speed_v_zero_ext.get(d) - char_speed_v_zero_int.get(d));
  }

  for (size_t d = 0; d < Dim; ++d) {
    phi_boundary_correction->get(d) =
        penalty_factor_ *
        (get(char_speed_v_plus_ext) - get(char_speed_v_plus_int));

    for (size_t i = 0; i < Dim; ++i) {
      phi_boundary_correction->get(d) +=
          penalty_factor_ *
          (char_speed_v_zero_ext.get(i) - char_speed_v_zero_int.get(i));
    }
  }
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID SimplePenalty<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(_, data) template class SimplePenalty<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace ScalarWave::BoundaryCorrections
