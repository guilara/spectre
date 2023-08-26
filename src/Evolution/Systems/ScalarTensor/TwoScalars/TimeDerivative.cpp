// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/TwoScalars/TimeDerivative.hpp"

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedContainers.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/StressEnergy.hpp"
#include "Evolution/Systems/ScalarTensor/TwoScalars/System.hpp"
#include "Evolution/Systems/ScalarTensor/TwoScalars/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarTensor::TwoScalars {
void TimeDerivative::apply(
    // GH dt variables
    gsl::not_null<tnsr::aa<DataVector, dim>*> dt_spacetime_metric,
    gsl::not_null<tnsr::aa<DataVector, dim>*> dt_pi,
    gsl::not_null<tnsr::iaa<DataVector, dim>*> dt_phi,
    // Scalar dt variables
    gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
    gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
    gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*> dt_phi_scalar,
    // Scalar dt variables (copy)
    gsl::not_null<Scalar<DataVector>*> dt_psi_scalar_2,
    gsl::not_null<Scalar<DataVector>*> dt_pi_scalar_2,
    gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*> dt_phi_scalar_2,

    // GH temporal variables
    gsl::not_null<Scalar<DataVector>*> temp_gamma1,
    gsl::not_null<Scalar<DataVector>*> temp_gamma2,
    gsl::not_null<tnsr::a<DataVector, dim>*> temp_gauge_function,
    gsl::not_null<tnsr::ab<DataVector, dim>*>
        temp_spacetime_deriv_gauge_function,
    gsl::not_null<Scalar<DataVector>*> gamma1gamma2,
    gsl::not_null<Scalar<DataVector>*> half_half_pi_two_normals,
    gsl::not_null<Scalar<DataVector>*> normal_dot_gauge_constraint,
    gsl::not_null<Scalar<DataVector>*> gamma1_plus_1,
    gsl::not_null<tnsr::a<DataVector, dim>*> pi_one_normal,
    gsl::not_null<tnsr::a<DataVector, dim>*> gauge_constraint,
    gsl::not_null<tnsr::i<DataVector, dim>*> half_phi_two_normals,
    gsl::not_null<tnsr::aa<DataVector, dim>*> shift_dot_three_index_constraint,
    gsl::not_null<tnsr::aa<DataVector, dim>*>
        mesh_velocity_dot_three_index_constraint,
    gsl::not_null<tnsr::ia<DataVector, dim>*> phi_one_normal,
    gsl::not_null<tnsr::aB<DataVector, dim>*> pi_2_up,
    gsl::not_null<tnsr::iaa<DataVector, dim>*> three_index_constraint,
    gsl::not_null<tnsr::Iaa<DataVector, dim>*> phi_1_up,
    gsl::not_null<tnsr::iaB<DataVector, dim>*> phi_3_up,
    gsl::not_null<tnsr::abC<DataVector, dim>*> christoffel_first_kind_3_up,
    gsl::not_null<Scalar<DataVector>*> lapse,
    gsl::not_null<tnsr::I<DataVector, dim>*> shift,
    gsl::not_null<tnsr::II<DataVector, dim>*> inverse_spatial_metric,
    gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
    gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
    gsl::not_null<tnsr::AA<DataVector, dim>*> inverse_spacetime_metric,
    gsl::not_null<tnsr::abb<DataVector, dim>*> christoffel_first_kind,
    gsl::not_null<tnsr::Abb<DataVector, dim>*> christoffel_second_kind,
    gsl::not_null<tnsr::a<DataVector, dim>*> trace_christoffel,
    gsl::not_null<tnsr::A<DataVector, dim>*> normal_spacetime_vector,

    // Scalar temporal variables
    gsl::not_null<Scalar<DataVector>*> result_gamma1_scalar,
    gsl::not_null<Scalar<DataVector>*> result_gamma2_scalar,
    // Scalar temporal variables (copy)
    gsl::not_null<Scalar<DataVector>*> result_gamma1_scalar_2,
    gsl::not_null<Scalar<DataVector>*> result_gamma2_scalar_2,

    // Extra temporal tags
    gsl::not_null<tnsr::aa<DataVector, dim>*> stress_energy,
    gsl::not_null<tnsr::aa<DataVector, dim>*> stress_energy_2,

    // GH spatial derivatives
    const tnsr::iaa<DataVector, dim>& d_spacetime_metric,
    const tnsr::iaa<DataVector, dim>& d_pi,
    const tnsr::ijaa<DataVector, dim>& d_phi,

    // scalar spatial derivatives
    const tnsr::i<DataVector, dim>& d_psi_scalar,
    const tnsr::i<DataVector, dim>& d_pi_scalar,
    const tnsr::ij<DataVector, dim>& d_phi_scalar,
    // scalar spatial derivatives
    const tnsr::i<DataVector, dim>& d_psi_scalar_2,
    const tnsr::i<DataVector, dim>& d_pi_scalar_2,
    const tnsr::ij<DataVector, dim>& d_phi_scalar_2,

    // GH argument variables
    const tnsr::aa<DataVector, dim>& spacetime_metric,
    const tnsr::aa<DataVector, dim>& pi, const tnsr::iaa<DataVector, dim>& phi,
    const Scalar<DataVector>& gamma0, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2,
    const gh::gauges::GaugeCondition& gauge_condition, const Mesh<dim>& mesh,
    double time,
    const tnsr::I<DataVector, dim, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const std::optional<tnsr::I<DataVector, dim, Frame::Inertial>>&
        mesh_velocity,

    // Scalar argument variables
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, dim>& phi_scalar,
    const Scalar<DataVector>& lapse_scalar,
    const tnsr::I<DataVector, dim>& shift_scalar,
    const tnsr::i<DataVector, dim>& deriv_lapse,
    const tnsr::iJ<DataVector, dim>& deriv_shift,
    const tnsr::II<DataVector, dim>& upper_spatial_metric,
    const tnsr::I<DataVector, dim>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1_scalar,
    const Scalar<DataVector>& gamma2_scalar,
    // Scalar argument variables (second scalar)
    const Scalar<DataVector>& pi_scalar_2,
    const tnsr::i<DataVector, dim>& phi_scalar_2,
    const Scalar<DataVector>& gamma1_scalar_2,
    const Scalar<DataVector>& gamma2_scalar_2,

    // Scalar sources
    const Scalar<DataVector>& scalar_source,
    const Scalar<DataVector>& scalar_source_2) {
  // Compute sourceless part of the RHS of the metric equations
  gh::TimeDerivative<dim>::apply(
      // GH dt variables
      dt_spacetime_metric, dt_pi, dt_phi,

      // GH temporal variables
      temp_gamma1, temp_gamma2, temp_gauge_function,
      temp_spacetime_deriv_gauge_function, gamma1gamma2,
      half_half_pi_two_normals, normal_dot_gauge_constraint, gamma1_plus_1,
      pi_one_normal, gauge_constraint, half_phi_two_normals,
      shift_dot_three_index_constraint,
      mesh_velocity_dot_three_index_constraint, phi_one_normal, pi_2_up,
      three_index_constraint, phi_1_up, phi_3_up, christoffel_first_kind_3_up,
      lapse, shift, inverse_spatial_metric, det_spatial_metric,
      sqrt_det_spatial_metric, inverse_spacetime_metric, christoffel_first_kind,
      christoffel_second_kind, trace_christoffel, normal_spacetime_vector,

      // GH argument variables
      d_spacetime_metric, d_pi, d_phi, spacetime_metric, pi, phi, gamma0,
      gamma1, gamma2, gauge_condition, mesh, time, inertial_coords,
      inverse_jacobian, mesh_velocity);

  // Compute sourceless part of the RHS of the scalar equation
  CurvedScalarWave::TimeDerivative<dim>::apply(
      // Scalar dt variables
      dt_psi_scalar, dt_pi_scalar, dt_phi_scalar,

      // Scalar temporal variables
      lapse, shift, inverse_spatial_metric,

      result_gamma1_scalar, result_gamma2_scalar,

      // Scalar argument variables
      d_psi_scalar, d_pi_scalar, d_phi_scalar, pi_scalar, phi_scalar,

      lapse_scalar, shift_scalar,

      deriv_lapse, deriv_shift, upper_spatial_metric, trace_spatial_christoffel,
      trace_extrinsic_curvature, gamma1_scalar, gamma2_scalar);

  // Compute sourceless part of the RHS of the scalar equation
  CurvedScalarWave::TimeDerivative<dim>::apply(
      // Scalar dt variables
      dt_psi_scalar_2, dt_pi_scalar_2, dt_phi_scalar_2,

      // Scalar temporal variables
      lapse, shift, inverse_spatial_metric,

      result_gamma1_scalar_2, result_gamma2_scalar_2,

      // Scalar argument variables
      d_psi_scalar_2, d_pi_scalar_2, d_phi_scalar_2, pi_scalar_2, phi_scalar_2,

      lapse_scalar, shift_scalar,

      deriv_lapse, deriv_shift, upper_spatial_metric, trace_spatial_christoffel,
      trace_extrinsic_curvature, gamma1_scalar_2, gamma2_scalar_2);

  // Compute the (trace-reversed) stress energy tensor here
  ScalarTensor::trace_reversed_stress_energy(stress_energy, pi_scalar,
                                             phi_scalar, lapse_scalar);

  ScalarTensor::trace_reversed_stress_energy(stress_energy_2, pi_scalar_2,
                                             phi_scalar_2, lapse_scalar);

  // Add stress tensors to the gh RHS
  ScalarTensor::add_stress_energy_term_to_dt_pi(dt_pi, *stress_energy,
                                                lapse_scalar);

  ScalarTensor::add_stress_energy_term_to_dt_pi(dt_pi, *stress_energy_2,
                                                lapse_scalar);

  // Add scalar sources to the scalar RHS
  ScalarTensor::add_scalar_source_to_dt_pi_scalar(dt_pi_scalar, scalar_source,
                                                  lapse_scalar);

  ScalarTensor::add_scalar_source_to_dt_pi_scalar(
      dt_pi_scalar_2, scalar_source_2, lapse_scalar);
}
}  // namespace ScalarTensor::TwoScalars
