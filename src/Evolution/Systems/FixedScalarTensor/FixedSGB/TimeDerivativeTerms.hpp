// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedContainers.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/TimeDerivative.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "Evolution/Systems/ScalarTensor/TimeDerivative.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
// Tag obtained from gh::TimeDerivative needs to be complete here to
// be used in TemporaryReference.
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace fe::sgb {
/*!
 * \brief Compute the RHS terms and flux values for the fixed system.
 *
 * \details See ScalarTensor for a similar implementation.
 */
struct TimeDerivativeTerms /*: public evolution::PassVariables*/ {
  // For now we do not package the arguments and temporaries in a Variables.
  // For this reason, this struct is not base of evolution::PassVariables
  // --see VolumeTermsImpl.tpp for details.
  using gh_dt_tags =
      db::wrap_tags_in<::Tags::dt,
                       typename ScalarTensor::System::variables_tag::tags_list>;
  using scalar_dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename fe::ScalarTensorDriver::System::variables_tag::tags_list>;
  using dt_tags = tmpl::append<gh_dt_tags, scalar_dt_tags>;
  using scalar_flux_tags = tmpl::transform<
      typename fe::ScalarTensorDriver::System::flux_variables,
      tmpl::bind<::Tags::Flux, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  using gh_temp_tags = tmpl::append<
      typename ScalarTensor::TimeDerivative::temporary_tags_no_extra_tags,
      tmpl::list<ScalarTensor::Tags::CompleteTraceReversedStressEnergy>>;
  using gh_gradient_tags = typename ScalarTensor::System::gradients_tags;
  using gh_arg_tags = typename ScalarTensor::TimeDerivative::argument_tags;
  using scalar_temp_tags =
      typename fe::ScalarTensorDriver::TimeDerivative::temporary_tags;
  using scalar_extra_temp_tags = tmpl::list<>;
  using scalar_gradient_tags =
      typename fe::ScalarTensorDriver::System::gradients_tags;
  using gradient_tags = tmpl::append<gh_gradient_tags, scalar_gradient_tags>;
  using scalar_arg_tags =
      typename fe::ScalarTensorDriver::TimeDerivative::argument_tags;
  using temporary_tags = tmpl::remove_duplicates<
      tmpl::append<gh_temp_tags, scalar_temp_tags, scalar_extra_temp_tags>>;
  using argument_tags = tmpl::remove_duplicates<
      tmpl::append<gh_arg_tags, scalar_arg_tags,
                   tmpl::list<ScalarTensor::Tags::TraceReversedStressEnergy<
                       DataVector, 3, ::Frame::Inertial>>>>;

  static void apply(
      // GH dt variables
      gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_pi,
      gsl::not_null<tnsr::iaa<DataVector, 3_st>*> dt_phi,
      // Scalar dt variables
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> dt_phi_scalar,
      // Scalar Tensor Driver dt variables
      gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_tensor_driver,
      gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_pi_tensor_driver,
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar_driver,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar_driver,

      // GH temporal variables
      gsl::not_null<Scalar<DataVector>*> temp_gamma1,
      gsl::not_null<Scalar<DataVector>*> temp_gamma2,
      gsl::not_null<tnsr::a<DataVector, 3_st>*> temp_gauge_function,
      gsl::not_null<tnsr::ab<DataVector, 3_st>*>
          temp_spacetime_deriv_gauge_function,
      gsl::not_null<Scalar<DataVector>*> gamma1gamma2,
      gsl::not_null<Scalar<DataVector>*> half_half_pi_two_normals,
      gsl::not_null<Scalar<DataVector>*> normal_dot_gauge_constraint,
      gsl::not_null<Scalar<DataVector>*> gamma1_plus_1,
      gsl::not_null<tnsr::a<DataVector, 3_st>*> pi_one_normal,
      gsl::not_null<tnsr::a<DataVector, 3_st>*> gauge_constraint,
      gsl::not_null<tnsr::i<DataVector, 3_st>*> half_phi_two_normals,
      gsl::not_null<tnsr::aa<DataVector, 3_st>*>
          shift_dot_three_index_constraint,
      gsl::not_null<tnsr::aa<DataVector, 3_st>*>
          mesh_velocity_dot_three_index_constraint,
      gsl::not_null<tnsr::ia<DataVector, 3_st>*> phi_one_normal,
      gsl::not_null<tnsr::aB<DataVector, 3_st>*> pi_2_up,
      gsl::not_null<tnsr::iaa<DataVector, 3_st>*> three_index_constraint,
      gsl::not_null<tnsr::Iaa<DataVector, 3_st>*> phi_1_up,
      gsl::not_null<tnsr::iaB<DataVector, 3_st>*> phi_3_up,
      gsl::not_null<tnsr::abC<DataVector, 3_st>*> christoffel_first_kind_3_up,
      gsl::not_null<Scalar<DataVector>*> lapse,
      gsl::not_null<tnsr::I<DataVector, 3_st>*> shift,
      gsl::not_null<tnsr::II<DataVector, 3_st>*> inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
      gsl::not_null<tnsr::AA<DataVector, 3_st>*> inverse_spacetime_metric,
      gsl::not_null<tnsr::abb<DataVector, 3_st>*> christoffel_first_kind,
      gsl::not_null<tnsr::Abb<DataVector, 3_st>*> christoffel_second_kind,
      gsl::not_null<tnsr::a<DataVector, 3_st>*> trace_christoffel,
      gsl::not_null<tnsr::A<DataVector, 3_st>*> normal_spacetime_vector,

      // Scalar temporal variables
      gsl::not_null<Scalar<DataVector>*> result_gamma1_scalar,
      gsl::not_null<Scalar<DataVector>*> result_gamma2_scalar,

      // Extra temporal tags
      // Avoid compute tags for deriv of lapse and shift by adding them here
      gsl::not_null<tnsr::aa<DataVector, 3_st>*>
          complete_trace_reversed_stress_energy,

      // Scalar Tensor Driver temporal variables

      // GH argument variables
      // GH spatial derivatives
      const tnsr::iaa<DataVector, 3_st>& d_spacetime_metric,
      const tnsr::iaa<DataVector, 3_st>& d_pi,
      const tnsr::ijaa<DataVector, 3_st>& d_phi,
      // Scalar spatial derivatives
      const tnsr::i<DataVector, 3_st>& d_psi_scalar,
      const tnsr::i<DataVector, 3_st>& d_pi_scalar,
      const tnsr::ij<DataVector, 3_st>& d_phi_scalar,
      // Scalar Tensor Driver spatial derivatives
      const tnsr::iaa<DataVector, 3_st>& d_tensor_driver,
      const tnsr::iaa<DataVector, 3_st>& d_pi_tensor_driver,
      const tnsr::i<DataVector, 3_st>& d_psi_scalar_driver,
      const tnsr::i<DataVector, 3_st>& d_pi_scalar_driver,

      const tnsr::aa<DataVector, 3_st>& spacetime_metric,
      const tnsr::aa<DataVector, 3_st>& pi,
      const tnsr::iaa<DataVector, 3_st>& phi, const Scalar<DataVector>& gamma0,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const gh::gauges::GaugeCondition& gauge_condition, const Mesh<3_st>& mesh,
      double time,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& inertial_coords,
      const InverseJacobian<DataVector, 3_st, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          mesh_velocity,
      // Scalar argument variables
      const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3_st>& phi_scalar,
      // These appear with the same name as temporals for the other system
      const Scalar<DataVector>& lapse_scalar,
      const tnsr::I<DataVector, 3_st>& shift_scalar,

      const tnsr::i<DataVector, 3_st>& deriv_lapse,
      const tnsr::iJ<DataVector, 3_st>& deriv_shift,
      const tnsr::II<DataVector, 3_st>& upper_spatial_metric,
      const tnsr::I<DataVector, 3_st>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& gamma1_scalar,
      const Scalar<DataVector>& gamma2_scalar,
      const Scalar<DataVector>& /* scalar_source*/,

      // Scalar Tensor Driver argument variables
      const tnsr::aa<DataVector, 3_st>& tensor_driver,
      const tnsr::aa<DataVector, 3_st>& pi_tensor_driver,

      const Scalar<DataVector>& psi_scalar_driver,
      const Scalar<DataVector>& pi_scalar_driver,

      // Note: we omit the repeated ::gr argument variables for the driver
      // system

      const tnsr::aa<DataVector, 3_st>& tensor_driver_source,
      const Scalar<DataVector>& scalar_driver_source,
      const Scalar<DataVector>& tau_parameter,
      const Scalar<DataVector>& sigma_parameter,
      const tnsr::aa<DataVector, 3_st>&
          canonical_trace_reversed_stress_energy) {
    // Note: Check that CurvedScalarWave does not update GH variables
    // to a different value. If it does, invert the order of application of the
    // corrections first, so that the GH update is applied at last

    // Compute sourceless part of the RHS of the metric equations
    gh::TimeDerivative<3>::apply(
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
        sqrt_det_spatial_metric, inverse_spacetime_metric,
        christoffel_first_kind, christoffel_second_kind, trace_christoffel,
        normal_spacetime_vector,

        // GH argument variables
        d_spacetime_metric, d_pi, d_phi, spacetime_metric, pi, phi, gamma0,
        gamma1, gamma2, gauge_condition, mesh, time, inertial_coords,
        inverse_jacobian, mesh_velocity);

    // Compute sourceless part of the RHS of the scalar equation
    CurvedScalarWave::TimeDerivative<3>::apply(
        // Scalar dt variables
        dt_psi_scalar, dt_pi_scalar, dt_phi_scalar,

        // Scalar temporal variables
        lapse, shift, inverse_spatial_metric,

        result_gamma1_scalar, result_gamma2_scalar,

        // Scalar argument variables
        d_psi_scalar, d_pi_scalar, d_phi_scalar, pi_scalar, phi_scalar,

        lapse_scalar, shift_scalar,

        deriv_lapse, deriv_shift, upper_spatial_metric,
        trace_spatial_christoffel, trace_extrinsic_curvature, gamma1_scalar,
        gamma2_scalar);

    // Compute the (trace-reversed) stress energy tensor here
    // ScalarTensor::trace_reversed_stress_energy(stress_energy, pi_scalar,
    //                                            phi_scalar, lapse_scalar);
    tenex::evaluate<ti::a, ti::b>(
        complete_trace_reversed_stress_energy,
        canonical_trace_reversed_stress_energy(ti::a, ti::b) +
            tensor_driver(ti::a, ti::b));

    ScalarTensor::add_stress_energy_term_to_dt_pi(
        dt_pi, *complete_trace_reversed_stress_energy, lapse_scalar);

    ScalarTensor::add_scalar_source_to_dt_pi_scalar(
        dt_pi_scalar,
        // We source the scalar equation with the driver
        psi_scalar_driver,  // scalar_source,
        lapse_scalar);

    // Call TimeDerivative for Scalar Tensor Driver
    fe::ScalarTensorDriver::TimeDerivative::apply(
        // Scalar Tensor Driver dt variables
        dt_tensor_driver, dt_pi_tensor_driver, dt_psi_scalar_driver,
        dt_pi_scalar_driver,

        // Scalar Tensor Driver gradients
        d_tensor_driver, d_pi_tensor_driver, d_psi_scalar_driver,
        d_pi_scalar_driver,

        // Scalar Tensor Driver argument variables
        tensor_driver, pi_tensor_driver,

        spacetime_metric,

        mesh, time, inertial_coords, inverse_jacobian, mesh_velocity,

        psi_scalar_driver, pi_scalar_driver, lapse_scalar, shift_scalar,

        tensor_driver_source, scalar_driver_source, tau_parameter,
        sigma_parameter);
  }
};

}  // namespace fe::sgb
