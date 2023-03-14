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
#include "Evolution/PassVariables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
// Add scalar time derivatives, stress-energy tensor, etc.
#include "Evolution/Systems/ScalarTensor/StressEnergy.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
//
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// Debug
#include "Parallel/Printf.hpp"

namespace ScalarTensor {
namespace detail {
// Wrap Scalar temporaries in this prefix to avoid data structure
// collisions. See doc in GRMHD.
template <typename Tag>
struct ScalarTempTag : db::SimpleTag, db::PrefixTag {
  using tag = Tag;
  using type = typename Tag::type;
};

} // namespace detail

/*!
 * \brief Compute the RHS terms and flux values for both the Generalized
 * Harmonic formulation of Einstein's equations and the scalar equations.
 *
 * \details The bulk of the computations in this class dispatch to
 * `GeneralizedHarmonic::TimeDerivative` and
 * `CurvedScalar::TimeDerivativeTerms` as a 'product system' -- each
 * independently operating on its own subset of the supplied variable
 * collections.
 * The additional step is taken to compute the trace-reversed stress energy
 * tensor associated with the scalar part of the system and add its contribution
 * to the \f$\partial_t \Pi_{a b}\f$ variable in the Generalized Harmonic
 * system, which is the only explicit coupling required to back-react the effect
 * of the scalar on the spacetime solution.
 */
struct TimeDerivativeTerms {


//   using scalar_dt_tags = db::wrap_tags_in<
//       ::Tags::dt,
//       typename CurvedScalarWave::System<3_st>::variables_tag::tags_list>;
//   using dt_tags = scalar_dt_tags;
//   using scalar_flux_tags = tmpl::transform<
//       typename CurvedScalarWave::System<3_st>::flux_variables,
//       tmpl::bind<::Tags::Flux, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
//                  tmpl::pin<Frame::Inertial>>>;
//   using scalar_temp_tags =
//       typename CurvedScalarWave::TimeDerivative<3_st>::temporary_tags;
//   using scalar_gradient_tags =
//       typename CurvedScalarWave::System<3_st>::gradients_tags;
//   using d_tags = scalar_gradient_tags;
//   using scalar_arg_tags =
//       typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags;
//   using temporary_tags = scalar_temp_tags;
//   using argument_tags = scalar_arg_tags;

//...
  using gh_dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename GeneralizedHarmonic::System<3_st>::variables_tag::tags_list>;
  using scalar_dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list>;
  using dt_tags = tmpl::append<gh_dt_tags, scalar_dt_tags>;
  using scalar_flux_tags = tmpl::transform<
      typename CurvedScalarWave::System<3_st>::flux_variables,
      tmpl::bind<::Tags::Flux, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  using gh_temp_tags =
      typename GeneralizedHarmonic::TimeDerivative<3_st>::temporary_tags;
  using gh_gradient_tags =
      typename GeneralizedHarmonic::System<3_st>::gradients_tags;
  using gh_arg_tags =
      typename GeneralizedHarmonic::TimeDerivative<3_st>::argument_tags;
  using scalar_temp_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::temporary_tags;
  using scalar_gradient_tags =
      typename CurvedScalarWave::System<3_st>::gradients_tags;
  using gradient_tags = tmpl::append<gh_gradient_tags, scalar_gradient_tags>;
  using scalar_arg_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags;
//   using temporary_tags = tmpl::append<gh_temp_tags, scalar_temp_tags>;
  using temporary_tags =
        tmpl::remove_duplicates<tmpl::append<gh_temp_tags,
        scalar_temp_tags>>;
  using argument_tags = tmpl::append<gh_arg_tags, scalar_arg_tags>;
  //...
  static void apply(
      // GH dt variables
      gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, 3_st>*> dt_pi,
      gsl::not_null<tnsr::iaa<DataVector, 3_st>*> dt_phi,
      // Scalar dt variables
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> dt_phi_scalar,

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
      gsl::not_null<tnsr::ii<DataVector, 3_st>*> spatial_metric,
      gsl::not_null<tnsr::II<DataVector, 3_st>*> inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
      gsl::not_null<tnsr::AA<DataVector, 3_st>*> inverse_spacetime_metric,
      gsl::not_null<tnsr::abb<DataVector, 3_st>*> christoffel_first_kind,
      gsl::not_null<tnsr::Abb<DataVector, 3_st>*> christoffel_second_kind,
      gsl::not_null<tnsr::a<DataVector, 3_st>*> trace_christoffel,
      gsl::not_null<tnsr::A<DataVector, 3_st>*> normal_spacetime_vector,
      gsl::not_null<tnsr::a<DataVector, 3_st>*> normal_spacetime_one_form,
      gsl::not_null<tnsr::abb<DataVector, 3_st>*> da_spacetime_metric,
      // Scalar temporal variables
      // These are duplicates
      //   gsl::not_null<Scalar<DataVector>*> result_lapse,
      //   gsl::not_null<tnsr::I<DataVector, 3_st>*> result_shift,
      //   gsl::not_null<tnsr::II<DataVector, 3_st>*>
      //   result_inverse_spatial_metric,
      //
      gsl::not_null<Scalar<DataVector>*> result_gamma1_scalar,
      gsl::not_null<Scalar<DataVector>*> result_gamma2_scalar,

      // GH argument variables
      // GH spatial derivatives
      const tnsr::iaa<DataVector, 3_st>& d_spacetime_metric,
      const tnsr::iaa<DataVector, 3_st>& d_pi,
      const tnsr::ijaa<DataVector, 3_st>& d_phi,
      // scalar spatial derivatives
      const tnsr::i<DataVector, 3_st>& d_psi_scalar,
      const tnsr::i<DataVector, 3_st>& d_pi_scalar,
      const tnsr::ij<DataVector, 3_st>& d_phi_scalar,
      //
      const tnsr::aa<DataVector, 3_st>& spacetime_metric,
      const tnsr::aa<DataVector, 3_st>& pi,
      const tnsr::iaa<DataVector, 3_st>& phi, const Scalar<DataVector>& gamma0,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const GeneralizedHarmonic::gauges::GaugeCondition& gauge_condition,
      const Mesh<3_st>& mesh, double time,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& inertial_coords,
      const InverseJacobian<DataVector, 3_st, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          mesh_velocity,
      // Scalar argument variables
      //   const tnsr::i<DataVector, 3_st>& d_psi_scalar,
      //   const tnsr::i<DataVector, 3_st>& d_pi_scalar,
      //   const tnsr::ij<DataVector, 3_st>& d_phi_scalar,
      const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3_st>& phi_scalar,
      // These appear with the same name as temporals for the other system
      const Scalar<DataVector>& lapse_scalar,
      const tnsr::I<DataVector, 3_st>& shift_scalar,
      //
      const tnsr::i<DataVector, 3_st>& deriv_lapse,
      const tnsr::iJ<DataVector, 3_st>& deriv_shift,
      const tnsr::II<DataVector, 3_st>& upper_spatial_metric,
      const tnsr::I<DataVector, 3_st>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& gamma1_scalar,
      const Scalar<DataVector>& gamma2_scalar) {
    // Call TimeDerivativeTerms for GH
    GeneralizedHarmonic::TimeDerivative<3_st>::apply(
        // Check for duplicates
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
        lapse, shift, spatial_metric, inverse_spatial_metric,
        det_spatial_metric, sqrt_det_spatial_metric, inverse_spacetime_metric,
        christoffel_first_kind, christoffel_second_kind, trace_christoffel,
        normal_spacetime_vector, normal_spacetime_one_form, da_spacetime_metric,

        // GH argument variables
        d_spacetime_metric, d_pi, d_phi, spacetime_metric, pi, phi, gamma0,
        gamma1, gamma2, gauge_condition, mesh, time, inertial_coords,
        inverse_jacobian, mesh_velocity);

    // Call TimeDerivativeTerms for scalar
    CurvedScalarWave::TimeDerivative<3_st>::apply(
        // Check for duplicates
        // Scalar dt variables
        dt_psi_scalar, dt_pi_scalar, dt_phi_scalar,

        // Scalar temporal variables
        // These are duplicates
        // result_lapse, result_shift, result_inverse_spatial_metric,
        // Use those of the other system:
        lapse, shift, inverse_spatial_metric,
        //
        result_gamma1_scalar, result_gamma2_scalar,

        // Scalar argument variables
        d_psi_scalar, d_pi_scalar, d_phi_scalar, pi_scalar, phi_scalar,
        // These appear with the same name as temporals for the other system
        lapse_scalar, shift_scalar,
        //
        deriv_lapse, deriv_shift, upper_spatial_metric,
        trace_spatial_christoffel, trace_extrinsic_curvature, gamma1_scalar,
        gamma2_scalar);

    // trace_reversed_stress_energy();
    // add_stress_energy_term_to_dt_pi();
  }
};
} // namespace ScalarTensor
