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

template <size_t Dim>
struct TimeDerivativeTermsImpl;

template <size_t Dim>
struct TimeDerivativeTermsImpl {
  static void apply(
      const gsl::not_null<Variables<tmpl::list<GhDtTags..., ScalarDtTags...>>*>
          dt_vars_ptr,
      const gsl::not_null<Variables<db::wrap_tags_in<
          ::Tags::Flux,
          // Should this change? Better remove?
          tmpl::size_t<3>, Frame::Inertial>>*>
          fluxes_ptr,
      //
      const gsl::not_null<Variables<TemporaryTagsList>*> temps_ptr,
      // Add gradient tags
      //   const Variables<tmpl::list<GhGradientTags...,
      //   ScalarGradientTags...>>&
      //       d_vars,

      const tuples::TaggedTuple<ExtraTags...>& arguments) {
    // Call TimeDerivativeTerms for GH
    // GeneralizedHarmonic::TimeDerivative<3_st>::apply(
    //     get<GhDtTags>(dt_vars_ptr)..., get<GhTempTags>(temps_ptr)...,
    //     d_spacetime_metric, d_pi, d_phi,
    //     get<Tags::detail::TemporaryReference<GhArgTags>>(arguments)...);

    // Call TimeDerivativeTerms for scalar
    CurvedScalarWave::TimeDerivative<3_st>::apply(
        dt_psi, dt_pi, dt_phi,

        result_lapse, result_shift, result_inverse_spatial_metric,
        result_gamma1, result_gamma2,

        d_psi, d_pi, d_phi, pi, phi, lapse, shift, deriv_lapse, deriv_shift,
        upper_spatial_metric, trace_spatial_christoffel,
        trace_extrinsic_curvature, gamma1, gamma2);

    // trace_reversed_stress_energy();
    // add_stress_energy_term_to_dt_pi();

  }
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
  using scalar_dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list>;
  using dt_tags = scalar_dt_tags;
  using scalar_flux_tags = tmpl::transform<
      typename CurvedScalarWave::System<3_st>::flux_variables,
      tmpl::bind<::Tags::Flux, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  using scalar_temp_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::temporary_tags;
  using scalar_gradient_tags =
      typename CurvedScalarWave::System<3_st>::gradients_tags;
  using d_tags = scalar_gradient_tags;
  using scalar_arg_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags;
  using temporary_tags = scalar_temp_tags;
  using argument_tags = scalar_arg_tags;

  static void apply(
      gsl::not_null<Scalar<DataVector>*> dt_psi,
      gsl::not_null<Scalar<DataVector>*> dt_pi,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> dt_phi,

      gsl::not_null<Scalar<DataVector>*> result_lapse,
      gsl::not_null<tnsr::I<DataVector, 3_st>*> result_shift,
      gsl::not_null<tnsr::II<DataVector, 3_st>*> result_inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> result_gamma1,
      gsl::not_null<Scalar<DataVector>*> result_gamma2,

      const tnsr::i<DataVector, 3_st>& d_psi,
      const tnsr::i<DataVector, 3_st>& d_pi,
      const tnsr::ij<DataVector, 3_st>& d_phi, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, 3_st>& phi, const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st>& shift,
      const tnsr::i<DataVector, 3_st>& deriv_lapse,
      const tnsr::iJ<DataVector, 3_st>& deriv_shift,
      const tnsr::II<DataVector, 3_st>& upper_spatial_metric,
      const tnsr::I<DataVector, 3_st>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2) {
    detail::TimeDerivativeTermsImpl<3_st>::apply(
        dt_psi, dt_pi, dt_phi,

        result_lapse, result_shift, result_inverse_spatial_metric,
        result_gamma1, result_gamma2,

        d_psi, d_pi, d_phi, pi, phi, lapse, shift, deriv_lapse, deriv_shift,
        upper_spatial_metric, trace_spatial_christoffel,
        trace_extrinsic_curvature, gamma1, gamma2);
  }
};
} // namespace ScalarTensor
