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

template <typename GhDtTagList, typename ScalarDtTagList,
          typename ScalarFluxTagList, typename GhTempTagList,
          typename ScalarTempTagList, typename GhGradientTagList,
          // Add this?
          typename ScalarGradientTagList,
          //
          typename GhArgTagList, typename ScalarArgTagList,
          // Should we keep this?
        //   typename ScalarTimeDerivativeArgTagList,
          //
          typename TraceReversedStressResultTagsList,
          typename TraceReversedStressArgumentTagsList>
struct TimeDerivativeTermsImpl;

template <typename... GhDtTags, typename... ScalarDtTags,
          typename... ScalarFluxTags, typename... GhTempTags,
          typename... ScalarTempTags, typename... GhGradientTags,
          // Add this?
          typename... ScalarGradientTags,
          //
          typename... GhArgTags,
          typename... ScalarArgTags,
          // Should we keep this?
          //   typename... ScalarTimeDerivativeArgTags,
          //
          typename... TraceReversedStressResultTags,
          typename... TraceReversedStressArgumentTags>
struct TimeDerivativeTermsImpl<
    tmpl::list<GhDtTags...>, tmpl::list<ScalarDtTags...>,
    tmpl::list<ScalarFluxTags...>, tmpl::list<GhTempTags...>,
    tmpl::list<ScalarTempTags...>, tmpl::list<GhGradientTags...>,
    // Add this?
    tmpl::list<ScalarGradientTags...>,
    //
    tmpl::list<GhArgTags...>, tmpl::list<ScalarArgTags...>,
    // Should we keep this?
    // tmpl::list<ScalarTimeDerivativeArgTags...>,
    //
    tmpl::list<TraceReversedStressResultTags...>,
    tmpl::list<TraceReversedStressArgumentTags...>> {
  template <typename TemporaryTagsList, typename... ExtraTags>
  static void apply(
      const gsl::not_null<Variables<tmpl::list<GhDtTags..., ScalarDtTags...>>*>
          dt_vars_ptr,
      const gsl::not_null<Variables<db::wrap_tags_in<
          ::Tags::Flux,
          // Should this change? Better remove?
          typename CurvedScalarWave::System<3_st>::flux_variables,
          tmpl::size_t<3>, Frame::Inertial>>*>
          fluxes_ptr,
      //
      const gsl::not_null<Variables<TemporaryTagsList>*> temps_ptr,
      // Add gradient tags
      //   const Variables<tmpl::list<GhGradientTags...,
      //   ScalarGradientTags...>>&
      //       d_vars,

      const tnsr::iaa<DataVector, 3>& d_spacetime_metric,
      const tnsr::iaa<DataVector, 3>& d_pi,
      const tnsr::ijaa<DataVector, 3>& d_phi,
      // Add by hand the scalar gradients here
      const tnsr::i<DataVector, 3>& d_psi_scalar,
      const tnsr::i<DataVector, 3>& d_pi_scalar,
      const tnsr::ij<DataVector, 3>& d_phi_scalar,
      //

      const tuples::TaggedTuple<ExtraTags...>& arguments) {
    // Call TimeDerivativeTerms for GH
    GeneralizedHarmonic::TimeDerivative<3_st>::apply(
        get<GhDtTags>(dt_vars_ptr)..., get<GhTempTags>(temps_ptr)...,
        // Add gradients from tags?
        // get<GhGradientTags>(d_vars)...,
        //
        d_spacetime_metric, d_pi, d_phi,
        get<Tags::detail::TemporaryReference<GhArgTags>>(arguments)...);

    // Additional computations needed here?
    // lapse derivative, phi, inv_spatial_metric, shift, deriv shift,
    // deriv spatial metric, extrinsic curvature, ...

    using extra_tags_list = tmpl::list<ExtraTags...>;

    // Call TimeDerivativeTerms for scalar
    // Check arguments in CurvedScalarWave::TimeDerivative
    CurvedScalarWave::TimeDerivative<3_st>::apply(
        // This argument list needs to be corrected
        get<ScalarDtTags>(dt_vars_ptr)...,
        //
        // get<ScalarFluxTags>(fluxes_ptr)...,
        //
        get<ScalarTempTags>(temps_ptr)...,
        // Add gradients from tags
        // get<ScalarGradientTags>(d_vars)...,
        //
        d_psi_scalar, d_pi_scalar, d_phi_scalar,
        // get<tmpl::conditional_t<
        //     tmpl::list_contains_v<extra_tags_list,
        //                           Tags::detail::TemporaryReference<
        //                               ScalarTimeDerivativeArgTags>>,
        //     Tags::detail::TemporaryReference<ScalarTimeDerivativeArgTags>,
        //     ScalarTimeDerivativeArgTags>>(arguments, *temps_ptr)...
        get<Tags::detail::TemporaryReference<ScalarArgTags>>(arguments)...);
    // Coupling terms between the two systems are added here, e.g.,
    // backreaction of the stress-energy tensor on the metric
    // As a first stage we will set the stress-energy tensor to zero

    // Defined in another StressEnergy files
    // trace_reversed_stress_energy(
    //     get<TraceReversedStressResultTags>(temps_ptr)...,
    //     get<tmpl::conditional_t<
    //         tmpl::list_contains_v<extra_tags_list,
    //                               Tags::detail::TemporaryReference<
    //                                   TraceReversedStressArgumentTags>>,
 //         Tags::detail::TemporaryReference<TraceReversedStressArgumentTags>,
    //         TraceReversedStressArgumentTags>>(*temps_ptr, arguments)...);
    trace_reversed_stress_energy(
        get<TraceReversedStressResultTags>(temps_ptr)...);

    // The addition to dt Pi is independent of the specific form of the stress
    // tensor.
    add_stress_energy_term_to_dt_pi(
        get<::Tags::dt<GeneralizedHarmonic::Tags::Pi<3>>>(dt_vars_ptr),
        get<ScalarTensor::Tags::TraceReversedStressEnergy>(
            *temps_ptr),
        get<gr::Tags::Lapse<DataVector>>(*temps_ptr));
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
struct TimeDerivativeTerms : evolution::PassVariables {

// Add using statements
  using gh_dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename GeneralizedHarmonic::System<3_st>::variables_tag::tags_list>;
  using scalar_dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list>;
  using dt_tags = tmpl::append<gh_dt_tags, scalar_dt_tags>;
  // We do not need this
  using scalar_flux_tags = tmpl::transform<
      typename CurvedScalarWave::System<3_st>::flux_variables,
      tmpl::bind<::Tags::Flux, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  //
  using gh_temp_tags =
      typename GeneralizedHarmonic::TimeDerivative<3_st>::temporary_tags;
  using gh_gradient_tags =
      typename GeneralizedHarmonic::System<3_st>::gradients_tags;
  using gh_arg_tags =
      typename GeneralizedHarmonic::TimeDerivative<3_st>::argument_tags;
  using scalar_temp_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::temporary_tags;
  // Add gradient tags
  using scalar_gradient_tags =
      typename CurvedScalarWave::System<3_st>::gradients_tags;
  using d_tags = tmpl::append<gh_gradient_tags, scalar_gradient_tags>;
  //

  // We do not need the following line. But keep it.
  // Additional temp tags are the derivatives of the metric since GH doesn't
  // explicitly calculate those.
//   using scalar_extra_temp_tags = tmpl::list<
//       // Remove this
//       ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
//                     Frame::Inertial>,
//       ::Tags::deriv<gr::Tags::Shift<3, Frame::Inertial, DataVector>,
//                     tmpl::size_t<3>, Frame::Inertial>,
//       ::Tags::deriv<gr::Tags::SpatialMetric<3, Frame::Inertial, DataVector>,
//                     tmpl::size_t<3>, Frame::Inertial>,
//       gr::Tags::ExtrinsicCurvature<3>
//       //
//       >;
//   using scalar_arg_tags = tmpl::list_difference<
//       typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags,
//       tmpl::append<gh_temp_tags, scalar_extra_temp_tags>>;
  using scalar_arg_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags;
  using trace_reversed_stress_result_tags =
      tmpl::list<Tags::TraceReversedStressEnergy
                 /* Removed FourVelocity, MagneticField */>;
  using trace_reversed_stress_argument_tags = tmpl::list<>;
  /* Add scalar and scalar gradient tags when needed */
  // Remove this
  //   gr::Tags::SpacetimeMetric<3>,
  //   gr::Tags::Shift<3_st, Frame::Inertial, DataVector>,
  //   gr::Tags::Lapse<DataVector>
  //
  //   >;
  //   using temporary_tags = tmpl::remove_duplicates<
  //       tmpl::append<gh_temp_tags, scalar_temp_tags, scalar_extra_temp_tags,
  //                    trace_reversed_stress_result_tags>>;
  using temporary_tags = tmpl::remove_duplicates<
      tmpl::append<gh_temp_tags, scalar_temp_tags,
                   trace_reversed_stress_result_tags>>;
  using argument_tags = tmpl::append<gh_arg_tags, scalar_arg_tags>;

  template <typename... Args>
//   template <typename... GradientVariables, typename... Args>
  static void apply(
      const gsl::not_null<Variables<dt_tags>*> dt_vars_ptr,
      const gsl::not_null<Variables<db::wrap_tags_in<
          ::Tags::Flux,
          // Should this change? Better remove?
          typename CurvedScalarWave::System<3_st>::flux_variables,
          tmpl::size_t<3>, Frame::Inertial>>*>
          fluxes_ptr,
          //
      const gsl::not_null<Variables<temporary_tags>*> temps_ptr,

      // Add gradients from tags
      // Bug, not expanding the pack here
    //   const GradientVariables&... d_vars,
      //
      const tnsr::iaa<DataVector, 3>& d_spacetime_metric,
      const tnsr::iaa<DataVector, 3>& d_pi,
      const tnsr::ijaa<DataVector, 3>& d_phi,

      const tnsr::i<DataVector, 3>& d_psi_scalar,
      const tnsr::i<DataVector, 3>& d_pi_scalar,
      const tnsr::ij<DataVector, 3>& d_phi_scalar,

      const Args&... args
      ) {
    // const tuples::tagged_tuple_from_typelist<
    //     db::wrap_tags_in<Tags::detail::TemporaryReference, argument_tags>>
    //     arguments{args...};
    // const tuples::tagged_tuple_from_typelist<
    //     db::wrap_tags_in<Tags::detail::TemporaryReference, d_tags>>
    //     d_variables{d_vars...};
    detail::TimeDerivativeTermsImpl<
        gh_dt_tags, scalar_dt_tags,
        //
        scalar_flux_tags,
        //
        gh_temp_tags, scalar_temp_tags, gh_gradient_tags,
        // Add scalar gradients from tags
        scalar_gradient_tags,
        //
        gh_arg_tags, scalar_arg_tags,
        // Should we keep this?
        // typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags,
        //
        trace_reversed_stress_result_tags,
        trace_reversed_stress_argument_tags>::apply(dt_vars_ptr, fluxes_ptr,
                                                    temps_ptr,
                                                    // d_vars...,
                                                    d_spacetime_metric, d_pi,
                                                    d_phi,

                                                    d_psi_scalar, d_pi_scalar,
                                                    d_phi_scalar,
                                                    args...);
  }
};
} // namespace ScalarTensor
