// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>
#include <random>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "Evolution/Systems/ScalarTensor/TimeDerivativeTerms.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/DerivSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace {
template <typename ComputeVolumeTimeDerivativeTerms, size_t Dim,
          typename EvolvedTagList, typename FluxTagList, typename TempTagList,
          typename GradientTagList, typename ArgTagList>
struct ComputeVolumeTimeDerivativeTermsHelper;

template <typename ComputeVolumeTimeDerivativeTerms, size_t Dim,
          typename... EvolvedTags, typename... FluxTags, typename... TempTags,
          typename... GradientTags, typename... ArgTags>
struct ComputeVolumeTimeDerivativeTermsHelper<
    ComputeVolumeTimeDerivativeTerms, Dim, tmpl::list<EvolvedTags...>,
    tmpl::list<FluxTags...>, tmpl::list<TempTags...>,
    tmpl::list<GradientTags...>, tmpl::list<ArgTags...>> {
  template <typename EvolvedVariables, typename FluxVariables,
            typename TemporaryVariables, typename GradientVariables,
            typename ArgumentVariables>
  static void apply(
      const gsl::not_null<EvolvedVariables*> dt_vars_ptr,
      [[maybe_unused]] const gsl::not_null<FluxVariables*> volume_fluxes,
      const gsl::not_null<TemporaryVariables*> temporaries,
      const GradientVariables& partial_derivs,
      const ArgumentVariables& time_derivative_args) {
    ComputeVolumeTimeDerivativeTerms::apply(
        make_not_null(&get<::Tags::dt<EvolvedTags>>(*dt_vars_ptr))...,
        make_not_null(&get<FluxTags>(*volume_fluxes))...,
        make_not_null(&get<TempTags>(*temporaries))...,
        get<GradientTags>(partial_derivs)..., [](const auto& t) -> const auto& {
          if constexpr (tt::is_a_v<std::unique_ptr,
                                   std::decay_t<decltype(t)>>) {
            return *t;
          } else {
            return t;
          }
        }(tuples::get<ArgTags>(time_derivative_args))...);
  }

  template <typename EvolvedVariables, typename FluxVariables,
            typename TemporaryVariables, typename GradientVariables,
            typename ArgumentVariables>
  static void apply_packed(
      const gsl::not_null<EvolvedVariables*> dt_vars_ptr,
      [[maybe_unused]] const gsl::not_null<FluxVariables*> volume_fluxes,
      const gsl::not_null<TemporaryVariables*> temporaries,
      const GradientVariables& partial_derivs,
      const ArgumentVariables& time_derivative_args) {
    ComputeVolumeTimeDerivativeTerms::apply(
        dt_vars_ptr, volume_fluxes, temporaries,
        get<GradientTags>(partial_derivs)..., [](const auto& t) -> const auto& {
          if constexpr (tt::is_a_v<std::unique_ptr,
                                   std::decay_t<decltype(t)>>) {
            return *t;
          } else {
            return t;
          }
        }(tuples::get<ArgTags>(time_derivative_args))...);
  }
};

} // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.ScalarTensor.TimeDerivativeTerms",
    "[Unit][Evolution]") {
  using scalar_variables_tags =
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list;
  using scalar_dt_variables_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_dt_tags;
  using dt_variables_type = Variables<scalar_dt_variables_tags>;
  using scalar_flux_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_flux_tags;
  using flux_variables_type = Variables<scalar_flux_tags>;
  using scalar_temp_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_temp_tags;
  using temp_variables_type = Variables<
      typename ScalarTensor::TimeDerivativeTerms::temporary_tags>;
  using scalar_gradient_tags = tmpl::transform<
      ScalarTensor::TimeDerivativeTerms::scalar_gradient_tags,
      tmpl::bind<::Tags::deriv, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  using gradient_variables_type = Variables<scalar_gradient_tags>;
  using scalar_arg_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_arg_tags;
  using all_scalar_arg_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags;
  using arg_variables_type =
      tuples::tagged_tuple_from_typelist<scalar_arg_tags>;

  const size_t element_size = 10_st;
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<> dist(0.1, 1.0);

  dt_variables_type expected_dt_variables{element_size};
  dt_variables_type dt_variables{element_size};

  flux_variables_type expected_flux_variables{element_size};
  flux_variables_type flux_variables{element_size};

  temp_variables_type temp_variables{element_size};
  temp_variables_type expected_temp_variables{element_size};

  const auto gradient_variables =
      make_with_random_values<gradient_variables_type>(
          make_not_null(&gen), make_not_null(&dist), element_size);
  arg_variables_type arg_variables;
  tmpl::for_each<scalar_arg_tags>([&gen, &dist, &arg_variables](auto tag_v) {
    using tag = typename decltype(tag_v)::type;
    if constexpr (std::is_same_v<
                      typename tag::type,
                      std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>>) {
      tuples::get<tag>(arg_variables) = make_with_random_values<
          typename tnsr::I<DataVector, 3, Frame::Inertial>>(
          make_not_null(&gen), make_not_null(&dist), DataVector{element_size});
    } else if constexpr (tt::is_a_v<Tensor, typename tag::type>) {
      tuples::get<tag>(arg_variables) =
          make_with_random_values<typename tag::type>(make_not_null(&gen),
                                                      make_not_null(&dist),
                                                      DataVector{element_size});
    }
  });

  // The logic of the test is the following
  // We compute use the individual time derivative functions for each system
  // and then compare the results with the time derivative function for the
  // combined system

  // The time derivative function for GeneralizedHarmonic is
  // ComputeVolumeTimeDerivativeTermsHelper<
  //     GeneralizedHarmonic::TimeDerivative<3_st>, 3_st, gh_variables_tags,
  //     gh_flux_tags, gh_temp_tags, gh_gradient_tags,
  //     gh_arg_tags>::apply(make_not_null(&expected_dt_variables),
  //                         make_not_null(&expected_flux_variables),
  //                         make_not_null(&expected_temp_variables),
  //                         gradient_variables, arg_variables);

// The time derivative function for CurvedScalarWave is
  ComputeVolumeTimeDerivativeTermsHelper<
      CurvedScalarWave::TimeDerivative<3_st>, 3_st,
      scalar_variables_tags, scalar_flux_tags, scalar_temp_tags,
      scalar_gradient_tags,
      scalar_arg_tags>::apply(make_not_null(&expected_dt_variables),
                                    make_not_null(&expected_flux_variables),
                                    make_not_null(&expected_temp_variables),
                                    gradient_variables, arg_variables);

// The time derivative function for the combined system is
  ComputeVolumeTimeDerivativeTermsHelper<
      ScalarTensor::TimeDerivativeTerms, 3_st, scalar_variables_tags,
      scalar_flux_tags, scalar_temp_tags, scalar_gradient_tags,
      scalar_arg_tags>::apply(make_not_null(&dt_variables),
                              make_not_null(&flux_variables),
                              make_not_null(&temp_variables),
                              gradient_variables, arg_variables);

// When we have backreaction we also need to compute and apply the correction
// to dt pi for the expected variables
//   ScalarTensor::trace_reversed_stress_energy(...);
//   ScalarTensor::add_stress_energy_term_to_dt_pi(...);

  // Finally we compare
  CHECK_VARIABLES_APPROX(dt_variables, expected_dt_variables);
  CHECK_VARIABLES_APPROX(flux_variables, expected_flux_variables);
  CHECK_VARIABLES_APPROX(temp_variables, expected_temp_variables);
}
