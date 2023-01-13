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
  using gh_variables_tags =
      typename GeneralizedHarmonic::System<3_st>::variables_tag::tags_list;
  using scalar_variables_tags =
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list;

  using gh_dt_variables_tags =
      ScalarTensor::TimeDerivativeTerms::gh_dt_tags;
  using scalar_dt_variables_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_dt_tags;
  using dt_variables_type =
      Variables<tmpl::append<gh_dt_variables_tags, scalar_dt_variables_tags>>;

  using gh_flux_tags = tmpl::list<>;
  using scalar_flux_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_flux_tags;
  using flux_variables_type =
      Variables<tmpl::append<gh_flux_tags, scalar_flux_tags>>;

  using gh_temp_tags =
      ScalarTensor::TimeDerivativeTerms::gh_temp_tags;
  using scalar_temp_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_temp_tags;
  using temp_variables_type = Variables<
      typename ScalarTensor::TimeDerivativeTerms::temporary_tags>;

  using gh_gradient_tags = tmpl::transform<
      ScalarTensor::TimeDerivativeTerms::gh_gradient_tags,
      tmpl::bind<::Tags::deriv, tmpl::_1, tmpl::pin<tmpl::size_t<3_st>>,
                 tmpl::pin<Frame::Inertial>>>;
  using scalar_gradient_tags = tmpl::list<>;
  using gradient_variables_type = Variables<gh_gradient_tags>;

  using gh_arg_tags =
      ScalarTensor::TimeDerivativeTerms::gh_arg_tags;
  using scalar_arg_tags =
      ScalarTensor::TimeDerivativeTerms::scalar_arg_tags;
  using all_scalar_arg_tags =
      typename CurvedScalarWave::TimeDerivative<3_st>::argument_tags;
  using arg_variables_type = tuples::tagged_tuple_from_typelist<
      tmpl::append<gh_arg_tags, scalar_arg_tags>>;
}
