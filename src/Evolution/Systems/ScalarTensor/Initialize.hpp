// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>  // IWYU pragma: keep
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/Norms.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
//
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
//
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintGammas.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/DerivSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
//
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
//
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalOneForm.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
/// \endcond

namespace ScalarTensor::Actions {
// template <size_t Dim>
struct InitializeScalarTensorAnd3Plus1Variables {
  static constexpr size_t Dim = 3_st;
  using frame = Frame::Inertial;
  using compute_tags = db::AddComputeTags<
      // Needed to compute the characteristic speeds for the AH finder
      gr::Tags::SpatialMetricCompute<DataVector, Dim, frame>,
      gr::Tags::DetAndInverseSpatialMetricCompute<DataVector, Dim, frame>,
      gr::Tags::ShiftCompute<DataVector, Dim, frame>,
      gr::Tags::LapseCompute<DataVector, Dim, frame>,

      gr::Tags::SpacetimeNormalVectorCompute<DataVector, Dim, frame>,
      gh::Tags::DerivLapseCompute<Dim, frame>,

      gr::Tags::InverseSpacetimeMetricCompute<DataVector, Dim, frame>,
      gh::Tags::DerivShiftCompute<Dim, frame>,

      gh::Tags::DerivSpatialMetricCompute<Dim, frame>,

      // Add trace compute tags for Christoffel and Extrinsic curvature

      gr::Tags::SpatialChristoffelFirstKindCompute<DataVector, Dim, frame>,
      gr::Tags::SpatialChristoffelSecondKindCompute<DataVector, Dim, frame>,
      gr::Tags::TraceSpatialChristoffelSecondKindCompute<DataVector, Dim,
                                                         frame>,
      gh::Tags::ExtrinsicCurvatureCompute<Dim, frame>,
      gh::Tags::TraceExtrinsicCurvatureCompute<Dim, frame>,

      // Compute constraint damping parameters.
      gh::ConstraintDamping::Tags::ConstraintGamma0Compute<Dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::ConstraintGamma1Compute<Dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::ConstraintGamma2Compute<Dim, Frame::Grid>,
      //   ScalarTensor::Sources::Tags::ScalarSourceCompute,
      ScalarTensor::Sources::Tags::ScalarCurvatureSourceCompute<Dim, frame,
                                                                DataVector>,

      // Extra tags for curvatures
      ::Tags::DerivTensorCompute<
          gr::Tags::ExtrinsicCurvature<DataVector, Dim, frame>,
          ::domain::Tags::InverseJacobian<Dim, ::Frame::ElementLogical,
                                          ::Frame::Inertial>>,
      gh::Tags::GradExtrinsicCurvatureCompute<Dim, frame>,
      ::Tags::DerivTensorCompute<
          gr::Tags::SpatialChristoffelSecondKind<DataVector, Dim, frame>,
          ::domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                          Frame::Inertial>>,

      gr::Tags::SpatialRicciCompute<DataVector, Dim, frame>,
      gr::Tags::SpatialRicciScalarCompute<DataVector, Dim, frame>,

      gr::Tags::WeylElectricCompute<DataVector, Dim, frame>,
      gr::Tags::WeylElectricScalarCompute<DataVector, Dim, frame>,

      gr::Tags::SqrtDetSpatialMetricCompute<DataVector, Dim, frame>,
      gr::Tags::WeylMagneticForGBCompute<DataVector, Dim, frame>,
      gr::Tags::WeylMagneticScalarCompute<DataVector, Dim, frame>>;

  using const_global_cache_tags = tmpl::list<
      gh::ConstraintDamping::Tags::DampingFunctionGamma0<Dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma1<Dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma2<Dim, Frame::Grid>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& /*box*/,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

// template <size_t Dim>
struct InitializeEvolvedScalarVariables {
 // using flat_variables_tag = typename ScalarWave::System<3_st>::variables_tag;
  // using curved_variables_tag =
  //     typename CurvedScalarWave::System<3_st>::variables_tag;
  using curved_variables_tag =
      typename ScalarTensor::System::variables_tag;
  using return_tags = tmpl::list<curved_variables_tag>;
  using argument_tags =
      tmpl::list<
      // ::Tags::Time, domain::Tags::Coordinates<3_st, Frame::Inertial>,
      //           //  ::Tags::AnalyticSolutionOrData,
                 gr::Tags::Lapse<DataVector>
                //  ,
      //            gr::Tags::Shift<DataVector, 3_st>
                 >;
  // template <typename AnalyticSolutionOrData>
  static void apply(
      const gsl::not_null<typename curved_variables_tag::type*> evolved_vars,
      // const double initial_time,
      // const tnsr::I<DataVector, 3_st>& inertial_coords,
      // // const AnalyticSolutionOrData& solution_or_data,
      [[maybe_unused]] const Scalar<DataVector>& lapse
      // ,
      // [[maybe_unused]] const tnsr::I<DataVector, 3_st>& shift
      ) {

    // Set variables to zero for now
    // Ideally we would like to
      get(get<CurvedScalarWave::Tags::Psi>(*evolved_vars)) = 0.0 * get(lapse);
      auto& scalar_phi =
                  get<CurvedScalarWave::Tags::Phi<3_st>>(*evolved_vars);
      for (size_t i = 0; i < 3_st; i++) {
        scalar_phi.get(i) = 0.0 * get(lapse);
      }
      get(get<CurvedScalarWave::Tags::Pi>(*evolved_vars)) =
                  - 1.0e-3 * (get(lapse) - 1.0);
    }
  };

}  // namespace ScalarTensor::Actions
