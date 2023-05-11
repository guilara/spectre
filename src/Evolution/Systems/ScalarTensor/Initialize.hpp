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
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintGammas.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/DerivSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
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
      gr::Tags::SpatialMetricCompute<Dim, frame, DataVector>,
      gr::Tags::DetAndInverseSpatialMetricCompute<Dim, frame, DataVector>,
      gr::Tags::ShiftCompute<Dim, frame, DataVector>,
      gr::Tags::LapseCompute<Dim, frame, DataVector>,

      gr::Tags::SpacetimeNormalVectorCompute<Dim, frame, DataVector>,
      GeneralizedHarmonic::Tags::DerivLapseCompute<Dim, frame>,

      gr::Tags::InverseSpacetimeMetricCompute<Dim, frame, DataVector>,
      GeneralizedHarmonic::Tags::DerivShiftCompute<Dim, frame>,

      GeneralizedHarmonic::Tags::DerivSpatialMetricCompute<Dim, frame>,

      // Add trace compute tags for Christoffel and Extrinsic curvature

      gr::Tags::SpatialChristoffelFirstKindCompute<Dim, frame, DataVector>,
      gr::Tags::SpatialChristoffelSecondKindCompute<Dim, frame, DataVector>,
      gr::Tags::TraceSpatialChristoffelSecondKindCompute<Dim, frame,
                                                         DataVector>,
      GeneralizedHarmonic::Tags::ExtrinsicCurvatureCompute<Dim, frame>,
      GeneralizedHarmonic::Tags::TraceExtrinsicCurvatureCompute<Dim, frame>,

      // Compute constraint damping parameters.
      GeneralizedHarmonic::ConstraintDamping::Tags::ConstraintGamma0Compute<
          Dim, Frame::Grid>,
      GeneralizedHarmonic::ConstraintDamping::Tags::ConstraintGamma1Compute<
          Dim, Frame::Grid>,
      GeneralizedHarmonic::ConstraintDamping::Tags::ConstraintGamma2Compute<
          Dim, Frame::Grid>,
      ScalarTensor::Sources::Tags::ScalarSourceCompute>;

  using const_global_cache_tags = tmpl::list<
      GeneralizedHarmonic::ConstraintDamping::Tags::DampingFunctionGamma0<
          Dim, Frame::Grid>,
      GeneralizedHarmonic::ConstraintDamping::Tags::DampingFunctionGamma1<
          Dim, Frame::Grid>,
      GeneralizedHarmonic::ConstraintDamping::Tags::DampingFunctionGamma2<
          Dim, Frame::Grid>>;

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
}  // namespace ScalarTensor::Actions
