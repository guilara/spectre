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
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
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

struct InitializeScalarTensorAnd3Plus1Variables {
  static constexpr size_t dim = 3_st;
  using frame = Frame::Inertial;
  using compute_tags = db::AddComputeTags<
      // Needed to compute the characteristic speeds for the AH finder
      gr::Tags::SpatialMetricCompute<DataVector, dim, frame>,
      gr::Tags::DetAndInverseSpatialMetricCompute<DataVector, dim, frame>,
      gr::Tags::ShiftCompute<DataVector, dim, frame>,
      gr::Tags::LapseCompute<DataVector, dim, frame>,

      gr::Tags::SpacetimeNormalVectorCompute<DataVector, dim, frame>,
      gh::Tags::DerivLapseCompute<dim, frame>,

      gr::Tags::InverseSpacetimeMetricCompute<DataVector, dim, frame>,
      gh::Tags::DerivShiftCompute<dim, frame>,

      gh::Tags::DerivSpatialMetricCompute<dim, frame>,

      // Compute tags for Trace of Christoffel and Extrinsic curvature
      gr::Tags::SpatialChristoffelFirstKindCompute<DataVector, dim, frame>,
      gr::Tags::SpatialChristoffelSecondKindCompute<DataVector, dim, frame>,
      gr::Tags::TraceSpatialChristoffelSecondKindCompute<DataVector, dim,
                                                         frame>,
      gh::Tags::ExtrinsicCurvatureCompute<dim, frame>,
      gh::Tags::TraceExtrinsicCurvatureCompute<dim, frame>,

      // Compute constraint damping parameters.
      gh::ConstraintDamping::Tags::ConstraintGamma0Compute<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::ConstraintGamma1Compute<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::ConstraintGamma2Compute<dim, Frame::Grid>,

      ScalarTensor::Tags::ScalarSourceCompute>;

  using const_global_cache_tags = tmpl::list<
      gh::ConstraintDamping::Tags::DampingFunctionGamma0<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma1<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma2<dim, Frame::Grid>>;

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

namespace ScalarTensor::Initialization {

/// \ingroup InitializationGroup
/// \brief Initialize the compute tags required by the ScalarTensor system
///
/// DataBox changes:
/// - Adds:
///   * compute_tags
/// - Removes: nothing
/// - Modifies: nothing
///
/// \details The compute tags required include those specified in
/// ::gh::Actions::InitializeGhAnd3Plus1Variables as well as the tags required
/// to compute spacetime quantities appearing in the scalar evolution equations.
/// Namely, we add the compute tags associated to the trace of the extrinsic
/// curvature and the trace of the spatial Christoffel symbol, as well as the
/// compute tag required to compute the source term of the scalar equation.
template <typename Metavariables>
struct ScalarTensor3Plus1Variables {
  static constexpr size_t dim = Metavariables::volume_dim;
  using frame = Frame::Inertial;

  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using simple_tags = tmpl::list<>;

  using compute_tags = tmpl::list<
      // Needed to compute the characteristic speeds for the AH finder
      gr::Tags::SpatialMetricCompute<DataVector, dim, frame>,
      gr::Tags::DetAndInverseSpatialMetricCompute<DataVector, dim, frame>,
      gr::Tags::ShiftCompute<DataVector, dim, frame>,
      gr::Tags::LapseCompute<DataVector, dim, frame>,

      gr::Tags::SpacetimeNormalVectorCompute<DataVector, dim, frame>,
      gh::Tags::DerivLapseCompute<dim, frame>,

      gr::Tags::InverseSpacetimeMetricCompute<DataVector, dim, frame>,
      gh::Tags::DerivShiftCompute<dim, frame>,

      gh::Tags::DerivSpatialMetricCompute<dim, frame>,

      // Compute tags for Trace of Christoffel and Extrinsic curvature
      gr::Tags::SpatialChristoffelFirstKindCompute<DataVector, dim, frame>,
      gr::Tags::SpatialChristoffelSecondKindCompute<DataVector, dim, frame>,
      gr::Tags::TraceSpatialChristoffelSecondKindCompute<DataVector, dim,
                                                         frame>,
      gh::Tags::ExtrinsicCurvatureCompute<dim, frame>,
      gh::Tags::TraceExtrinsicCurvatureCompute<dim, frame>,

      // Compute constraint damping parameters.
      gh::ConstraintDamping::Tags::ConstraintGamma0Compute<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::ConstraintGamma1Compute<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::ConstraintGamma2Compute<dim, Frame::Grid>,

      ScalarTensor::Tags::ScalarSourceCompute>;

  using const_global_cache_tags = tmpl::list<
      gh::ConstraintDamping::Tags::DampingFunctionGamma0<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma1<dim, Frame::Grid>,
      gh::ConstraintDamping::Tags::DampingFunctionGamma2<dim, Frame::Grid>>;

  using argument_tags = tmpl::list<>;
  using return_tags = tmpl::list<>;

  static void apply() {}
};
}  // namespace ScalarTensor::Initialization
