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

namespace ScalarTensor::Initialization {

/// \brief List of compute tags to be initialized in the ScalarTensor system
///
/// \details The compute tags required include those specified in
/// ::gh::Actions::InitializeGhAnd3Plus1Variables as well as the tags required
/// to compute spacetime quantities appearing in the scalar evolution equations.
/// Namely, we include the compute tags associated to the trace of the extrinsic
/// curvature and the trace of the spatial Christoffel symbol, as well as the
/// compute tag required to calculate the source term of the scalar equation.
template <size_t Dim, typename Fr = Frame::Inertial>
using scalar_tensor_3plus1_compute_tags = tmpl::list<
    // Needed to compute the characteristic speeds for the AH finder
    gr::Tags::SpatialMetricCompute<DataVector, Dim, Fr>,
    gr::Tags::DetAndInverseSpatialMetricCompute<DataVector, Dim, Fr>,
    gr::Tags::ShiftCompute<DataVector, Dim, Fr>,
    gr::Tags::LapseCompute<DataVector, Dim, Fr>,

    gr::Tags::SpacetimeNormalVectorCompute<DataVector, Dim, Fr>,
    gh::Tags::DerivLapseCompute<Dim, Fr>,

    gr::Tags::InverseSpacetimeMetricCompute<DataVector, Dim, Fr>,
    gh::Tags::DerivShiftCompute<Dim, Fr>,

    gh::Tags::DerivSpatialMetricCompute<Dim, Fr>,

    // Compute tags for Trace of Christoffel and Extrinsic curvature
    gr::Tags::SpatialChristoffelFirstKindCompute<DataVector, Dim, Fr>,
    gr::Tags::SpatialChristoffelSecondKindCompute<DataVector, Dim, Fr>,
    gr::Tags::TraceSpatialChristoffelSecondKindCompute<DataVector, Dim, Fr>,
    gh::Tags::ExtrinsicCurvatureCompute<Dim, Fr>,
    gh::Tags::TraceExtrinsicCurvatureCompute<Dim, Fr>,

    // Compute constraint damping parameters.
    gh::ConstraintDamping::Tags::ConstraintGamma0Compute<Dim, Frame::Grid>,
    gh::ConstraintDamping::Tags::ConstraintGamma1Compute<Dim, Frame::Grid>,
    gh::ConstraintDamping::Tags::ConstraintGamma2Compute<Dim, Frame::Grid>,

    ScalarTensor::Tags::ScalarSourceCompute>;

}  // namespace ScalarTensor::Initialization

namespace ScalarTensor::Actions {

struct InitializeEvolvedScalarVariables {
  using curved_variables_tag = typename ScalarTensor::System::variables_tag;
  using return_tags = tmpl::list<curved_variables_tag>;
  using argument_tags = tmpl::list<gr::Tags::Lapse<DataVector>>;

  static void apply(
      const gsl::not_null<typename curved_variables_tag::type*> evolved_vars,
      [[maybe_unused]] const Scalar<DataVector>& lapse) {
    // Set variables to zero for now
    get(get<CurvedScalarWave::Tags::Psi>(*evolved_vars)) = 0.0 * get(lapse);
    auto& scalar_phi = get<CurvedScalarWave::Tags::Phi<3_st>>(*evolved_vars);
    for (size_t i = 0; i < 3_st; i++) {
      scalar_phi.get(i) = 0.0 * get(lapse);
    }
    get(get<CurvedScalarWave::Tags::Pi>(*evolved_vars)) =
        0.0 * (get(lapse) - 1.0);
  }
};

}  // namespace ScalarTensor::Actions
