// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/Protocols/Mutator.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/ConstraintDamping/ConstraintGammas.hpp"
#include "Evolution/Systems/ScalarTensor/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintGammas.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/DerivSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalOneForm.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/ProtocolHelpers.hpp"
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

    ScalarTensor::ConstraintDamping::Tags::ConstraintGamma1Compute<Dim,
                                                                   Frame::Grid>,
    ScalarTensor::ConstraintDamping::Tags::ConstraintGamma2Compute<Dim,
                                                                   Frame::Grid>,

    // ScalarTensor::Tags::ScalarSourceCompute>;
    ScalarTensor::Tags::ScalarCurvatureSourceCompute<DataVector, Dim, Fr>,

    // Extra tags for curvatures
    ::Tags::DerivTensorCompute<
        gr::Tags::ExtrinsicCurvature<DataVector, Dim, Fr>,
        ::domain::Tags::InverseJacobian<Dim, ::Frame::ElementLogical,
                                        ::Frame::Inertial>,
        ::domain::Tags::Mesh<Dim>>,
    gh::Tags::GradExtrinsicCurvatureCompute<Dim, Fr>,
    ::Tags::DerivTensorCompute<
        gr::Tags::SpatialChristoffelSecondKind<DataVector, Dim, Fr>,
        ::domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                        Frame::Inertial>,
        ::domain::Tags::Mesh<Dim>>,

    gr::Tags::SpatialRicciCompute<DataVector, Dim, Fr>,
    gr::Tags::SpatialRicciScalarCompute<DataVector, Dim, Fr>,

    gr::Tags::WeylElectricCompute<DataVector, Dim, Fr>,
    gr::Tags::WeylElectricScalarCompute<DataVector, Dim, Fr>,

    gr::Tags::SqrtDetSpatialMetricCompute<DataVector, Dim, Fr>,
    gr::Tags::WeylMagneticForGBCompute<DataVector, Dim, Fr>,
    gr::Tags::WeylMagneticScalarCompute<DataVector, Dim, Fr>>;

struct InitializeEvolvedScalarVariables
    : tt::ConformsTo<db::protocols::Mutator> {
  using curved_variables_tag = typename ScalarTensor::System::variables_tag;
  using return_tags = tmpl::list<curved_variables_tag>;
  using argument_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>,
                 //  CurvedScalarWave::Tags::ConstraintGamma2,
                 domain::Tags::Coordinates<3, Frame::Inertial>>;
  static void apply(
      const gsl::not_null<typename curved_variables_tag::type*> evolved_vars,
      [[maybe_unused]] const Scalar<DataVector>& lapse,
      // [[maybe_unused]] const Scalar<DataVector>& gamma2_scalar,
      const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords) {
    get(get<CurvedScalarWave::Tags::Psi>(*evolved_vars)) = 0.0 * get(lapse);
    auto& scalar_phi = get<CurvedScalarWave::Tags::Phi<3>>(*evolved_vars);
    for (size_t i = 0; i < 3; i++) {
      scalar_phi.get(i) = 0.0 * get(lapse);
    }

    const auto ones_scalar = make_with_value<Scalar<DataVector>>(lapse, 1.0);

    const double AhA_x = 9.0;
    const double AhB_x = -9.0;
    const double beta = 1.0;
    const double normalization_factor = beta / std::sqrt(2.0 * M_PI);

    const double AhA_sign = 1.0;
    const double AhB_sign = -1.0;

    auto argument_gaussian_A = make_with_value<Scalar<DataVector>>(lapse, 0.0);
    auto argument_gaussian_B = make_with_value<Scalar<DataVector>>(lapse, 0.0);
    get(argument_gaussian_A) = -square(
        (0.5 * beta) * (get<0>(inertial_coords) - AhA_x * get(ones_scalar)));
    get(argument_gaussian_B) = -square(
        (0.5 * beta) * (get<0>(inertial_coords) - AhB_x * get(ones_scalar)));
    for (size_t i = 1; i < 3; i++) {
      get(argument_gaussian_A) += -square(0.5 * beta * inertial_coords.get(i));
      get(argument_gaussian_B) += -square(0.5 * beta * inertial_coords.get(i));
    }

    auto& scalar_pi = get<CurvedScalarWave::Tags::Pi>(*evolved_vars);
    const size_t num_points = get_size(get(lapse));
    for (size_t i = 0; i < num_points; i++) {
      get(scalar_pi)[i] =
          AhA_sign * normalization_factor * exp(get(argument_gaussian_A)[i]) +
          AhB_sign * normalization_factor * exp(get(argument_gaussian_B)[i]);
    }

    // get(scalar_pi) = -1.0e-7 * get<0>(inertial_coords) * get(gamma2_scalar);
  }
};

}  // namespace ScalarTensor::Initialization
