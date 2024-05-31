// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/Protocols/Mutator.hpp"
#include "DataStructures/Tensor/EagerMath/Norms.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/ConstraintDamping/ConstraintGammas.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Sources.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/ConstraintDamping/ConstraintGammas.hpp"
#include "Evolution/Systems/ScalarTensor/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/StressEnergy.hpp"
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
#include "PointwiseFunctions/ScalarTensor/GBScalarSourceTerms.hpp"
#include "PointwiseFunctions/ScalarTensor/GBTensorSourceTerms.hpp"
#include "PointwiseFunctions/ScalarTensor/WeylElectricRicciPart.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::sgb::Initialization {

/// \brief List of compute tags to be initialized in the fixed ScalarTensor
/// system
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
    // fe::sgb::ConstraintDamping::Tags::ConstraintGamma1Compute<Dim,
    // Frame::Grid>,
    // fe::sgb::ConstraintDamping::Tags::ConstraintGamma2Compute<Dim,
    // Frame::Grid>,

    // Scalar Tensor Driver parameters
    fe::sgb::ConstraintDamping::Tags::SigmaParameterCompute<Dim, Frame::Grid>,
    fe::sgb::ConstraintDamping::Tags::TauParameterCompute<Dim, Frame::Grid>,

    // ScalarTensor::Tags::ScalarSourceCompute>;
    ScalarTensor::Tags::TraceReversedStressEnergyCompute,
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
    gr::Tags::WeylMagneticScalarCompute<DataVector, Dim, Fr>,

    // BR
    ScalarTensor::Tags::WeylElectricFullCompute<DataVector, Dim, Fr>,
    ScalarTensor::Tags::WeylElectricFullScalarCompute<DataVector, Dim, Fr>,
    ScalarTensor::Tags::OrderReducedGBScalarCompute<Fr>,
    ScalarTensor::Tags::RhsPsiCompute, ScalarTensor::Tags::RhsPiCompute,
    ScalarTensor::Tags::RhsPhiCompute,
    ScalarTensor::Tags::SpacetimeDerivScalarCompute<Fr>,
    ScalarTensor::Tags::nnDDKGCompute<Fr>,
    ScalarTensor::Tags::nsDDKGCompute<Fr>,
    ScalarTensor::Tags::ssDDKGCompute<Fr>,
    // ScalarTensor::Tags::OrderReducednnHCompute<Fr>,
    // ScalarTensor::Tags::SCrossBCompute<Fr>,
    // ScalarTensor::Tags::JCrossBCompute<Fr>,
    // ScalarTensor::Tags::OrderReducednsHCompute<Fr>,
    // ScalarTensor::Tags::OrderReducedssHCompute<Fr>,
    // ScalarTensor::Tags::OrderReducedHTensorCompute<Fr>,
    // ScalarTensor::Tags::OrderReducedQTensorCompute<Fr>,
    ScalarTensor::Tags::DDKGTensorCompute<Fr>,
    // ScalarTensor::Tags::OrderReducedHTensorRicciPartCompute<Fr>,
    // ScalarTensor::Tags::DDFPsiTensorCompute<Fr>,
    // ScalarTensor::Tags::nnDDFPsiCompute<Fr>,
    // ScalarTensor::Tags::nsDDFPsiCompute<Fr>,
    // ScalarTensor::Tags::ssDDFPsiCompute<Fr>,
    // ScalarTensor::Tags::OrderReducedTraceReversedStressEnergyCompute<Fr>,
    // Fixing diagnostics
    // ::Tags::PointwiseL2NormCompute<
    //     ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, Dim, Fr>>,
    // ::Tags::PointwiseL2NormCompute<
    //     ScalarTensor::Tags::OrderReducedTraceReversedStressEnergy>,

    // Tags for the scalar driver
    fe::ScalarTensorDriver::Tags::TargetTensorCompute<Fr, DataVector>,
    fe::ScalarTensorDriver::Tags::TargetScalarCompute<Fr, DataVector>,
    fe::ScalarTensorDriver::Tags::TensorDriverSourceCompute<Fr, DataVector>,
    fe::ScalarTensorDriver::Tags::ScalarDriverSourceCompute<Fr, DataVector>>;

struct InitializeEvolvedScalarVariables
    : tt::ConformsTo<db::protocols::Mutator> {
  using curved_variables_tag = typename fe::sgb::System::variables_tag;
  using return_tags = tmpl::list<curved_variables_tag>;
  using argument_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>,
                 domain::Tags::Coordinates<3, Frame::Inertial>>;
  static void apply(
      const gsl::not_null<typename curved_variables_tag::type*> evolved_vars,
      [[maybe_unused]] const Scalar<DataVector>& lapse,
      [[maybe_unused]] const tnsr::I<DataVector, 3, Frame::Inertial>&
          inertial_coords) {
    get(get<CurvedScalarWave::Tags::Psi>(*evolved_vars)) = 0.0 * get(lapse);
    auto& scalar_phi = get<CurvedScalarWave::Tags::Phi<3>>(*evolved_vars);
    for (size_t i = 0; i < 3; i++) {
      scalar_phi.get(i) = 0.0 * get(lapse);
    }
    const auto ones_scalar = make_with_value<Scalar<DataVector>>(lapse, 1.0);
    // get(get<CurvedScalarWave::Tags::Pi>(*evolved_vars)) =
    //     -1.0e-4 * get<0>(inertial_coords) * (get(lapse) - get(ones_scalar));
    get(get<CurvedScalarWave::Tags::Pi>(*evolved_vars)) =
        -1.0e-3 * (get(lapse) - get(ones_scalar));

    // Driver values
    get(get<fe::ScalarTensorDriver::Tags::Psi>(*evolved_vars)) =
        0.0 * get(lapse);
    get(get<fe::ScalarTensorDriver::Tags::PiScalar>(*evolved_vars)) =
        0.0 * get(lapse);

    auto& tensor_driver =
        get<fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, 3>>(
            *evolved_vars);
    auto& pi_tensor_driver =
        get<fe::ScalarTensorDriver::Tags::Pi<DataVector, 3>>(*evolved_vars);
    for (size_t a = 0; a < 4; a++) {
      for (size_t b = a; b < 4; b++) {
        tensor_driver.get(a, b) = 0.0 * get(lapse);
        pi_tensor_driver.get(a, b) = 0.0 * get(lapse);
      }
    }
  }
};

}  // namespace fe::sgb::Initialization
