// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
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
    gr::Tags::WeylMagneticScalarCompute<DataVector, Dim, Fr>,

    // Tags for the scalar driver
    fe::ScalarDriver::Tags::TargetPsiCompute<Fr, DataVector>,
    fe::ScalarDriver::Tags::ScalarDriverSourceCompute<Fr, DataVector>>;

}  // namespace ScalarTensor::Initialization
