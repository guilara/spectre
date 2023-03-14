// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"

// Initialization::Actions::AddComputeTags<
//     StepChoosers::step_chooser_compute_tags<
// TemplateBase, local_time_stepping>,
//     GeneralizedHarmonic::Tags::ExtrinsicCurvatureCompute<volume_dim,
//                                                          ::Frame::Inertial>,
//     GeneralizedHarmonic::Tags::TraceExtrinsicCurvatureCompute<
//         volume_dim, ::Frame::Inertial>,
//     ::Tags::DerivTensorCompute<
//         gr::Tags::Lapse<DataVector>,
//         ::domain::Tags::InverseJacobian<volume_dim, ::Frame::ElementLogical,
//                                         ::Frame::Inertial>>,
//     ::Tags::DerivTensorCompute<
//         gr::Tags::Shift<volume_dim, Frame::Inertial, DataVector>,
//         ::domain::Tags::InverseJacobian<volume_dim, ::Frame::ElementLogical,
//                                         ::Frame::Inertial>>>;
