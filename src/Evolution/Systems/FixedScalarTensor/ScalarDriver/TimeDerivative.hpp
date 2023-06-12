// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
#include "Evolution/Systems/FixedScalarTensor/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <typename X, typename Symm, typename IndexList>
class Tensor;
/// \endcond

namespace fe::ScalarDriver {

template <size_t Dim>
struct TimeDerivative {
 public:
  using temporary_tags = CurvedScalarWave::TimeDerivative::temporary_tags;

  using argument_tags =
      tmpl::list<Tags::Pi, Tags::Phi<Dim>, gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, Dim>,
                 ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<Dim>,
                               Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::Shift<DataVector, Dim>,
                               tmpl::size_t<Dim>, Frame::Inertial>,
                 gr::Tags::InverseSpatialMetric<DataVector, Dim>,
                 gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, Dim>,
                 gr::Tags::TraceExtrinsicCurvature<DataVector>,
                 Tags::ConstraintGamma1, Tags::ConstraintGamma2>;

  static void apply(
      gsl::not_null<Scalar<DataVector>*> dt_psi,
      gsl::not_null<Scalar<DataVector>*> dt_pi,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*> dt_phi,

      gsl::not_null<Scalar<DataVector>*> result_lapse,
      gsl::not_null<tnsr::I<DataVector, Dim>*> result_shift,
      gsl::not_null<tnsr::II<DataVector, Dim>*> result_inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> result_gamma1,
      gsl::not_null<Scalar<DataVector>*> result_gamma2,

      const tnsr::i<DataVector, Dim>& d_psi,
      const tnsr::i<DataVector, Dim>& d_pi,
      const tnsr::ij<DataVector, Dim>& d_phi, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, Dim>& phi, const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim>& shift,
      const tnsr::i<DataVector, Dim>& deriv_lapse,
      const tnsr::iJ<DataVector, Dim>& deriv_shift,
      const tnsr::II<DataVector, Dim>& upper_spatial_metric,
      const tnsr::I<DataVector, Dim>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2);
};

}  // namespace fe::ScalarDriver
