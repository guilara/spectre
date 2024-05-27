// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
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

struct TimeDerivative {
 public:
  using temporary_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3_st>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3_st>,
                 Tags::ConstraintGamma1, Tags::ConstraintGamma2>;

  using argument_tags =
      tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3_st>,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3_st>,
                 ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3_st>,
                               Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::Shift<DataVector, 3_st>,
                               tmpl::size_t<3_st>, Frame::Inertial>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3_st>,
                 gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, 3_st>,
                 gr::Tags::TraceExtrinsicCurvature<DataVector>,
                 Tags::ConstraintGamma1, Tags::ConstraintGamma2,
                 fe::ScalarDriver::Tags::ScalarDriverSource,
                 fe::ScalarDriver::Tags::ScalarTauParameter,
                 fe::ScalarDriver::Tags::ScalarSigmaParameter,
                 domain::Tags::MeshVelocity<3_st, Frame::Inertial>>;

  static void apply(
      gsl::not_null<Scalar<DataVector>*> dt_psi,
      gsl::not_null<Scalar<DataVector>*> dt_pi,
      gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> dt_phi,

      gsl::not_null<Scalar<DataVector>*> result_lapse,
      gsl::not_null<tnsr::I<DataVector, 3_st>*> result_shift,
      gsl::not_null<tnsr::II<DataVector, 3_st>*> result_inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> result_gamma1,
      gsl::not_null<Scalar<DataVector>*> result_gamma2,

      const tnsr::i<DataVector, 3_st>& d_psi,
      const tnsr::i<DataVector, 3_st>& d_pi,
      const tnsr::ij<DataVector, 3_st>& d_phi,
      const Scalar<DataVector>& /*psi*/, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, 3_st>& phi, const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st>& shift,
      const tnsr::i<DataVector, 3_st>& deriv_lapse,
      const tnsr::iJ<DataVector, 3_st>& deriv_shift,
      const tnsr::II<DataVector, 3_st>& upper_spatial_metric,
      const tnsr::I<DataVector, 3_st>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const Scalar<DataVector>& scalar_driver_source,
      const Scalar<DataVector>& scalar_tau_parameter,
      const Scalar<DataVector>& scalar_sigma_parameter,
      const std::optional<tnsr::I<DataVector, 3_st, Frame::Inertial>>&
          mesh_velocity) {
    // Use the definition from the CurvedScalarWave system
    // CurvedScalarWave::TimeDerivative<3_st>::apply(
    //     dt_psi, dt_pi, dt_phi,

    //     result_lapse, result_shift, result_inverse_spatial_metric,
    //     result_gamma1, result_gamma2,

    //     d_psi, d_pi, d_phi, pi, phi, lapse, shift, deriv_lapse, deriv_shift,
    //     upper_spatial_metric, trace_spatial_christoffel,
    //     trace_extrinsic_curvature, gamma1, gamma2);

    *result_lapse = lapse;
    *result_shift = shift;
    *result_inverse_spatial_metric = upper_spatial_metric;
    *result_gamma1 = gamma1;
    *result_gamma2 = gamma2;

    if (mesh_velocity.has_value()) {
      // Psi equation
      tenex::evaluate(dt_psi, -lapse() * pi() -
                                  mesh_velocity.value()(ti::I) * d_psi(ti::i));
    } else {
      // Psi equation
      tenex::evaluate(dt_psi, -lapse() * pi());
    }
    // Pi equation
    tenex::evaluate(dt_pi, shift(ti::I) * d_pi(ti::i));

    // Phi equation. Not needed so set to zero.
    for (size_t index = 0; index < 3_st; ++index) {
      dt_phi->get(index) = 0.0 * get(lapse) * phi.get(index);
    }

    // Add extra terms to the Klein-Gordon equation
    // Make sure all variables called here are in the arguments of apply
    // and in the DataBox
    Sources::add_scalar_driver_friction_term_to_dt_pi_scalar(
        dt_pi, pi, lapse, shift, scalar_tau_parameter, scalar_sigma_parameter);

    Sources::add_scalar_driver_source_to_dt_pi_scalar(
        dt_pi, scalar_driver_source, lapse);
  }
};

}  // namespace fe::ScalarDriver
