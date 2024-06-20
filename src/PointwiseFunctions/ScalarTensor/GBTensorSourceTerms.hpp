// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/StressEnergy.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

void spacetime_derivative_scalar(const gsl::not_null<tnsr::a<DataVector, 3>*>
                                     spacetime_derivative_scalar_result,
                                 const Scalar<DataVector>& lapse,
                                 const tnsr::I<DataVector, 3>& shift,
                                 const Scalar<DataVector>& pi_scalar,
                                 const tnsr::i<DataVector, 3>& phi_scalar);

void DDKG_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> DDKG_normal_normal_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    // Scalar quantities
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_pi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const Scalar<DataVector>& dt_pi_scalar,

    const tnsr::i<DataVector, 3>& d_lapse);

void DDKG_normal_spatial_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> DDKG_normal_spatial_result,

    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,

    // Scalar quantities
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_pi_scalar);

void DDKG_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> DDKG_spatial_spatial_result,

    // Metric quantities

    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const tnsr::Ijj<DataVector, 3>& spatial_christoffel_second_kind,

    // Scalar quantities
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::ij<DataVector, 3>& d_phi_scalar

    // Provide them with RHS compute tags or from dt<> prefixes
);

void trace_of_DDKG_spatial_spatial_projection_diagnostic(
    const gsl::not_null<Scalar<DataVector>*> result,
    const tnsr::ii<DataVector, 3>& ssDDKG,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    const tnsr::ij<DataVector, 3>& d_phi_scalar);

void DDKG_tensor_from_projections(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDKG_tensor_result,
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const Scalar<DataVector>& nnDDKG, const tnsr::i<DataVector, 3>& nsDDKG,
    const tnsr::ii<DataVector, 3>& ssDDKG);

// template <typename Frame>
void DDKG_tensor_from_projections(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDKG_tensor_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const tnsr::Ijj<DataVector, 3>& spatial_christoffel_second_kind,

    // Scalar quantities
    const Scalar<DataVector>& psi_scalar, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_psi_scalar,
    const tnsr::i<DataVector, 3>& d_pi_scalar,
    const tnsr::ij<DataVector, 3>& d_phi_scalar,

    // Provide them with RHS compute tags or from dt<> prefixes
    const Scalar<DataVector>& dt_psi_scalar,
    const Scalar<DataVector>& dt_pi_scalar,
    const tnsr::i<DataVector, 3>& dt_phi_scalar);

void DDKG_trace_minus_eom(
    const gsl::not_null<Scalar<DataVector>*> diagnostic,
    const tnsr::aa<DataVector, 3>& DDKG,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
    const Scalar<DataVector>& dt_pi_scalar,
    const Scalar<DataVector>& scalar_driver);

void DDKG_trace_minus_eom(
    const gsl::not_null<Scalar<DataVector>*> diagnostic,
    const tnsr::aa<DataVector, 3>& DDKG, const Scalar<DataVector>& nnDDKG,
    const tnsr::ii<DataVector, 3>& ssDDKG,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& dt_pi_scalar,
    const Scalar<DataVector>& scalar_driver);

// TODO: Compute projections as well
void DDFPsi_tensor_from_DDKG_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDFPsi_tensor_result,
    const tnsr::aa<DataVector, 3>& DDKG,
    const tnsr::a<DataVector, 3>& spacetime_derivative_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi);

void raise_indices_DDFPsi(
    const gsl::not_null<tnsr::AA<DataVector, 3>*> DDFPsiUpUp_tensor_result,
    const tnsr::aa<DataVector, 3>& DDFPsi,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric);

void DDFPsi_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> DDFPsi_normal_normal_result,
    // Scalar quantities
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    // DDKG projections
    const Scalar<DataVector>& DDKG_normal_normal_projection,
    // Coupling function parameters
    const double first_coupling_psi, const double second_coupling_psi);

void DDFPsi_spatial_normal_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> DDFPsi_spatial_normal_result,
    // Scalar quantities
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    // DDKG projections
    const tnsr::i<DataVector, 3>& DDKG_spatial_normal_projection,
    // Coupling function parameters
    const double first_coupling_psi, const double second_coupling_psi);

void DDFPsi_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> DDFPsi_spatial_spatial_result,
    // Scalar quantities
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    // DDKG projections
    const tnsr::ii<DataVector, 3>& DDKG_spatial_spatial_projection,
    // Coupling function parameters
    const double first_coupling_psi, const double second_coupling_psi);

void order_reduced_gb_H_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> nnH_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const tnsr::ii<DataVector, 3>& ssDDKG);

void compute_S_cross_B(
    const gsl::not_null<tnsr::i<DataVector, 3>*> S_cross_B_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& ssDDKG);

void compute_j_cross_B(
    const gsl::not_null<tnsr::ij<DataVector, 3>*> j_cross_B_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_magnetic,
    const tnsr::i<DataVector, 3>& nsDDKG);

void order_reduced_gb_H_normal_spatial_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> nsH_result,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const tnsr::i<DataVector, 3>& nsDDKG,
    const tnsr::i<DataVector, 3>& S_cross_B);

void order_reduced_gb_H_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> ssH_result,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3>& weyl_electric,
    const Scalar<DataVector>& nnDDKG, const tnsr::ii<DataVector, 3>& ssDDKG,
    const tnsr::ij<DataVector, 3>& j_cross_B, const Scalar<DataVector>& nnH);

void order_reduced_gb_H_tensor_weyl_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const Scalar<DataVector>& nnH, const tnsr::i<DataVector, 3>& nsH,
    const tnsr::ii<DataVector, 3>& ssH);

void order_reduced_Q_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> Q_tensor_result,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar);

// void order_reduced_gb_H_tensor_ricci_part(
//     const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
//     const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
//     const tnsr::aa<DataVector, 3>& T, const tnsr::aa<DataVector, 3>& DDKG);

void order_reduced_gb_H_tensor_ricci_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
    const tnsr::aa<DataVector, 3>& T, const Scalar<DataVector>& trace_T,
    const tnsr::AA<DataVector, 3>& DDKGUpUp);

void order_reduced_gb_H_tensor_ricci_part(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
    const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
    const tnsr::aa<DataVector, 3>& T, const Scalar<DataVector>& trace_T,
    const tnsr::AA<DataVector, 3>& DDKGUpUp,
    const tnsr::aa<DataVector, 3>& Tdriv,
    const Scalar<DataVector>& trace_Tdriv);

void order_reduced_trace_reversed_stress_energy(
    const gsl::not_null<tnsr::aa<DataVector, 3>*>
        order_reduced_trace_reversed_stress_energy_result,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
    const tnsr::aa<DataVector, 3>& gb_H_tensor_ricci_part,
    const tnsr::aa<DataVector, 3>& gb_H_tensor_weyl_part,
    const tnsr::aa<DataVector, 3>& trace_reversed_canonical_stress_energy);

/*
void order_reduced_gb_H_tensor(
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> result,
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric const
        tnsr::aa<DataVector, Dim>& spacetime_metric,
    const tnsr::A<DataVector, 3>& normal_spacetime_vector,
    const tnsr::a<DataVector, 3>& normal_spacetime_one_form,
    const tnsr::aa<DataVector, 3>& dd_coupling_function);
*/

namespace Tags {

/*!
 * \brief Compute the spacetime derivative of the scalar.
 */
template <typename Frame>
struct SpacetimeDerivScalarCompute : SpacetimeDerivScalar, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, 3, Frame>,
                 CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3>>;
  using return_type = tnsr::a<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::a<DataVector, 3>*>
          spacetime_derivative_scalar_result,
      const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
      const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3>& d_pi_scalar) = &spacetime_derivative_scalar;
  using base = SpacetimeDerivScalar;
};

/*!
 * \brief Compute tag for normal-normal projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct nnDDKGCompute : nnDDKG, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3, Frame>,
      gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
      CurvedScalarWave::Tags::Phi<3>,
      ::Tags::deriv<CurvedScalarWave::Tags::Pi, tmpl::size_t<3>, Frame>,
      ScalarTensor::Tags::RhsPi,
      ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>, Frame>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> DDKG_normal_normal_result,

      // Metric quantities
      const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
      const tnsr::II<DataVector, 3>& inverse_spatial_metric,
      // Scalar quantities
      const tnsr::i<DataVector, 3>& phi_scalar,

      // Scalar gradients
      const tnsr::i<DataVector, 3>& d_pi_scalar,

      // Provide them with RHS compute tags or from dt<> prefixes
      const Scalar<DataVector>& dt_pi_scalar,

      const tnsr::i<DataVector, 3>& d_lapse) = &DDKG_normal_normal_projection;
  using base = nnDDKG;
};

/*!
 * \brief Compute tag for normal-spatial projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct nsDDKGCompute : nsDDKG, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
      gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame>,
      CurvedScalarWave::Tags::Phi<3>,
      ::Tags::deriv<CurvedScalarWave::Tags::Pi, tmpl::size_t<3>, Frame>>;
  using return_type = tnsr::i<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> DDKG_normal_spatial_result,

      const tnsr::II<DataVector, 3>& inverse_spatial_metric,
      const tnsr::ii<DataVector, 3>& extrinsic_curvature,

      // Scalar quantities
      const tnsr::i<DataVector, 3>& phi_scalar,

      // Scalar gradients
      const tnsr::i<DataVector, 3>& d_pi_scalar) =
      &DDKG_normal_spatial_projection;
  using base = nsDDKG;
};

/*!
 * \brief Compute tag for normal-spatial projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct ssDDKGCompute : ssDDKG, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame>,
      gr::Tags::SpatialChristoffelSecondKind<DataVector, 3, Frame>,
      CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3>,
      ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>, Frame>>;
  using return_type = tnsr::ii<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ii<DataVector, 3>*> result,
      const tnsr::ii<DataVector, 3>&, const tnsr::Ijj<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const tnsr::ij<DataVector, 3>&) = &DDKG_spatial_spatial_projection;
  using base = ssDDKG;
};

/*!
 * \brief Compute diagnostic for DDKG tensor.
 */
template <typename Frame>
struct TraceOfSsDDKGTensorDiagnosticCompute : TraceOfSsDDKGTensorDiagnostic,
                                              db::ComputeTag {
  using argument_tags = tmpl::list<
      ScalarTensor::Tags::ssDDKG,
      gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
      gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, 3, Frame>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>, CurvedScalarWave::Tags::Pi,
      CurvedScalarWave::Tags::Phi<3>,
      ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>, Frame>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> result,
      const tnsr::ii<DataVector, 3>& ssDDKG,
      const tnsr::II<DataVector, 3>& inverse_spatial_metric,
      const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3>& phi_scalar,
      const tnsr::ij<DataVector, 3>& d_phi_scalar) =
      &trace_of_DDKG_spatial_spatial_projection_diagnostic;
  using base = TraceOfSsDDKGTensorDiagnostic;
};

/*!
 * \brief Compute tag for second covariant derivative of the scalar.
 */
template <typename Frame>
struct DDKGTensorCompute : DDKGTensor, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nnDDKG, ScalarTensor::Tags::nsDDKG,
                 ScalarTensor::Tags::ssDDKG>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> result,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>& shift,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) = &DDKG_tensor_from_projections;
  using base = DDKGTensor;
};

// /*!
//  * \brief Compute diagnostic for DDKG tensor.
//  */
// template <typename Frame>
// struct EomFromDDKGTensorDiagnosticCompute : EomFromDDKGTensorDiagnostic,
//                                             db::ComputeTag {
//   using argument_tags =
//       tmpl::list<ScalarTensor::Tags::DDKGTensor,
//                  gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>,
//                  ScalarTensor::Tags::RhsPi,
//                  fe::ScalarTensorDriver::Tags::Psi>;
//   using return_type = Scalar<DataVector>;
//   static constexpr void (*function)(
//       const gsl::not_null<Scalar<DataVector>*> diagnostic,
//       const tnsr::aa<DataVector, 3>& DDKG,
//       const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
//       const Scalar<DataVector>& dt_pi_scalar,
//       const Scalar<DataVector>& scalar_driver) = &DDKG_trace_minus_eom;
//   using base = EomFromDDKGTensorDiagnostic;
// };

/*!
 * \brief Compute diagnostic for DDKG tensor.
 */
template <typename Frame>
struct EomFromDDKGTensorDiagnosticCompute : EomFromDDKGTensorDiagnostic,
                                            db::ComputeTag {
  using argument_tags =
      tmpl::list<ScalarTensor::Tags::DDKGTensor, ScalarTensor::Tags::nnDDKG,
                 ScalarTensor::Tags::ssDDKG,
                 gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 ScalarTensor::Tags::RhsPi, fe::ScalarTensorDriver::Tags::Psi>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> diagnostic,
      const tnsr::aa<DataVector, 3>& DDKG, const Scalar<DataVector>& nnDDKG,
      const tnsr::ii<DataVector, 3>& ssDDKG,
      const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
      const tnsr::II<DataVector, 3>& inverse_spatial_metric,
      const Scalar<DataVector>& dt_pi_scalar,
      const Scalar<DataVector>& scalar_driver) = &DDKG_trace_minus_eom;
  using base = EomFromDDKGTensorDiagnostic;
};

/*!
 * \brief Compute tag for normal-normal projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct nnDDFPsiCompute : nnDDFPsi, db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<3>, ScalarTensor::Tags::nnDDKG,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> DDFPsi_normal_normal_result,
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3>& phi_scalar,
      const Scalar<DataVector>& DDKG_normal_normal_projection,
      const double first_coupling_psi,
      const double second_coupling_psi) = &DDFPsi_normal_normal_projection;
  using base = nnDDFPsi;
};

/*!
 * \brief Compute tag for spatial-normal projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct nsDDFPsiCompute : nsDDFPsi, db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<3>, ScalarTensor::Tags::nsDDKG,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter>;
  using return_type = tnsr::i<DataVector, 3>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> DDFPsi_spatial_normal_result,
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3>& phi_scalar,
      const tnsr::i<DataVector, 3>& DDKG_spatial_normal_projection,
      const double first_coupling_psi,
      const double second_coupling_psi) = &DDFPsi_spatial_normal_projection;
  using base = nsDDFPsi;
};

/*!
 * \brief Compute tag for spatial-normal projection of the second covariant
 * derivative of the scalar.
 */
template <typename Frame>
struct ssDDFPsiCompute : ssDDFPsi, db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<3>, ScalarTensor::Tags::ssDDKG,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter>;
  using return_type = tnsr::ii<DataVector, 3>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ii<DataVector, 3>*>
          DDFPsi_spatial_spatial_result,
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, 3>& phi_scalar,
      const tnsr::ii<DataVector, 3>& DDKG_spatial_spatial_projection,
      const double first_coupling_psi,
      const double second_coupling_psi) = &DDFPsi_spatial_spatial_projection;
  using base = ssDDFPsi;
};

/*!
 * \brief Compute tag for second covariant derivative of the scalar.
 */
template <typename Frame>
struct DDFPsiTensorCompute : DDFPsiTensor, db::ComputeTag {
  using argument_tags =
      tmpl::list<ScalarTensor::Tags::DDKGTensor,
                 ScalarTensor::Tags::SpacetimeDerivScalar,
                 CurvedScalarWave::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> DDFPsi_tensor_result,
      const tnsr::aa<DataVector, 3>& DDKG,
      const tnsr::a<DataVector, 3>& spacetime_derivative_scalar,
      const Scalar<DataVector>& psi, const double first_coupling_psi,
      const double second_coupling_psi) = &DDFPsi_tensor_from_DDKG_tensor;
  using base = DDFPsiTensor;
};

/*!
 * \brief Compute tag for second covariant derivative of the scalar.
 */
template <typename Frame>
struct DDFPsiUpUpTensorCompute : DDFPsiUpUpTensor, db::ComputeTag {
  using argument_tags =
      tmpl::list<ScalarTensor::Tags::DDFPsiTensor,
                 gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>>;
  using return_type = tnsr::AA<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::AA<DataVector, 3>*> DDFPsiUpUp_tensor_result,
      const tnsr::aa<DataVector, 3>& DDFPsi,
      const tnsr::AA<DataVector, 3>& inverse_spacetime_metric) =
      &raise_indices_DDFPsi;
  using base = DDFPsiUpUpTensor;
};

/*!
 * \brief Compute tag for normal-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducednnHCompute : OrderReducednnH, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 ScalarTensor::Tags::WeylElectricFull<DataVector, 3, Frame>,
                 ScalarTensor::Tags::ssDDFPsi>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> result,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) =
      &order_reduced_gb_H_normal_normal_projection;
  using base = OrderReducednnH;
};

/*!
 * \brief Compute S cross B.
 */
template <typename Frame>
struct SCrossBCompute : SCrossB, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::WeylMagnetic<DataVector, 3, Frame>,
                 ScalarTensor::Tags::ssDDFPsi>;
  using return_type = tnsr::i<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> result,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) = &compute_S_cross_B;
  using base = SCrossB;
};

/*!
 * \brief Compute S cross B.
 */
template <typename Frame>
struct JCrossBCompute : JCrossB, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::WeylMagnetic<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nsDDFPsi>;
  using return_type = tnsr::ij<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ij<DataVector, 3>*> result,
      const tnsr::II<DataVector, 3>&, const tnsr::ii<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&) = &compute_j_cross_B;
  using base = JCrossB;
};

/*!
 * \brief Compute tag for normal-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducednsHCompute : OrderReducednsH, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::SqrtDetSpatialMetric<DataVector>,
                 ScalarTensor::Tags::WeylElectricFull<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nsDDFPsi, ScalarTensor::Tags::SCrossB>;
  using return_type = tnsr::i<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::i<DataVector, 3>*> result,
      const tnsr::II<DataVector, 3>&, const Scalar<DataVector>&,
      const tnsr::ii<DataVector, 3>&, const tnsr::i<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&) =
      &order_reduced_gb_H_normal_spatial_projection;
  using base = OrderReducednsH;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedssHCompute : OrderReducedssH, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>,
                 gr::Tags::SqrtDetSpatialMetric<DataVector>,
                 ScalarTensor::Tags::WeylElectricFull<DataVector, 3, Frame>,
                 ScalarTensor::Tags::nnDDFPsi, ScalarTensor::Tags::ssDDFPsi,
                 ScalarTensor::Tags::JCrossB,
                 ScalarTensor::Tags::OrderReducednnH>;
  using return_type = tnsr::ii<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ii<DataVector, 3>*> result,
      const tnsr::ii<DataVector, 3>& spatial_metric,
      const tnsr::II<DataVector, 3>& inverse_spatial_metric,
      const Scalar<DataVector>& sqrt_det_spatial_metric,
      const tnsr::ii<DataVector, 3>& weyl_electric,
      const Scalar<DataVector>& nnDDKG, const tnsr::ii<DataVector, 3>& ssDDKG,
      const tnsr::ij<DataVector, 3>& j_cross_B, const Scalar<DataVector>& nnH) =
      &order_reduced_gb_H_spatial_spatial_projection;
  using base = OrderReducedssH;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedHTensorCompute : OrderReducedHTensor, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3, Frame>,
      ScalarTensor::Tags::OrderReducednnH, ScalarTensor::Tags::OrderReducednsH,
      ScalarTensor::Tags::OrderReducedssH>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> result,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>& shift,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const tnsr::ii<DataVector, 3>&) = &order_reduced_gb_H_tensor_weyl_part;
  using base = OrderReducedHTensor;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedQTensorCompute : OrderReducedQTensor, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3, Frame>,
                 ScalarTensor::Tags::WeylElectricFullScalar<DataVector>,
                 gr::Tags::WeylMagneticScalar<DataVector>>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> Q_tensor_result,
      const tnsr::aa<DataVector, 3>&, const Scalar<DataVector>&,
      const Scalar<DataVector>&) = &order_reduced_Q_tensor;
  using base = OrderReducedQTensor;
};

// /*!
//  * \brief Compute tag for spatial-spatial projection of the order reduced H
//  * tensor.
//  */
// template <typename Frame>
// struct OrderReducedHTensorRicciPartCompute : OrderReducedHTensorRicciPart,
//                                              db::ComputeTag {
//   using argument_tags = tmpl::list<
//       gr::Tags::SpacetimeMetric<DataVector, 3, Frame>,
//       gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>,
//       ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, 3, Frame>,
//       ScalarTensor::Tags::DDFPsiTensor>;
//   using return_type = tnsr::aa<DataVector, 3, Frame>;
//   static constexpr void (*function)(
//       const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
//       const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
//       const tnsr::aa<DataVector, 3>& T, const tnsr::aa<DataVector, 3>& DDKG)
//       = &order_reduced_gb_H_tensor_ricci_part;
//   using base = OrderReducedHTensorRicciPart;
// };

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedHTensorRicciPartCompute : OrderReducedHTensorRicciPart,
                                             db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::SpacetimeMetric<DataVector, 3, Frame>,
      gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>,
      ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, 3, Frame>,
      ScalarTensor::Tags::TraceOfTraceReversedStressEnergy<DataVector>,
      ScalarTensor::Tags::DDFPsiUpUpTensor,
      fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, 3, Frame>,
      fe::ScalarTensorDriver::Tags::TensorDriverTrace>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> gb_H_tensor_result,
      const tnsr::aa<DataVector, 3>& g, const tnsr::AA<DataVector, 3>& inv_g,
      const tnsr::aa<DataVector, 3>& T, const Scalar<DataVector>& trace_T,
      const tnsr::AA<DataVector, 3>& DDKGUpUp,
      const tnsr::aa<DataVector, 3>& Tdriv,
      const Scalar<DataVector>& trace_Tdriv) =
      &order_reduced_gb_H_tensor_ricci_part;
  using base = OrderReducedHTensorRicciPart;
};

/*!
 * \brief Compute tag for spatial-spatial projection of the order reduced H
 * tensor.
 */
template <typename Frame>
struct OrderReducedTraceReversedStressEnergyCompute
    : OrderReducedTraceReversedStressEnergy,
      db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::SpacetimeMetric<DataVector, 3, Frame>,
      gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>,
      ScalarTensor::Tags::OrderReducedHTensorRicciPart,
      ScalarTensor::Tags::OrderReducedHTensor,
      ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, 3, Frame>>;
  using return_type = tnsr::aa<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*>
          order_reduced_trace_reversed_stress_energy_result,
      const tnsr::aa<DataVector, 3>& spacetime_metric,
      const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
      const tnsr::aa<DataVector, 3>& gb_H_tensor_ricci_part,
      const tnsr::aa<DataVector, 3>& gb_H_tensor_weyl_part,
      const tnsr::aa<DataVector, 3>& trace_reversed_canonical_stress_energy) =
      &order_reduced_trace_reversed_stress_energy;
  using base = OrderReducedTraceReversedStressEnergy;
};

}  // namespace Tags

}  // namespace ScalarTensor

namespace fe::ScalarTensorDriver {

void tensor_driver_normal_normal_projection(
    const gsl::not_null<Scalar<DataVector>*> tensor_driver_nn,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::A<DataVector, 3>& spacetime_normal_vector);

void tensor_driver_spatial_normal_projection(
    const gsl::not_null<tnsr::i<DataVector, 3>*> tensor_driver_sn,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::A<DataVector, 3>& spacetime_normal_vector,
    const tnsr::a<DataVector, 3>& spacetime_normal_covector);

void tensor_driver_spatial_spatial_projection(
    const gsl::not_null<tnsr::ii<DataVector, 3>*> tensor_driver_ss,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::A<DataVector, 3>& spacetime_normal_vector,
    const tnsr::a<DataVector, 3>& spacetime_normal_covector);

namespace Tags {
// gr::Tags::SpacetimeNormalVector<DataVector, Dim>

/*!
 * \brief Compute tag for normal-normal projection of the tensor driver.
 */
template <typename Frame>
struct TensorDriverSpatialProjectionCompute : TensorDriverSpatialProjection,
                                              db::ComputeTag {
  using argument_tags = tmpl::list<
      fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, 3, Frame>,
      gr::Tags::SpacetimeNormalVector<DataVector, 3, Frame>,
      gr::Tags::SpacetimeNormalOneForm<DataVector, 3, Frame>>;
  using return_type = tnsr::ii<DataVector, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::ii<DataVector, 3>*> tensor_driver_ss,
      const tnsr::aa<DataVector, 3>& tensor_driver,
      const tnsr::A<DataVector, 3>& spacetime_normal_vector,
      const tnsr::a<DataVector, 3>& spacetime_normal_covector) =
      &tensor_driver_spatial_spatial_projection;
  using base = TensorDriverSpatialProjection;
};

/*!
 * \brief Compute tag for normal-normal projection of the tensor driver.
 */
template <typename Frame>
struct TensorDriverTraceCompute : TensorDriverTrace, db::ComputeTag {
  using argument_tags = tmpl::list<
      fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, 3, Frame>,
      gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> trace,
      const tnsr::aa<DataVector, 3>& tensor_driver,
      const tnsr::AA<DataVector, 3>& inverse_spacetime_metric) =
      &::ScalarTensor::trace_of_trace_reversed_stress_energy;
  using base = TensorDriverTrace;
};

}  // namespace Tags

}  // namespace fe::ScalarTensorDriver
