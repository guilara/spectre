// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

/*!
 * \brief Add in the source term to the \f$\Pi\f$
 * evolved variable of the ::CurvedScalarWave system.
 *
 * \details The only source term in the wave equation
 * \f[
 *  \Box \Psi = \mathcal{S} ~,
 *  \f]
 *
 * is in the equation for \f$\Pi\f$:
 * \f[
 *  \partial_t \Pi + \text{\{spatial derivative
 * terms\}} = \alpha \mathcal{S}
 * ~,
 * \f]
 *
 * where \f$\mathcal{S}\f$ is the source term (e. g. in the Klein-Gordon
 * equation, the source term is the derivative of the scalar potential
 * \f$\mathcal{S} \equiv \partial V / \partial \Psi \f$.)
 *
 * This function adds that contribution to the existing value of `dt_pi_scalar`.
 * The wave equation terms in the scalar equation should be computed before
 * passing the `dt_pi_scalar` to this function for updating.
 *
 * \param dt_pi_scalar Time derivative terms of $\Pi$. The sourceless part
 * should be computed before with ::CurvedScalarWave::TimeDerivative.
 * \param scalar_source Source term $\mathcal{S}$ for the scalar equation.
 * \param lapse Lapse $\alpha$.
 *
 * \see `CurvedScalarWave::TimeDerivative` for details about the source-less
 * part of the time derivative calculation.
 */
void add_scalar_source_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
    const Scalar<DataVector>& scalar_source, const Scalar<DataVector>& lapse);

/*!
 * \brief Computes the source term given by the mass of the scalar.
 *
 * \details For a scalar field with mass parameter \f$ m_\Psi \f$,
 * the wave equation takes the form
 * \f[
 *   \Box \Psi = \mathcal{S} ~,
 * \f]
 *
 * where the source is given by
 * \f[
 *   \mathcal{S} \equiv m^2_\Psi \Psi~.
 * \f]
 *
 * Here the mass parameter value is an option that needs to be specified in the
 * input file.
 *
 * \param scalar_source Source term $\mathcal{S}$ for the scalar equation.
 * \param psi Scalar field $\Psi$.
 * \param mass_psi Mass of the scalar field $m_\Psi$.
 *
 * \see `ScalarTensor::Tags::ScalarMass` for details about the mass.
 */
void mass_source(gsl::not_null<Scalar<DataVector>*> scalar_source,
                 const Scalar<DataVector>& psi, const double mass_psi);

void compute_scalar_curvature_source(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi);

Scalar<DataVector> compute_scalar_curvature_source(
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi);

void multiply_by_coupling_function_prime_quartic(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi);

// Extra functions for debugging
void compute_coupling_function_derivative(
    gsl::not_null<Scalar<DataVector>*> result, const Scalar<DataVector>& psi,
    const double first_coupling_psi, const double second_coupling_psi);

void compute_gb_scalar(gsl::not_null<Scalar<DataVector>*> gb_scalar,
                       const Scalar<DataVector>& weyl_electric_scalar,
                       const Scalar<DataVector>& weyl_magnetic_scalar);

void compute_rhs_psi(
    gsl::not_null<Scalar<DataVector>*> dt_psi,
    const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
    const tnsr::ij<DataVector, 3>& d_phi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3>& phi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2);

void compute_rhs_pi(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
    const tnsr::ij<DataVector, 3>& d_phi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3>& phi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& scalar_source);

void compute_rhs_phi(
    gsl::not_null<tnsr::i<DataVector, 3>*> dt_phi,
    const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
    const tnsr::ij<DataVector, 3>& d_phi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3>& phi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2);

namespace Tags {

/*!
 * \brief Compute tag for the scalar source.
 *
 * \details Compute the scalar source from data box items using
 * `mass_source`.
 */
struct ScalarSourceCompute : ScalarSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>&,
                                    const double) = &mass_source;
  using base = ScalarSource;
};

/*!
 * \brief Compute tag for the scalar source given by the background curvature.
 *
 * \details Call compute_scalar_curvature_source. Needs that WeylElectric is in
 * data box.
 */
template <typename DataType, size_t SpatialDim, typename Frame>
struct ScalarCurvatureSourceCompute : ScalarSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::WeylElectricScalar<DataType>,
                 gr::Tags::WeylMagneticScalar<DataType>,
                 CurvedScalarWave::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter,
                 ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
      const Scalar<DataVector>&, const Scalar<DataVector>&, const double,
      const double, const double) = &compute_scalar_curvature_source;
  using base = ScalarSource;
};

// Extra compute tags for debugging
/*!
 * \brief Compute tag for the coupling function.
 *
 * \details Call ....
 */
template <typename DataType>
struct CouplingFunctionDerivativeCompute : CouplingFunctionDerivative,
                                           db::ComputeTag {
  //   using argument_tags = tmpl::list<fe::ScalarDriver::Psi>;
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
      const double, const double) = &compute_coupling_function_derivative;
  using base = CouplingFunctionDerivative;
};

/*!
 * \brief Compute tag for the GB scalar.
 *
 * \details Call compute_gb_scalar.
 */
template <typename DataType>
struct GBScalarCompute : GBScalar, db::ComputeTag {
  //   using argument_tags = tmpl::list<fe::ScalarDriver::Psi>;
  using argument_tags = tmpl::list<gr::Tags::WeylElectricScalar<DataType>,
                                   gr::Tags::WeylMagneticScalar<DataType>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
      const Scalar<DataVector>&) = &compute_gb_scalar;
  using base = GBScalar;
};

/*!
 * \brief Compute the Rhs of the Psi equation.
 *
 * \details Call ....
 */
// template <typename DataType>
struct RhsPsiCompute : RhsPsi, db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<CurvedScalarWave::Tags::Pi, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3>,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>,
                 ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::Shift<DataVector, 3>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3>,
                 gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, 3>,
                 gr::Tags::TraceExtrinsicCurvature<DataVector>,
                 CurvedScalarWave::Tags::ConstraintGamma1,
                 CurvedScalarWave::Tags::ConstraintGamma2>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const tnsr::i<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::ij<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::iJ<DataVector, 3>&,
      const tnsr::II<DataVector, 3>&, const tnsr::I<DataVector, 3>&,
      const Scalar<DataVector>&, const Scalar<DataVector>&,
      const Scalar<DataVector>&) = &compute_rhs_psi;
  using base = RhsPsi;
};

/*!
 * \brief Compute the Rhs of the Pi equation.
 *
 * \details Call ....
 */
// template <typename DataType>
struct RhsPiCompute : RhsPi, db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<CurvedScalarWave::Tags::Pi, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3>,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>,
                 ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::Shift<DataVector, 3>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3>,
                 gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, 3>,
                 gr::Tags::TraceExtrinsicCurvature<DataVector>,
                 CurvedScalarWave::Tags::ConstraintGamma1,
                 CurvedScalarWave::Tags::ConstraintGamma2,
                 ScalarTensor::Tags::ScalarSource>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const tnsr::i<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::ij<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::iJ<DataVector, 3>&,
      const tnsr::II<DataVector, 3>&, const tnsr::I<DataVector, 3>&,
      const Scalar<DataVector>&, const Scalar<DataVector>&,
      const Scalar<DataVector>&, const Scalar<DataVector>&) = &compute_rhs_pi;
  using base = RhsPi;
};

/*!
 * \brief Compute the Rhs of the Phi equation.
 *
 * \details Call ....
 */
// template <typename DataType>
struct RhsPhiCompute : RhsPhi, db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<CurvedScalarWave::Tags::Pi, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<CurvedScalarWave::Tags::Phi<3>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3>,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>,
                 ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::Shift<DataVector, 3>, tmpl::size_t<3>,
                               Frame::Inertial>,
                 gr::Tags::InverseSpatialMetric<DataVector, 3>,
                 gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, 3>,
                 gr::Tags::TraceExtrinsicCurvature<DataVector>,
                 CurvedScalarWave::Tags::ConstraintGamma1,
                 CurvedScalarWave::Tags::ConstraintGamma2>;
  using return_type = tnsr::i<DataVector, 3>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const tnsr::i<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::ij<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const Scalar<DataVector>&, const tnsr::I<DataVector, 3>&,
      const tnsr::i<DataVector, 3>&, const tnsr::iJ<DataVector, 3>&,
      const tnsr::II<DataVector, 3>&, const tnsr::I<DataVector, 3>&,
      const Scalar<DataVector>&, const Scalar<DataVector>&,
      const Scalar<DataVector>&) = &compute_rhs_phi;
  using base = RhsPhi;
};

}  // namespace Tags
}  // namespace ScalarTensor
