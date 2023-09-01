// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

namespace ScalarTensor {

template <typename Frame>
void gb_scalar(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector> weyl_electric_scalar,
    const Scalar<DataVector> weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3> trace_reversed_stress_energy,
    const tnsr::AA<DataType, SpatialDim, Frame> inverse_spacetime_metric);

/*!
 * \brief Compute the source term of the scalar equation with the reduction of
 * order scheme.
 *
 * \begin{align}
 *     \mathcal{G}_\text{GB}
 *          &= R_{abcd}R^{abcd} - 4 * R_{ab} * R^{ab} + R^2 \\
 *           &= 8 \left(E_{ab}E^{ab} - B_{ab} B^{ab} \right)
 *              -2 R_{ab}R^{ab} + \dfrac{2}{3} R^2
 *           &=  8 \left(E_{ab}E^{ab} - B_{ab} B^{ab} \right) \\
 *            &  -2 \kappa^2 T^{(\Psi, \mathrm{TR})}_{ab}
 *               T^{(\Psi, \mathrm{TR})\, ab}
 *               + \dfrac{2}{3} \kappa^2
 *               \left(T^{(\Psi, \mathrm{TR})}\right)^2
 *                + \mathcal{O}(\epsilon)~,
 * \end{align}
 *
 * where \f$ T^{(\Psi, \mathrm{TR})}_{ab}\f$ is the canonical trace reversed
 * stress energy tensor.
 */
template <typename Frame>
void gb_scalar_with_tenex(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector> weyl_electric_scalar,
    const Scalar<DataVector> weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3, Frame> trace_reversed_stress_energy,
    const tnsr::AA<DataType, 3, Frame> inverse_spacetime_metric);

namespace Tags {

/*!
 * \brief Compute tag for the trace reversed stress energy tensor.
 *
 * \details Call trace_reversed_stress_energy.
 */
template <typename Frame = Frame::Inertial>
struct GBScalarCompute : GBScalar, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::WeylElectricScalar<DataVector>,
                 gr::Tags::WeylMagneticScalar<DataVector>,
                 ScalarTensor::Tags::TraceReversedStressEnergy<Frame>,
                 gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> result,
      const Scalar<DataVector>&, const Scalar<DataVector>&,
      const tnsr::aa<DataVector, 3, Frame>&,
      const tnsr::AA<DataType, 3, Frame>&) = &gb_scalar_with_tenex<Frame>;
  using base = GBScalar;
};
}  // namespace Tags

}  // namespace ScalarTensor
