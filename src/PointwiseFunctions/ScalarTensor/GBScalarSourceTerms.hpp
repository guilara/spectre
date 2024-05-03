// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

/*!
 * \brief Compute the source term of the scalar equation with the reduction of
 * order scheme.
 *
 * \details The computation is done as follows:
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
void order_reduced_gb_scalar_with_tenex(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3, Frame>& trace_reversed_stress_energy,
    const tnsr::AA<DataVector, 3, Frame>& inverse_spacetime_metric,
    const Scalar<DataVector>& weyl_electric_scalar_complement);

namespace Tags {

/*!
 * \brief Compute tag for the GB scalar with reduction of order.
 *
 * \details Call gb_scalar_with_tenex.
 *
 * The computation is done as follows:
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
 * \note Should replace with a function including the couplings.
 */
template <typename Frame>
struct OrderReducedGBScalarCompute : OrderReducedGBScalar, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::WeylElectricScalar<DataVector>,
      gr::Tags::WeylMagneticScalar<DataVector>,
      ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, 3, Frame>,
      gr::Tags::InverseSpacetimeMetric<DataVector, 3, Frame>,
      ScalarTensor::Tags::WeylElectricRicciScalarComplement<DataVector>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> result,
      const Scalar<DataVector>&, const Scalar<DataVector>&,
      const tnsr::aa<DataVector, 3, Frame>&,
      const tnsr::AA<DataVector, 3, Frame>&,
      const Scalar<DataVector>&) = &order_reduced_gb_scalar_with_tenex;
  using base = OrderReducedGBScalar;
};
}  // namespace Tags

}  // namespace ScalarTensor
