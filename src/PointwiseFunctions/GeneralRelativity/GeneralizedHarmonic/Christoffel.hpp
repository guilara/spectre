// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace GeneralizedHarmonic {
/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Computes spatial Christoffel symbol of the 2nd kind from the
 * the generalized harmonic spatial derivative variable and the
 * inverse spatial metric.
 *
 * \details
 * If \f$ \Phi_{kab} \f$ is the generalized harmonic spatial derivative
 * variable \f$ \Phi_{kab} = \partial_k g_{ab}\f$ and \f$\gamma^{ij}\f$ is the
 * inverse spatial metric, the Christoffel symbols are
 * \f[
 *   \Gamma^m_{ij} = \frac{1}{2}\gamma^{mk}(\Phi_{ijk}+\Phi_{jik}-\Phi_{kij}).
 * \f]
 *
 * In the not_null version, no memory allocations are performed if the
 * output tensor already has the correct size.
 *
 */
template <size_t SpatialDim, typename Frame, typename DataType>
void christoffel_second_kind(
    const gsl::not_null<tnsr::Ijj<DataType, SpatialDim, Frame>*> christoffel,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi,
    const tnsr::II<DataType, SpatialDim, Frame>& inv_metric);

template <size_t SpatialDim, typename Frame, typename DataType>
auto christoffel_second_kind(
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi,
    const tnsr::II<DataType, SpatialDim, Frame>& inv_metric)
    -> tnsr::Ijj<DataType, SpatialDim, Frame>;
/// @}

/// @{
/*!
 * \brief Compute \f$\Gamma_a\f$ from the generalized harmonic evolved
 * variables.
 *
 * Starting from Eq. (40) of \cite Lindblom2005qh we get
 *
 * \f{align*}{
 *  \Gamma_a &= \gamma^{ij}\Phi_{ija} + n^b \Pi_{ab} -
 *   \frac{1}{2}\gamma^{i}{}_ag^{bc}\Phi_{ibc} - \frac{1}{2}n_ag^{bc}\Pi_{bc}
 *   \\
 *   \Gamma_a &= \gamma^{ij}\Phi_{ija} + n^b \Pi_{ab} -\frac{1}{2}
 *               n_a g^{bc} \left(n^i \Phi_{ibc} + \Pi_{bc}\right)
 *               - \frac{1}{2} \delta^i_a g^{bc} \Phi_{ibc}
 * \f}
 */
template <size_t SpatialDim, typename Frame, typename DataType>
tnsr::a<DataType, SpatialDim, Frame> trace_christoffel(
    const tnsr::a<DataType, SpatialDim, Frame>& spacetime_normal_one_form,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::AA<DataType, SpatialDim, Frame>& inverse_spacetime_metric,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi);

template <size_t SpatialDim, typename Frame, typename DataType>
void trace_christoffel(
    gsl::not_null<tnsr::a<DataType, SpatialDim, Frame>*> trace,
    const tnsr::a<DataType, SpatialDim, Frame>& spacetime_normal_one_form,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::AA<DataType, SpatialDim, Frame>& inverse_spacetime_metric,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi);
/// @}

namespace Tags {
template <size_t SpatialDim, typename Frame>
struct TraceSpatialChristoffelSecondKindCompute
    : gr::Tags::TraceSpatialChristoffelSecondKind<SpatialDim, Frame,
                                                  DataVector>,
      db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::SpacetimeNormalOneForm<SpatialDim, Frame, DataVector>,
      gr::Tags::SpacetimeNormalVector<SpatialDim, Frame, DataVector>,
      gr::Tags::InverseSpatialMetric<SpatialDim, Frame, DataVector>,
      gr::Tags::InverseSpacetimeMetric<SpatialDim, Frame, DataVector>,
      Pi<SpatialDim, Frame>, Phi<SpatialDim, Frame>>;

  using return_type = tnsr::I<DataVector, SpatialDim, Frame>;

  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<tnsr::a<DataVector, SpatialDim, Frame>*>,
      const tnsr::a<DataVector, SpatialDim, Frame>&,
      const tnsr::A<DataVector, SpatialDim, Frame>&,
      const tnsr::II<DataVector, SpatialDim, Frame>&,
      const tnsr::AA<DataVector, SpatialDim, Frame>&,
      const tnsr::aa<DataVector, SpatialDim, Frame>&,
      const tnsr::iaa<DataVector, SpatialDim, Frame>&)>(
      &trace_christoffel<SpatialDim, Frame, DataVector>);

  using base = gr::Tags::TraceSpatialChristoffelSecondKind<SpatialDim, Frame,
                                                           DataVector>;
};
} // namespace Tags

}  // namespace GeneralizedHarmonic
