// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Trace.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace Tags
}  // namespace domain
class DataVector;
template <typename X, typename Symm, typename IndexList>
class Tensor;
/// \endcond

namespace gh {
/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Computes extrinsic curvature from generalized harmonic variables
 *        and the spacetime normal vector.
 *
 * \details If \f$ \Pi_{ab} \f$ and \f$ \Phi_{iab} \f$ are the generalized
 * harmonic conjugate momentum and spatial derivative variables, and if
 * \f$n^a\f$ is the spacetime normal vector, then the extrinsic curvature
 * is computed as
 * \f{align}
 *     K_{ij} &= \frac{1}{2} \Pi_{ij} + \Phi_{(ij)a} n^a
 * \f}
 */
template <typename DataType, size_t SpatialDim, typename Frame>
void extrinsic_curvature(
    gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*> ex_curv,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi);

template <typename DataType, size_t SpatialDim, typename Frame>
tnsr::ii<DataType, SpatialDim, Frame> extrinsic_curvature(
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi);
/// @}

template <typename DataType, size_t SpatialDim, typename Frame>
void grad_extrinsic_curvature(
    const gsl::not_null<tnsr::ijj<DataType, SpatialDim, Frame>*> grad_ex_curv,
    const tnsr::ijj<DataType, SpatialDim, Frame>& d_ex_curv,
    const tnsr::ii<DataType, SpatialDim, Frame>& ex_curv,
    const tnsr::Ijj<DataType, SpatialDim, Frame>&
        spatial_christoffel_second_kind);

namespace Tags {
/*!
 * \brief Compute item to get extrinsic curvature from generalized harmonic
 * variables and the spacetime normal vector.
 *
 * \details See `extrinsic_curvature()`. Can be retrieved using
 * `gr::Tags::ExtrinsicCurvature`.
 */
template <size_t SpatialDim, typename Frame>
struct ExtrinsicCurvatureCompute
    : gr::Tags::ExtrinsicCurvature<DataVector, SpatialDim, Frame>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpacetimeNormalVector<DataVector, SpatialDim, Frame>,
                 Pi<DataVector, SpatialDim, Frame>,
                 Phi<DataVector, SpatialDim, Frame>>;

  using return_type = tnsr::ii<DataVector, SpatialDim, Frame>;

  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<tnsr::ii<DataVector, SpatialDim, Frame>*>,
      const tnsr::A<DataVector, SpatialDim, Frame>&,
      const tnsr::aa<DataVector, SpatialDim, Frame>&,
      const tnsr::iaa<DataVector, SpatialDim, Frame>&)>(
      &extrinsic_curvature<DataVector, SpatialDim, Frame>);

  using base = gr::Tags::ExtrinsicCurvature<DataVector, SpatialDim, Frame>;
};

/*!
 * \brief Compute item to get the trace of extrinsic curvature from generalized
 * harmonic variables and the spacetime normal vector.
 *
 * \details See `extrinsic_curvature()` for how the extrinsic curvature
 * \f$ K_{ij}\f$ is computed. Its trace is taken as
 * \f{align}
 *     tr(K) &= g^{ij} K_{ij}.
 * \f}
 *
 * Can be retrieved using `gr::Tags::TraceExtrinsicCurvature`.
 */
template <size_t SpatialDim, typename Frame>
struct TraceExtrinsicCurvatureCompute
    : gr::Tags::TraceExtrinsicCurvature<DataVector>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::ExtrinsicCurvature<DataVector, SpatialDim, Frame>,
                 gr::Tags::InverseSpatialMetric<DataVector, SpatialDim, Frame>>;

  using return_type = Scalar<DataVector>;

  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<DataVector>*>,
                           const tnsr::ii<DataVector, SpatialDim, Frame>&,
                           const tnsr::II<DataVector, SpatialDim, Frame>&)>(
          &trace);

  using base = gr::Tags::TraceExtrinsicCurvature<DataVector>;
};

template <size_t SpatialDim, typename Frame>
struct GradExtrinsicCurvatureCompute
    : gr::Tags::GradExtrinsicCurvature<DataVector, SpatialDim, Frame>,
      db::ComputeTag {
  using argument_tags = tmpl::list<
     ::Tags::deriv<gr::Tags::ExtrinsicCurvature<DataVector, SpatialDim, Frame>,
                    tmpl::size_t<SpatialDim>, Frame>,
        gr::Tags::ExtrinsicCurvature<DataVector, SpatialDim, Frame>,
        gr::Tags::SpatialChristoffelSecondKind<DataVector, SpatialDim, Frame>>;

  using return_type = tnsr::ijj<DataVector, SpatialDim, Frame>;

  static constexpr auto function = static_cast<void (*)(
    const gsl::not_null<tnsr::ijj<DataVector, SpatialDim, Frame>*>,
    const tnsr::ijj<DataVector, SpatialDim, Frame>&,
    const tnsr::ii<DataVector, SpatialDim, Frame>&,
    const tnsr::Ijj<DataVector, SpatialDim, Frame>&
        spatial_christoffel_second_kind)>(
      &grad_extrinsic_curvature<DataVector, SpatialDim, Frame>);

  using base = gr::Tags::GradExtrinsicCurvature<DataVector, SpatialDim, Frame>;
};

}  // namespace Tags
}  // namespace gh
