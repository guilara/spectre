// Distributed under the MIT License.
// See LICENSE.txt for details.

///\file
/// Declares function templates to calculate the Riemann tensor

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"

/// \cond
namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl
/// \endcond

namespace GeneralizedHarmonic {
/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Compute spatial Riemann tensor using the equations of motion for the
 * evolved variables.
 * \details As a start, in Kerr-Shild coordinates the metric is time
 * independent, i.e. \f$ \partial_t u = 0 \f$.
 *
 * We compute the lower index 4-Riemann as
 * \f[
 * R_{abcd} =
 * \f]
 */
template <size_t SpatialDim, typename Frame, typename DataType>
void spacetime_riemann_tensor(
    gsl::not_null<tnsr::abcc<DataType, VolumeDim, Frame>*> riemann,
    const tnsr::aa<DataType, SpatialDim, Frame> pi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& phi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& deriv_pi,
    const tnsr::ijaa<DataType, VolumeDim, Frame>& deriv_phi,
    const tnsr::II<DataType, VolumeDim, Frame>& inverse_spatial_metric);

template <size_t SpatialDim, typename Frame, typename DataType>
tnsr::abcc<DataType, VolumeDim, Frame> spacetime_riemann_tensor(
    const tnsr::aa<DataType, SpatialDim, Frame> pi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& phi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& deriv_pi,
    const tnsr::ijaa<DataType, VolumeDim, Frame>& deriv_phi,
    const tnsr::II<DataType, VolumeDim, Frame>& inverse_spatial_metric);
/// @}
}  // namespace GeneralizedHarmonic
