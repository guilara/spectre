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
 */
template <size_t SpatialDim, typename Frame, typename DataType>
void spacetime_riemann_tensor();

template <size_t SpatialDim, typename Frame, typename DataType>
tnsr::abcc<DataType, VolumeDim, Frame> spacetime_riemann_tensor();
/// @}
}  // namespace GeneralizedHarmonic
