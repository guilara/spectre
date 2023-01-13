// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup ScalarTensorGroup
 * \brief Items related to scalar fields with backreaction in the metric.
 */
namespace ScalarTensor {

template <size_t Dim>
struct System {};

} // namespace ScalarTensor
