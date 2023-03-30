// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiZeroScalar.hpp"
#include "Utilities/TMPL.hpp"

namespace GeneralizedHarmonic::Solutions {
/// ScalarTensor solutions wrapped for GH
namespace ScalarTensor {
/// \brief List of all analytic solutions
using all_solutions = tmpl::list<GeneralizedHarmonic::Solutions::WrappedGr<
    ::ScalarTensor::Solutions::MinkowskiZeroScalar>>;
}  // namespace ScalarTensor
}  // namespace GeneralizedHarmonic::Solutions
